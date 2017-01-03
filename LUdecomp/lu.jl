module lu

import Base: dot, deepcopy
using consts
using various
import various: tofloat, deepcopy, mcopy
export LUdata

type LUdata{T<:Real}

luMat::Array{T,2}
perm::Array{Int64,1}
scales::Array{T,1}
end

mcopy(lu::LUdata)=LUdata(mcopy(lu.luMat),deepcopy(lu.perm),mcopy(lu.scales))
mcopy(o::LUdata,lu::LUdata)=(mcopy(o.luMat,lu.luMat); o.perm=deepcopy(lu.perm); mcopy(o.scales,lu.scales); o)
tofloat(lu::LUdata)=LUdata{Float64}(convert(Array{Float64,2},lu.luMat),lu.perm,convert(Array{Float64,1},lu.scales))
LUdata{T<:Real}(A::Array{T,2})=LUdata!(mcopy(A))
LUdata(O::Array{BigFloat,2},A::Array{BigFloat,2})=LUdata!(mcopy(O,A))

function LUdataBACKUP!(A::Array{BigFloat,2})

        tmp=BigFloat(0)
        n=stride(A,2)


        perm=[1:n]


        for k=1:n-1


            #p=findfirst(x->abs(x)>LU_FUDGE,A[k:n,k]) #find pivot - pivoting to get rid of zeros
			#pivot, p=findmax(map(abs,A[k:n,k])) #find pivot: partial pivoting for num stability
			# find pivot: scaled partial pivoting for num stability
			si=[findmax(map(abs,A[i,1:n]))[1]::BigFloat for i=k:n]                                     
			p=findmax(map(abs,A[k:n,k]./si))[2]
			pivot=abs(A[k:n,k][p])   #
			##########################
			
			
            if pivot<LU_FUDGE error("Singular matrix") end
			#if p==0 error("Singular matrix") end
            p=p+k-1

            if p!=k
                (perm[k],perm[p])=(perm[p],perm[k]) #swap
                (A[k,:],A[p,:])=(A[p,:],A[k,:]) #swap rows
            end
            if abs(A[k,k])<10*LU_FUDGE println("Small pivot: $(A[k,k])") end #DEBUGs
            for i=k+1:n
                mdiv(A[i,k],A[k,k])
                #A[i,k]=A[i,k]/A[k,k]
                for j=k+1:n
                    mmult(tmp,A[i,k],A[k,j])
                    msub(A[i,j],tmp)
                    #A[i,j]-=A[i,k]*A[k,j]
                end
            end
        end


        return LUdata(A,perm)
end

#Added rescalings
function LUdata!(A::Array{BigFloat,2})

        tmp=BigFloat(0)
        n=stride(A,2)
		

        perm=collect(1:n)
		scales=[mcopy(maximum(map(abs,A[i,:])))::BigFloat for i=1:n]                                     
		
		for i=1:n
			for j=1:n
			mdiv(A[i,j],scales[i])
			end
		end
		
		
        for k=1:n-1

			#si=scales[perm]
            #p=findfirst(x->abs(x)>LU_FUDGE,A[k:n,k]) #find pivot - pivoting to get rid of zeros
			#pivot, p=findmax(map(abs,A[k:n,k])) #find pivot: partial pivoting for num stability
			# find pivot: scaled partial pivoting for num stability
			#scales=[findmax(map(abs,A[i,1:n]))[1]::BigFloat for i=k:n]                                     
			p=findmax(map(abs,A[k:n,k]))[2]   #partial pivoting
			pivot=abs(A[k:n,k][p])   #
			##########################
			
			
            if pivot<LU_FUDGE error("Singular matrix") end
			#if p==0 error("Singular matrix") end
            p=p+k-1

            if p!=k
                (perm[k],perm[p])=(perm[p],perm[k]) #swap
                (A[k,:],A[p,:])=(A[p,:],A[k,:]) #swap rows				
            end
            if abs(A[k,k])<10*LU_FUDGE println("Small pivot: $(A[k,k])") end #DEBUGs
            for i=k+1:n
                mdiv(A[i,k],A[k,k])
                #A[i,k]=A[i,k]/A[k,k]
                for j=k+1:n
                    mmult(tmp,A[i,k],A[k,j])
                    msub(A[i,j],tmp)
                    #A[i,j]-=A[i,k]*A[k,j]
                end
            end
        end


        return LUdata(A,perm,scales)
end





######################################################################
#
#   Methods
#
######################################################################

# We have PA=LU; where P permutes columns of A. Hence:
# A.x=b => LU x=P b; solve L y= P.b first for y; then Ux=y for x


function dot{T,S}(lu::LUdata{T},b::Array{S,1})  # this should be understood as (A^-1).v !

        A=lu.luMat
        n=stride(A,2)
        if n!=length(b) error("Incompatible dimensions") end
                #Solve Ly=P b
				si=lu.scales[lu.perm]
                Pb=b[lu.perm]./si
                y=Array(T,n)

                for i=1:n
                    y[i]=Pb[i]-sum([(A[i,j]y[j])::T for j=1:i-1])
                end


                #Solve for x
                x=Array(T,n)

                for i=n:-1:1
                    x[i]=(y[i]-sum([(A[i,j]*x[j])::T for j in i+1:n]))/A[i,i]
                end

                return x

end


dot(lu::LUdata{BigFloat},b::Array{BigFloat,1})=(o=[BigFloat(0) for i=1:length(b)]; dot(o,lu,b))

function dot(o::Array{BigFloat,1},lu::LUdata{BigFloat},b::Array{BigFloat,1})  # this should be understood as o=(A^-1).v !

        A=lu.luMat
        n=stride(A,2)
        tmp1=mcopy(zerobf); tmp2=mcopy(zerobf); #two temporary variables initialized at zero

        if n!=length(b) error("Incompatible dimensions") end

        #Solve Ly=P b
		si=lu.scales[lu.perm]
        Pb=b[lu.perm]./si


        for i=1:n
            #y[i]=Pb[i]-sum([(A[i,j]y[j])::T for j=1:i-1])
            mcopy(tmp2,zerobf);
            for j=1:i-1
               mmult(tmp1,A[i,j],o[j])
               mplus(tmp2,tmp1)
            end
            msub(o[i],Pb[i],tmp2)

        end


                #Solve for x


        for i=n:-1:1
            #x[i]=(y[i]-sum([(A[i,j]*x[j])::T for j in i+1:n]))/A[i,i]
            mcopy(tmp2,zerobf)
            for j=i+1:n
                    mmult(tmp1,A[i,j],o[j])
                    mplus(tmp2,tmp1)
            end
            msub(o[i],tmp2)
            mdiv(o[i],A[i,i])
        end

        return o

end




# We have PA=LU; where P permutes columns of A. Hence:
# x.A=b want to solve for x.  Then x.(P^-1 LU)=b. Solve y.U=b first for y. Then x(P^-1).L)=y for x


dot(b::Array{BigFloat,1},lu::LUdata{BigFloat})=(o=[BigFloat(0) for i=1:length(b)]; dot(o,b,lu))

function dot(o::Array{BigFloat,1},b::Array{BigFloat,1},lu::LUdata{BigFloat})  # this should be understood as v.(A^-1) !

        A=lu.luMat
        n=stride(A,2)
		
        tmp1=mcopy(zerobf); tmp2=mcopy(zerobf); #two temporary variables initialized at zero
        if n!=length(b) error("Incompatible dimensions") end

                for i=1:n
                    #y[i]=(b[i]-sum([(A[j,i]y[j])::T for j=1:i-1]))/A[i,i]
                    mcopy(tmp2,zerobf);
                    for j=1:i-1
                        mmult(tmp1,A[j,i],o[j])
                        mplus(tmp2,tmp1)
                    end
                    msub(o[i],b[i],tmp2)
                    mdiv(o[i],A[i,i])
                end


                #Solve for x
                #x=Array(T,n)

                for i=n:-1:1
                    #x[i]=(y[i]-sum([(A[j,i]*x[j])::T for j in i+1:n]))
                    mcopy(tmp2,zerobf)
                    for j=i+1:n
                        mmult(tmp1,A[j,i],o[j])
                        mplus(tmp2,tmp1)
                    end
                    msub(o[i],tmp2)
                end

                ipermute!(o,lu.perm)    #inverse permutation
				for i=1:n
				mdiv(o[i],lu.scales[i])
				end
				o

end








end
