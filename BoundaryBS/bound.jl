module bound

using juliboots
using various
using qfunc
using cb
using table
using main
using LP
using consts
using main
using JLD


function LPload_B(file::String;reduced=true,tablefile="None")

	lp=load(file,"lp")
	if !reduced return lp end

	tab= tablefile=="None" ? tab=lp.extra[1] : tablefile
	tmp_prob=setupLP_B(lp.extra[2],tab,ders=lp.extra[3]) #could add vectortypes too
	for (i,lpf) in enumerate(lp.lpFunctions)
		lpf.vecfunc=tmp_prob.lpFunctions[i].vecfunc
	end
	lp.extra=(tab,lp.extra[2],lp.extra[3]);
	lp
end


bf=BigFloat
setupLP_B(sig::T,file::String; ders="all") where {T<:Real}=(tab=table.loadTable(file); setupLP_B(tab,convert(BigFloat,sig),ders=ders,file=file))
setupLP_B(sigs::Array{T,1},file::String;ders="all") where {T<:Real}=setupLP_B(convert(Array{BigFloat,1},sigs),file,ders=ders,file=file)
setupLP_B(sigs::Array{BigFloat,1},file::String;ders="all")=(tab=table.loadTable(file); [setupLP_B(tab,s,ders=ders,file=file)::LP.LinearProgram{BigFloat} for s in sigs])

function setupLP_B(tab::table.Table,sigma::BigFloat; ders="all",file="N/A")

        println("Setting up LP...")
        Lmax=tab.Lmax

        eps=tab.eps
        vecfuncs=convTable_B(sigma,tab.table) # Convolve

        #----- allows one to choose which derivatives to work with
        #
        if ders!="all" vecfuncs=[vf[ders] for vf in vecfuncs] end
        #
        #---------------------------------------------------------

        dim0(L::Int)= L==0 ? maximum([eps+FUDGE,zerobf]) : L+2*eps    # unitarity bound
        zerop=qfunc.Polynomial([bf(0)])     # This is the cost associated with the vectors: zero

        spins= [0,1]

        #let's create a linear program based on this data
        trgt=-value(vecfuncs[2].vec,bf(0))    # the target: bulk identity vector
		
		
        lpVectorFuncs=[LP.LPVectorFunction(
        [dim0(0),DELTAMAX],vecfuncs[i].vec,zerop,i==1 ? "Boundary" : "Bulk")::LP.LPVectorFunction{BigFloat} for (i,L) in enumerate(spins)]


        dim=convert(Float64,2*eps+2)
        prob=LP.initLP(lpVectorFuncs,trgt,"Bootstrap Bound\nD = $dim\tsigma=$(convert(Float64,sigma))\t m = $(tab.mmax)",extra=(file,mcopy(sigma),ders))
        
		tmp=main.cullpoles(prob,CULLPOLES)
		
		bndid=LP.LPVector(tmp.lpFunctions[1],bf(0))    # boundary identity, must be included. We cannot get it directly because of an annoying order of limits issue.
		bndid.vector=-[(-1)^k*pochhammer(-sigma/2,k) for k=0:(length(bndid.vector)-1)];
		bndid.label=(bf(0),"Boundary - Id")
		push!(tmp.lpVectors,bndid)
        println("Done")
        return tmp
end

############################################################################
#
# Convolution
#
############################################################################

function xi_convCoeffs(i::Int64,sigma::Real,tup::Tuple{Int64,Int64}) #gives a dictionnary with the coefficients
                                                     # necessary to convolve with to obtain the
                                                     # derivatives of u^sigma G
        m=tup[1] #rho derivatives       
		n=tup[2] #this is just 0.
		m=convert(BigInt,m)
		n=convert(BigInt,n)
        d=sigma

		#derivatives of \xi^\Delta_sigma\equiv (\Delta_1+\Delta_2)/2
		
        if i==1  
			return Dict((k,l)=>-(-1)^(m-k)*binomial(m,k)*pochhammer(-d/2,m-k) for k=0:m, l=0:n)  #derivatives of xi^{\Delta/2}. Yes, sign of d is correct.
		else
			return Dict( (k,l)=>(-1)^(m-k)*binomial(m,k)*pochhammer(d/2,m-k) for k=0:m, l=0:n)
		end
end


function convTable_B(sigma::BigFloat,tab::Array{CBVec_Q{BigFloat},1})

        dict=tab[1].dict
        ders=sort(collect(keys(dict)))

        #all coefficients required

        mmax=maximum([ders[i][1] for i in 1:length(ders)])
        nmax=maximum([ders[i][2] for i in 1:length(ders)])
        maxders=(mmax,nmax)


        convtab=[ConvVec_Q(tab[i],sigma,"All")::ConvVec_Q for i=1:length(tab)] 
		mcopy(convtab[2].rho,sqrt(convtab[2].rho))  # bulk blocks come with a (4 rho)^\Delta/2 prefactor, not (4rho)^\Delta

        ct=1;

        for (m,n) in ders

                #if m%2==0 && sign==-1 continue end #skip even m, since they lead to zeros for components
                #if m%2==1 && sign==1 continue end #skip odd m, since they lead to zeros for components

                for i=1:length(tab)
                        convtab[i].dict[(m,n)]=ct
                end
                ct+=1

                
                comps=[(k,l) for k=0:m, l=0:n]

                for i=1:2 
                    cbvec=tab[i] #vector of derivatives, spin i
					tmp=mcopy(cbvec[1]) #to hold temporary results
					coeffs=xi_convCoeffs(i,sigma,(m,n))  #i=1, boundary, i=2, bulk. We convolve both, xi^{-\Delta/2}(1+fbulk)-xi^{\Delta/2}(fboundary)=0
                    for (k,der) in enumerate(comps)   # do convolution: go through derivatives, multiply them by appropriate coefficients, add them up, put them into convtab.
                        getindex(tmp,cbvec,der) #get component der from cbvec, place it in tmp                                                
						mmult(tmp,coeffs[der])  #sum rule is 1+bulk guys-bound guys=0
                        if k==1
                            push!(convtab[i].vec,mcopy(tmp.func))                            
                        else
                            mplus(convtab[i].vec[end],tmp.func)
                        end
                    end
                end
        end

        return convtab
end


end



