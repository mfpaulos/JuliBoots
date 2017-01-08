module main

using various
using consts
using JLD

import bb
import qfunc
import cb
import Base: show
using LP
using table
export chooseTable, setupLP, bissect, bissect!,bissect_pair, value, dropOdd!, changeTarget!,dropEven!,opemax, avgSpec
export filter,filter!,iterate!,status,cost,solution,makeVector,resume_opemax # LP routines, this makes them accessible when 'using main'
export saveresults


bf=BigFloat
#push!(LOAD_PATH,"C:\\Users\\Miguel_Paulos\\Google\ Drive\ N\\Julia\ Code")

######## TYPES #############################

type bissect_pair{T<:Real}
	last_sol::LinearProgram{T}
	last_func::LinearProgram{T}
	upper::T
	lower::T
	criteria::LP.LabelF
end
	
show(io::IO,bp::bissect_pair)=println(io,"Bissect Pair: $(bp.criteria), ($(tofloat(bp.upper)),$(tofloat(bp.lower))) - accuracy: $(tofloat(bp.upper-bp.lower))")

	
############################################




function chooseTable()
    a=map(chomp,readlines(`ls ../Tables`))
    println("Files in \'./Tables\':")
    for i=1:length(a)
            println("[$i]\t - \t $(a[i]) ")
    end
    println("\nChoose a table to load:")
    nr=int(chomp(readline(STDIN)))
    return "./Tables/$(a[nr])"
end



function cullpoles{T<:Real}(lp0::LP.LinearProgram{T},cutoff::T)
	lp=mcopy(lp0)
	cullpoles(lp.lpFunctions,cutoff)
	return lp
end

function cullpoles{T<:Real}(lpfs::Array{LP.LPVectorFunction{T},1},cutoff::T) #Absolute cutoff, would be nice to implement relative

    
    for lpf in lpfs
        smallguys=Array(Array{Int64,1},0)
            for qf in lpf.vecfunc
                push!(smallguys,find(x->abs(x.coeff)<cutoff,qf.poles))
            end

        common=reduce(intersect,smallguys)

        for qf in lpf.vecfunc
                bigguys=setdiff(1:length(qf.poles),common)
                qf.poles=[qf.poles[i]::qfunc.Pole{BigFloat} for i in bigguys]

        end
    end

    return lpfs
end

############################################################
#                                                          #
#  Setting up Linear Programs                              #
#                                                          #
############################################################


# ======= Setting up a basic LP (no global symmetry) ===========


# Utility for dropping odd spins from a Linear Program
function dropOdd!(prob::LinearProgram)
        ll=length(prob.lpFunctions)
        prob.lpFunctions=[prob.lpFunctions[2i-1] for i=1:(floor(ll/2)+1)]
        prob
end

function dropEven!(prob::LinearProgram)
        ll=length(prob.lpFunctions)
        prob.lpFunctions=[prob.lpFunctions[2i] for i=1:floor(ll/2)]
        prob
end


# Setup LP functions

setupLP{T<:Real}(sig::T,file::String; ders="all")=(tab=table.loadTable(file); setupLP(tab,convert(BigFloat,sig),file=file,ders=ders))
setupLP{T<:Real}(sigs::Array{T,1},file::String;ders="all")=setupLP(convert(Array{BigFloat,1},sigs),file=file,ders=ders)
setupLP(sigs::Array{BigFloat,1},file::String;ders="all")=(tab=table.loadTable(file); [setupLP(tab,s,ders=ders,file=file)::LP.LinearProgram{BigFloat} for s in sigs])

function setupLP()

       file=chooseTable()
       println("\nInitial Sigma:")
       s0=BigFloat(chomp(readline()))
       println("Final Sigma:")
       sf=BigFloat(chomp(readline()))
       println("How many points:")
       nr=int(chomp(readline()))

       sigmas= nr==1 ? [s0] : [s0:(abs(sf-s0)/(nr-1)):sf]
       return setupLP(sigmas,file)

end

function setupLP(tab::table.Table,sigma::BigFloat;file="N/A", ders="all")

        println("Setting up LP...")
        Lmax=tab.Lmax

        eps=tab.eps
        vecfuncs=cb.convTable(sigma,tab.table,-1) # Convolve

        #----- allows one to choose which derivatives to work with
        #
        if ders!="all" vecfuncs=[vf[ders] for vf in vecfuncs] end
        #
        #---------------------------------------------------------

        dim0(L::Int)= L==0 ? maximum([eps+FUDGE,zerobf]) : L+2*eps    # unitarity bound
        zerop=qfunc.Polynomial([bf(0)])     # This is the cost associated with the vectors: zero

        spins= tab.OddL ? collect(0:1:Lmax) : collect(0:2:Lmax)

        #let's create a linear program based on this data
        trgt=-value(vecfuncs[1].vec,bf(0))    # the target: identity vector

        lpVectorFuncs=[LP.LPVectorFunction(
        [dim0(L),DELTAMAX],vecfuncs[i].vec,zerop,"L=$L")::LP.LPVectorFunction{BigFloat} for (i,L) in enumerate(spins)]


        dim=convert(Float64,2*eps+2)
        prob=LP.initLP(lpVectorFuncs,trgt,"Basic Bound\nD = $dim\tsigma=$(convert(Float64,sigma))\t(m,n) = $((tab.mmax,tab.nmax))\tLmax = $(tab.Lmax)\tOdd spins: $(tab.OddL)",extra=(file,mcopy(sigma)))        
        tmp=cullpoles(prob,CULLPOLES)

        # The following normalizes the components by dividing by the scalar with \Delta=1
        #----------
#        val=main.value(tmp.lpFunctions[1].vecfunc,BigFloat(1))
#        val=[1/v for v in val]

#        main.mmult(tmp.target,val)
#        for lpf in tmp.lpFunctions
#                main.mmult(lpf.vecfunc,val)
#        end
#        for vec in tmp.lpVectors
#                main.mmult(vec.vector,val)
#        end
#        for vec in tmp.solVecs
#                main.mmult(vec.vector,val)
#        end
#        tmp.invA=LP.LPInverse(tmp.invA,LP.getA(tmp))
#
#        LP.updateCoeffs!(tmp)
#        LP.updateFunctional!(tmp)
        #----------

        println("Done")
        return tmp
end

#------ Routine for updating the target of a Linear Program (has to be done when solution has only auxiliary vectors)

changeTarget!{T<:Real}(lp::LinearProgram{T},targetvec::LP.LPVector{T})=changeTarget!(lp,targetvec.vector)
function changeTarget!{T<:Real}(lp::LinearProgram{T},targetvec::Array{T,1})

        labels=[v.label[2] for v in lp.solVecs]
        for l in labels
            if l!="AUX"
                     println("Can only act on Linear Programs with purely AUX type vectors.")
                     return
            end
        end

        mcopy(lp.target,targetvec)

        for (i,comp) in enumerate(targetvec)
                lp.solVecs[i].vector[i]= comp<0 ? -BigFloat(1) : BigFloat(1)
                mcopy(lp.lpVectors[i].vector,lp.solVecs[i].vector)
        end
        updateInverse!(lp)
        updateCoeffs!(lp)
end




#----- GLOBAL SYMMETRY STUFF

function buildVector(vtype,ZVec,FVec,HVec)

        funcarray=Array(qfunc.QFunc{BigFloat},0)
        for i=1:length(vtype)-2
            if vtype[i][2]=="Z" v=ZVec
            elseif vtype[i][2]=="F" v=FVec
            elseif vtype[i][2]=="H" v=HVec
            end
            v=BigFloat(vtype[i][1])*v
            for el in v push!(funcarray,mcopy(el)) end
        end
        return funcarray
end



setupLP{T<:Real}(sig::T,file::String,vectortypes)=(tab=table.loadTable(file); setupLP(tab,convert(BigFloat,sig),vectortypes))
setupLP{T<:Real}(sigs::Array{T,1},file::String,vectortypes)=setupLP(convert(Array{BigFloat,1},sigs),file,vectortypes)
setupLP(sigs::Array{BigFloat,1},file::String,vectortypes)=(tab=table.loadTable(file); [setupLP(tab,s,vectortypes,file=file)::LP.LinearProgram for s in sigs])


function setupLP(tab::table.Table,sigma::BigFloat, vectortypes,file="N/A")

        println("Setting up LP...")

        # General data
        Lmax=tab.Lmax        
        zerop=qfunc.Polynomial([bf(0)])     # This is the cost associated with the vectors: zero
        spins= tab.OddL ? [0:1:Lmax] : [0:2:Lmax]
        eps=tab.eps
        dim0(L::Int)= L==0 ? maximum([eps+FUDGE,zerobf]) : L+2*eps    # unitarity bound

        # Convolve
        Fvecs=cb.convTable(sigma,tab.table,-1)
        Hvecs=cb.convTable(sigma,tab.table,1)
        Zvecs=[(i*BigFloat(0))::cb.ConvVec_Q{BigFloat} for i in Fvecs]


        #Fill out vector functions..
        lpVectorFuncs=Array(LP.LPVectorFunction{BigFloat},0)
        for (i,L) in enumerate(spins)  #... for every spin                
            for vtype in vectortypes   #... and for every vector type

                    if vtype[end-1]=="even" && mod(L,2)==0
                            v=buildVector(vtype,Zvecs[i].vec,Fvecs[i].vec,Hvecs[i].vec)
                            push!(lpVectorFuncs,LP.LPVectorFunction(
                                    (dim0(L),DELTAMAX),v,zerop,"L=$L - $(vtype[end])")
                            )
                    elseif vtype[end-1]=="odd" && mod(L,2)==1
                            v=buildVector(vtype,Zvecs[i].vec,Fvecs[i].vec,Hvecs[i].vec)
                            push!(lpVectorFuncs,LP.LPVectorFunction(
                                    (dim0(L),DELTAMAX),v,zerop,"L=$L - $(vtype[end])")
                            )
                    elseif vtype[end-1]=="all"
                            v=buildVector(vtype,Zvecs[i].vec,Fvecs[i].vec,Hvecs[i].vec)
                            push!(lpVectorFuncs,LP.LPVectorFunction(
                                    (dim0(L),DELTAMAX),v,zerop,"L=$L - $(vtype[end])")
                            )
                    end
            end
        end

        trgt=-value(lpVectorFuncs[1],bf(0))   # the target: identity vector. We add an extra bit because the identity vector has a bunch of zeros otherwise.
                                                # We must remember to disregard this vector in the final solution

        dim=convert(Float64,2*eps+2)
        prob=LP.initLP(lpVectorFuncs,trgt,"Basic Bound\nD = $dim\tsigma=$(convert(Float64,sigma))\t(m,n) = $((tab.mmax,tab.nmax))\tLmax = $(tab.Lmax)\tOdd spins: $(tab.OddL)",extra=(file,mcopy(sigma)))
        tmp=cullpoles(prob,CULLPOLES)
        println("Done")
        return tmp
end




function saveresults(file::String,prob::LinearProgram)

        f=open(file,"w")
        sol=solution(prob)
        write(f,"{\n")
        for (i,s) in enumerate(sol)
                # dimension, OPE, type
                write(f,"{$(s[1][1]),$(s[2]),\"$(s[1][2])\"}")
                if i<length(sol) write(f,",\n") else write(f,"\n}") end
        end
        close(f)
end


bissect(lp::LinearProgram{BigFloat},top::Real, bot::Real, acc::Real,criteria::LP.LabelF; method="mcv",quiet=false,bak_file="NoBak",bak_iters=100,log_file=LP.LPFILE)=bissect(lp,BigFloat(top),BigFloat(bot),BigFloat(acc),criteria, method=method,quiet=quiet,bak_file=bak_file,bak_iters=bak_iters,log_file=log_file)
# No longer allow generic functions for criteria
function bissect(lp::LinearProgram{BigFloat},top::BigFloat, bot::BigFloat, acc::BigFloat,criteria::LP.LabelF; method="mcv", quiet=false,bak_file="NoBak",bak_iters=100,log_file=LP.LPFILE)

        upper=maximum([bot,top])
        bottom=minimum([bot,top])
		criteriaf=x->x==criteria

        lastfunctional=mcopy(lp)
        lastsol=mcopy(lp)
        lp1=mcopy(lp)
        tmp=mcopy(lp)
		bis_pair=bissect_pair(lastsol,lastfunctional,upper,bottom,criteria)
        
		while (upper-bottom)>acc
                x=1/2*(upper+bottom)
                if !quiet println("x= $x") end
                tmp=mcopy(lp)
                filter!(tmp,[-BigFloat(1),x],criteriaf)

                #----  Hotstart  ----
                tmp.solVecs=lp1.solVecs
                tmp.invA=lp1.invA
                tmp.coeffs=lp1.coeffs
                lp1=tmp

                for v in lp1.solVecs
                        if criteriaf(v.label[2])
                                mcopy(v.cost, v.label[1]<x ? onebf : zerobf)
                        end
                end
                updateFunctional!(lp1)

                 #------------ iterations with backup
					
				if bak_file=="NoBak" 
					iter_bak_file="NoBak"
				else
					iter_bak_file="$(bak_file[1:end-4])_lp1.jld";
				end
				iterate!(lp1,LP_ITERMAX,method=method,quiet=quiet,bak_file=iter_bak_file,bak_iters=bak_iters)			
				#----- end iterations
                
				if cost(lp1)==zerobf
                        bottom=x; lastsol=mcopy(lp1)
                else
                        upper=x;  lastfunctional=mcopy(lp1)
                end
				# backup bissection
				if bak_file!="NoBak"
					mcopy(bis_pair.last_sol,lastsol)
					mcopy(bis_pair.last_func,lastfunctional)
					mcopy(bis_pair.upper,upper)
					mcopy(bis_pair.lower,bottom)
					save(bak_file,"bis_pair",bis_pair)
				end				
        end

        return bissect_pair(lastsol,lastfunctional,upper,bottom,criteria)
end

function bissect!(bis_pair::bissect_pair{BigFloat}, acc::Real; method="mcv", quiet=false, bak_file="NoBak",bak_iters=100,log_file=LP.LPFILE)

        upper=bis_pair.upper
        bottom=bis_pair.lower

		criteria=bis_pair.criteria
		criteriaf=x->x==criteria
        lastfunctional=mcopy(bis_pair.last_func)
        lastsol=mcopy(bis_pair.last_sol)
        lp1=mcopy(lastsol)
        tmp=mcopy(lastsol)

        while (upper-bottom)>acc
                x=1/2*(upper+bottom)
                if !quiet println("x= $x") end
                tmp=mcopy(lastsol)
                filter!(tmp,[-BigFloat(1),x],criteriaf)

                #----  Hotstart  ----
                tmp.solVecs=lp1.solVecs
                tmp.invA=lp1.invA
                tmp.coeffs=lp1.coeffs
                lp1=tmp

                for v in lp1.solVecs
                        if criteriaf(v.label[2])
                                mcopy(v.cost, v.label[1]<x ? onebf : zerobf)
                        end
                end
                updateFunctional!(lp1)
 #------------ iterations with backup
					
				if bak_file=="NoBak" 
					iter_bak_file="NoBak"
				else
					iter_bak_file="$(bak_file[1:end-4])_lp1.jld";
				end
				iterate!(lp1,LP_ITERMAX,method=method,quiet=quiet,bak_file=iter_bak_file,bak_iters=bak_iters,log_file=log_file)			
				#----- end iterations
                
				if cost(lp1)==zerobf
                        bottom=x; lastsol=mcopy(lp1)
                else
                        upper=x;  lastfunctional=mcopy(lp1)
                end
				# backup bissection
				if bak_file!="NoBak"
					mcopy(bis_pair.last_sol,lastsol)
					mcopy(bis_pair.last_func,lastfunctional)
					mcopy(bis_pair.upper,upper)
					mcopy(bis_pair.lower,bottom)
					save(bak_file,"bis_pair",bis_pair)
				end				
        end

		mcopy(bis_pair.last_sol,lastsol)
		mcopy(bis_pair.last_func,lastfunctional)
		mcopy(bis_pair.upper,upper)
		mcopy(bis_pair.lower,bottom)
		bis_pair.criteria=criteria
        return bis_pair
end




function opemax{T<:Real}(lp::LinearProgram{T},confdim::Real,label::LP.LabelF;itermax=LP_ITERMAX,bak_file="NoBak",bak_iters=100,log_file=LP.LPFILE)

        lp2=mcopy(lp)
        x=convert(T,confdim)
		vector=makeVector(lp2.lpFunctions[label],x)
		push!(lp2.lpVectors,vector)
		
        iterate!(lp2,itermax,bak_file=bak_file,bak_iters=bak_iters,log_file=log_file)
        if cost(lp2)!=0 println("Feasible solution not found!"); return lp2 end
        println("Feasible solution found, maximizing OPE...")

		lp2.lpVectors[end].cost=-convert(T,1) #give negative cost to inserted vector
		for v in lp2.lpVectors
			if v.label[2]=="AUX" mcopy(v.cost,BigFloat(1e10)) end #get rid of auxiliary vectors
		end
		
        #push!(lp2.lpVectors,vector)
        iterate!(lp2,itermax,bak_file=bak_file,bak_iters=bak_iters,log_file=log_file)

        return lp2
end

function resume_opemax{T<:Real}(lp::LinearProgram{T};itermax=LP_ITERMAX,bak_file="NoBak",bak_iters=100,log_file=LP.LPFILE)
		
		lp2=mcopy(lp)
		if cost(lp2)>0
			iterate!(lp2,itermax,bak_file=bak_file,bak_iters=bak_iters,log_file=log_file)
			if cost(lp2)!=0 println("Feasible solution not found!"); return lp2 end
		    println("Feasible solution found, maximizing OPE...")
			lp2.lpVectors[end].cost=-convert(T,1) #give negative cost to inserted vector
			for v in lp2.lpVectors
				if v.label[2]=="AUX" mcopy(v.cost,BigFloat(1e10)) end #get rid of auxiliary vectors
			end
		end
        #push!(lp2.lpVectors,vector)
        iterate!(lp2,itermax,bak_file=bak_file,bak_iters=bak_iters,log_file=log_file)

        return lp2
end

#############################
#
# Average spectrum
#
#############################

function avgSpec(lp::LP.LinearProgram;cutoff=1e-6)

	sol=sort(solution(lp));
	Ds=[s[1][1] for s in sol];
	Ls=[s[1][2]::LP.LabelF for s in sol];
	Cs=[s[2] for s in sol];
	Ranges=Dict([lpf.label => lpf.range for lpf in lp.lpFunctions])

    distinct_Ls=[Ls[1]]
    for i=2:length(Ls)
        if findfirst(distinct_Ls,Ls[i])>0 continue end
        push!(distinct_Ls,Ls[i])
    end
    Ls_pos=[find(x->x==l,Ls) for l in distinct_Ls]

	# Need to average

	i=1;
	avDs=Array(BigFloat,0)
	avCs=Array(BigFloat,0)
	doubled=Array(Int64,0)
	labels=Array(LP.LabelF,0)

    Ranges=Dict([lpf.label => lpf.range for lpf in lp.lpFunctions])
    fixed=Array(Int,0);

	ct=1;
    
    for (m,L) in enumerate(distinct_Ls)
        
        if L=="AUX" continue end      
        Ds_L=Ds[Ls_pos[m]]
        Cs_L=Cs[Ls_pos[m]]
        
        i=1
        while i<=length(Ds_L)
            dim=Ds_L[i]
		    ope=Cs_L[i]
		    push!(labels,L)
            if i==length(Ds_L)

			    push!(avDs,dim);
			    push!(avCs,ope);
                ct+=1
			    break
		    end

		    eps=Ds_L[i+1]-Ds_L[i]

		    k=1
		    while abs(eps)<cutoff
        	    push!(doubled,ct)
			    dim=(ope*dim+Cs_L[i+k]*Ds_L[i+k])/(ope+Cs_L[i+k])
			    ope=ope+Cs_L[i+k];
			    k+=1
			    if i+k>length(Ds_L) break; end #eps=2*cutoff; break end
			    eps=Ds_L[i+k]-dim;        
		    end
		    i=i+k;
		    push!(avDs,dim);
		    push!(avCs,ope);    
		    ct+=1;
	    end
    end
#get fixed guys


	for (i,d) in enumerate(avDs)
		if abs(d-Ranges[labels[i]][1])<=1e-10 || abs(d-Ranges[labels[i]][2])<=1e-10 push!(fixed,i) end
	end

	#get singles
	singles=setdiff(collect(1:length(avDs)),[fixed;doubled])
	return (avDs,avCs,fixed,singles,doubled,labels)
end





end
