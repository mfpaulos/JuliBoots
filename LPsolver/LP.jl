###################################################

# Class Linear Problem

####################################################

module LP

using JLD
#using ParallelDataTransfer
using various
import various: value, tofloat, mcopy, mmult, mplus,mdiv, derivative

using consts
import consts: LPFILE
using LPlinks
import LPlinks: LPFindMinimum, LPInverse, VecFunc, Func, CostFunction, Inverse

import Base: getindex, length, show,filter,filter!
#import PyPlot

export LinearProgram, iterate!,cost,updateFunctional!,
        updateInverse!,updateCoeffs!,LabelF,solution,status,makeVector,filter,filter!,
		LPsave
export findMRC
export findBVar


#export LPFunction, LPVector, LPVectorFunction

bf=BigFloat

########## TYPES #########

# All possible kinds. Choice is made when setting up the LP



#-- own types

typealias LabelF Any   # how to label functions
typealias LabelV Tuple{Real,LabelF}   # how to label vectors (conf dimension, other labels)


type LPVectorFunction{T<:Real}
    range::Array{T,1}
    vecfunc::VecFunc
    cost::CostFunction
    label::LabelF
end

mcopy(v::LPVectorFunction{BigFloat})=LPVectorFunction(mcopy(v.range),mcopy(v.vecfunc),mcopy(v.cost),deepcopy(v.label))
mcopy(o::LPVectorFunction{BigFloat},v::LPVectorFunction{BigFloat})=(mcopy(o.range,v.range); mcopy(o.vecfunc,v.vecfunc);
                                                                   mcopy(o.cost,v.cost); o.label=deepcopy(v.label); o)
mcopy(v::Array{LPVectorFunction{BigFloat},1})=(o=Array(LPVectorFunction{BigFloat},0); for vf in v push!(o,mcopy(vf)) end; o)
mcopy(lv::LabelV)=(mcopy(lv[1]),deepcopy(lv[2]))

function mcopy(o::Array{LPVectorFunction{BigFloat},1},v::Array{LPVectorFunction{BigFloat},1})

        if length(o)!=length(v) println("wrong dims") end
        for i=1:minimum([length(o),length(v)])
                mcopy(o[i],v[i])
        end
      #  if length(o)>length(v) splice!(o,(length(v)+1):length(o)) end
      #  if length(o)<length(v)
      #          for i=(length(o)+1):length(v)
      #                  push!(o,mcopy(v[i]))
      #          end
      #  end
      #  return o
end




type LPVector{T<:Real}
    vector::Array{T,1}
    cost::T
    label::Tuple{T,LabelF}
end

mcopy(v::LPVector{BigFloat})=LPVector(mcopy(v.vector),mcopy(v.cost),(mcopy(v.label[1]),deepcopy(v.label[2])))
mcopy(o::LPVector{BigFloat},v::LPVector{BigFloat})=(mcopy(o.vector,v.vector); mcopy(o.cost,v.cost); mcopy(o.label[1],v.label[1]); o.label=(o.label[1],deepcopy(v.label[2])); o)
mcopy(v::Array{LPVector{BigFloat},1})=(o=Array(LPVector{BigFloat},0); for vf in v push!(o,mcopy(vf)) end; o)

function mcopy(o::Array{LPVector{BigFloat},1},v::Array{LPVector{BigFloat},1})

        if length(o)!=length(v) println("Incorrect dimensions!") end
        for i=1:minimum([length(o),length(v)])
                mcopy(o[i],v[i])
        end
        #if length(o)>length(v) splice!(o,(length(v)+1):length(o)) end
        #if length(o)<length(v)
        #        for i=(length(o)+1):length(v)
        #                push!(o,mcopy(v[i]))
        #        end
        #end
        return o
end



show(io::IO, z::LPVector)=println(io,"LPVector - ",z.label,", cost - ", z.cost)
show(io::IO, z::LPVectorFunction)=println(io,"LPVectorFunc - ",z.label,", range - ",z.range)


LPVector{T<:Real}(lpf::LPVectorFunction{T},xx::Real)=(x=convert(T,xx);
                    LPVector(value(lpf,x),value(lpf.cost,x),(x,lpf.label)))  #constructing a vector by evaluating
                                                                   #a function at a particular point
#shorthand:

makeVector{T<:Real}(lpf::LPVectorFunction{T},xx::Real)=LPVector(lpf,xx)


type LinearProgram{T<:Real}

    lpFunctions::Array{LPVectorFunction{T},1}  # vector functions
    lpVectors::Array{LPVector{T},1}      # discrete vectors

    target::Array{T,1}             # typically the identity vector
    solVecs::Array{LPVector{T},1}      # LP vectors in the current solution

    functional::Array{T,1}        # the current linear functional, \simeq costs.A^{-1}.
    invA::Inverse                # this is basically the set of vectors in the solution, inverted

    label::String               # A description of this linear problem
    coeffs::Array{T,1}          # the coefficients in the current solution
    status::String
	extra::Any					# Any extra useful information - for bootstrap this would be the table file used, the convolution parameter, and the derivatives used.
end


mcopy(lp::LinearProgram{BigFloat})=LinearProgram(mcopy(lp.lpFunctions),mcopy(lp.lpVectors),mcopy(lp.target),
                                        mcopy(lp.solVecs),mcopy(lp.functional),mcopy(lp.invA),deepcopy(lp.label),mcopy(lp.coeffs),deepcopy(lp.status),deepcopy(lp.extra))

mcopy(o::LinearProgram{BigFloat},lp::LinearProgram{BigFloat})=(mcopy(o.lpFunctions,lp.lpFunctions); mcopy(o.lpVectors,lp.lpVectors);
                                                            mcopy(o.target,lp.target); mcopy(o.solVecs,lp.solVecs); mcopy(o.functional,lp.functional);
                                                            mcopy(o.invA,lp.invA); o.label=deepcopy(lp.label); mcopy(o.coeffs,lp.coeffs); o.status=deepcopy(lp.status);o.extra=deepcopy(lp.extra); o)



function show(io::IO,lp::LinearProgram)
        print(io,"Linear Problem:")
        print(io,lp.label)
        print(io,"\t Status: $(lp.status)")
        #for f in lp.lpFunctions show(io,f) end
        #for v in lp.lpVectors show(io,v) end
end

function LPsave{T<:Real}(file::String,lp::LinearProgram{T};reduced=true)

	lp2=mcopy(lp)
	for lpf in lp2.lpFunctions
		lpf.vecfunc=LPlinks.qfunc.QFunc{BigFloat}[]  #Get rid of data
	end
	reduced ? save(file,"lp",lp2) : save(file,"lp",lp) 		#assumes filename ends in JLD
end










##################################################################################################################
#
#                                   METHODS
#
#################################################################################################################

#---- Algebra




#---- Various
 	
derivative(lv::LPVectorFunction,i::Int64)=LPVectorFunction(lv.range,derivative(lv.vecfunc,i),lv.cost,lv.label)

getindex{T<:Real}(lpa::Array{LPVectorFunction{T},1},label::String)=(lpa[find(x->x.label==label,lpa)[1]])
getindex(lv::LPVector,i::Int64)=lv.vector[i]
length(lv::LPVector)=length(lv.vector)
length(lv::LPVectorFunction)=length(lv.vecfunc)

value{T<:Real}(lv::LPVectorFunction{T},x::Real)=value(lv.vecfunc,convert(T,x))

tofloat(lf::LPVectorFunction)=LPVectorFunction([tofloat(lf.range[1])::Float64,tofloat(lf.range[2])::Float64],tofloat(lf.vecfunc),tofloat(lf.cost),lf.label)
tofloat(lv::LPVector)=LPVector{Float64}(tofloat(lv.vector)::Array{Float64,1},tofloat(lv.cost)::Float64,(tofloat(lv.label[1]),lv.label[2]))


tofloat(lp::LinearProgram)=LinearProgram([tofloat(i)::LPVectorFunction for i in lp.lpFunctions],
                                         [tofloat(i)::LPVector{Float64} for i in lp.lpVectors],
                                          tofloat(lp.target),
                                         [tofloat(i)::LPVector{Float64} for i in lp.solVecs],
                                         tofloat(lp.functional),
                                         tofloat(lp.invA),
                                         lp.label,
                                         lp.status
                                         )

#-------- Converts a problem back to big float precision. Takes a Float64 and a BigFloat linear problems, where the first has been
#-------- already improved via simplex --- WARNING: BUGGY AT THIS POINT 4/15/2014

function tobigfloat(lp1::LinearProgram,lp0::LinearProgram)

        res=deepcopy(lp0)
        global newvec
        for (i,v) in enumerate(lp1.solVecs)
            dim=BigFloat(v.label[1])
            lbl=v.label[2]

            for lvf in lp0.lpFunctions
                if lvf.label==lbl
                        newvec=LPVector(lvf,dim); break
                end
            end

            for lv in lp0.lpVectors
                if lv.label[2]==v.label[2] && isequal(tofloat(lv.label[1]),v.label[1])
                   newvec=lv; break
                end
            end
            res.solVecs[i]=newvec
        end

        res.invA=LPInverse(getA(res))
        updateFunctional!(res)

        return res
end






#--------------------------------------------------------------------


#------------------------Simplex helpers


cost(lp::LinearProgram)=(dot(lp.coeffs,[vec.cost for vec in lp.solVecs]))

getA(lp::LinearProgram)=(n=length(lp.solVecs[1]); [lp.solVecs[j][i] for i=1:n, j=1:n])


solution(lp::LinearProgram)=[(lp.solVecs[i].label,lp.coeffs[i]) for i=1:length(lp.coeffs)]


function updateFunctional!{T<:Real}(lp::LinearProgram{T})

        costvec=[v.cost::T for v in lp.solVecs]
        dot(lp.functional,costvec,lp.invA)
end

updateCoeffs!{T<:Real}(lp::LinearProgram{T})=dot(lp.coeffs,lp.invA,lp.target)


updateInverse!(lp::LinearProgram)=(lp.invA=LPInverse(lp.invA,getA(lp)); return)


#####################################################################################
#
#
#                   SIMPLEX METHOD
#
#####################################################################################



#--------- Find MRC --------------------------------------------

function findAllRC{T<:Real}(lp::LinearProgram{T}; minMethod="bbLocal")   #assumes functional has been computed. returns a tuple LPVector,mrc


        htime=0.
        mrcs=Array(Tuple{LPVector{T},T},0)
        for i=1:length(lp.lpFunctions)
            t1=@elapsed tmp=findMRC(lp.lpFunctions[i],lp.functional; minMethod=minMethod) #this returns a list of minima (eventually with a single element)
			if VERBOSE println("In findAllRC: findMRC took $t1\n") end
            for min in tmp push!(mrcs,min) end    #populates mrcs with the minima (which take the form (LPVector, mrc) )
			htime+=t1
            
        end
        if VERBOSE println("Vec Funcs findMRC total time: $(htime)") end

        t=@elapsed for vec in lp.lpVectors push!(mrcs,findMRC(vec,lp.functional)::Tuple{LPVector{T},T}) end
        if VERBOSE println("Vectors findMRC time: $t") end

        return mrcs
end








function findMRC(lp::LinearProgram)

        min=(0.,Inf)
        cands=findAllRC(lp::LinearProgram)

        for cand in cands
            if cand[2]<min[2] min=cand end
        end
        return min
end

findMRC{T<:Real}(lpvec::LPVector{T},functional::Array{T,1})=(lpvec,(lpvec.cost-dot(functional,lpvec.vector)))

function findMRC{T<:Real}(lpf::LPVectorFunction{T},functional::Array{T,1}; minMethod ="bbGlobal",verbose=VERBOSE)

       if verbose println("Starting findMRC.") end
       stfindmrc=time()

       t1=@elapsed dottedfunc=dot(functional, lpf.vecfunc);    #there have to be preexistent ways to dot; these are effectively loaded at LPlinks
                                                   #dottedfunc will have type Func            

       t2=@elapsed minima=LPFindMinimum(lpf.range,dottedfunc,lpf.cost, minMethod=minMethod) #finds global or local minima over the range lpf.range.
                                                                     #LPFindMinimum function is defined at LPLinks
                                                                     #dottedfunc and cost are sent separately since they have different
                                                                  #types.
	   
	   negative_minima=minima[find(x->x[2]<0,minima)] #only keep negative mrcs
       t3=@elapsed output=[(LPVector(lpf,x)::LPVector{T},mrc::T) for (x,mrc) in negative_minima]

       if verbose println("Inside findMRC: $(time()-stfindmrc); dotting took $t1 ; minima took $t2 output took $t3 ;  ") end
       return output
end


#--------- Find basic variable  --------------------------------------------

findBVar{T<:Real}(lp::LinearProgram{T},nb::LPVector{T})=(l=findBVar(lp,[nb]); l[1])

function findBVar{T<:Real}(lp::LinearProgram{T},nba::Array{LPVector{T},1})

        ivec=dot(lp.invA,lp.target)
        minx_bvar=Array((T,Int64),0)
        icol=[BigFloat(0) for i=1:length(ivec)]

        for nbv in nba
            dot(icol,lp.invA,nbv.vector)
            xvals=[icol[i]> zero(T) ? ivec[i]/icol[i] : typemax(T) for i=1:length(ivec)]

            (minx,bvar)=findmin(xvals)
            #if minx==inf(T) println("Problem unbounded"); return "unbounded" end
            push!(minx_bvar,(minx,bvar))
        end
        return minx_bvar
end

function findBVar{T<:Real}(invA::Inverse,coeffs::Array{T,1},nba::Array{LPVector{T},1})

        ivec=coeffs
        minx_bvar=Array((T,Int64),0)
        icol=[BigFloat(0) for i=1:length(ivec)]

        for nbv in nba
            icol=invA*nbv.vector
            xvals=[icol[i]> zero(T) ? ivec[i]/icol[i] : typemax(T) for i=1:length(ivec)]

            (minx,bvar)=findmin(xvals)
            #if minx==inf(T) println("Problem unbounded"); return "unbounded" end
            push!(minx_bvar,(minx,bvar))
        end
        return minx_bvar
end


findBVar(lp::LinearProgram{BigFloat},nb::LPVector{BigFloat})=(l=findBVar(lp,[nb]); l[1])


function findBVar(invA::Inverse,coeffs::Array{BigFloat,1},nba::Array{LPVector{BigFloat},1})

        ivec=coeffs
        minx_bvar=Array(Tuple{BigFloat,Int64},0)
        icol=[BigFloat(0) for i=1:length(ivec)]
        infs=[BigFloat(Inf) for i=1:length(ivec)]
        xvals=mcopy(infs)

        for nbv in nba
            mcopy(xvals,infs)
            dot(icol,invA,nbv.vector)
            #xvals=[icol[i]> zero(T) ? ivec[i]/icol[i] : inf(T) for i=1:length(ivec)]
            for i=1:length(ivec)
                if icol[i]>zerobf mdiv(xvals[i],abs(ivec[i]),icol[i]) end      #abs(..) shouldn't affect anything, since lp.coeffs are always >0
            end

            (minx,bvar)=findmin(xvals)
           # if minx==infs[1] println("Problem unbounded"); return "unbounded" end
            push!(minx_bvar,(mcopy(minx),bvar))
        end

        return minx_bvar
end

function findBVar(lp::LinearProgram{BigFloat},nba::Array{LPVector{BigFloat},1})

        ivec=lp.coeffs
        minx_bvar=Array(Tuple{BigFloat,Int64},0)
        icol=[BigFloat(0) for i=1:length(ivec)]
        infs=[BigFloat(Inf) for i=1:length(ivec)]
        xvals=mcopy(infs)

        for nbv in nba
            mcopy(xvals,infs)
            dot(icol,lp.invA,nbv.vector)
            #xvals=[icol[i]> zero(T) ? ivec[i]/icol[i] : inf(T) for i=1:length(ivec)]
            for i=1:length(ivec)
                if icol[i]>zerobf mdiv(xvals[i],abs(ivec[i]),icol[i]) end      #abs(..) shouldn't affect anything, since lp.coeffs are always >0
            end

            (minx,bvar)=findmin(xvals)
           # if minx==infs[1] println("Problem unbounded"); return "unbounded" end
            push!(minx_bvar,(mcopy(minx),bvar))
        end

        return minx_bvar
end

########## Do find_AllRC and findBVar together in one go. This is useful for paralellization


function find_swap{T<:Real}(lpf::LPVectorFunction{T},functional::Array{T,1},invA::Inverse,coeffs::Array{T,1};minMethod="bbLocal")
		nb_rc=findMRC(lpf,functional,minMethod=minMethod)
		allrc=[nb_rc[i][2] for i=1:length(nb_rc)]
		pos_neg=find(x->x<0,allrc)
		
		neg_nb_rc=nb_rc[find(x->x<0,allrc)]
		neg_nb=[_[1]::LPVector{T} for _ in neg_nb_rc]     #all the LPVector candidates
    
		minx_bvar=findBVar(invA,coeffs,neg_nb)	
			
		#return: minx, bvar,vector,reduced cost
		return [(minx_bvar[i][1],minx_bvar[i][2],neg_nb_rc[i][1],neg_nb_rc[i][2]) for i=1:length(neg_nb_rc)]
end

# This function assumes that at least part of the lp is known for each worker if initWorkers is false.

function find_swaps{T<:Real}(lp::LinearProgram{T};minMethod="bbLocal",initWorkers=false,useWorkers="all")
	
	func=lp.functional
	invA=lp.invA
	lpfs=lp.lpFunctions
	lpvs=lp.lpVectors
	coeffs=lp.coeffs
	
	if useWorkers=="all" whichWorkers=workers() else whichWorkers=useWorkers end
	
	if initWorkers
		for i in whichWorkers
			sendto(i,lpfs=lpfs,lpvs=lp.lpVectors)
		end
	end
	t=0.
	for i in whichWorkers
		t+=@elapsed sendto(i,func=func,invA=invA,coeffs=coeffs)
	end
	#println(t)
	
	#otherwise, all workers should already have lpfs, func and invA (with these exact names) defined locally.
	
	#Create some containers for the results, namely triplets of minx,bvar,and 
		
	#getfrom(1,:func)
	
	for p in whichWorkers
			doat(p,:(res_lpf=Tuple{$T,Int64,LP.LPVector{$T},$T}[]))
			doat(p,:(res_lpv=Tuple{$T,Int64,LP.LPVector{$T},$T}[]))
	end
	
	# Go through lpfunctions
    i = 1
	n=length(lpfs)
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p in whichWorkers
            @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
						cmd=:(res=LP.find_swap(lpfs[$idx],func,invA,coeffs,minMethod=$minMethod); res_lpf=[res_lpf;res])
						fetch(doat(p,cmd))
						#println("$idx done at $p")
						#fetch(doat(p,:(sleep(2))))
						#println("done sleeping")
                    end
                end
            #end
        end
    end
	
			
	# Go through lpvectors (unfinished)
    
	i = 1
	n=length(lpvs)
    @sync begin
        for p=1:0
            #if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end 
						cmd=:(res=LP.find_swap(lpvs[$idx],func,invA,coeffs,minMethod=$minMethod); res_lpv=[res_lpv;res]);
						doat(p,cmd)					
                    end
                end
           # end
        end
	end
    
	#collect results
	
	mm=whichWorkers[1]
	res_lpf=getfrom(mm,:res_lpf)
	res_lpv=getfrom(mm,:res_lpv)
	for p in whichWorkers[2:end]
		res_lpf=[res_lpf;getfrom(p,:res_lpf)]
		res_lpv=[res_lpv;getfrom(p,:res_lpv)]
	end
	
	return total_res=[res_lpf;res_lpv]
end



#-------- Swap basic with non-basic ------------------------------------------------------------

function swapBNB!(lp::LinearProgram, nb::LPVector, bvar::Int64)
        #lp.solVecs[bvar]=nb
        mcopy(lp.solVecs[bvar],nb)
        return
end



#-------- Basic iteration of simplex method

function iterate!{T<:Real}(lp::LinearProgram{T},n::Int64; minMethod="bbLocal", method="mcv", quiet=false,bak_file="NoBak",bak_iters=100,log_file=LPFILE,log_iters=500)


        #Initializations
        log=open(log_file,"w")		
		close(log)
        tmp=BigFloat(0)         #temporary variable
        stopiters=0
        starttime=time()
        if !quiet println("Started at: $(strtime(time()))\t\t Initial Cost: $(convert(Float64,cost(lp)))") end

        for i=1:n
			i%log_iters==0 ? log=open(log_file,"w") : log=open(log_file,"a")		
            t0=time()
            currentcost=cost(lp)            
            write(log,"$i - $(strtime(time())) - Iteration: $i\n")
            write(log,"$i - $(strtime(time())) - Current cost: $(currentcost)\n")
            if i%100==0 && !quiet println("Iteration $i\nCurrent cost: $(convert(Float64,currentcost))\t\t Elapsed: $(time()-starttime)") end
			
			if i%bak_iters==0 && bak_file!="NoBak"
				LPsave(bak_file,lp)
			end

            #---- Find a vector to bring in -------

            
            if method=="mrc"                #in this case simplex uses the vector with smallest minimum reduced cost
					t=@elapsed nb_rc=findAllRC(lp,minMethod=minMethod)     #this is a list of LPVectors and associated reduced costs
					write(log,"$i - $(strtime(time())) - mrc done in $t\n")
					allrc=[nb_rc[i][2] for i=1:length(nb_rc)]       #all reduced costs: one per vector, and a set of local minima for each lpFunction
					
					(mrc,posmin)=findmin(allrc)
					if mrc>=zero(mrc)
						if !quiet println("Min cost achieved"); lp.status="Minimized" end
						break
					end
                  nb=nb_rc[posmin][1]

                  #---- Swapping ---------------
                  t=@elapsed minx, bvar=findBVar(lp,nb)
                  swapped=lp.solVecs[bvar].label
                  write(log,"$i - $(strtime(time())) - BVar done in $t\n")
                  t=@elapsed swapBNB!(lp,nb,bvar)
                  write(log,"$i - $(strtime(time())) - Swapping done in $t \n")
                  #-----------------------------
            end

            if method=="mcv"       # maximum cost variation: simplex brings in the vector leading to the greatest cost decrease
                t=@elapsed nb_rc=findAllRC(lp,minMethod=minMethod)     #this is a list of LPVectors and associated reduced costs
				write(log,"$i - $(strtime(time())) - mrc done in $t\n")
				allrc=[nb_rc[i][2] for i=1:length(nb_rc)]       #all reduced costs: one per vector, and a set of local minima for each lpFunction			  
				allnb=[nb_rc[i][1]::LPVector{T} for i=1:length(nb_rc)]     #all the LPVector candidates
                write(log,"$i - $(strtime(time())) - Candidate minima: $(length(allnb))\n")
				#only check those vectors whose reduced costs are negative
				pos_neg=find(x->x<0,allrc)
				negnb=allnb[find(x->x<0,allrc)]
				negrc=allrc[pos_neg]
				write(log,"$i - $(strtime(time())) - Negative minima: $(length(negrc))\n")
				minx_bvar=findBVar(lp,negnb)
				write(log,"$i - $(strtime(time())) - Finished findBVar\n")

                
                costvars=[minx_bvar[i][1]*negrc[i] for i=1:length(minx_bvar)] # all cost variations
				write(log,"$i - $(strtime(time())) - Finished costvars\n")
				if length(costvars)==0
					if !quiet println("Min cost achieved") end
					lp.status="Minimized"
					break
				end
                (costvar,pos)=findmin(costvars)
				write(log,"$i - $(strtime(time())) - Finished findmin costvars\n")
                if costvar==-typemax(BigFloat) println("Problem unbounded"); lp.status="Unbounded"; break end
                nb=negnb[pos]
                mrc=negrc[pos]
				if mrc>=zero(mrc) 
				  	if !quiet println("Min cost achieved") end
					lp.status="Minimized"
				break
				end
                minx=minx_bvar[pos][1]
                bvar=minx_bvar[pos][2]
                swapped=lp.solVecs[bvar].label
				write(log,"$i - $(strtime(time())) - Finished mcopy.\n")
                t=@elapsed swapBNB!(lp,nb,bvar)
                write(log,"$i - $(strtime(time())) - Swapping done in $t \n")                  
            end
			
			if method=="mcvParallel"       # Same as mcv, but suitable for parallelization - TESTING
				###	THIS IS REPLACED find_swaps
				#nb_rc=findAllRC(lp,minMethod=minMethod)     #this is a list of LPVectors and associated reduced costs
				#allrc=[nb_rc[i][2] for i=1:length(nb_rc)]       #all reduced costs
				#allnb=[nb_rc[i][1]::LPVector{T} for i=1:length(nb_rc)]     #all the LPVector candidates
    
				
				#pos_neg=find(x->x<0,allrc)
				#negnb=allnb[find(x->x<0,allrc)]
				#negrc=allrc[pos_neg]
	
				#minx_bvar=findBVar(lp,negnb)
				#### THIS IS REPLACED BY find_swaps
				
				minxs,bvars,negnb,negrc=find_swaps(lp,minMethod=minMethod)           #Find minx, bvar, and negative reduced costs for all lpFunctions and lpVectors
				####
                
                costvars=[minx[i]*negrc[i] for i=1:length(minx)] # all cost variations
				write(log,"$i - $(strtime(time())) - Finished costvars\n")
				if length(costvars)==0
					if !quiet println("Min cost achieved") end
					lp.status="Minimized"
					break
				end
				
                (costvar,pos)=findmin(costvars)
				write(log,"$i - $(strtime(time())) - Finished findmin costvars\n")
                if costvar==-typemax(BigFloat) println("Problem unbounded"); lp.status="Unbounded"; break end
                nb=negnb[pos]
                mrc=negrc[pos]
				if mrc>=zero(mrc) 
				  	if !quiet println("Min cost achieved") end
					lp.status="Minimized"
				break
				end
                minx=minxs[pos]
                bvar=bvars[pos]
                swapped=lp.solVecs[bvar].label
				write(log,"$i - $(strtime(time())) - Finished mcopy.\n")
                t=@elapsed swapBNB!(lp,nb,bvar)
                write(log,"$i - $(strtime(time())) - Swapping done in $t \n")                  
            end	
			
			
			


            #updateInverse
            t=@elapsed updateInverse!(lp) #this updates the inverse
            write(log,"$i - $(strtime(time())) - Inverting done in $t \n")
            t=@elapsed updateFunctional!(lp)
            t+=@elapsed updateCoeffs!(lp)
            write(log,"$i - $(strtime(time())) - Updates done in $t \n")

            t=@elapsed cc=cost(lp)
            write(log,"$i - $(strtime(time())) - Cost computed in $t \n")
            write(log,"$i - $(strtime(time())) - Total time: $(time()-t0)\n\n")


            write(log,"$i - $(strtime(time())) - mrc: $mrc\n")
            write(log,"$i - $(strtime(time())) - nb: $(nb.label)\n")
            write(log,"$i - $(strtime(time())) - minx: $(minx)\n")
            write(log,"$i - $(strtime(time())) - bvar: $(swapped)\n")
            write(log,"$i - $(strtime(time())) - newcost: $(cc)\n")
            write(log,"$i - $(strtime(time())) - cost var: $(cc-currentcost)\n")
            write(log,"$i - $(strtime(time())) - minx*mrc: $(minx*mrc)\n")
            write(log,"\n--------------------------------------------\n")


            if abs(cc-currentcost)/maximum([abs(currentcost),FUDGE])<= LP_STOPGOAL stopiters+=1 else stopiters=0 end
            if stopiters== LP_STOPITERS 
				if !quiet println("Relative cost variation too slow -- Not improving any more."); lp.status="SlowVar" end
				break
			end
            if cc>currentcost println("Cost increased at iteration $i") end
			close(log)
        end
		
		if isopen(log) close(log) end
        
        #println(cost(lp))
        return lp
end


#########################################################################################################
#
#       Setting up problems
#
#########################################################################################################

#Returns a linear problem object based on a target and an array of LPvector functions; initial solution is obtained using auxiliary vectors

auxLPVector{T<:Real}(i::Int64,n::Int64,cost::T)=LPVector([j==i ? one(cost) : zero(cost) for j=1:n],cost,(i*one(cost),"AUX"))
auxLPMatrix{T<:Real}(n::Int64,cost::T)=[auxLPVector(i,n,cost)::LPVector{T} for i=1:n]


function initLP{T<:Real}(lpf::Array{LPVectorFunction{T},1},t::Array{T,1},description::String;extra="N/A")

        n=length(lpf[1])
        solVecs=auxLPMatrix(n,one(T))

        for i=1:n
                if t[i]<= zero(T) solVecs[i].vector=solVecs[i].vector*(-1) end
        end

        A=[solVecs[j][i] for i=1:n, j=1:n]
        coeffs=[BigFloat(0) for i=1:n]

        res=LinearProgram(lpf,
                        solVecs,                     
                        t,
                        mcopy(solVecs),
                        [zero(T) for k=1:n],
                        LPInverse(A),
                        description,
                        coeffs,
                        "Initialized",
						extra
                        )

        updateFunctional!(res)
        updateCoeffs!(res)        
        return res
end



# excise the range


filter(lp::LinearProgram{BigFloat},x::Real,criteria::LabelF)=(y=BigFloat(x); filter(lp,[BigFloat(-100),y],z->z==criteria))
filter!(lp::LinearProgram{BigFloat},x::Real,criteria::LabelF)=(y=BigFloat(x); filter!(lp,[BigFloat(-100),y],z->z==criteria))
filter!(lp::LinearProgram{BigFloat},x::Real,criteria::Function)=(y=BigFloat(x); filter!(lp,[BigFloat(-100),y],criteria) )
filter(lp::LinearProgram{BigFloat},x::Real,criteria::Function)=(y=BigFloat(x); filter(lp,[BigFloat(-100),y],criteria) )
filter{T<:Real}(lp::LinearProgram{BigFloat},range::Array{T,1},criteria::Function)=(lp0=mcopy(lp); filter!(lp0,range,criteria))

function filter!{T<:Real}(lp::LinearProgram{BigFloat},range::Array{T,1},criteria::Function)

        lpfuncs=Array(LPVectorFunction{BigFloat},0)
        l=range[1]; u=range[2];

        for vecfunc in lp.lpFunctions
                if !(criteria(vecfunc.label)) push!(lpfuncs,mcopy(vecfunc)); continue end #check if criterium is satisfied

                newfunc=vecfunc
                vl=vecfunc.range[1]; vu=vecfunc.range[2]

                if u>vl && u<vu && l<vl newfunc.range=[u,vu] end
                if u>vu && l<vl continue end
                if l>vl && l<vu && u>vu newfunc.range=[vl,l] end
                if l>vl && u<vu
                    newfunc2=mcopy(newfunc);
                    newfunc.range=(vl,l); newfunc2.range=[u,vu];
                    push!(lpfuncs,mcopy(newfunc2))
                end

                push!(lpfuncs,mcopy(newfunc))

        end

        lpvecs=Array(LPVector{BigFloat},0)
        for vector in lp.lpVectors
                if !(criteria(vector.label[2])) push!(lpvecs,mcopy(vector)); continue end #check if criterium is satisfied

                if l<vector.label[1]<u continue end
                push!(lpvecs,mcopy(vector))
        end

        lp.lpVectors=lpvecs
        lp.lpFunctions=lpfuncs
        return lp
end



#######################################################
#
#  Outputs
#
#######################################################

function status(lp::LinearProgram)

        sol=sort(solution(lp))

        for si in sol
                println(
                si[1][2]," ",convert(Float64,si[1][1])," - OPE = ",convert(Float64,si[2])
                )
        end
        println("Cost: $(cost(lp))")
        println("Status: $(lp.status)")
        return
end

functionalCurves(lp::LinearProgram)=[dot(lp.functional,lpf.vecfunc)::Func for lpf in lp.lpFunctions]

function plotFunctional(lp::LinearProgram,i::Int64;nrpoints=100,range=(NaN,NaN),logplot=true)

        PyPlot.figure()

        lpf=lp.lpFunctions[i]
        curve=dot(lp.functional,lpf.vecfunc)

        (d0,df) = isequal(range,(NaN,NaN)) ? lpf.range : convert(Tuple{BigFloat,BigFloat},range)

        xs=[x::Real for x in d0:((df-d0)/nrpoints):df]

        values=convert(Array{Float64,1},map(y->value(curve,y),xs))
        xs=convert(Array{Float64,1},xs)

        if !logplot
                PyPlot.plot(xs,values)
                return (xs,values)
        else
            signs=map(sign,values)
            xs1=Array(Float64,0)
            values1=Array(Float64,0)
            xs2=Array(Float64,0)
            values2=Array(Float64,0)
            for i=1:length(signs)
                if signs[i]>0
                        push!(values1,log(abs(values[i])))
                        push!(xs1,xs[i])
                else
                        push!(values2,log(abs(values[i])))
                        push!(xs2,xs[i])
                end
            end

            PyPlot.plot(xs1,values1)
            PyPlot.plot(xs2,values2)
            return (xs,values)
        end
end

end


















