###################################################

# Class Linear Problem

####################################################

module LP

using various
import various: value, tofloat, mcopy, mmult, mplus,mdiv, derivative

using consts
import consts: LPFILE
using LPlinks
import LPlinks: LPFindMinimum, LPInverse, VecFunc, Func, CostFunction, Inverse

import Base: getindex, length, show
#import PyPlot

export LinearProgram, iterate!,filter,filter!,cost,updateFunctional!,
        updateInverse!,updateCoeffs!,LabelF,solution,status,makeVector
export findMRC
export findBVar


#export LPFunction, LPVector, LPVectorFunction

bf=BigFloat

########## TYPES #########

# All possible kinds. Choice is made when setting up the LP



#-- own types

typealias LabelF Any   # how to label functions
typealias LabelV (Real,LabelF)   # how to label vectors (conf dimension, other labels)


type LPVectorFunction{T<:Real}
    range::(T,T)
    vecfunc::VecFunc
    cost::CostFunction
    label::LabelF
end

mcopy(v::LPVectorFunction{BigFloat})=LPVectorFunction(mcopy(v.range),mcopy(v.vecfunc),mcopy(v.cost),deepcopy(v.label))
mcopy(o::LPVectorFunction{BigFloat},v::LPVectorFunction{BigFloat})=(mcopy(o.range,v.range); mcopy(o.vecfunc,v.vecfunc);
                                                                   mcopy(o.cost,v.cost); o.label=deepcopy(v.label); o)
mcopy(v::Array{LPVectorFunction{BigFloat},1})=(o=Array(LPVectorFunction{BigFloat},0); for vf in v push!(o,mcopy(vf)) end; o)

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
    label::(T,LabelF)
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

end


mcopy(lp::LinearProgram{BigFloat})=LinearProgram(mcopy(lp.lpFunctions),mcopy(lp.lpVectors),mcopy(lp.target),
                                        mcopy(lp.solVecs),mcopy(lp.functional),mcopy(lp.invA),deepcopy(lp.label),mcopy(lp.coeffs))

mcopy(o::LinearProgram{BigFloat},lp::LinearProgram{BigFloat})=(mcopy(o.lpFunctions,lp.lpFunctions); mcopy(o.lpVectors,lp.lpVectors);
                                                            mcopy(o.target,lp.target); mcopy(o.solVecs,lp.solVecs); mcopy(o.functional,lp.functional);
                                                            mcopy(o.invA,lp.invA); o.label=deepcopy(lp.label); mcopy(o.coeffs,lp.coeffs); o)



function show(io::IO,lp::LinearProgram)
        print(io,"Linear Problem:")
        print(io,lp.label)
        #for f in lp.lpFunctions show(io,f) end
        #for v in lp.lpVectors show(io,v) end
end










##################################################################################################################
#
#                                   METHODS
#
#################################################################################################################

#---- Algebra




#---- Various
 	
derivative(lv::LPVectorFunction,i::Int64)=LPVectorFunction(lv.range,derivative(lv.vecfunc,i),lv.cost,lv.label)

getindex{T<:Real}(lpa::Array{LPVectorFunction{T},1},label::ASCIIString)=(lpa[find(x->x.label==label,lpa)[1]])
getindex(lv::LPVector,i::Int64)=lv.vector[i]
length(lv::LPVector)=length(lv.vector)
length(lv::LPVectorFunction)=length(lv.vecfunc)

value{T<:Real}(lv::LPVectorFunction{T},x::Real)=value(lv.vecfunc,convert(T,x))

tofloat(lf::LPVectorFunction)=LPVectorFunction((tofloat(lf.range[1])::Float64,tofloat(lf.range[2])::Float64),tofloat(lf.vecfunc),tofloat(lf.cost),lf.label)
tofloat(lv::LPVector)=LPVector{Float64}(tofloat(lv.vector)::Array{Float64,1},tofloat(lv.cost)::Float64,(tofloat(lv.label[1]),lv.label[2]))


tofloat(lp::LinearProgram)=LinearProgram([tofloat(i)::LPVectorFunction for i in lp.lpFunctions],
                                         [tofloat(i)::LPVector{Float64} for i in lp.lpVectors],
                                          tofloat(lp.target),
                                         [tofloat(i)::LPVector{Float64} for i in lp.solVecs],
                                         tofloat(lp.functional),
                                         tofloat(lp.invA),
                                         lp.label
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

function findAllRC{T<:Real}(lp::LinearProgram{T}; minMethod="bbGlobal")   #assumes functional has been computed. returns a tuple LPVector,mrc


        htime=0.
        mrcs=Array((LPVector{T},T),0)
        for i=1:length(lp.lpFunctions)
            t=@elapsed tmp=findMRC(lp.lpFunctions[i],lp.functional; minMethod=minMethod) #this returns a list of minima (eventually with a single element)
            t+=@elapsed for min in tmp push!(mrcs,min) end    #populates mrcs with the minima (which take the form (LPVector, mrc) )
            htime+=t
            if VERBOSE println("In findAllRC: $t\n") end
        end

        if VERBOSE println("Vec Funcs findMRC time: $(htime)") end
        #---
        # to parallelize it's as simple as this:
        # tmp=[(@spawnat i findMRC(lp.lpFunctions[i],lp.functional)) for i=1:length(lp.lpFunctions)]
        # mrc1=[fetch(t) for t in tmp]
        #--
        t=@elapsed for vec in lp.lpVectors push!(mrcs,findMRC(vec,lp.functional)::(LPVector{T},T)) end
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

function findMRC{T<:Real}(lpf::LPVectorFunction{T},functional::Array{T,1}; minMethod ="bbGlobal")

       if VERBOSE println("Starting findMRC.") end
       stfindmrc=time()

       t1=@elapsed dottedfunc=dot(functional, lpf.vecfunc);    #there have to be preexistent ways to dot; these are effectively loaded at LPlinks
                                                   #dottedfunc will have type Func            

       minima=LPFindMinimum(lpf.range,dottedfunc,lpf.cost, minMethod=minMethod) #finds global or local minima over the range lpf.range.
                                                                     #LPFindMinimum function is defined at LPLinks
                                                                     #dottedfunc and cost are sent separately since they have different
                                                                  #types.

       t2=@elapsed output=[(LPVector(lpf,x)::LPVector{T},mrc::T) for (x,mrc) in minima]

       if VERBOSE println("Inside findMRC: $(time()-stfindmrc); dotting took $t1 ; output took $t2 ;  ") end
       return output
end


#--------- Find basic variable  --------------------------------------------

findBVar{T<:Real}(lp::LinearProgram{T},nb::LPVector{T})=(l=findBVar(lp,[nb]); l[1])

function findBVar{T<:Real}(lp::LinearProgram{T},nba::Array{LPVector{T},1})

        ivec=dot(lp.invA,lp.target)
        minx_bvar=Array((T,Int64),0)
        icol=[BigFloat(0) for i=1:length(ivec)]

        for nbv in nba
            icol=dot(lp.invA,nbv.vector)
            xvals=[icol[i]> zero(T) ? ivec[i]/icol[i] : inf(T) for i=1:length(ivec)]

            (minx,bvar)=findmin(xvals)
            #if minx==inf(T) println("Problem unbounded"); return "unbounded" end
            push!(minx_bvar,(minx,bvar))
        end
        return minx_bvar
end



findBVar(lp::LinearProgram{BigFloat},nb::LPVector{BigFloat})=(l=findBVar(lp,[nb]); l[1])

function findBVar(lp::LinearProgram{BigFloat},nba::Array{LPVector{BigFloat},1})

        ivec=lp.coeffs
        minx_bvar=Array((BigFloat,Int64),0)
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




#-------- Swap basic with non-basic ------------------------------------------------------------

function swapBNB!(lp::LinearProgram, nb::LPVector, bvar::Int64)
        #lp.solVecs[bvar]=nb
        mcopy(lp.solVecs[bvar],nb)
        return
end



#-------- Basic iteration of simplex method

function iterate!{T<:Real}(lp::LinearProgram{T},n::Int64; minMethod="bbLocal", method="mcv", quiet=false)


        #Initializations
        log=open(LPFILE,"w")
        tmp=BigFloat(0)         #temporary variable
        stopiters=0
        starttime=time()
        if !quiet println("Started at: $(strftime(time()))\t\t Initial Cost: $(convert(Float64,cost(lp)))") end

        for i=1:n            
            t0=time()
            currentcost=cost(lp)            
            write(log,"$i - $(strftime(time())) - Iteration: $i\n")
            write(log,"$i - $(strftime(time())) - Current cost: $(currentcost)\n")
            if i%100==0 && !quiet println("Iteration $i\nCurrent cost: $(convert(Float64,currentcost))\t\t Elapsed: $(time()-starttime)") end


            #---- Find a vector to bring in -------

            t=@elapsed nb_rc=findAllRC(lp,minMethod=minMethod)     #this is a list of LPVectors and associated reduced costs
            write(log,"$i - $(strftime(time())) - mrc done in $t\n")


            allrc=[nb_rc[i][2] for i=1:length(nb_rc)]       #all reduced costs: one per vector, and a set of local minima for each lpFunction

            if method=="mrc"                #in this case simplex uses the vector with smallest minimum reduced cost
					(mrc,posmin)=findmin(allrc)
					if mrc>=zero(mrc)
						if !quiet println("Min cost achieved") end
						break
					end
                  nb=nb_rc[posmin][1]

                  #---- Swapping ---------------
                  t=@elapsed minx, bvar=findBVar(lp,nb)
                  swapped=(mcopy(lp.solVecs[bvar].label[1]),deepcopy(lp.solVecs[bvar].label[2]));
                  write(log,"$i - $(strftime(time())) - BVar done in $t\n")
                  t=@elapsed swapBNB!(lp,nb,bvar)
                  write(log,"$i - $(strftime(time())) - Swapping done in $t \n")
                  #-----------------------------
            end

            if method=="mcv"       # maximum cost variation: simplex brings in the vector leading to the greatest cost decrease
                  allnb=[nb_rc[i][1]::LPVector{T} for i=1:length(nb_rc)]     #all the LPVector candidates
                  write(log,"$i - $(strftime(time())) - Candidate minima: $(length(allnb))\n")
                  minx_bvar=findBVar(lp,allnb)


                  # More efficient in principle, but less readable
                  #(costvar,pos)=(BigFloat(Inf),0)
                  #for i=1:length(minx_bvar)
                  #      mmult(tmp,minx_bvar[i][1],nb_rc[i][2])
                  #      if tmp<costvar mcopy(costvar,tmp); pos=i; end
                  #end

                  costvars=[minx_bvar[i][1]*nb_rc[i][2] for i=1:length(minx_bvar)] # all cost variations
                  (costvar,pos)=findmin(costvars)
                  if costvar==-inf(BigFloat) println("Problem unbounded"); close(log); return lp end
                  nb=nb_rc[pos][1]
                  mrc=nb_rc[pos][2]
				  if mrc>=zero(mrc) 
				  	if !quiet println("Min cost achieved") end
					break
				  end
                  minx=minx_bvar[pos][1]
                  bvar=minx_bvar[pos][2]
                  swapped=(mcopy(lp.solVecs[bvar].label[1]),deepcopy(lp.solVecs[bvar].label[2]));
                  t=@elapsed swapBNB!(lp,nb,bvar)
                  write(log,"$i - $(strftime(time())) - Swapping done in $t \n")                  
            end


            #updateInverse
            t=@elapsed updateInverse!(lp) #this updates the inverse



            write(log,"$i - $(strftime(time())) - Inverting done in $t \n")
            t=@elapsed updateFunctional!(lp)
            t+=@elapsed updateCoeffs!(lp)
            write(log,"$i - $(strftime(time())) - Updates done in $t \n")

            t=@elapsed cc=cost(lp)
            write(log,"$i - $(strftime(time())) - Cost computed in $t \n")
            write(log,"$i - $(strftime(time())) - Total time: $(time()-t0)\n\n")


            write(log,"$i - $(strftime(time())) - mrc: $mrc\n")
            write(log,"$i - $(strftime(time())) - nb: $(nb.label)\n")
            write(log,"$i - $(strftime(time())) - minx: $(minx)\n")
            write(log,"$i - $(strftime(time())) - bvar: $(swapped)\n")
            write(log,"$i - $(strftime(time())) - newcost: $(cc)\n")
            write(log,"$i - $(strftime(time())) - cost var: $(cc-currentcost)\n")
            write(log,"$i - $(strftime(time())) - minx*mrc: $(minx*mrc)\n")
            write(log,"\n--------------------------------------------\n")


            if abs(cc-currentcost)/maximum([abs(currentcost),FUDGE])<= LP_STOPGOAL stopiters+=1 else stopiters=0 end
            if stopiters== LP_STOPITERS 
				if !quiet println("Not improving any more.") end
				break
			end
            if cc>currentcost println("Cost increased at iteration $i") end
        end


        close(log)
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


function initLP{T<:Real}(lpf::Array{LPVectorFunction{T},1},t::Array{T,1},description::String)

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
                        coeffs
                        )

        updateFunctional!(res)
        updateCoeffs!(res)
        return res
end


# excise the range


filter(lp::LinearProgram{BigFloat},x::Real,criteria::LabelF)=(y=BigFloat(x); filter(lp,(BigFloat(-100),y),z->z==criteria))
filter(lp::LinearProgram{BigFloat},x::Real,criteria::Function)=(y=BigFloat(x); filter(lp,(BigFloat(-100),y),criteria) )
filter(lp::LinearProgram{BigFloat},range::(BigFloat,BigFloat),criteria::Function)=(lp0=mcopy(lp); filter!(lp0,range,criteria))

function filter!(lp::LinearProgram{BigFloat},range::(BigFloat,BigFloat),criteria::Function)

        lpfuncs=Array(LPVectorFunction{BigFloat},0)
        l=range[1]; u=range[2];

        for vecfunc in lp.lpFunctions
                if !(criteria(vecfunc.label)) push!(lpfuncs,mcopy(vecfunc)); continue end #check if criterium is satisfied

                newfunc=vecfunc
                vl=vecfunc.range[1]; vu=vecfunc.range[2]

                if u>vl && u<vu && l<vl newfunc.range=(u,vu) end
                if u>vu && l<vl continue end
                if l>vl && l<vu && u>vu newfunc.range=(vl,l) end
                if l>vl && u<vu
                    newfunc2=mcopy(newfunc);
                    newfunc.range=(vl,l); newfunc2.range=(u,vu);
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
        return
end

functionalCurves(lp::LinearProgram)=[dot(lp.functional,lpf.vecfunc)::Func for lpf in lp.lpFunctions]

function plotFunctional(lp::LinearProgram,i::Int64;nrpoints=100,range=(NaN,NaN),logplot=true)

        PyPlot.figure()

        lpf=lp.lpFunctions[i]
        curve=dot(lp.functional,lpf.vecfunc)

        (d0,df)= isequal(range,(NaN,NaN)) ? lpf.range::(BigFloat,BigFloat) : convert((BigFloat,BigFloat),range)

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


















