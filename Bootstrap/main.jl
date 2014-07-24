module main

using various
using consts

import bb
import qfunc
import cb
#include("LP.jl")
using LP
import table
export chooseTable, convTable, setupLP, bissect, value


bf=BigFloat
#push!(LOAD_PATH,"C:\\Users\\Miguel_Paulos\\Google\ Drive\ N\\Julia\ Code")

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

function convTable(sig::Real,tab::table.CBDerTable_Q,sign::Int64)
        sigma=BigFloat(sig)
        convtab=cb.convTable(sigma,tab.table,sign)
        ders=sort(collect(keys(convtab[1].dict)))
        table.CBDerTable_Q(convtab,sigma,tab.eps,tab.binprec,tab.nmax,tab.mmax,tab.Lmax,tab.OddL,ders)
end



function cullpoles{T<:Real}(lp0::LP.LinearProblem{T},cutoff::T)

    lp=mcopy(lp0)
    for lpf in lp.lpFunctions
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

    return lp
end


# Setting up a basic LP (no global symmetry)

setupLP{T<:Real}(sig::T,file::String; ders="all")=(tab=table.loadTable(file); setupLP(tab,convert(BigFloat,sig),ders=ders))
setupLP{T<:Real}(sigs::Array{T,1},file::String,ders="all")=setupLP(convert(Array{BigFloat,1},sigs),file,ders=ders)
setupLP(sigs::Array{BigFloat,1},file::String,ders="all")=(tab=table.loadTable(file); [setupLP(tab,s,ders=ders)::LP.LinearProblem for s in sigs])

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

function setupLP(tab::table.Table,sigma::BigFloat; ders="all")

        println("Setting up LP...")
        Lmax=tab.Lmax

        eps=tab.eps
        vecfuncs=cb.convTable(sigma,tab.table,-1) # Convolve

        #----- New ----- 19-07-14

        if ders!="all" vecfuncs=[vf[ders] for vf in vecfuncs] end

        #--------

        dim0(L::Int)= L==0 ? maximum([eps+FUDGE,zerobf]) : L+2*eps    # unitarity bound
        zerop=qfunc.Polynomial([bf(0)])     # This is the cost associated with the vectors: zero

        spins= tab.OddL ? [0:1:Lmax] : [0:2:Lmax]

        #let's create a linear problem based on this data
        trgt=-value(vecfuncs[1].vec,bf(0))    # the target: identity vector

        lpVectorFuncs=[LP.LPVectorFunction(
        (dim0(L),DELTAMAX),vecfuncs[i].vec,zerop,"L=$L")::LP.LPVectorFunction{BigFloat} for (i,L) in enumerate(spins)]


        dim=convert(Float64,2*eps+2)
        prob=LP.initLP(lpVectorFuncs,trgt,"Basic Bound\nD = $dim\tsigma=$(convert(Float64,sigma))\t(m,n) = $((tab.mmax,tab.nmax))\tLmax = $(tab.Lmax)\tOdd spins: $(tab.OddL)")        
        tmp=cullpoles(prob,CULLPOLES)

        # The following normalizes the components by dividing by the scalar with \Delta=1
        val=1/main.value(tmp.lpFunctions[1].vecfunc,BigFloat(2))

        main.mmult(tmp.target,val)
        for lpf in tmp.lpFunctions
                main.mmult(lpf.vecfunc,val)
        end
        for vec in tmp.lpVectors
                main.mmult(vec.vector,val)
        end
        for vec in tmp.solVecs
                main.mmult(vec.vector,val)
        end
        tmp.invA=LP.LPInverse(tmp.invA,LP.getA(tmp))

        LP.updateCoeffs!(tmp)
        #----------

        println("Done")
        return tmp
end



#--- GLOBAL SYMMETRY STUFF

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




function setupLP(tab::table.Table,sigma::BigFloat, vectortypes)

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
        prob=LP.initLP(lpVectorFuncs,trgt,"Basic Bound\nD = $dim\tsigma=$(convert(Float64,sigma))\t(m,n) = $((tab.mmax,tab.nmax))\tLmax = $(tab.Lmax)\tOdd spins: $(tab.OddL)")
        tmp=cullpoles(prob,CULLPOLES)
        println("Done")
        return tmp
end



bissect(lp::LinearProblem{BigFloat},top::Real, bot::Real, acc::Real,criteria::LP.LabelF)=bissect(lp,BigFloat(top),BigFloat(bot),BigFloat(acc),x->x==criteria)
function bissect(lp::LinearProblem{BigFloat},top::BigFloat, bot::BigFloat, acc::BigFloat,criteria::Function)

        upper=maximum([bot,top])
        bottom=minimum([bot,top])


        lastfunctional=mcopy(lp)
        lastsol=mcopy(lp)
        lp1=mcopy(lp)
        tmp=mcopy(lp)

        while (upper-bottom)>acc
                x=1/2*(upper+bottom)
                println("x= $x")
              #  println(upper)
              #  println(bottom)
                tmp=mcopy(lp)
                filter!(tmp,(-BigFloat(1),x),criteria)
                #Hotstart----
                tmp.solVecs=lp1.solVecs
                tmp.invA=lp1.invA
                tmp.coeffs=lp1.coeffs
                lp1=tmp

                for v in lp1.solVecs
                        if criteria(v.label[2])
                                mcopy(v.cost, v.label[1]<x ? onebf : zerobf)
                        end
                end
                updateFunctional!(lp1)

                #------------
                iterate!(lp1,LP_ITERMAX,method="mrc")


                if cost(lp1)==zerobf
                        bottom=x; lastsol=mcopy(lp1)
                else
                        upper=x;  lastfunctional=mcopy(lp1)
                end

        end

        return (lastsol,lastfunctional)
end


#prob=setupLP(0.5182,"./Tables/eps0.5n3m1.txt")

end
