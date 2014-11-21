function saveresults(file::String,prob::LinearProblem)

        f=open(file,"w")
        sol=solution(prob)
        write(f,"{\n")
        for (i,s) in enumerate(sol)
                # dimension, OPE, type
                write(f,"{$(convert(Float64,s[1][1])),$(convert(Float64,s[2])),\"$(s[1][2])\"}")
                if i<length(sol) write(f,",\n") else write(f,"\n}") end
        end
        close(f)
end


bissect(lp::LinearProblem{BigFloat},top::Real, bot::Real, acc::Real,criteria::LP.LabelF; method="mcv")=bissect(lp,BigFloat(top),BigFloat(bot),BigFloat(acc),x->x==criteria, method=method)
function bissect(lp::LinearProblem{BigFloat},top::BigFloat, bot::BigFloat, acc::BigFloat,criteria::Function; method="mcv")

        upper=maximum([bot,top])
        bottom=minimum([bot,top])


        lastfunctional=mcopy(lp)
        lastsol=mcopy(lp)
        lp1=mcopy(lp)
        tmp=mcopy(lp)

        while (upper-bottom)>acc
                x=1/2*(upper+bottom)
                println("x= $x")
                tmp=mcopy(lp)
                filter!(tmp,(-BigFloat(1),x),criteria)

                #----  Hotstart  ----
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

                iterate!(lp1,LP_ITERMAX,method=method)


                if cost(lp1)==zerobf
                        bottom=x; lastsol=mcopy(lp1)
                else
                        upper=x;  lastfunctional=mcopy(lp1)
                end

        end

        return (lastsol,lastfunctional)
end

function opemax{T<:Real}(lp::LinearProblem{T},confdim::Real,label::LP.LabelF;itermax=LP_ITERMAX)

        lp2=mcopy(lp)
        x=convert(T,confdim)
        iterate!(lp2,itermax)
        if cost(lp2)!=0 println("Feasible solution not found!"); return lp2 end
        println("Feasible solution found, maximizing OPE...")

        local vector
        for i=1:length(lp2.lpFunctions)
                lpf=lp2.lpFunctions[i]
                if lpf.label==label && x<= lpf.range[2] && x>=lpf.range[1]
                        vector=makeVector(lpf,x); vector.cost=-convert(T,1)
                end
        end

        push!(lp2.lpVectors,vector)
        iterate!(lp2,itermax)

        return lp2
end



