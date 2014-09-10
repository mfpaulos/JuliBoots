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
