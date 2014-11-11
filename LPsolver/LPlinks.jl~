####################################################################
#
#    Defines CostFunctions types
#
###################################################################

module LPlinks

using various, consts
import qfunc, cb, bb, lu
export LPFindMinimum, LPInverse, VecFunc, Func, CostFunction, Inverse

CostFunction=Union(Function,qfunc.Qpiece,cb.CB)
VecFunc=Union(Array{Function,1}, cb.ConvVec_Q, Array{qfunc.QFunc{BigFloat},1},Array{qfunc.QFunc{Float64},1})
Func=Union(Function,cb.Conv_Q,qfunc.QFunc)
Inverse=Union(lu.LUdata{BigFloat})

#--------- Inverse matrix. Only Doolittle method so far

function LPInverse{T<:Real}(mat::Array{T,2}; method="doolittle")

        if method=="doolittle"
                return lu.LUdata(mat)
        end
end

function LPInverse{T<:Real}(out::Inverse,mat::Array{T,2}; method="doolittle")

        if method=="doolittle"
                return lu.LUdata(out.luMat,mat)
        end
end



#--------- Find minimum. Right now, only Branch and Bound is implemented

function LPFindMinimum(range::(Real,Real),funcc::Func, cost::CostFunction; minMethod="bbGlobal")



    if minMethod=="bbGlobal" || minMethod=="bbLocal"
        t0=time()

        #func=qfunc.trim(funcc,value(funcc,0.5*(range[1]+range[2]))/BigFloat("1e30")) #HACK to make it faster. cut out poles which affect the result to less than one part in 10^-15
        func=funcc


        #totalfunc=cost-func #there have to be methods defined for this
        dcost=derivative(cost,1)
        d2cost=derivative(dcost,1)
        d3cost=derivative(d2cost,1)
        dfunc=derivative(func,1)
        d2func=derivative(dfunc,1)
        d3func=derivative(d2func,1)
        f(x::BigFloat)=value(cost,x)-value(func,x)
        df(x::BigFloat)=value(dcost,x)-value(dfunc,x)
        d2f(x::BigFloat)=value(d2cost,x)-value(d2func,x)
        d3f(x::BigFloat)=value(d3cost,x)-value(d3func,x)



        if VERBOSE println("Starting find min") end
        if minMethod=="bbGlobal"                    #finds global minimum

            res=[bb.FindMinimum(bb.MinFunction(range,f,df,d2f,d3f))]

            if VERBOSE println("Overall time in FindMinimum: $(time()-t0)") end
            return res
        end

        if minMethod=="bbLocal"

            res=bb.FindLocalMinima(bb.MinFunction(range,f,df,d2f,d3f))
            if VERBOSE println("Overall time in FindMinimum: $(time()-t0)") end
            return res
        end

    end

    error("No such minimization method defined")   # if we got here, user put in an undefined method

end

end
