
#Hook up MinFunctions in bb to QFuncs
import qfunc
function MinFunction(range::Array{<:Real,1},func::qfunc.Qpiece) # corrected issue with syntax ::Real replaced with <:Real
        dfunc=derivative(func,1)
        d2func=derivative(dfunc,1)
        d3func=derivative(d2func,1)
        f(x::Real)=value(func,x)
        df(x::Real)=value(dfunc,x)
        d2f(x::Real)=value(d2func,x)
        d3f(x::Real)=value(d3func,x)
        MinFunction(range,f,df,d2f,d3f)
end
