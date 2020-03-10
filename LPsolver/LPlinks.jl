####################################################################
#
#    Defines CostFunctions types
#
###################################################################

module LPlinks

using various, consts
import qfunc, cb, bb, lu
export LPFindMinimum, LPInverse, VecFunc, Func, CostFunction, Inverse

CostFunction=Union{Function,qfunc.Qpiece,cb.CB,qfunc.Power}
VecFunc=Union{Array{Function,1}, cb.ConvVec_Q, Array{qfunc.QFunc{BigFloat},1},Array{qfunc.QFunc{Float64},1}}
Func=Union{Function,cb.Conv_Q,qfunc.QFunc}
Inverse=Union{lu.LUdata{BigFloat}}

#--------- Inverse matrix. Only Doolittle method so far

function LPInverse(mat::Array{T,2}; method="doolittle")  where {T<:Real}

        if method=="doolittle"
                return lu.LUdata(mat)
        end
end

function LPInverse(out::Inverse,mat::Array{T,2}; method="doolittle")  where {T<:Real}

        if method=="doolittle"
                return lu.LUdata(out.luMat,mat)
        end
end



#--------- Find minimum. Right now, only Branch and Bound is implemented



function LPFindMinimum(range::Array{T,1},funcc::Func, cost::CostFunction; minMethod="bbGlobal")  where {T<:Real}



    if minMethod=="bbGlobal" || minMethod=="bbLocal" 
        t0=time()


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
		
		
		t1=time()
        if VERBOSE println("Starting find min") end
        if minMethod=="bbGlobal"                    #finds global minimum

            res=[bb.FindMinimum(bb.MinFunction(range,f,df,d2f,d3f))]

            if VERBOSE println("Overall time in FindMinimum: $(time()-t0)") end
            return res
        end

        if minMethod=="bbLocal"

            res=bb.FindLocalMinima(bb.MinFunction(range,f,df,d2f,d3f))
            if VERBOSE println("Overall time in FindMinimum: $(time()-t0), of which $(t1-t0) seting up functions.") end
            return res
        end

    end

	if  minMethod=="bbLocalFloat" 
        t0=time()
		rrange=various.tofloat(range)
        func=various.tofloat(funcc)
		ccost=various.tofloat(cost)
		

        #totalfunc=cost-func #there have to be methods defined for this
        dcost=derivative(ccost,1)
        d2cost=derivative(dcost,1)
        d3cost=derivative(d2cost,1)
        dfunc=derivative(func,1)
        d2func=derivative(dfunc,1)
        d3func=derivative(d2func,1)
        #f(x::BigFloat)=value(cost,x)-value(func,x)
        #df(x::BigFloat)=value(dcost,x)-value(dfunc,x)
        #d2f(x::BigFloat)=value(d2cost,x)-value(d2func,x)
        #d3f(x::BigFloat)=value(d3cost,x)-value(d3func,x)
		ff(x::Float64)=value(ccost,x)-value(func,x)
        dff(x::Float64)=value(dcost,x)-value(dfunc,x)
        d2ff(x::Float64)=value(d2cost,x)-value(d2func,x)
        d3ff(x::Float64)=value(d3cost,x)-value(d3func,x)
		
		
		t1=time()
        if VERBOSE println("Starting find min") end
      
        if minMethod=="bbLocalFloat"

            res=bb.FindLocalMinima(bb.MinFunction(rrange,ff,dff,d2ff,d3ff))
            if VERBOSE println("Overall time in FindMinimum: $(time()-t0), of which $(t1-t0) seting up functions.") end
			if T==Float64 return res end #If the computations was done with Floats, we're done.
			# Otherwise refine it, essentially by calling Newton; only do it for negative guys
			xs=convert(Array{BigFloat,1},[r[1] for r in res])
			pos=find(x->x<0,[r[2] for r in res])
			xs=xs[pos]                #Do this later, after testing
			
			func=funcc
		    dcost=derivative(cost,1)
			d2cost=derivative(dcost,1)
			d3cost=derivative(d2cost,1)
			dfunc=derivative(func,1)
			d2func=derivative(dfunc,1)
			d3func=derivative(d2func,1)
			nf(x::BigFloat)=value(cost,x)-value(func,x)
			ndf(x::BigFloat)=value(dcost,x)-value(dfunc,x)
			nd2f(x::BigFloat)=value(d2cost,x)-value(d2func,x)
			nd3f(x::BigFloat)=value(d3cost,x)-value(d3func,x)
			
			bf_res=Array(Tuple{BigFloat,BigFloat},0)
			eps=parse(BigFloat,"1e-3")   #HACK! This should be a tunable parameter
			eps2=parse(BigFloat,"1e-15")   #HACK! This should be a tunable parameter
			for x in xs
				
				xl=x-eps;xr=x+eps;
				xl=max(xl,range[1]); xr=min(xr,range[2])
				mf=bb.MinFunction([xl,xr],nf,ndf,nd2f,nd3f)
				if abs(x-range[1])<eps2
					push!(bf_res,(range[1],nf(range[1])))
				elseif abs(x-range[2])<eps2
					push!(bf_res,(range[2],nf(range[2])))
				else				
					push!(bf_res,bb.Newton(mf,xl,xr))
				end
			end
			
            return bf_res
        end

    end
	
	
    error("No such minimization method defined")   # if we got here, user put in an undefined method

end

end
