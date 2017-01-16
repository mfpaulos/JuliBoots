##### branch and bound
#####
#####

module bb

using consts
import various: mcopy
export MinFunction, FindMinimum, FindLocalMinima
using various

#bf(x)=convert(Float64,x)
#bf=BigFloat


type MinFunction{T<:Real}

        range::Array{T,1}
        f::Function
        df::Function
        d2f::Function
        d3f::Function
end

type Interval{T<:Real}

        xL::T
        xR::T
        label::String
        #fL::T
        #fR::T
        dfL::T
        dfR::T
        d2fL::T
        d2fR::T
        d3fL::T
        d3fR::T
end

mcopy(i::Interval{BigFloat},j::Interval{BigFloat})=(mcopy(i.xL,j.xL);mcopy(i.xR,j.xR);
                                                    mcopy(i.dfL,j.dfL);mcopy(i.dfR,j.dfR);
                                                    mcopy(i.d2fL,j.d2fL);mcopy(i.d2fR,j.d2fR);
                                                    mcopy(i.d3fL,j.d3fL);mcopy(i.d3fR,j.d3fR);
                                                    i.label=copy(j.label); i)

type GoodInterval{T<:Real}

        xL::T
        xR::T
        label::String

end

function Interval{T<:Real}(mf::MinFunction{T},label::String)

        xL=mf.range[1]; xR=mf.range[2];
        #fL=mf.f(xL); fR=mf.f(xR);
        dfL=mf.df(xL); dfR=mf.df(xR);
        d2fL=mf.d2f(xL); d2fR=mf.d2f(xR);
        d3fL=mf.d3f(xL); d3fR=mf.d3f(xR);
        Interval(xL,xR,label,dfL,dfR,d2fL,d2fR,d3fL,d3fR)
end


type BBproblem{T<:Real}

        mf::MinFunction{T}
        range::Array{T,1}
        badList::Array{Interval{T},1}
        goodList::Array{GoodInterval{T},1}
end

BBproblem{T<:Real}(mf::MinFunction{T})=BBproblem(mf,mf.range,[Interval(mf,"Bad")],Array(GoodInterval{T},0))


########################################################################################################
#
#                   Minimization Functions
#
#########################################################################################################

# For reals
function divide!{T<:Real}(bb::BBproblem{T})

       ct=0
       htime=0.
       t0=time()
       
       while length(bb.badList)!=0 && ct<=5000
                i=splice!(bb.badList,endof(bb.badList)) #last interval on the list

                xl=i.xL; xr=i.xR;
                dx=xr-xl
				
                (dfL,d2fL,d3fL)=(i.dfL,i.d2fL,i.d3fL)
                (dfR,d2fR,d3fR)=(i.dfR,i.d2fR,i.d3fR)


				taylor_dfL=dfL+dx*d2fL/2+1/8*dx*dx*d3fL
                taylor_dfR=dfR-dx*d2fR/2+1/8*dx*dx*d3fR

                               
                xc=1/2*(xl+xr)

                htime+= @elapsed (dfC,d2fC,d3fC)=(bb.mf.df(xc),bb.mf.d2f(xc),bb.mf.d3f(xc))
                # alloc1+= @allocated (bb.mf.df(xc),bb.mf.d2f(xc),bb.mf.d3f(xc))
                ct+=3




                q=  dfC==0 ? 1. : dfC    # this is to take into account the case where derivative itself is zero
                val=false

                #could be optimized
                if abs((dfC-taylor_dfL)/q) < BB_ISGOOD && abs((dfC-taylor_dfR)/q) < BB_ISGOOD
                            val=true

                       #     if dfR>0 && dfL>0 && d3fC <0
                       #         dxmin=-d2fC/d3fC
                       #         if abs(dxmin)<dx/2
                       #                 val= dfC+dxmin*d2fC/2 > 1/2*minimum([dfL,dfR]) ? true : false
                       #         end
                       #     end
                       #     if dfR<0 && dfL<0 && d3fC> 0
                       #         dxmin=-d2fC/d3fC
                       #         if abs(dxmin)<dx/2
                       #                 val= dfC+dxmin*d2fC/2 > 1/2*minimum([dfL,dfR]) ? true : false
                       #         end
                       #     end
                end

                if val                        
                       push!(bb.goodList,GoodInterval(xl,xr,labelInt(dfL,dfR)))
                else

                        # update lists                                               
                        newint1=Interval(xl,xc,"Bad",dfL,dfC,d2fL,d2fC,d3fL,d3fC); # recall that i is the last interval on the bad list
                        newint2=Interval(xc,xr,"Bad",dfC,dfR,d2fC,d2fR,d3fC,d3fR);
                        push!(bb.badList,newint1);
                        push!(bb.badList,newint2);
                end

       end

       if VERBOSE println("Inside Divide\nFunction valuations: $ct\nTime: $htime / $(time()-t0)") end

end


function divide!(bb::BBproblem{BigFloat})

       ct=0
       htime=0.
       t0=time()
       tmp1=BigFloat(0);     #temporary bigfloats
       tmp2=BigFloat(0);
       tmp3=BigFloat(0);
       tmp4=BigFloat(0);
       tmp5=BigFloat(0);


       while length(bb.badList)!=0 && ct<=5000
                i=splice!(bb.badList,endof(bb.badList)) #last interval on the list

                xl=i.xL; xr=i.xR;
                msub(tmp1,xr,xl) #dx=xr-xl
                dx=tmp1          #--------

                (dfL,d2fL,d3fL)=(i.dfL,i.d2fL,i.d3fL)
                (dfR,d2fR,d3fR)=(i.dfR,i.d2fR,i.d3fR)



                #taylor_dfL=dfL+dx*d2fL/2+1/8*dx*dx*d3fL ----

                mcopy(tmp2,dfL); mcopy(tmp3,tmp2);
                mmult(tmp2,dx,d2fL); mmult(tmp2,0.5); mplus(tmp3,tmp2);
                mmult(tmp2,dx,d3fL); mmult(tmp2,dx); mmult(tmp2,0.125); mplus(tmp3,tmp2);
                mcopy(tmp4,tmp3);
                taylor_dfL=tmp4


                #taylor_dfR=dfR-dx*d2fR/2+1/8*dx*dx*d3fR

                mcopy(tmp2,dfR); mcopy(tmp3,tmp2);
                mmult(tmp2,dx,d2fR); mmult(tmp2,0.5); msub(tmp3,tmp2);
                mmult(tmp2,dx,d3fR); mmult(tmp2,dx); mmult(tmp2,0.125); mplus(tmp3,tmp2);
                mcopy(tmp5,tmp3);

                taylor_dfR=tmp5


                #--------------------------------------------

                #xc=1/2*(xl+xr)
                mplus(tmp1,xl,xr);
                xc=tmp1*0.5 #  Need to create a new xc

                #dfC=tmp1; d2fC=tmp2; d3fC=tmp3;
                #Need to creat new bigfloats anyway, they will be stored - maybe this needs to change...
                htime+= @elapsed (dfC,d2fC,d3fC)=(bb.mf.df(xc),bb.mf.d2f(xc),bb.mf.d3f(xc))
                # alloc1+= @allocated (bb.mf.df(xc),bb.mf.d2f(xc),bb.mf.d3f(xc))
                ct+=3




                q=  dfC==zerobf ? mcopy(tmp2,onebf) : dfC    # this is to take into account the case where derivative itself is zero
                val=false

                #could be optimized
                if abs((dfC-taylor_dfL)/q) < BB_ISGOOD && abs((dfC-taylor_dfR)/q) < BB_ISGOOD
                            val=true

                       #     if dfR>0 && dfL>0 && d3fC <0
                       #         dxmin=-d2fC/d3fC
                       #         if abs(dxmin)<dx/2
                       #                 val= dfC+dxmin*d2fC/2 > 1/2*minimum([dfL,dfR]) ? true : false
                       #         end
                       #     end
                       #     if dfR<0 && dfL<0 && d3fC> 0
                       #         dxmin=-d2fC/d3fC
                       #         if abs(dxmin)<dx/2
                       #                 val= dfC+dxmin*d2fC/2 > 1/2*minimum([dfL,dfR]) ? true : false
                       #         end
                       #     end
                end

                if val                        
                       push!(bb.goodList,GoodInterval(xl,xr,labelInt(dfL,dfR)))
                else

                        # update lists                                               
                        newint1=Interval(xl,xc,"Bad",dfL,dfC,d2fL,d2fC,d3fL,d3fC); # recall that i is the last interval on the bad list
                        newint2=Interval(xc,xr,"Bad",dfC,dfR,d2fC,d2fR,d3fC,d3fR);
                        push!(bb.badList,newint1);
                        push!(bb.badList,newint2);
                end

       end

       if VERBOSE println("Inside Divide\nFunction valuations: $ct\nTime: $htime / $(time()-t0)") end

end




#auxiliary function
function labelInt{T<:Real}(dfL::T,dfR::T)

        if dfL==zero(T) && dfR==zero(T) return "const" end
        if dfL<=zero(T) && dfR<=zero(T) return "down" end
        if dfL>=zero(T) && dfR>=zero(T) return "up" end
        if dfL>=zero(T) && dfR<=zero(T) return "max" end
        if dfL<=zero(T) && dfR>=zero(T) return "min" end

end

# In practice, T is Float64
function Newton{T<:Real}(mf::MinFunction{T},xl::T,xr::T;verbose=false)

        
        x0=1/2*(xl+xr)  #initial guess
        diff=zero(T)
        df0=mf.df(x0)

        x=x0
        iter=0

        tmp=mf.d2f(x0)^2-2*mf.df(x0)*mf.d3f(x0)
        
        if tmp <0
            #println("domain error Newton") #DEBUG
            #diff=BigFloat(0)
        else
            diff=-(sqrt(tmp)-mf.d2f(x0))/mf.d3f(x0)            
        end
		
        goal=max(BB_NEWTON_GOAL,1e-12) #HACK! No use going beyond machine precision
		while abs(diff)>goal && iter<=BB_ITERMAX
		
                iter+=1
                xnew=x-diff
                

                if xl<xnew && xnew<xr
                    x=xnew;
                else
                    if mf.df(x0)<0
                            xl=x
                    else
                            xr=x
                    end
                    x=1/2*(xl+xr)
                end
                diff=mf.df(x)/mf.d2f(x)				
        end
        if iter>=BB_ITERMAX && verbose println("Newton saturated") end

        return (x,mf.f(x))
end


function Newton(mf::MinFunction{BigFloat},xl::BigFloat,xr::BigFloat;verbose=false)

        tmp1=BigFloat(0)
        tmp2=BigFloat(0)
        diff=BigFloat(0)

        x0=1/2*(xl+xr)  #initial guess
        diff=zero(BigFloat)
        df0=mf.df(x0)

        x=x0
        iter=0

        #tmp=mf.d2f(x0)^2-2*mf.df(x0)*mf.d3f(x0)
        tmp1=mf.d3f(x0); mmult(tmp1,df0); mmult(tmp1,-2);
        mpow(tmp2,mf.d2f(x0),2); mplus(tmp2,tmp1)



        if tmp2 <0
            #println("domain error Newton") #DEBUG
            #diff=BigFloat(0)
        else
            #diff=-(sqrt(tmp)-mf.d2f(x0))/mf.d3f(x0)
            tmp1=mf.d2f(x0); mpow(tmp2,BigFloat(1/2)); msub(tmp2,tmp1); mdiv(diff,tmp2,-1*mf.d3f(x0));
        end

        while abs(diff)>BB_NEWTON_GOAL && iter<=BB_ITERMAX
                iter+=1
                #xnew=x-diff
                xnew=tmp1; msub(xnew,x,diff);

                if xl<xnew && xnew<xr
                    x=xnew;
                else
                    if mf.df(x0)<0
                            xl=x
                    else
                            xr=x
                    end
                    mplus(x,xl,xr); mmult(x,0.5);
                end
                mdiv(diff,mf.df(x),mf.d2f(x))
        end
        if iter>=BB_ITERMAX println("Newton saturated") end

        return (x,mf.f(x))
end

FindMinimum{T<:Real}(mf::MinFunction{T})=(
                            minList=FindLocalMinima(mf); (min,pos)=findmin([m[2]::T for m in minList]);
                            return minList[pos])

function FindLocalMinima{T<:Real}(mf::MinFunction{T})

        bb=BBproblem(mf)
        t0=time()
        divide!(bb)

        if VERBOSE println("Overall time spent on divide: $(time()-t0)") end
        minList=Array(Tuple{T,T},0)    # xmin, min
        t0=time()

        for i in bb.goodList

                if i.label=="min" push!(minList,Newton(mf,i.xL,i.xR)) end
                if (i.label=="up" || i.label=="max") && i.xL==bb.range[1]  push!(minList,(i.xL,mf.f(i.xL) ) ) end
                if (i.label=="down" || i.label=="min") && i.xR==bb.range[2] push!(minList,(i.xR,mf.f(i.xR))) end
                if i.label=="const" push!(minList,((i.xL+i.xR)/2, mf.f((i.xL+i.xR)/2))) end
        end

        if VERBOSE println("Time spent on Newton: $(time()-t0)") end
        return minList
end


###### LINKS

include("bblinks.jl")



##### tests


#f=x->x+BigFloat(10)*sin(BigFloat(10)x)
#df=x->1+10*BigFloat(10)*cos(BigFloat(10)*x)
#ddf=x->-100*BigFloat(10)*sin(BigFloat(10)x)
#dddf=x->-1000*BigFloat(10)*cos(BigFloat(10)x)

#f=x->0*x
#df=x->0*x
#ddf=x->0*x
#dddf=x->0*x
#bf=BigFloat
#mf=MinFunction([bf(1),bf(15)],f,df,ddf,dddf)


#bbp=BBproblem(mf)
#tot=FindMinimum(mf)



end


