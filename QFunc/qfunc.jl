#####################################################################
#
#   Rational Functions
#
#####################################################################


module qfunc

using various                          #with this, can use functions in various without various.something
import various: derivative, padL, padR, value, tofloat, msum, msub,mmult, mdiv,mplus,mcopy   #we will extend the definitions of derivative, padL, padR;
                                        # to do this, Julia requires explicit import

import Base: convert, promote_rule, promote, isequal, getindex, setindex!, length, +, -, *, /,==, show, dot, endof
export Polynomial, QFunc, pochhammer, value, fastvalue, trim!,Pole,residue, Qpiece, invert, shift_arg, derivative


bf=BigFloat



###########################################################################################
#
#               General methods and types
#
###########################################################################################

#--------------- General types

abstract Qpiece  #a piece of a rational function

########################################################################
#
#   Polynomials
#
#########################################################################


# Types and constructors  ----------------------------------------------------------------

type Polynomial{T<:Real} <:Qpiece

    coeffs::Array{T,1}   # the polynomial coefficients

end

mcopy(p::Polynomial{BigFloat})=Polynomial(mcopy(p.coeffs))
function mcopy(p::Polynomial{BigFloat},q::Polynomial{BigFloat})

       if length(p)>length(q) splice!(p.coeffs,(length(q)+1):length(p)) end
       if length(p)<length(q) append!(p.coeffs,[BigFloat(0) for i=1:(length(q)-length(p))]) end

        mcopy(p.coeffs,q.coeffs); return p
end

function show(io::IO, r::Polynomial)
   pm(x,i)= i==0 ? "$x" : (x < 0 ? " - $(-x) X^$i" :  " + $x X^$i")

   zz=convert(Array{Float64,1},r.coeffs)
   nonzeropos=findn(zz)[1]
   z=[zz[i] for i in nonzeropos]
   ppm(i)=(if i!=1 return pm(z[i],nonzeropos[i]-1) end;
           return (nonzeropos[1]-1)==0 ? "$(z[1])" : "$(z[1]) X^$(nonzeropos[i]-1)"
          )

   if length(z)<=5
         for i=1:length(z) print(io,ppm(i)) end
   else
           print(io, ppm(1), ppm(2),ppm(3),ppm(4), "...",ppm(length(z)))
   end
end




# Methods  ------------------------------------------------------------

==(a::Polynomial,b::Polynomial)= (a.coeffs==b.coeffs)

#--------- ALGEBRA

function +(a::Polynomial,b::Polynomial)

        (l,s)= length(a)>=length(b) ? (a,b) : (b,a) #check which one is longest
        ss=padR(s,length(l)-length(s))
        return Polynomial(l.coeffs+ss.coeffs)
end

# Operations in place for bigfloat polynomials

for fJ in (:mplus,:msub)
    @eval begin
        ($fJ)(a::Polynomial{bf},b::Polynomial{bf})=($fJ)(a,a,b)
        function($fJ)(o::Polynomial{bf},a::Polynomial{bf},b::Polynomial{bf})

            #pad if necessary
            (l,s)= length(a)>=length(b) ? (a,b) : (b,a)
            if length(o)<length(l)
                    for i=1:(length(l)-length(o))
                        push!(o.coeffs,BigFloat(0))
                    end
            end

            for i=1:length(s)
                ($fJ)(o.coeffs[i],s.coeffs[i],l.coeffs[i])
            end
            for i=(length(s)+1):length(l)
                mcopy(o.coeffs[i],l.coeffs[i])
            end
            return o
        end

    end
end
#---


# -----------

*(x::Real,a::Polynomial)= Polynomial(a.coeffs*x)
*(a::Polynomial,x::Real)= Polynomial(a.coeffs*x)

mmult(a::Polynomial{bf},x::Real)=mmult(a,a,x)
function mmult(o::Polynomial{bf},a::Polynomial{bf},x::Real)
         mcopy(o,a)
         for i=1:length(o)
                 mmult(o.coeffs[i],a[i],x)
         end
         return o
end



-(a::Polynomial)=-1*a
-(a::Polynomial,b::Polynomial)=a+(-1)*b


function trim!(p::Polynomial) #remove trailing zeros

        nz=findn(p.coeffs)
       # println("at trim, ",nz)
       # return
        if length(first(nz))==0
            p.coeffs=[zero(p[1])]    # so that it's a zero of the same type (bf, Real, etc)
            return p
        end
        lastnonzero=last(first(nz)) #for compatibility with v0.3.0
        p.coeffs=p[1:lastnonzero]
        return p
end


# ----- Poly multiplication using divide and conquer algorithm ---

function *(pp::Polynomial, qq::Polynomial)

        #pad shortest
        (p,q)= length(pp)>=length(qq) ? (pp,padR(qq,length(pp)-length(qq))) : (qq,padR(pp,length(qq)-length(pp)))

        #divide and conquer algorithm
        n=length(p)
        m=iceil(n/2)


        if n==1 return Polynomial([p[1]*q[1]]) end

        a=Polynomial(p[m+1:n])
        b=Polynomial(p[1:m])
        c=Polynomial(q[m+1:n])
        d=Polynomial(q[1:m])

        tmp1=(a+b)*(c+d)
        tmp2=a*c
        tmp3=b*d

        return trim!(padL(tmp2,2m)+padL(tmp1-tmp2-tmp3,m)+tmp3) # trimming is to remove annoying zeros
end

#----- Polynomial division --------

function /(n::Polynomial, d::Polynomial) #returns the quotient and remainder

        (q,r)=(Polynomial([zero(n[1])]),n)

        while length(r) >= length(d) && endof(d) != 0

                t=endof(r)/endof(d)
                pt=padL(Polynomial([t]),length(r)-length(d))                
                (q,r)= (q+pt, r-(pt*d))                
                trim!(q)                               
                trim!(r)                                
        end

        return (q,r)
end


#------ VARIOUS -------------

shift_arg(p::Polynomial,x::Real)=Polynomial([value(derivative(p,i),x)/factorial(i) for i=0:length(p)-1]) #shifts argument of the polynomial

derivative(p::Polynomial,n::Int64) = n>=length(p) ?
                                            Polynomial([zero(p[1])]) :
                                            Polynomial([p[i+n]*(-1)^(n)*pochhammer((-i-n+1),n) for i=1:length(p)-n])
# REMOVED A BIGFLOAT INSIDE POCHHAMMER 13/11/14
#---- to treat polynomials as arrays...

padL(p::Polynomial,n::Int64)=(q=deepcopy(p); q.coeffs=padL(q.coeffs,n); return q)
padR(p::Polynomial,n::Int64)=(q=deepcopy(p); q.coeffs=padR(q.coeffs,n); return q)


getindex(v::Polynomial,i::Int64)=v.coeffs[i]
getindex(v::Polynomial,i::Range1)=getindex(v.coeffs,i)
setindex!(v::Polynomial,value::Real,i::Int64)=setindex!(v.coeffs,value,i)
endof(x::Polynomial)=x.coeffs[end]
length(p::Polynomial)=length(p.coeffs)



#----- Valuations, using Horner's method

#-- non bigfloat
function value(f::Polynomial,x::Real)
        order=length(f.coeffs)
        output=f.coeffs[order]
        if order==1 return output; end

        for i=0:order-2
            output=x*output+f.coeffs[order-i-1]
        end
        output
end

#-- bigfloat
value(f::Polynomial{BigFloat},x::Real)=(output=zero(BigFloat); return value(f,x,output))
function value(f::Polynomial{BigFloat},x::Real,output::BigFloat)

        order=length(f.coeffs)
        mcopy(output,f.coeffs[order])
        if order==1 return output; end


        for i=0:order-2
            mmult(output,x)
            mplus(output,f.coeffs[order-i-1])
        end

        return output
end





#------ Convert to floats

tofloat(p::Polynomial)=Polynomial(tofloat(p.coeffs))

#------

########################################################################
#
#   Poles
#
#########################################################################

#----------- Types

type Pole{T<:Real} <: Qpiece

        order::Int64
        pole::T
        coeff::T
end

mcopy(p::Pole{BigFloat})=Pole(copy(p.order),mcopy(p.pole),mcopy(p.coeff))
mcopy(q::Pole{BigFloat},p::Pole{BigFloat})=(q.order=copy(p.order); mcopy(q.pole,p.pole); mcopy(q.coeff,p.coeff); q)


function show(io::IO, p::Pole)
    pm(x)= x < 0 ? "+ $(-x)" : "- $x"
    c=convert(Float64,p.coeff)
    po=convert(Float64,p.pole)
    if p.order!=1
        print(io,"$c / (x $(pm(po)))^$(p.order)")
    else
        print(io,"$c / (x $(pm(po)))")
    end
end


#----------- Algebra

*(p::Pole,x::Real)=Pole(p.order,p.pole,p.coeff*x)
*(x::Real,p::Pole)=p*x
mmult(p::Pole,x::Real)=mmult(p,p,x)
mmult(o::Pole,p::Pole,x::Real)=(mmult(o.coeff,p.coeff,x); mcopy(o.pole,p.pole); o.order=copy(p.order); o)

*{T<:Real}(x::Real,p::Array{Pole{T},1})=[(i*x)::Pole for i in p]
*{T<:Real}(p::Array{Pole{T},1},x::Real)=x*p

mmult(p::Array{Pole{bf},1},x::Real)=mmult(p,p,x)
function mmult(o::Array{Pole{bf},1},p::Array{Pole{bf},1},x::Real)

        if length(o)!=length(p) println("Incompatible dimensions! at array p, array p"); return o end
        for i=1:length(o)
                mmult(o[i],p[i],x)
        end
        return o
end





#----------- various Methods

derivative(p::Pole,n::Int64)=(n==0 ? Pole(p.order,p.pole,p.coeff) : (-1)^n*pochhammer(p.order,n)*Pole(p.order+n,p.pole,p.coeff))
derivative(p::Array{Pole,1},n::Int64)=[derivative(pp,n)::Pole for pp in p]


residue(p::Polynomial,q::Pole)=value(derivative(p,q.order-1),q.pole)*q.coeff

isequal(p1::Pole, p2::Pole)= (p1.order==p2.order && p1.pole-p2.pole==zero(p1.pole))

invert(p1::Pole)= 1/(p1.coeff) * Polynomial([(binomial(p1.order,i)*(-1*p1.pole)^(p1.order-i))::typeof(p1.coeff) for i=0:p1.order])


value{T<:Real}(p::Pole{T},x::Real)=p.coeff/((x-p.pole)^(p.order))
value(p::Pole{BigFloat},x::Real)=(res=BigFloat(1); value(p,BigFloat(x),res))
value(p::Pole{BigFloat},x::Real,res::BigFloat)=(msub(res,x,p.pole);mpow(res,p.order); mdiv(res,p.coeff,res))


shift_arg(q::Pole, x::Real)=Pole(q.order,q.pole-x,q.coeff)
shift_arg(q::Array{Pole,1}, x::Real)=[shift_arg(q[i],x)::Pole for i=1:length(q)]

tofloat(p::Pole)=Pole(p.order,tofloat(p.pole),tofloat(p.coeff))
tofloat(p::Array{Pole,1})=[tofloat(i)::Pole for i in p]



########################################################################
#
#   Rational Functions (QFunc)
#
#########################################################################

type QFunc{T<:Real} <: Qpiece

    poly::Polynomial{T}
    poles::Array{Pole{T},1}

end

mcopy(qf::QFunc{BigFloat})=QFunc(mcopy(qf.poly),[mcopy(p)::Pole{bf} for p in qf.poles])

function mcopy(o::QFunc{BigFloat},qf::QFunc{BigFloat})

        if length(qf.poles)!=length(o.poles) println("Wrong dimensions in mcopy(qf,qf)"); return o end
        mcopy(o.poly,qf.poly); for i=1:length(qf.poles) mcopy(o.poles[i],qf.poles[i]) end

        return o
end

mcopy(qf::Array{QFunc{BigFloat},1})=[mcopy(t)::QFunc{BigFloat} for t in qf]
mcopy(qfo::Array{QFunc{BigFloat},1},qf::Array{QFunc{BigFloat},1})=(for i=1:length(qf) mcopy(qfo[i],qf[i]) end; qfo)

show(io::IO,qf::Array{QFunc{BigFloat},1})=print(io,typeof(qf))



#------ promotion and conversion

promote_rule{T<:Real,S<:Real}(::Type{QFunc{T}},::Type{Polynomial{S}})= QFunc{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{QFunc{T}},::Type{Pole{S}})= QFunc{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{Polynomial{T}},::Type{Pole{S}})= QFunc{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{Polynomial{S}},::Type{T})=Polynomial{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{Pole{T}},::Type{S})= QFunc{promote_type(T,S)}
promote_rule{T<:Real,S<:Real}(::Type{QFunc{T}},::Type{S})= QFunc{promote_type(T,S)}


convert{T<:Real}(::Type{QFunc{T}},x::Polynomial{T})=QFunc(x,Array(Pole{T},0))
convert{T<:Real}(::Type{QFunc},x::Polynomial{T})=QFunc(x,Array(Pole{T},0))
convert{T<:Real}(::Type{QFunc{T}},x::Pole{T})=QFunc(Polynomial([zero(x.coeff)]),[x::Pole])
convert{T<:Real}(::Type{QFunc},x::Pole{T})=QFunc(Polynomial([zero(x.coeff)]),[x::Pole])
convert{T<:Real,S<:Real}(::Type{Polynomial{S}},x::T) = Polynomial([x])
convert{T<:Real}(::Type{Polynomial},x::T) = Polynomial([x])
convert{S<:Real,T<:Real}(::Type{QFunc{S}},x::T) = convert(QFunc{T},convert(Polynomial,x))
convert{T<:Real}(::Type{QFunc},x::T) = convert(QFunc{T},convert(Polynomial{T},x))


#----- sum -----------------------------------


function +{T<:Real,S<:Real}(a::QFunc{T},b::QFunc{S})

    polelist1=a.poles
    polelist2=b.poles

    #fast sum:
    pls1=[polelist1[i].pole for i=1:length(polelist1)]
    pld1=[polelist1[i].order for i=1:length(polelist1)]
    pls2=[polelist2[i].pole for i=1:length(polelist2)]
    pld2=[polelist2[i].order for i=1:length(polelist2)]

    if pls1==pls2 && pld1==pld2
            return QFunc(a.poly+b.poly,[Pole(pld1[i],pls1[i],polelist1[i].coeff+polelist2[i].coeff) for i=1:length(polelist1)])
    end

    polelist_res= T==BigFloat ? mcopy(polelist1) : deepcopy(polelist1)



    for i=1:length(polelist2)

            tmp=polelist2[i]
            pos=find(x->isequal(x,tmp),polelist1)

            if length(pos)==0
                    push!(polelist_res,Pole(tmp.order,tmp.pole,tmp.coeff)); continue
            else
                    polelist_res[pos[1]].coeff+=polelist2[i].coeff
            end
    end

    ples=polelist_res
    if length(ples)==0 return convert(QFunc,a.poly+b.poly) end
    return QFunc(a.poly+b.poly,ples)
end


#now we can add up poles and polynomials

+(x::Qpiece,y::Qpiece)= +(promote(x,y)...)
-(x::Qpiece,y::Qpiece)= -(promote(x,y)...)
+{T<:Real}(x::Qpiece,y::T)= +(promote(x,y)...)
+{T<:Real}(y::T,x::Qpiece)= x+y
+(x::Pole, y::Pole)= +(convert(QFunc,x),convert(QFunc,y))

#---- In place -----

#ASSUMES POLE LISTS ARE EQUAL!
for fJ in (:mplus,:msub)
    @eval begin
        ($fJ)(a::QFunc{bf},b::QFunc{bf})=($fJ)(a,a,b)
        function($fJ)(o::QFunc{bf},a::QFunc{bf},b::QFunc{bf})
            if length(o.poles)!=length(a.poles) ||length(o.poles)!=length(b.poles) println("Incompatible dimensions - at qf,qf"); return o end
            ($fJ)(o.poly,a.poly,b.poly)
            for i=1:length(o.poles)
                ($fJ)(o.poles[i].coeff,a.poles[i].coeff,b.poles[i].coeff)
            end

            return o
        end

    end
end
#---



#-------- Multiplication ------------------------------


function *(p::Polynomial,q::Pole)

        if q.order==0 return p*q.coeff end
        if length(p)==1 return p[1]*q end
        #multiplying a polynomial by some pole. result is a poly + sum of poles of various degrees
        #poles part
        poles_part=sum(i-> value(derivative(p,q.order-i),q.pole)/factorial(q.order-i) * Pole(i,q.pole,q.coeff), 1:q.order)
        # poly piece is the quotient of polynomial division
        poly= length(p) > q.order ? (p/invert(q))[1] : Polynomial([zero(p[1])])
        return poly+poles_part
end

*(q::Pole,p::Polynomial)=p*q

function *(q1::Pole, q2::Pole)

        if q1.order==0 && q2.order!=0 return q1.coeff*q2 end
        if q2.order==0 && q1.order!=0 return q2.coeff*q1 end
        if q1.order==0 && q2.order==0 return q1.coeff*q2.coeff end

        if q1.pole==q2.pole
              return Pole(q1.order+q2.order,q1.pole,q1.coeff*q2.coeff)
        end

        a=q1.coeff
        b=q2.coeff
        p1=q1.pole
        p2=q2.pole

        r1= q1.order-1 == 0 ? one(a) : Pole(q1.order-1,p1,one(a))
        r2= q2.order-1 == 0 ? one(a) : Pole(q2.order-1,p2,one(a))

        return ( ((p1-p2))^(-1)*b*q1+((p2-p1))^(-1)*a*q2) * r1*r2
end

function *(qf1::QFunc, qf2::QFunc)

        p1=qf1.poly
        p2=qf2.poly
        q1=qf1.poles
        q2=qf2.poles

        res=p1*p2

        res+=sum(i->p1*q2[i],1:length(q2))
        res+=sum(i->p2*q1[i],1:length(q1))

        for i=1:length(q1)
            res+=sum(j->q1[i]*q2[j],1:length(q2))
        end

        return res
end



#------------- Algebra

*(x::Real,a::QFunc)= QFunc(x*a.poly,x*a.poles)
*(a::QFunc,x::Real)= QFunc(x*a.poly,x*a.poles)
mmult(a::QFunc{BigFloat},x::Real)=mmult(a,a,x)
mmult(o::QFunc{BigFloat},a::QFunc{BigFloat},x::Real)=(mmult(o.poly,a.poly,x); mmult(o.poles,a.poles,x); return o)

mmult(o::QFunc{BigFloat},a::QFunc{BigFloat},x::Real)=(mmult(o.poly,a.poly,x); mmult(o.poles,a.poles,x); return o)


mmult{T<:Real}(a::Array{QFunc{BigFloat},1},x::Array{T,1})=mmult(a,a,x)
mmult{T<:Real}(o::Array{QFunc{BigFloat},1},a::Array{QFunc{BigFloat},1},x::Array{T,1})=(for i=1:length(o) mmult(o[i],a[i],x[i]) end; return o)


-(a::QFunc)=-1*a
-(a::QFunc,b::QFunc)=a+(-1)*b
*(a::Qpiece,b::Qpiece)=*(promote(a,b)...)

+{T<:Real}(a::Array{QFunc{T},1},b::Array{QFunc{T},1})=[(a[i]+b[i])::QFunc{T} for i=1:length(a)]
-{T<:Real}(a::Array{QFunc{T},1},b::Array{QFunc{T},1})=[(a[i]-b[i])::QFunc{T} for i=1:length(a)]
*{T<:Real}(a::Array{QFunc{T},1},x::T)=[(x*qf)::QFunc{T} for qf in a]
*{T<:Real}(x::Real,a::Array{QFunc{T},1})=a*x
*{T<:Real,S<:Real}(x::Array{T,1},a::Array{QFunc{S},1})=[(a[i]*x[i])::QFunc{S} for i=1:length(a)]
*{T<:Real,S<:Real}(a::Array{QFunc{S},1},x::Array{T,1})=x*a


dot{T<:Real}(a::Array{T,1},b::Array{QFunc{T},1})=sum([(a[i]*b[i])::QFunc for i=1:length(b)])
dot{T<:Real}(b::Array{QFunc{T},1},a::Array{T,1})=dot(a,b)

dot{T<:Real}(b::Array{QFunc{bf},1},a::Array{T,1})=(o=mcopy(b[1]); dot(o,a,b))
dot{T<:Real}(a::Array{T,1},b::Array{QFunc{bf},1})=dot(b,a)
function dot{T<:Real}(o::QFunc{bf},a::Array{T,1},b::Array{QFunc{bf},1})
        tmp=mcopy(o)
        mmult(o,b[1],a[1])
        for i=2:length(b)
                mmult(tmp,b[i],a[i])
                mplus(o,tmp)
        end
        return o
end

#------------ VARIOUS METHODS



value(f::QFunc{BigFloat},x::Real)=(output=zero(BigFloat); value(f,x,output))

function value(f::QFunc{BigFloat},x::Real,output::BigFloat)
    value(f.poly,x,output)
    tmp=zero(BigFloat)
    for p in f.poles
            value(p,x,tmp)
            mplus(output,tmp)
    end
    return output
end

function value(f::QFunc,x::Real)
    o=value(f.poly,x)
    o+=sum([value(p,x) for p in f.poles])
    o
end




value{T<:Real}(fvec::Array{QFunc{T},1},x::T)=[value(fi,x)::T for fi in fvec]

#This is supposed to be faster; gives values quicker when we have an array of QFuncs which share the same poles.
#NOTE: Actually it doesn't seem to lead to any particular improvement. Not currently used.

function fastvalue(fvec::Array{QFunc{BigFloat},1},x::Real)

        res=[value(v.poly,x) for v in fvec]
        tmp=BigFloat(1);
        tmp2=BigFloat(1);
        for i=1:length(fvec[1].poles)
                msub(tmp,x,fvec[1].poles[i].pole)
                mpow(tmp,-fvec[1].poles[i].order)
                for j=1:length(fvec)
                        mmult(tmp2,tmp,fvec[j].poles[i].coeff)
                        mplus(res[j],tmp2)
                end

        end

        return res
end

#----------------------------------------------------------------



shift_arg(qf::QFunc,x::Real)=QFunc(shift_arg(qf.poly,x),shift_arg(qf.poles,x))

function trim{T<:Real}(qff::QFunc{T},cutoff::T)

        qf=deepcopy(qff)
        ples=find(x->abs(x)>=cutoff, [p.coeff::T for p in qf.poles])
        newples=[qf.poles[i]::Pole{T} for i in ples]
        qf.poles=newples
        return qf
end

function trim{T<:Real}(qff::QFunc{BigFloat},cutoff::T)

        qf=mcopy(qff)
        ples=find(x->abs(x)>=cutoff, [p.coeff::T for p in qf.poles])
        newples=[qf.poles[i]::Pole{T} for i in ples]
        qf.poles=newples
        return qf
end




trim{T<:Real}(qf::Array{QFunc{T},1},cutoff::T)=[trim(i,cutoff)::QFunc{T} for i in qf]

#----- Derivative

derivative(qf::QFunc,i::Int64)=QFunc(derivative(qf.poly,i),[derivative(qf.poles[j],i)::Pole for j=1:length(qf.poles)])

#----- To Float --- not currently used

tofloat(q::QFunc)=QFunc(tofloat(q.poly),tofloat(q.poles))
tofloat{T<:Real}(q::Array{QFunc{T},1})=[tofloat(i)::QFunc{T} for i in q]




end


