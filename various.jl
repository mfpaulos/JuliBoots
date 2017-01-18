########################################################################
#
# Various
#
#
# Derivative function base is defined here too, so that modules can communicate with each other.
#
#########################################################################


module various

#Add directories to Julia's default load path

push!(LOAD_PATH,"$(pwd())/CBlock")
push!(LOAD_PATH,"$(pwd())/LPsolver")
push!(LOAD_PATH,"$(pwd())/LUdecomp")
push!(LOAD_PATH,"$(pwd())/Minimizer")
push!(LOAD_PATH,"$(pwd())/QFunc")
push!(LOAD_PATH,"$(pwd())/Bootstrap")
#

import consts
import Base: deepcopy,dot
export derivative, pochhammer, padL, padR, value, tofloat,msum,mplus,msub,mdiv,mpow,mmult,MPFR_clear, mcopy, onebf,zerobf,strtime
export sendto, getfrom, passobj,doat
strtime=Libc.strftime

const ROUNDING_MODE = [0]
const DEFAULT_PRECISION = [consts.PRECISION]

# Methods

#---- Seem to remember there was a problem with using gamma function ratios in Julia for bigfloats..
pochhammer(x::Number,n::Integer)= n>0 ? (x+n-1)*pochhammer(x,n-1) : 1
#pochhammer(x::Number,y::Number)=gamma(x+y)/gamma(x)
#------------------


padL{T}(p::Array{T,1},n::Int64)=vcat(zeros(T,n), p)      # padleft and padright with zeros
padR{T}(p::Array{T,1},n::Int64)=vcat(p,zeros(T,n))
tofloat(x::Number)=convert(Float64,x)
tofloat{T<:Number}(a::Array{T,1})=[tofloat(i)::Float64 for i in a]


derivative(x::Number,i::Int64)=zero(x)

#-----


derivative{T}(a::Array{T,1},i::Int64)=[derivative(aa,i)::T for aa in a]

#-----

value(f::Function,x::Number)=f(x)



#--------- MPFR stuff

onebf=one(BigFloat)
zerobf=zero(BigFloat)


for (fJ, fC) in ((:mplus,:add), (:msub,:sub), (:mmult,:mul), (:mdiv,:div), (:mpow, :pow))
    @eval begin
        function ($fJ)(x::BigFloat, y::BigFloat)
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &x, &x, &y, ROUNDING_MODE[end])
            return x
        end
    end
end

for (fJ, fC) in ((:mplus,:add), (:msub,:sub), (:mmult,:mul), (:mdiv,:div), (:mpow, :pow))
    @eval begin
        function ($fJ)(z::BigFloat,x::BigFloat, y::BigFloat)
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Ptr{BigFloat}, Int32), &z, &x, &y, ROUNDING_MODE[end])
            return z
        end
    end
end


for f in (:mplus,:msub,:mmult,:mdiv,:mpow)
        @eval begin
            ($f)(o::BigFloat,x::BigFloat,y::Real)=($f)(o,x,convert(BigFloat,y))
            ($f)(o::BigFloat,y::Real,x::BigFloat)=($f)(o,convert(BigFloat,y),x)
            ($f)(x::BigFloat,y::Real)=($f)(x,x,y)
        end
end

function mpow(x::BigFloat, y::Signed)
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Clong, Int32), &x, &x, y, ROUNDING_MODE[end])
    return x
end

function mpow(z::BigFloat,x::BigFloat, y::Signed)
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Clong, Int32), &z, &x, y, ROUNDING_MODE[end])
    return z
end


# Signed, Float multiplication

mmult(x::BigFloat, c::Signed)=mmult(x,x,c)
function mmult(output::BigFloat, x::BigFloat, c::Signed)
    ccall((:mpfr_mul_si, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Clong, Int32), &output, &x, convert(Clong, c), ROUNDING_MODE[end])
    return x
end

mmult(x::BigFloat, c::Float64)=mmult(x,x,c)
function mmult(output::BigFloat, x::BigFloat, c::Float64)
    ccall((:mpfr_mul_d, :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Float64, Int32), &output, &x, c, ROUNDING_MODE[end])
    return x
end





#--------- Array operations

msum(a::Array{BigFloat,1})=(res=BigFloat(0); for i in a mplus(res,i) end; res)      #foldl(mplus,BigFloat(0),a) for later versions of julia
msum(a::Array{Float64,1})=reduce(+,0.,a)



for fJ in (:mplus,:msub, :mmult, :mdiv, :mpow)
    @eval begin
        ($fJ)(a::Array{BigFloat},b::Array{BigFloat})=($fJ)(a,a,b)
        function($fJ)(o::Array{BigFloat},a::Array{BigFloat},b::Array{BigFloat})

            if size(a)!=size(b) || size(a)!=size(o) println("Incompatible dimensions"); return o end
            for i=1:length(o)
                ($fJ)(o[i],a[i],b[i])
            end
            return o
        end

    end
end



#Copying

function mcopy(a::BigFloat,b::BigFloat)        
        ccall((:mpfr_set, :libmpfr),Int32,(Ptr{BigFloat}, Ptr{BigFloat}, Int32), &a, &b, ROUNDING_MODE[end])
        return a
end

function mcopy(a::Array{BigFloat},b::Array{BigFloat})

        if size(a)!=size(b) println("Incompatible dimensions"); return a; end
        for i=1:length(a)
            mcopy(a[i],b[i])
        end
        return a
end


mcopy(b::BigFloat)=(a=BigFloat(1); mcopy(a,b))
#mcopy(b::Array{BigFloat})=(a=Array(BigFloat,0); for i=1:length(b) push!(a,mcopy(b[i])) end; return reshape(a,size(b)))
mcopy{T<:Any}(b::Array{T})=(a=Array(T,0); for i=1:length(b) push!(a,mcopy(b[i])) end; return reshape(a,size(b)))
mcopy(o::Tuple{BigFloat,BigFloat},a::Tuple{BigFloat,BigFloat})=(mcopy(o[1],a[1]),mcopy(o[2],a[2]))::Tuple{BigFloat,BigFloat}
mcopy(a::Tuple{BigFloat,BigFloat})=(mcopy(a[1]),mcopy(a[2]))::Tuple{BigFloat,BigFloat}


#----- Dotting

#function dot(a::Array{BigFloat,1},b::Array{BigFloat,1})

#         if size(a)!=size(b) println("In various.dot: Different dimensions!") end
#         res=BigFloat(0)
#         tmp=BigFloat(1)

#         for i=1:length(a)
#             mmult(tmp,a[i],b[i])
#             mplus(res,tmp)
#         end
#         return res
#end

#------ PARALLEL DATA TRANSFER --- Taken from ParallelDataTransfer Package, see also
# http://stackoverflow.com/questions/27677399/julia-how-to-copy-data-to-another-processor-in-julia
#

function doat(p::Int,expr) #this one is mine
	@spawnat(p,eval(Main,expr))
end

function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end


function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

getfrom(p::Int, nm::Symbol; mod=Main) = fetch(@spawnat(p, getfield(mod, nm)))

function passobj(src::Int, target::Vector{Int}, nm::Symbol;
                 from_mod=Main, to_mod=Main)
    r = RemoteRef(src)
    @spawnat(src, put!(r, getfield(from_mod, nm)))
    for to in target
        @spawnat(to, eval(to_mod, Expr(:(=), nm, fetch(r))))
    end
    nothing
end


function passobj(src::Int, target::Int, nm::Symbol; from_mod=Main, to_mod=Main)
    passobj(src, [target], nm; from_mod=from_mod, to_mod=to_mod)
end


function passobj(src::Int, target, nms::Vector{Symbol};
                 from_mod=Main, to_mod=Main)
    for nm in nms
        passobj(src, target, nm; from_mod=from_mod, to_mod=to_mod)
    end
end









end
