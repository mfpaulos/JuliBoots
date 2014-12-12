##############################################################################

#               Conformal Blocks and Convolutions

# This module provide types and methods for holding conformal blocks data.

# dd/mm/yy
#
# 06/03/14 - added support for polynomial representation of
# conformal blocks at z=zb=1/2 and their derivatives

##############################################################################

module cb

using various, qfunc

import various: derivative, value, tofloat, mcopy, mplus,mmult,msub                          #we will change all these
import qfunc: shift_arg                                             #
import Base: convert, promote_rule, promote,isequal, getindex,      #
        setindex!, length, dot, +, *, -, /, show, merge

export CB_Q, CBVec_Q, Conv_Q, ConvVec_Q, DerivativeVec
export convBlock

#####################################################################################
#
#                  General Methods and Types
#
######################################################################################

abstract Comp    # means we're dealing with a component
abstract Vec     # means we're dealign with a vector type


abstract CB <: Comp       # abstract type for conformal blocks
                  # all conformal block types should have methods for adding, multiplying by rational functions
                  # a value method that gives the value of CB at a given conf dimension

abstract CBVec <: Vec   # Same thing, but now for a vector of CB's. It is convenient to use a separate type
                  # instead of Array{CB,1}

abstract Convolved <: Comp       # abstract type for convolved conformal blocks
abstract ConvolvedVec <: Vec


#----- Methods

length(v::Vec)=length(v.vec)
value{T<:Real}(v::Vec,x::T)=[value(v[i],x)::T for i=1:length(v)]



######################################################################################
#                                                                                    #
#       CB's with rational representation                                            #
#                                                                                    #
######################################################################################



#-----------------------------------------------------------------------------------#
#------- Types ---------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#


# Vanilla CB's - equal outside scalars, N=0 SUSY


type CB_Q{T<:Real} <: CB    # stands for some component of a CB with a rational (Q) representation

    rho::T     #the value of a CB is the value of the rational function times 4^delta rho^delta
    func::QFunc{T}
    spin::Int64
    label::(Int64,Int64,String)     #label is the derivative (m,n) it corresponds to. String is some name

end

function show(io::IO, cb::CB_Q)
        print(io,"CB_Q - Label = $(cb.label[3]) - \tSpin = $(cb.spin)\t (m,n) = ($(cb.label[1]),$(cb.label[2]))")
end


mcopy(cb::CB_Q{BigFloat})=CB_Q(mcopy(cb.rho),mcopy(cb.func),copy(cb.spin),deepcopy(cb.label))
mcopy(cbo::CB_Q{BigFloat},cb::CB_Q{BigFloat})=(mcopy(cbo.rho,cb.rho); mcopy(cbo.func,cb.func); cbo.spin=copy(cb.spin); cbo=deepcopy(cb.label); cbo)



type CBVec_Q{T<:Real} <:CBVec

        rho::T
        vec::Array{QFunc{T},1}  # a set of CB's representing different derivatives
        spin::Int64
        dict::Dict{(Int64,Int64),Int64}  #derivatives, and corresponding index
        label::String
end

mcopy(cb::CBVec_Q{BigFloat})=CBVec_Q(mcopy(cb.rho),mcopy(cb.vec),copy(cb.spin),deepcopy(cb.dict),deepcopy(cb.label))
mcopy(cbo::CBVec_Q{BigFloat},cb::CBVec_Q{BigFloat})=(mcopy(cbo.rho,cb.rho); mcopy(cbo.vec,cb.vec); cbo.spin=copy(cb.spin); co.dict=deepcopy(cb.dict); cbo=deepcopy(cb.label); cbo)



#---- Convolved types

type Conv_Q{T<:Real} <: Convolved

        rho::T
        sigma::T
        func::QFunc{T}
        spin::Int64
        label::(Int64,Int64,String) #label is the (m,n) derivative it corresponds to, and the string
                                    # tells us if it's convolved in a symmetric or antisymmetric fashion
end

mcopy(cb::Conv_Q{BigFloat})=Conv_Q(mcopy(cb.rho),mcopy(cb.sigma),mcopy(cb.func),copy(cb.spin),deepcopy(cb.label))
mcopy(cbo::Conv_Q{BigFloat},cb::Conv_Q{BigFloat})=
        (mcopy(cbo.rho,cb.rho); mcopy(cbo.sigma,cb.sigma); mcopy(cbo.func,cb.func); cbo.spin=copy(cb.spin); cbo.label=deepcopy(cb.label); cbo)


type ConvVec_Q{T<:Real} <: ConvolvedVec

        rho::T
        sigma::T
        vec::Array{QFunc{T},1}
        spin::Int64
        dict::Dict{(Int64,Int64),Int64}
        label::String
end

mcopy(cb::ConvVec_Q{BigFloat})=ConvVec_Q(mcopy(cb.rho),mcopy(cb.sigma),mcopy(cb.vec),copy(cb.spin),deepcopy(cb.dict),copy(cb.label))
mcopy(cbo::ConvVec_Q{BigFloat},cb::ConvVec_Q{BigFloat})=
        (mcopy(cbo.rho,cb.rho); mcopy(cbo.sigma,cb.sigma); mcopy(cbo.vec,cb.vec); cbo.spin=copy(cb.spin); cbo.dict=deepcopy(cb.dict); cbo.label=copy(cb.label); cbo)

#constructor

ConvVec_Q{T<:Real}(cbvec::CBVec_Q{T},sigma::Real,label::String)=ConvVec_Q(cbvec.rho,sigma,Array(QFunc{T},0),cbvec.spin,Dict{(Int64,Int64),Int64}(),label)


# convenient to define derivative vec type
DerivativeVec=Union(CBVec_Q{BigFloat},CBVec_Q{Float64},ConvVec_Q{BigFloat},ConvVec_Q{Float64})   #all vectors of derivatives
CBDerVec=Union(CBVec_Q{BigFloat},CBVec_Q{Float64})     # vectors of derivatives of conformal blocks


#-----------------------------------------------------------------------------------#
#------- Methods ---------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#


# Get derivatives out of vectors using their dictionaries;
getindex(dv::DerivativeVec,tup::(Int64,Int64))=dv[dv.dict[tup]]
getindex(o::Comp,dv::DerivativeVec,tup::(Int64,Int64))=getindex(o,dv,dv.dict[tup])

# Join different vectors as a single one
function mjoin(vecs::Vec...)

        o=mcopy(vecs[1])

        for v in vecs[2:end]
                o.vec=union(o.vec,v.vec)
                tmp=v.dict
                for k in keys(tmp)
                    tmp[k]=tmp[k]+length(o.dict)
                end
                println(tmp)
                println(o.dict)
                o.dict=merge(o.dict,tmp)
                println(o.dict)
        end
        return o
end



# The following methods can effectively act as constructors

# -------- Get Index ---------

#utility
function orderedkeys(dic::Dict{(Int64,Int64),Int64})
        o=Array((Int64,Int64),0)
        for i=1:length(dic)
            derpos=findfirst([j for j in values(dic)],i)
            push!(o,getindex([k for k in keys(dic)],derpos))
        end
        o
end


#In place
function getindex(o::Conv_Q,v::ConvVec_Q,i::Int64)
        der=orderedkeys(v.dict)[i]
        label=(der[1],der[2],v.label)
        mcopy(o.rho,v.rho); mcopy(o.sigma,v.sigma); mcopy(o.func,v.vec[i]);
        o.spin=v.spin; o.label=label

        return o
end

function getindex(v::ConvVec_Q,i::Int64)

        derpos=findfirst([j for j in values(v.dict)],i)
        der=getindex([i for i in keys(v.dict)],derpos)
        label=(der[1],der[2],v.label)
        Conv_Q(mcopy(v.rho),mcopy(v.sigma),mcopy(v.vec[i]),v.spin,label)
end
#
function getindex(o::CB_Q,v::CBVec_Q,i::Int64)

        derpos=findfirst([j for j in values(v.dict)],i)
        der=getindex([i for i in keys(v.dict)],derpos)
        label=(der[1],der[2],v.label)
        mcopy(o.rho,v.rho);mcopy(o.func,v.vec[i]);
        o.spin=v.spin; o.label=label
        return o
end

function getindex(v::CBVec_Q,i::Int64)


        derpos=findfirst([j for j in values(v.dict)],i)
        der=getindex([i for i in keys(v.dict)],derpos)
        label=(der[1],der[2],v.label)        
        CB_Q(mcopy(v.rho),mcopy(v.vec[i]),v.spin,label)
end

#-------- testing, 19-07-14
function getindex{T<:Real}(v::ConvVec_Q{T},is::Array{Int64,1})

        # To construct new dictionary
        derspos=[findfirst([j for j in values(v.dict)],i) for i in is]
        ders=[getindex([i for i in keys(v.dict)],dp) for dp in derspos]
        #
        newvec=[v.vec[i]::QFunc{T} for i in is]
        newdict=Dict{(Int64,Int64),Int64}()
        for (i,d) in enumerate(ders) push!(newdict,d,i) end

        ConvVec_Q{T}(v.rho,v.sigma,newvec,v.spin,newdict,v.label)
end

getindex{T<:Real}(v::ConvVec_Q{T},ders::Array{(Int64,Int64),1})=getindex(v,[v.dict[d]::Int64 for d in ders])

#------- end testing

#--------------------------

# Dotting up

dot{T}(a::Array{T,1},b::ConvVec_Q)=(c=getindex(b,1); c.func=dot(a,b.vec); c)
dot{T}(b::ConvVec_Q,a::Array{T,1})=dot(a,b)
dot{T}(a::Array{T,1},b::CBVec_Q)=(c=getindex(b,1); c.func=dot(a,b.vec); c)


# Valuations -- vector valuations have already been defined in terms of those of components

value(c::CB_Q,x::Real)=value(c.func,x)*(4*c.rho)^x
value(c::Conv_Q,x::Real)=value(c.func,x)*(4*c.rho)^x

# Shift arguments (makes f(y)->f(y+x) )

shift_arg(cb::CB_Q,x::Real)=CB_Q(cb.rho,(4*cb.rho)^x*shift_arg(cb.func,x),cb.spin,cb.label)

# Algebra

+(a::Comp,b::Comp)=(c=mcopy(a); c.func=a.func+b.func; return c)
+(a::Vec,b::Vec)=(c=mcopy(a); c.vec=a.vec+b.vec; return c)

-(a::Comp,b::Comp)=(c=mcopy(a); c.func=a.func-b.func; return c)
-(a::Vec,b::Vec)=(c=mcopy(a); c.vec=a.vec-b.vec; return c)

*(a::Comp,x::Real)=(c=mcopy(a); c.func=x*c.func; return c)
*(x::Real,a::Comp)=a*x

*(a::Vec,x::Real)=(c=mcopy(a); c.vec=x*c.vec; return c)
*(x::Real,a::Vec)=a*x

#Algebra - in place

mmult(o::Comp,a::Comp,x::Real)=(mcopy(o,a); mmult(o.func,a.func,x); o)
mmult(a::Comp,x::Real)=mmult(a,a,x)
mplus(o::Comp,a::Comp,b::Comp)=(mcopy(o,a); mplus(o.func,a.func,b.func); o)
mplus(a::Comp,x::Real)=mplus(a,a,x)
msub(o::Comp,a::Comp,b::Comp) =(mcopy(o,a); msub(o.func,a.func,b.func); o)


# Derivatives

function derivative(c::Conv_Q,i::Int64)

        if i==0 return mcopy(c) end

        res=mcopy(c)
        x=log(4*c.rho)

        res.func=binomial(i,0)*derivative(res.func,i)
        for k=1:i
            res.func+=binomial(i,k)*(x^k)*derivative(res.func,i-k)
        end


        return res
end


# To Float


tofloat(cb::ConvVec_Q)=ConvVec_Q(tofloat(cb.rho),tofloat(cb.sigma),tofloat(cb.vec),cb.spin,cb.dict,cb.label)
tofloat(cb::Conv_Q)=Conv_Q(tofloat(cb.rho),tofloat(cb.sigma),tofloat(cb.func),cb.spin,cb.label)
tofloat(cb::CB_Q)=CB_Q(tofloat(cb.rho),tofloat(cb.func),cb.spin,cb.label)
tofloat(cb::CBVec_Q)=CBVec_Q(tofloat(cb.rho),tofloat(cb.vec),cb.spin,cb.dict,cb.label)

###############################################################################################
#
# Convolution methods
#
###############################################################################################

function v_convCoeffs(sigma::Real,tup::(Int64,Int64)) #gives a dictionnary with the coefficients
                                                     # necessary to convolve with to obtain the
                                                     # derivatives of u^sigma G
        m=tup[1] #a derivatives
        n=tup[2] #b derivatives
        d=sigma

        #comment: extra powers of 2 relative to my Mathematica code is due to the way the derivatives are defined by the Pythonists

        return [ (k,l)=> 1/(2^(k+2*l))*binomial(n,n-l)*binomial(m,m-k)*(4)^(-d)*pochhammer(-2d+2*(n-l),m-k)*pochhammer(-d, n-l) for k=0:m, l=0:n]
end

function u_convCoeffs(sigma::Real,tup::(Int64,Int64)) #gives a dictionnary with the coefficients
                                                     # necessary to convolve with to obtain the
                                                     # derivative (m,n) (in the tuple) of u^sigma G
        m=tup[1] #a derivatives
        n=tup[2] #b derivatives
        d=sigma

        return [ (k,l)=> 1/(2^(k+2*l))*(-1)^(m-k)*binomial(n,n-l)*binomial(m,m-k)*(4)^(-d)*pochhammer(-2d+2*(n-l),m-k)*pochhammer(-d, n-l) for k=0:m, l=0:n]
end


convBlockAS(sigma::Real,cbvec::CBVec_Q)=convBlock(sigma,cbvec,-1)
convBlockS(sigma::Real,cbvec::CBVec_Q)=convBlock(sigma,cbvec,1)


function convBlock{T<:Real}(sigma::T,cbvec::CBVec_Q{T},sign::Int64)

        dict=cbvec.dict
        ders=sort(collect(keys(dict)))

        #all coefficients required

        mmax=maximum([ders[i][1] for i in 1:length(ders)])
        nmax=maximum([ders[i][2] for i in 1:length(ders)])
        maxders=(mmax,nmax)


        convblock= sign==1 ? ConvVec_Q(cbvec,sigma,"S") : ConvVec_Q(cbvec,sigma,"AS")


        ct=1;

        tmp=cbvec[1]; #to hold temporary results

        for (m,n) in ders

                if m%2==0 && sign==-1 continue end #skip even m, since they lead to zeros for components
                if m%2==1 && sign==1 continue end #skip odd m, since they lead to zeros for components

                push!(convblock.dict,(m,n),ct); ct+=1

                vCoeffs=v_convCoeffs(sigma,(m,n))
                comps=[(k,l) for k=0:m, l=0:n]

                for (i,der) in enumerate(comps)
                            getindex(tmp,cbvec,der)
                            mmult(tmp,2*vCoeffs[der])
                        if i==1
                            push!(convblock.vec,mcopy(tmp.func))
                        else
                            mplus(convblock.vec[end],tmp.func)
                        end
                end

        end
        #println("htime: $htime\notime: $otime") #DEBUG
        return convblock
end




function convTable(sigma::BigFloat,tab::Array{CBVec_Q{BigFloat},1},sign::Int64)

        dict=tab[1].dict
        ders=sort(collect(keys(dict)))

        #all coefficients required

        mmax=maximum([ders[i][1] for i in 1:length(ders)])
        nmax=maximum([ders[i][2] for i in 1:length(ders)])
        maxders=(mmax,nmax)


        convtab= sign==1 ? [ConvVec_Q(tab[i],sigma,"S")::ConvVec_Q for i=1:length(tab)] : [ConvVec_Q(tab[i],sigma,"AS")::ConvVec_Q for i=1:length(tab)]


        ct=1;

        for (m,n) in ders

                if m%2==0 && sign==-1 continue end #skip even m, since they lead to zeros for components
                if m%2==1 && sign==1 continue end #skip odd m, since they lead to zeros for components

                for i=1:length(tab)
                        push!(convtab[i].dict,(m,n),ct)
                end
                ct+=1

                vCoeffs=v_convCoeffs(sigma,(m,n))
                comps=[(k,l) for k=0:m, l=0:n]

                for i=1:length(tab)
                    cbvec=tab[i]
                    tmp=mcopy(cbvec[1]) #to hold temporary results
                    for (k,der) in enumerate(comps)
                        getindex(tmp,cbvec,der)
                        #println(value(tmp,BigFloat("0.7"))) #DeBUG
                        mmult(tmp,2*vCoeffs[der])
                        if k==1
                            push!(convtab[i].vec,mcopy(tmp.func))
                            #println(value(tmp,BigFloat("0.7"))) #DeBUG
                        else
                            mplus(convtab[i].vec[end],tmp.func)
                        end
                    end
                end
        end

        return convtab
end







end
