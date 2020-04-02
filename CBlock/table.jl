module table


#import qfunc, cb
using qfunc, cb
import cb
import Base: show, getindex
export Table, loadTable, convTable

abstract type Table end

mutable struct CBDerTable_Q{T<:cb.DerivativeVec}<:Table
        table::Array{T}
        sigma::BigFloat
        eps::BigFloat
        binprec::Real
        nmax::Int64
        mmax::Int64
        Lmax::Int64
        OddL::Bool
        ders::Array{Tuple{Int64,Int64}}
end

function show(io::IO, t::CBDerTable_Q)
    println("CB derivatives table - D = $(convert(Float64,2+2*t.eps))\n")
    if t.sigma!=BigFloat(NaN) println("Convolved with sigma = $(convert(Float64,t.sigma))") end
    println("(nmax,mmax) = ($(t.nmax),$(t.mmax))\t $(length(t.table[1])) components")
    println("Lmax = $(t.Lmax)\t Odd Spins - $(t.OddL)")
    println("Precision: $(t.binprec)")
end

getindex(t::Table,i::Int64)=getindex(t.table,i)



myread(file::IO) = parse(BigFloat,chomp(readline(file)))


function loadTable(file::String; label="Vanilla N=0")

    f=open(file)
    #label=chomp(readline(file))
    eps=myread(f)
    nmax=convert(Int64,myread(f))
    mmax=convert(Int64,myread(f))

    #dictionary telling me where the derivatives are located. Convention dependent!
    dict=Dict{Tuple{Int64,Int64},Int64}();
    ct=1

    for n=0:nmax

            for m=0:mmax+2(nmax-n)
            dict[(m,n)]=ct #push!(dict::Dict,(m,n),ct)
            ct+=1
            end
    end

    #--------------------------------

    Lmax=convert(Int64,myread(f))

    oddL = myread(f)==1 ? true : false

    binprec=convert(Int64,myread(f))
    plen=convert(Int64,myread(f))         # nr of components


    #set_bigfloat_precision(binprec)     # set precision to be that of the tables
    rho=BigFloat(3-2*sqrt(BigFloat(2)))       # set rho; this corresponds to z=zb=1/2

    nrspins=oddL ? 1+Lmax : div(Lmax,2)+1
    table=Array{CBVec_Q{BigFloat}}(undef,nrspins);
    

    for lct=1:nrspins
        spin=oddL ? lct-1 : 2*(lct-1)                 # current spin; this holds if we only have even spins in the table
        delta0=myread(f)                      # currently this is always zero in any table. It's a historical artifact

        maxPlength=convert(Int64,myread(f))     # the maximum degree of the polynomial piece, across components

        polycoeffs=Array{BigFloat}(undef,plen,maxPlength)

        singleSpin=CBVec_Q(rho, Array{QFunc{BigFloat}}(undef,plen), spin,dict,label)  #create CBVec object. dictionary tells us what contents are

        # read polynomials
		
        for pct=1:plen
                for Pct=1:maxPlength
                    polycoeffs[pct,Pct]=myread(f)
                end
        end

        # read single poles
        nrpoles1=convert(Int,myread(f))
        polelist1=Array{BigFloat}(undef,nrpoles1)
        polecoeffs1=Array{BigFloat}(undef,plen,nrpoles1)

        for ct=1:nrpoles1
                polelist1[ct]=myread(f)
        end
        for pct=1:plen
                for ct=1:nrpoles1
                        polecoeffs1[pct,ct]=myread(f)
                end
        end


        # read double poles
        tmp=myread(f)
        #nrpoles2=convert(Int,myread(f))
        nrpoles2=convert(Int,tmp)

        polelist2=Array{BigFloat}(undef,nrpoles2)
        polecoeffs2=Array{BigFloat}(undef,plen,nrpoles2)
        for ct=1:nrpoles2
                polelist2[ct]=myread(f)
        end
        for pct=1:plen
                for ct=1:nrpoles2
                        polecoeffs2[pct,ct]=myread(f)
                end
        end

        #the hack a[i,:][1:end] picks out column i. Couldn't figure out better way to do it
        # I'm assuming all poles are degree 1. But this can (and will) be changed in the future.

		
        for i=1:plen
            polearray=Array{Pole{BigFloat}}(undef,0)
            for j=1:length(polelist1) push!(polearray,Pole(1,polelist1[j],polecoeffs1[i,j])) end
            for j=1:length(polelist2) push!(polearray,Pole(2,polelist2[j],polecoeffs2[i,j])) end
            qf=QFunc(trim!(Polynomial(polycoeffs[i,:][1:end])),polearray)
            singleSpin.vec[i]=qf
        end
    table[lct]=singleSpin
    end
    #return table
    tab=CBDerTable_Q(table,BigFloat(NaN),eps,binprec,nmax,mmax,Lmax,oddL,cb.orderedkeys(dict))
    #PYTHON STUFF - TEMPORARY
     #for (i,cbb) in enumerate(tab.table)
     #       for k=1:length(cbb)
     #               cbb.vec[k].poly=qfunc.shift_arg(cbb.vec[k].poly,-BigFloat(2*0.00001+(i-1)))
     #       end
    #end
    return tab
end



#---- Convolution -----
# This method allows one to convolve a table of type table.CBDerTable_Q.
# The method with the same name in cb.jl does not return a nice object, but rather
# it works directly with an array of conformal blocks.
# The presently defined method allows one to have a nice table object

function convTable(sig::Real,tab::table.CBDerTable_Q,sign::Int64;delsigma=false)
        sigma=BigFloat(sig)
        convtab=cb.convTable(sigma,tab.table,sign,delsigma=delsigma)
        ders=sort(collect(keys(convtab[1].dict)))
        table.CBDerTable_Q(convtab,sigma,tab.eps,tab.binprec,tab.nmax,tab.mmax,tab.Lmax,tab.OddL,ders)
end





end #module


######


