module qcb

using juliboots
using LP
using various
using qfunc


myread(file::IO) = parse(BigFloat,chomp(readline(file)))

function load_qcbtable(file::String)

    f=open(file)
    #label=chomp(readline(file))
    eps=myread(f)
    #nmax=convert(Int64,myread(f))
    #mmax=convert(Int64,myread(f))
	Lmax=convert(Int64,myread(f))
	oddL = myread(f)==1 ? true : false
    binprec=convert(Int64,myread(f))
    plen=convert(Int64,myread(f))         # nr of 	
	
	rhovals=Array{BigFloat}(undef,plen)
	etavals=Array{BigFloat}(undef,plen)
	
	for i=1:plen
		rhovals[i]=myread(f)
	end
	for i=1:plen
		etavals[i]=myread(f)
	end
	
    nrspins=oddL ? 1+Lmax : div(Lmax,2)+1
    table=Array{Array{Power{BigFloat}}}(undef,nrspins);
	tableSimp=Array{Array{Power{BigFloat}}}(undef,nrspins);
    

    for lct=1:nrspins
        spin=oddL ? lct-1 : 2*(lct-1)                 # current spin; this holds if we only have even spins in the table
        #delta0=myread(f)                      # currently this is always zero in any table. It's a historical artifact

        maxPlength=convert(Int64,myread(f))     # the maximum degree of the polynomial piece, across components

        polycoeffs=Array{BigFloat}(undef,plen,maxPlength)
		

       # singleSpin=CBVec_Q(rho, Array{QFunc{BigFloat}}(undef,plen), spin,dict,label)  #create CBVec object. dictionary tells us what contents are
		singleSpin=Array{Power{BigFloat}}(undef,plen)
		singleSpinSimp=Array{Power{BigFloat}}(undef,plen)
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
        #tmp=myread(f)
        #nrpoles2=convert(Int,myread(f))
        #nrpoles2=convert(Int,tmp)

        #polelist2=Array{BigFloat}(undef,nrpoles2)
        #polecoeffs2=Array{BigFloat}(undef,plen,nrpoles2)
        #for ct=1:nrpoles2
        #        polelist2[ct]=myread(f)
        #end
        #for pct=1:plen
        #        for ct=1:nrpoles2
        #                polecoeffs2[pct,ct]=myread(f)
        #        end
        #end

        #the hack a[i,:][1:end] picks out column i. Couldn't figure out better way to do it
        
		
        for i=1:plen
            polearray=Array{Pole{BigFloat}}(undef,0)
            for j=1:length(polelist1) push!(polearray,Pole(1,polelist1[j],polecoeffs1[i,j])) end
        #    for j=1:length(polelist2) push!(polearray,Pole(2,polelist2[j],polecoeffs2[i,j])) end
            qf=QFunc(trim!(Polynomial(polycoeffs[i,:][1:end])),polearray)
            singleSpin[i]=Power(rhovals[i]/(3-2sqrt(BigFloat(2))),qf)
			singleSpinSimp[i]=Power(rhovals[i]/(3-2sqrt(BigFloat(2))),Polynomial([polycoeffs[i]]))
        end
    table[lct]=singleSpin
	tableSimp[lct]=singleSpinSimp
    end
    #return table
    #tab=CBDerTable_Q(table,BigFloat(NaN),eps,binprec,nmax,mmax,Lmax,oddL,cb.orderedkeys(dict))

    #return tab
	return [table,tableSimp,rhovals,etavals]
end




end #module

