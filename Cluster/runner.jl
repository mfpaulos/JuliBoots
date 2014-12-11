module runner

batchdir=ARGS[1]
jobdir="$(batchdir)/$(ARGS[2])"
include("$(batchdir)/specs.jl")
push!(LOAD_PATH,maindir)

push!(LOAD_PATH,codedir)
push!(LOAD_PATH,"$(codedir)/CBlock")
push!(LOAD_PATH,"$(codedir)/LPsolver")
push!(LOAD_PATH,"$(codedir)/LUdecomp")
push!(LOAD_PATH,"$(codedir)/Minimizer")
push!(LOAD_PATH,"$(codedir)/QFunc")
push!(LOAD_PATH,"$(codedir)/Bootstrap")

using consts
using various
using main


function dopoint()

	sigfile="$(jobdir)/sig.txt"
	sf=open(sigfile,"r")
	sigma=BigFloat(chomp(readline(sf)))
	prob=setupLP(sigma,tableFile)

	result=bissect(prob,top,bottom,delta,name)       	
	output=open("$(jobdir)/results.txt","w")
	eps=convert(Float64,result[1].lpFunctions[1].range[1])
	sig=convert(Float64,sigma)

	write(output,"$sig\n$eps\n")
	close(output)
	main.saveresults("$(jobdir)/spectrum.txt",result[1])
end

dopoint()

end



