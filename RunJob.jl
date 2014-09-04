#
# THIS IS A SCRIPT FOR RUNNING THE JOB, DESCRIBED IN SPECS.JL, IN PARALLEL, AND SAVE THE RESULT.
#


#Initialize threads
addprocs(threads)

#Load code
@everywhere using various,consts,LP,main
@everywhere include("$(pwd())/specs.jl")  #pwd() gives current working directory
@everywhere outputFile="$(outputFile)_$(strftime("%F_%H-%M",time()))"

#Create/reset output file
output=open(outputFile,"w")
    write(output,"$version\n")
    write(output,"$tableFile\n")
    write(output,"$sigInit\n"); write(output,"$sigFinal\n"); write(output,"$sigStep\n");
    write(output,"$top\n");write(output,"$bottom\n");write(output,"$delta\n");
    write(output,"$PRECISION\n")
    write(output,"$BB_NEWTON_GOAL\n")
    write(output,"$BB_ISGOOD\n")
    write(output,"$BB_ITERMAX\n")
    write(output,"$DELTAMAX\n")
    write(output,"$FUDGE\n")
    write(output,"$LP_STOPGOAL\n")
    write(output,"$LP_STOPITERS\n")
    write(output,"$LP_ITERMAX\n")
    write(output,"$CULLPOLES\n")
    write(output,"$LU_FUDGE\n")
close(output)


#Set-up Linear Problems.
problist=setupLP(sigs,tableFile)

#This function does most of the work

@everywhere function runandsave(sigma,prob,bottom,top,delta,outputFile)


        res=bissect(prob,bottom,top,delta,"L=0");
        eps=convert(Float64,res[1].lpFunctions[1].range[1])
        output=open(outputFile,"a")
        write(output,"$sigma\n$eps\n")
        close(output)
end

@everywhere rs(x)=runandsave(x[1],x[2],bottom,top,delta,outputFile)


#Do the Run

#ll=[(sigs[i],problist[i]) for i=1:length(sigs)]
#pmap(rs,ll)







