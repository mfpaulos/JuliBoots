module batchfuncs

export makebatch, localbatch, cernbatch, testrun, relaunchbatch, cernbatchfermion

function incl(file)
    eval(parse(open(readall, file)))
end


function makebatch()

	incl("./specs.jl")
	batchdir="$(rundir)/$(strftime("%F_%H-%M",time()))"
	run(`mkdir $(batchdir)`)
	run(`cp $(maindir)/consts.jl $(batchdir)`)
	run(`cp $(maindir)/specs.jl $(batchdir)`)

	for (i,s) in enumerate(sigs)
		run(`mkdir $(batchdir)/$(i)`)
		run(`echo $(s)` |> "$(batchdir)/$(i)/sig.txt")
	end
end


function get_batch_dir(rundir)
    a=map(chomp,readlines(`ls $(rundir)`))
    println("Batches:")
    for i=1:length(a)
            println("[$i]\t - \t $(a[i]) ")
    end
    println("\nChoose batch: (ENTER for most recent)")
    chr=chomp(readline(STDIN))
    nr= chr=="" ? length(a) : int(chr)
    batchdir="$(rundir)/$(a[nr])"
    return batchdir 
end

function get_batch_dir(rundir,i::Int64)
	a=map(chomp,readlines(`ls $(rundir)`))	
	batchdir= i<0 ? "$(rundir)/$(a[end])" : "$(rundir)/$(a[i])"
	return batchdir
end


function localbatch()
        file="$(pwd())/specs.jl"	
	incl(file)
 	batchdir=get_batch_dir(rundir)
	for (i,s) in enumerate(sigs)
		spawn(`julia runner.jl $(batchdir) $i`)
	end
end

function testrun()
 	batchdir=get_batch_dir(rundir)
        file="$(pwd())/specs.jl"	
	incl(file)
	run(`julia runner.jl $(batchdir) 1`)
end


function cernbatch()
        file="$(pwd())/specs.jl"	
	incl(file)
	batchdir=get_batch_dir(rundir)
	for (i,s) in enumerate(sigs)
		### COMMAND TO BE LAUNCHED ON EACH NODE
		spawn(`bsub -o $(batchdir)/$(i)/output.out -q $(queue) julia $(maindir)/runner.jl $(batchdir) $i`)
		###
	end
end

function cernbatch(i::Int64)
        file="$(pwd())/specs.jl"	
	incl(file)
	batchdir=get_batch_dir(rundir,i)
	for (i,s) in enumerate(sigs)
		### COMMAND TO BE LAUNCHED ON EACH NODE. The Julia part of the command itself should be kept as is, but the job submission
		### (here adapted for the CERN cluster) should be modified.
		spawn(`bsub -o $(batchdir)/$(i)/output.out -q $(queue) julia $(maindir)/runner.jl $(batchdir) $i`)
	end
end
 	
function relaunchbatch(dir)
	file="$(pwd())/$(dir)/specs.jl"
	incl(file)
	batchdir="$(pwd())/$(dir)"
	for (i,s) in enumerate(sigs)
		if !isfile("$(batchdir)/$i/spectrum.txt")
			spawn(`bsub -o $(batchdir)/$(i)/output.out -q $(queue) julia $(maindir)/runner.jl $(batchdir) $i`)
		end
	end
end		



end



