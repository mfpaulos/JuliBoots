module restart

using batchfuncs
jobdir="$(ARGS[1])"
relaunchbatch(jobdir)

end
