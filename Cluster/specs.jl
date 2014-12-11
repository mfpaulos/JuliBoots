# specs.jl
# 
# This file specifies all the information required for a particular run.
#

begin


maindir="/afs/cern.ch/user/m/mfpaulos/private/CERN_Julia"   # Directory where cluster code (including this file) is located.
codedir="/afs/cern.ch/user/m/mfpaulos/private/JuliBootS"    # Directory where Julia code is located
tabledir="$(maindir)/Tables"                                # Directory where table files are located
rundir="$(maindir)/proj/FracDim"							   # Directory where run data will be stored



tableFile="$(tabledir)/n7/eps-0.25n7m1_allspins.txt"        # Table file to be used in the run.

# Specifies values of \Delta_\sigma to be considered.
# Each point will require one thread!

s0=0.2501 
s1=0.6
nrpts=50   						

sigs=[s0:(s1-s0)/(nrpts-1):s1]                            
#--------------------------

#---- Bissection parameters

bottom=1.5001
top=5.
delta=1e-3            
name="L=0"            # label of vector function to be bissected


#-- other cluster-specific information that might be used in the runner.jl file when launching the job ==


end
