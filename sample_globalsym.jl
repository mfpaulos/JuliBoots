module sample
using various
using main

sigma=BigFloat(0.518)

# Define types of vectors and their makeup
N=2
v1=[(1/2,"F"), (1,"F"), (1,"H"), "even","S"]
v2=[(3/2-1/N,"F"), (1-2/N,"F"), (-(1+2/N),"H"),"even","T"]
v3=[(-1/2,"F"), (1,"F"), (-1,"H"), "odd","A"]
vectortypes=(v1,v2,v3)

prob=setupLP(sigma,"./Tables/eps0.5n3m1_allspins_prec80.txt", vectortypes)

result=bissect(prob,3,1,0.1,"L=0 - S")


end

