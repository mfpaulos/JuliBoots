module sample
using juliboots

sigma=BigFloat(0.52)

# Define types of vectors and their makeup
N=3
v1=[(1/2,"F"), (1,"F"), (1,"H"), "even","S"]
v2=[(3/2-1/N,"F"), (1-2/N,"F"), (-(1+2/N),"H"),"even","T"]
v3=[(-1/2,"F"), (1,"F"), (-1,"H"), "odd","A"]
vectortypes=(v1,v2,v3)

prob=setupLP(sigma,"./Tables/D3_n3m1_L20_allspins.txt", vectortypes)

result=bissect(prob,3,1,0.1,"L=0 - S")


end

