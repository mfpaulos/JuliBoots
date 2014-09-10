module sample
using various
using main

sigma=BigFloat(0.518)

# Define types of vectors and their makeup
v1=[(1,"Z") (1,"F") (1,"H") "Scalar"]
v2=[(1,"F") (1/3,"F") (-5/3,"H") "Tensor"]
v3=[(-1,"F") (1,"F") (-1,"H") "Vector"]
vectortypes=[v1,v2,v3]

prob=setupLP(sigma,"./Tables/eps0.5n3m1_allspins_prec80.txt", vectortypes)
prob2=main.mcopy(prob)
prob2=dropOdd!(prob2)

result=bissect(prob2,3,1,0.1,"L=0")


end
