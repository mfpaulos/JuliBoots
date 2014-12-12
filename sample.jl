module sample
using juliboots


sigma=BigFloat(0.50001)
prob=setupLP(sigma,"./Tables/D3_n3m1_L20.txt")

result=bissect(prob2,3,1,1e-3,"L=0",method="mcv")


end
