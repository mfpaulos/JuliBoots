module tests

push!(LOAD_PATH,"$(pwd())/..")
using various
using qfunc

bf=BigFloat
# Polynomials

#c=[1.,2.,3.]
c=[bf(1),bf(2),bf(3)]


# Division
p=Polynomial(c)


p2=p*p
(q,r)=p2/p
if q==p println("check.") end



poly=Polynomial([bf(4),bf(6)])
pole=Pole(1,BigFloat(2),BigFloat(3))
q=poly+pole


end
