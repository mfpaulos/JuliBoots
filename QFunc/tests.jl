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



q=Polynomial([4.,6.])





end
