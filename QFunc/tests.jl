module tests

push!(LOAD_PATH,"$(pwd())/..")
using various
using qfunc

# Polynomials

p=Polynomial([1.,2.,3.])
p2=p*p

(q,r)=p2/p


end
