module sample
using various
using main

sigma=BigFloat(0.518)
prob=setupLP(sigma,"./Tables/eps0.5n3m1_allspins_prec80.txt")


prob2=main.mcopy(prob)
function dropOdd!(prob)
        ll=length(prob.lpFunctions)
        prob.lpFunctions=[prob.lpFunctions[2i-1] for i=1:floor(ll/2)]
        prob
end

prob2=dropOdd!(prob2)
end
