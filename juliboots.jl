module juliboots
println("JuliBootS: a package for numerical bootstrap computations.")
println("Version: 1.0")

using various
using main
using cb
using LP
using bb
using qfunc
using table


export chooseTable, setupLP, bissect, value, dropOdd!, changeTarget!,dropEven!,opemax, avgSpec #main functions
export mcopy
export filter,filter!,iterate!,status,cost,solution,makeVector # LP routines, this makes them accessible when 'using main'
export QFunc, Polynomial, Pole
export MinFunction, FindMinimum, FindLocalMinima
export loadTable, convTable
export convBlock
export findMRC, findBVar

end

