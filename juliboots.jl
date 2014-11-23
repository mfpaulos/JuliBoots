module juliboots
println("JuliBootS: a package for numerical bootstrap computations.")
println("Version: 1.0")

using various
using main

export chooseTable, setupLP, bissect, value, dropOdd!, changeTarget!,dropEven!,opemax #main functions
export mcopy
export filter,iterate!,status,cost,solution,makeVector # LP routines, this makes them accessible when 'using main'

end

