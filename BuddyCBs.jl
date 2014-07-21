using table,cb

#this loads stuff
sigma=BigFloat("0.125")
tab=table.loadTable("../Tables/migN2_eps0.00001n11m1.txt",oddL=true)
convolved=cb.convTable(sigma,tab.table,-1) # Convolve; the minus sign means we want v^sigma G-u^sigma G; the plus is required for global sym stuff


#convolved is a table of vector functions, one for each spin.

convolved[1].spin

cb.value(convolved[1],BigFloat(0.3))         #evaluate L=0 at 0.3, gives all components

cb.value(convolved[1][(3,4)],BigFloat(0.3))         #evaluate L=0 at 0.3, (3,4) derivative

#or accessing components directly

cb.value(convolved[1][17],BigFloat(0.3))

#the dictionary is included in the object
convolved.dict



