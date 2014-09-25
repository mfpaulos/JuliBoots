# Setting up an LP

tableFile="/Tables/eps0.5n3m1.txt"
oddSpins=false
sigma=0.6
#vectortypes=[(1,"Z") (1,"F") (1,"H") "Scalar"; (1,"F") (1/3,"F") (-5/3,"H") "Tensor"; (-1,"F") (1,"F") (-1,"H") "Vector"]




sigs=[0.51:0.005:0.51]

#Bissection

threads=1
bottom=0.6
top=3.
delta=0.01



outputFile="res.txt"
