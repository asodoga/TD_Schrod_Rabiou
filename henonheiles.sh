#!/bin/bash
make clean
#make
#rm -rf    *.x  *.txt
#rm -rf    results*  


nbb=("2500" "4900" "8100" "10000" "14400" "22500")
nb=("50" "70" "90" "100" "120" "150")

#for i in  0 1 2 3 4 5
#do
#./calc-anharminic-std.sh    ${nb[i]}    ${nbb[i]} 
#done
#./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-std.sh    ${nb[2]}    ${nbb[2]} 


#./calc_diff.sh    ${nb[0]}    ${nb[1]} 
#./calc_diff.sh    ${nb[1]}    ${nb[2]}
#./calc_diff.sh    ${nb[2]}    ${nb[3]} 
#./calc_diff.sh    ${nb[3]}    ${nb[4]} 
#./calc_diff.sh    ${nb[4]}    ${nb[5]}  

#./calc-anharminic-h.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-h.sh    ${nb[2]}    ${nbb[2]} 
#./calc-anharminic-h.sh    ${nb[3]}    ${nbb[3]} 
#./calc-anharminic-h.sh    ${nb[4]}    ${nbb[4]} 
#./calc-anharminic-h.sh    ${nb[5]}    ${nbb[5]}

nbb2=("400" "900" "1600" "2500" )
nb2=("20" "30" "40" "50" )


./calc-anharminic-h.sh    ${nb2[0]}    ${nbb2[0]} 
./calc-anharminic-h.sh    ${nb2[1]}    ${nbb2[1]} 
./calc-anharminic-h.sh    ${nb2[2]}    ${nbb2[2]}
./calc-anharminic-h.sh    ${nb2[3]}    ${nbb2[3]} 



#./calc-anharminic-std.sh    ${nb2[0]}    ${nbb2[0]} 
#./calc-anharminic-std.sh    ${nb2[1]}    ${nbb2[1]} 
./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
./calc-anharminic-h.sh    ${nb[1]}    ${nbb[1]}


./calc_diff.sh    ${nb[1]}    ${nb2[0]} 
./calc_diff.sh    ${nb[1]}    ${nb2[1]} 
./calc_diff.sh    ${nb[1]}    ${nb2[2]} 
./calc_diff.sh    ${nb[1]}    ${nb2[3]}
./calc_diff.sh    ${nb[1]}    ${nb[1]}

gnuplot norm-diff-std-hag.gp 


#./calc_spectrum_H.sh    ${nbb2[0]} 
#./calc_spectrum_H.sh    ${nbb2[1]} 
#./calc_spectrum_H.sh    ${nbb2[2]} 
#
#./calc_spectrum_std.sh    ${nbb2[0]} 
#./calc_spectrum_std.sh    ${nbb2[1]} 
#./calc_spectrum_std.sh    ${nbb2[2]} 




