#!/bin/bash
make clean
#make
#rm -rf    *.x  *.txt
#rm -rf    results*  


nbb=("2500" "4900" "8100" "12100")
nb=("50" "70" "90" "110" )

nbb2=("25" "100" "400" "2500" "4900" )
nb2=("5" "10" "20" "50" "70" )



./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
./calc-anharminic-std.sh    5 25
./calc-anharminic-std.sh    10 100
./calc-anharminic-std.sh    20 400

./calc-anharminic-h.sh    5 25
./calc-anharminic-h.sh    10 100
./calc-anharminic-h.sh    20 400
#for i in  0 1 2# 3 
#do
#./calc-anharminic-std.sh    ${nb2[i]}    ${nbb2[i]}
#done

#for i in  0 1 2 #3 4  #0 1 2  
#do
#./calc-anharminic-h.sh    ${nb2[i]}    ${nbb2[i]}
#done


./calc_spectrum_H.sh    ${nbb2[0]} 
./calc_spectrum_H.sh    ${nbb2[1]} 
./calc_spectrum_H.sh    ${nbb2[2]} 

./calc_spectrum_std.sh    ${nbb2[0]} 
./calc_spectrum_std.sh    ${nbb2[1]} 
./calc_spectrum_std.sh    ${nbb2[2]}
./calc_spectrum_std.sh    ${nbb[1]} 

#mv  *.txt    ./henon-spec/doc-std/.
#
#cd plot-specdt01
#gnuplot spec-b=f-p=t-renorm=f.gp
#gnuplot spec-b=f-p=t-renorm=t.gp
#gnuplot spec-b=t-p=t-renorm=f.gp
#gnuplot spec-b=t-p=t-renorm=t.gp
#gnuplot spec-std.gp
#gnuplot spec-std-hag-b=t-p=t-renorm=f.gp
#gnuplot spec-std-hag-b=t-p=t-renorm=t.gp
#gnuplot spec-std-hag-b=f-p=t-renorm=f.gp
#gnuplot spec-std-hag-b=f-p=t-renorm=t.gp
#
#mkdir -p ./img-henon-spec
#mv  *.png ./img-henon-spec/


#
#mkdir -p ./henon-spec/doc-b=f-p=t-renorm=t
#mv  *.txt   ./henon-spec/doc-b=f-p=t-renorm=t/.
#cp -rf  results_std_nb4900  ./henon-spec/doc-b=f-p=t-renorm=f/.
#mv -f  results_H_nb*        ./henon-spec/doc-b=f-p=t-renorm=f/.
#mv -f  results_std_nb25     ./henon-spec/doc-b=f-p=t-renorm=f/.
#mv -f  results_std_nb100     ./henon-spec/doc-b=f-p=t-renorm=f/.
#mv -f  results_std_nb400     ./henon-spec/doc-b=f-p=t-renorm=f/.
###
#./calc_diff.sh    ${nb[1]}    ${nb2[0]} 
#./calc_diff.sh    ${nb[1]}    ${nb2[1]}
#./calc_diff.sh    ${nb[1]}    ${nb2[2]}
#./calc_diff.sh    ${nb[1]}    ${nb2[3]}
#./calc_diff.sh    ${nb[1]}    ${nb2[4]} 
#
#
#mv  diff_std*    ./henon-doc/doc-b=t-p=t-renorm=t/.
#cp -rf  results_std_nb4900    ./henon-doc/doc-b=t-p=t-renorm=t/.
#mv -f  results_H_nb*    ./henon-doc/doc-b=t-p=t-renorm=t/.

# ----------------  norm & norm-diff dt = 0.1-----------------------------
#cd plot-normdt=01-diff
#
#gnuplot norm-diff-b=f-p=t-renorm=f.gp
#gnuplot norm-diff-b=f-p=t-renorm=t.gp
#
#gnuplot norm-diff-b=f-p=f-renorm=f.gp
#gnuplot norm-diff-b=f-p=f-renorm=t.gp
#
#gnuplot norm-diff-b=t-p=t-renorm=f.gp
#gnuplot norm-diff-b=t-p=t-renorm=t.gp
#
#
#gnuplot  norm-diff-b=f-p=f-renorm=tnb10.gp
#gnuplot  norm-diff-b=f-p=f-renorm=tnb5.gp
#gnuplot  norm-diff-b=f-p=f-renorm=fnb5.gp
#gnuplot  norm-diff-b=f-p=f-renorm=fnb10.gp
#
#
#mkdir -p ./img-henon-norm-diff
#mv  *.png ./img-henon-norm-diff/.
#
#cd ..
#
#cd plot-norm13dt=01
#gnuplot norm13-b=f-p=t-renorm=f.gp
#gnuplot norm13-b=f-p=t-renorm=t.gp
#
#gnuplot norm13-b=f-p=f-renorm=f.gp
#gnuplot norm13-b=f-p=f-renorm=t.gp
#
#gnuplot norm13-b=t-p=t-renorm=f.gp
#gnuplot norm13-b=t-p=t-renorm=t.gp
#
#gnuplot  norm13-b=f-p=f-renorm=tnb10.gp
#gnuplot  norm13-b=f-p=f-renorm=tnb5.gp
#gnuplot  norm13-b=f-p=f-renorm=fnb5.gp
#gnuplot  norm13-b=f-p=f-renorm=fnb10.gp
#
#mkdir -p ./img-henon-norm13
#mv  *.png ./img-henon-norm13/.
#
#cd ..
##---------energy dt = 0.1------------------
#
#cd plot-energydt=01
#gnuplot energy-b=f-p=t-renorm=f.gp
#gnuplot energy-b=f-p=t-renorm=t.gp
##
#gnuplot energy-b=f-p=f-renorm=f.gp
#gnuplot energy-b=f-p=f-renorm=t.gp
##
#gnuplot energy-b=t-p=t-renorm=f.gp
#gnuplot energy-b=t-p=t-renorm=t.gp
#
#gnuplot  energy-b=f-p=f-renorm=tnb10.gp
#gnuplot  energy-b=f-p=f-renorm=tnb5.gp
#gnuplot  energy-b=f-p=f-renorm=fnb5.gp
#gnuplot  energy-b=f-p=f-renorm=fnb10.gp
##
#mkdir -p ./img-henon-energy
#mv  *.png ./img-henon-energy/.
#
#cd ..
##---------norm3 dt = 0.1------------------
#
#
#cd plot-norm3dt=01
#gnuplot norm-b=f-p=t-renorm=f.gp
#gnuplot norm-b=f-p=t-renorm=t.gp
##
#gnuplot norm-b=f-p=f-renorm=f.gp
#gnuplot norm-b=f-p=f-renorm=t.gp
##
#gnuplot norm-b=t-p=t-renorm=f.gp
#gnuplot norm-b=t-p=t-renorm=t.gp
##
#gnuplot  norm-b=f-p=f-renorm=tnb10.gp
#gnuplot  norm-b=f-p=f-renorm=tnb5.gp
#gnuplot  norm-b=f-p=f-renorm=fnb5.gp
#gnuplot  norm-b=f-p=f-renorm=fnb10.gp
#
#mkdir -p ./img-henon-norm3
#mv  *.png ./img-henon-norm3/.
#


#rm -r img-henon-norm13
#rm -r  img-henon-norm-diff 
#mkdir -p ./img-henon-norm-diff 
#mkdir -p ./img-henon-norm13
#cd ./plot-norm-diff
#gnuplot  norm-diff-b=t-p=f-renorm=tnb10.gp
#gnuplot  norm-diff-b=t-p=f-renorm=tnb5.gp
#gnuplot  norm-diff-b=t-p=f-renorm=fnb5.gp
#gnuplot  norm-diff-b=t-p=f-renorm=fnb10.gp
#
#gnuplot  norm-diff-b=f-p=f-renorm=tnb10.gp
#gnuplot  norm-diff-b=f-p=f-renorm=tnb5.gp
#gnuplot  norm-diff-b=f-p=f-renorm=fnb5.gp
#gnuplot  norm-diff-b=f-p=f-renorm=fnb10.gp
#
#gnuplot  norm-diff-b=t-p=f-renorm=t.gp
#gnuplot  norm-diff-b=t-p=f-renorm=f.gp
#
#gnuplot  norm-diff-b=t-p=t-renorm=t.gp
#gnuplot  norm-diff-b=t-p=t-renorm=f.gp
#
#gnuplot  norm-diff-b=f-p=f-renorm=t.gp
#gnuplot  norm-diff-b=f-p=f-renorm=f.gp
#
#gnuplot  norm-diff-b=f-p=t-renorm=t.gp
#gnuplot  norm-diff-b=f-p=t-renorm=f.gp
#
#mv  *.png ../img-henon-norm-diff/. 
#cd ..
#cd ./plot-norm13
#
#gnuplot  norm13-b=f-p=f-renorm=f.gp
#gnuplot  norm13-b=f-p=f-renorm=t.gp
#gnuplot  norm13-b=t-p=t-renorm=f.gp
#gnuplot  norm13-b=t-p=t-renorm=t.gp
#
#gnuplot  norm13-b=t-p=f-renorm=f.gp
#gnuplot  norm13-b=t-p=f-renorm=t.gp
#gnuplot  norm13-b=f-p=t-renorm=f.gp
#gnuplot  norm13-b=f-p=t-renorm=t.gp
#
#gnuplot  norm13-b=f-p=f-renorm=fnb5.gp
#gnuplot  norm13-b=f-p=f-renorm=tnb5.gp
#gnuplot  norm13-b=t-p=f-renorm=fnb5.gp
#gnuplot  norm13-b=t-p=f-renorm=tnb5.gp
#
#gnuplot  norm13-b=f-p=f-renorm=fnb10.gp
#gnuplot  norm13-b=f-p=f-renorm=tnb10.gp
#gnuplot  norm13-b=t-p=f-renorm=fnb10.gp
#gnuplot  norm13-b=t-p=f-renorm=tnb10.gp
#
#mv  *.png ../img-henon-norm13/. 
#cd ..
#cp -fr plot-norm13 plot-e13
#cd ./plot-e13
#gnuplot  norm13-b=f-p=f-renorm=f.gp
#gnuplot  norm13-b=f-p=f-renorm=t.gp
#gnuplot  norm13-b=t-p=t-renorm=f.gp
#gnuplot  norm13-b=t-p=t-renorm=t.gp
##
#gnuplot  norm13-b=t-p=f-renorm=f.gp
#gnuplot  norm13-b=t-p=f-renorm=t.gp
#gnuplot  norm13-b=f-p=t-renorm=f.gp
#gnuplot  norm13-b=f-p=t-renorm=t.gp
#
#gnuplot  e3-b=f-p=f-renorm=fnb5.gp
#gnuplot  e3-b=f-p=f-renorm=tnb5.gp
#gnuplot  e3-b=t-p=f-renorm=fnb5.gp
#gnuplot  e3-b=t-p=f-renorm=tnb5.gp
##
#gnuplot  e3-b=f-p=f-renorm=fnb10.gp
#gnuplot  e3-b=f-p=f-renorm=tnb10.gp
#gnuplot  e3-b=t-p=f-renorm=fnb10.gp
#gnuplot  e3-b=t-p=f-renorm=tnb10.gp



#./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-std.sh    ${nb[2]}    ${nbb[2]} 


#./calc_diff.sh    ${nb[0]}    ${nb[1]} 
#./calc_diff.sh    ${nb[1]}    ${nb[2]}
#./calc_diff.sh    ${nb[2]}    ${nb[3]} 
#./calc_diff.sh    ${nbb[3]}    ${nbb[3]} 
#./calc_diff.sh    ${nb[4]}    ${nb[5]}  

#./calc-anharminic-h.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-h.sh    ${nb[2]}    ${nbb[2]} 
#./calc-anharminic-h.sh    ${nb[3]}    ${nbb[3]} 
#./calc-anharminic-h.sh    ${nb[4]}    ${nbb[4]} 
#./calc-anharminic-h.sh    ${nb[5]}    ${nbb[5]}

#nbb2=("400" "900" "1600" "2500" )
#nb2=("20" "30" "40" "50" )
#
#
#./calc-anharminic-h.sh    ${nb2[0]}    ${nbb2[0]} 
#./calc-anharminic-h.sh    ${nb2[1]}    ${nbb2[1]} 
#./calc-anharminic-h.sh    ${nb2[2]}    ${nbb2[2]}
#./calc-anharminic-h.sh    ${nb2[3]}    ${nbb2[3]} 
#
#
#
##./calc-anharminic-std.sh    ${nb2[0]}    ${nbb2[0]} 
##./calc-anharminic-std.sh    ${nb2[1]}    ${nbb2[1]} 
#./calc-anharminic-std.sh    ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-h.sh    ${nb[1]}    ${nbb[1]}
#
#
#./calc_diff.sh    ${nb[1]}    ${nb2[0]} 
#./calc_diff.sh    ${nb[1]}    ${nb2[1]} 
#./calc_diff.sh    ${nb[1]}    ${nb2[2]} 
#./calc_diff.sh    ${nb[1]}    ${nb2[3]}
#./calc_diff.sh    ${nb[1]}    ${nb[1]}

#gnuplot norm-diff-std-hag.gp 


#./calc_spectrum_H.sh    ${nbb2[0]} 
#./calc_spectrum_H.sh    ${nbb2[1]} 
#./calc_spectrum_H.sh    ${nbb2[2]} 
#
#./calc_spectrum_std.sh    ${nbb2[0]} 
#./calc_spectrum_std.sh    ${nbb2[1]} 
#./calc_spectrum_std.sh    ${nbb2[2]} 




