#!/bin/bash
#make clean

nbb=("2500" "4900" "8100" "10000" "14400" "22500")
nb=("50" "70" "90" "100" "120" "150")

nbb2=( "25" "100" "400" "900" "1600" "2500" )
nb2=("5" "10" "20" "30" "40" "50" )


#./calc-anharminic-std.sh   ${nb[1]}    ${nbb[1]} 
#./calc-anharminic-h.sh    ${nb2[0]}    ${nbb2[0]} 
#./calc-anharminic-h.sh    ${nb2[1]}    ${nbb2[1]} 
#./calc-anharminic-h.sh    ${nb2[2]}    ${nbb2[2]}
#
#./calc-anharminic-std.sh    ${nb2[0]}    ${nbb2[0]} 
#./calc-anharminic-std.sh    ${nb2[1]}    ${nbb2[1]} 
#./calc-anharminic-std.sh    ${nb2[2]}    ${nbb2[2]} 
#
#
#./calc_spectrum_H.sh    ${nbb2[0]} 
#./calc_spectrum_H.sh    ${nbb2[1]} 
#./calc_spectrum_H.sh    ${nbb2[2]} 
#
#./calc_spectrum_std.sh    ${nbb2[0]} 
#./calc_spectrum_std.sh    ${nbb2[1]} 
#./calc_spectrum_std.sh    ${nbb2[2]}
#
#./calc_spectrum_std.sh    ${nbb[1]} 


#cat >  spectrum-hag-renorm=t.gp<< EOF
#reset
#set term pdf size 40,30 font "Arial, 90"
#
#set output 'spectrum-hag-renorm=t.pdf'
#set border linewidth 18.0
#
#set ylabel "spectrum " font "Arial,90"
#set xlabel "Time (ua)" font "Arial,90"
#set  title "Hag dt=0.01 renorm=T"
#
#
#plot [:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb25.txt" every 1 u 1:3 w l lw 10  t"{nb_k=5}",\
#"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb100.txt" every 1 u 1:3 w lp lw 10 pt 7 ps 1  t"{nb_k=10}", \
#"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb400.txt" every 1 u 1:3 w lp lw 10 pt 12 ps 1  t"{nb_k=20}"
#EOF
#
#
#cat >  spectrum-std-renorm=t.gp<< EOF
#reset
#set term pdf size 40,30 font "Arial, 90"
#
#set output 'spectrum-std-renorm=t.pdf'
#set border linewidth 18.0
#
#set ylabel "spectrum " font "Arial,90"
#set xlabel "Time (ua)" font "Arial,90"
#set  title " std dt=0.01 renorm=T"
#
#
#plot [:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_stdnb25.txt" every 1 u 1:3 w l lw 10  t"{nb_k=5}",\
#"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_stdnb100.txt" every 1 u 1:3 w lp lw 10 pt 7 ps 1  t"{nb_k=10}", \
#"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_stdnb400.txt" every 1 u 1:3 w lp lw 10 pt 12 ps 1  t"{nb_k=20}"
#
#EOF
#
#cat >  spectrum-std-hag-renorm=t.gp<< EOF
#reset
#set term pdf size 40,30 font "Arial, 90"
#
#set output 'spectrum-std-hag.pdf'
#set border linewidth 18.0

#set ylabel "spectrum " font "Arial,90"
#set xlabel "Time (ua)" font "Arial,90"
#
#plot [:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_stdnb4900.txt" every 1 u 1:3 w l lw 10  t"{nb_k=70}",\
#"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb100.txt" every 1 u 1:3 w lp lw 10 pt 7 ps 1  t"{nb_k=10}"
#
#EOF  



#gnuplot spectrum-hag-renorm=f.gp
#gnuplot spectrum-std-renorm=f.gp
#gnuplot spectrum-std-hag-renorm=f.gp

gnuplot spectrum-hag-renorm=t.gp
gnuplot spectrum-std-renorm=t.gp
gnuplot spectrum-std-hag-renorm=f.gp
gnuplot spectrum-std-hag-renorm=t.gp