 reset
set term pdf size 40,30 font "Arial, 90"

set output 'norm-diff.pdf'
set border linewidth 18.0

set ylabel "{/Symbol D}N" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"



plot [:][:]"/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb4900_H-nb4900" every 20 u 1:2 w lp lw  10 pt 7 ps 3 t"{/Symbol D}t=0.1",\
#'/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/diff_std-nb4900_H-nb4900' every 15 u 1:2 w lp lw 10 pt 8 ps 3 t"{/Symbol D}t=0.1",\
