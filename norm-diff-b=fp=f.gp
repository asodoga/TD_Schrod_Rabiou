
 reset
set term pdf size 40,30 font "Arial, 90"

set output 'norm-diff-std-hag-b=f-p=f.pdf'
set border linewidth 18.0
 set  title "dt=0.01 B=F P=F"
 
set ylabel "{/Symbol D}N(psi-std-psi-hag)" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set  logscale y
set key right bottom

plot[:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/henonheiles-b=f-p=f/diff_std-nb4900_H-nb400" every 1 u 1:2 w lp lw  5 pt 7 ps 1 t"nb=70-20" , "/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/henonheiles-b=f-p=f/diff_std-nb4900_H-nb900" every 10 u 1:2 w lp lw 5 pt 10 ps 2 t "nb=70-30", "/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/henonheiles-b=f-p=f/diff_std-nb4900_H-nb1600" every 10 u 1:2 w lp lw  5 pt 12 ps 1 t "nb=70-40", "/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/henonheiles-b=f-p=f/diff_std-nb4900_H-nb2500" every 10 u 1:2 w lp lw  5 pt 7 ps 1 t "nb=70-50",  "/Users/issa/Developpements/TD_Schrod_Rabiou/henonheiles-doc/henonheiles-b=f-p=f/diff_std-nb4900_H-nb4900" every 15 u 1:2 w lp lw 5 pt 5 ps 1 t "nb=70-70"

