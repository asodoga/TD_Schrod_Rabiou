reset
set term pdf size 40,30 font "Arial, 90"

set output 'con-wp-std.pdf'
set border linewidth 18.0

set ylabel "{/Symbol D}" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"


plot [-1:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb2500_std-nb4900" every 5 u 1:2 w l lw 5  t"{nb_k=50-70}","/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb4900_std-nb8100" every 6 u 1:2 w lp lw 5 pt 7 ps 1.5  t"{nb_k=70-90}", "/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb14400_std-nb22500" every 20 u 1:2 w lp lw 10 pt 12 ps 3  t"{nb_k=120-150}"




#'/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb8100_std-nb10000' every 10 u 1:2 w lp lw  10 pt 7 ps 3 t"{nb_k=90-100}",\
#'/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb10000_std-nb14400' every 15 u 1:2 w lp lw 10 pt 3 ps 3 t"{nb_k=100-120}",\
#'/Users/issa/Developpements/TD_Schrod_Rabiou/diff_std-nb14400_std-nb22500' every 20 u 1:2 w lp lw 10 pt 12 ps 3  t"{nb_k=120-150}"
   
