 
 set terminal pngcairo
set output 'diff_std-nb70_std-nb100.png'
 set title "norm-diff"
set xlabel "temps"
set ylabel "norm-diff"
set boxwidth 2


plot "diff_std-nb70_std-nb100"  w l lw 4 
