set terminal pngcairo
set output 'auto_cor_stdnb70Hnb70.png'
set title "auto_cor_stdnb70Hnb70"
set xlabel "temps"
set ylabel "auto_cor"
set boxwidth 2
 plot '~/rabiou.issa/TD_Schrod_Rabiou/results_H_nb70/auto_cor_hagedorn_taylor.txt'  w l lw 2 t 'H nb =70',\
   '~/rabiou.issa/TD_Schrod_Rabiou/results_std_nb70/auto_cor_non_hagedorn_taylor.txt'  w p pt 2 t 'std nb =70' 