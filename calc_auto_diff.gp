 
set terminal pngcairo
set output 'diff-auto.png'
 set title "auto-func"
set xlabel "temps"
set ylabel "auto-func"
set boxwidth 2


#plot "<paste ~/rabiou.issa/TD_Schrod_Rabiou/results_H_nb70/auto_cor_hagedorn_taylor.txt \
#~/rabiou.issa/TD_Schrod_Rabiou/results_H_nb70/auto_cor_hagedorn_taylor.txt" u 1:($2-$4) w l lw 4  t 'norm diff'
#plot  "aut_diff.txt" u 1:($2-$3) w l lw 3 t 'diff aut'
plot  "aut_diff.txt" u 1:2 w l lw 3 t 'H nb=70', "aut_diff.txt" u 1:3 w p pt 3  t 'std nb=70'