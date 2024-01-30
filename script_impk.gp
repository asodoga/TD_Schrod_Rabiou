 
 set terminal pngcairo
set output 'Imp_k.png'
 set title "Impk"
set xlabel "temps"
set ylabel "Impk"
set boxwidth 2
#plot "<paste SQt_non_hagedorn_taylor.txt  SQt_hagedorn_taylor.txt" u 1:($2-$5) w l lw 4  t 'scalQ diff '
p "Imp_k_non_hagedorn_taylor.txt" u 1:2 w l lw 2 t '1D std' , "Imp_k_hagedorn_taylor.txt" u 1:2 w l lw 2 t '1D H'
