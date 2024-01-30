 
 set terminal pngcairo
set output 'SQdiff.png'
 set title "scaleQ diff"
set xlabel "temps"
set ylabel "SQt_std-SQt_H"
set boxwidth 2
plot "<paste SQt_non_hagedorn_taylor.txt  SQt_hagedorn_taylor.txt" u 1:($2-$5) w l lw 4  t 'scalQ diff '
 #p "SQt_non_hagedorn_taylor.txt" u 1:2 w l lw 2 t '1D std' , "SQt_hagedorn_taylor.txt" u 1:2 w l lw 2 t '1D H'
