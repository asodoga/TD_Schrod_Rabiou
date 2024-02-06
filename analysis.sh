#!/bin/bash
./calc_spectrum_std.sh 70

 cat >  plot_spectrum_stdnb70.gp<< EOF

set terminal pngcairo
set output 'spectrum_stdnb70.png'
set title "spectrum_std"
set xlabel "E[ua]"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb70.txt'  u 1:3 w l lw 2 t 'spectrum_stdnb70 op=3'
   
EOF

gnuplot   plot_spectrum_stdnb70.gp
/usr/bin/rm -r   plot_spectrum_stdnb70.gp



for i in  10 15 25 35

do
./calc_spectrum_H.sh  $i
done


for i in  10 15 25 35
 do

 cat >  plot_spectrum_Hnb$i.gp<< EOF

set terminal pngcairo
set output 'spectrum_Hnb$i.png'
set title "spectrum_H"
set xlabel "E[ua]"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_Hnb$i.txt'  u 1:3 w l lw 2 t 'spectrum_Hnb$i op=3'
   
EOF

gnuplot   plot_spectrum_Hnb$i.gp
/usr/bin/rm -r  plot_spectrum_Hnb$i.gp

   done



for  i in  10 15 25 35
 do

 cat >  plot_spectrum_stdnb70_Hnb$i.gp<< EOF

set terminal pngcairo
set output 'spectrum_stdnb70_Hnb$i.png'
set title "spectrum_diff"
set xlabel "E[ua]"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb70.txt'  w l lw 1 t 'file_spectrum_stdnb70.txt op=3',\
 'file_spectrum_Hnb$i.txt'  u 1:3 w p pt 7 ps 0.6  t'spectrum_H$i op=3'
   
EOF

gnuplot   plot_spectrum_stdnb70_Hnb$i.gp
/usr/bin/rm -r  plot_spectrum_stdnb70_Hnb$i.gp

   done

for i in  10 15 25 35
do
paste  file_spectrum_stdnb70.txt file_spectrum_Hnb$i.txt |awk '{print $1,"    ", $3-$8}' > file_spectrum_stdnb70_Hnb$i.txt
done


for i in 10 15 25 35
do
cat > plot_spectrum_diff_stdnb70_Hnb$i.gp<< EOF

set terminal pngcairo
set output 'spectrum_diff_stdnb70_Hnb$i.png'
set title "spectrum_diff"
set xlabel "E[ua]"
set ylabel "spectrum"
set boxwidth 2
plot 'file_spectrum_stdnb70_Hnb$i.txt'  w l lw 1 t 'spectrum_diff_stdnb70_Hnb$i op=3'

EOF

gnuplot plot_spectrum_diff_stdnb70_Hnb$i.gp
/usr/bin/rm -r plot_spectrum_diff_stdnb70_Hnb$i.gp
/usr/bin/rm -r file_spectrum_stdnb70_Hnb$i.txt

done



/usr/bin/rm -r file_spectrum_stdnb70.txt

for i in 10 15 25 35
do

/usr/bin/rm -r file_spectrum_Hnb$i.txt

done
