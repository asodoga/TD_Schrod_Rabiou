#!/bin/bash

rm -fr IMG-A-VAR 
#rm  *.png 
make clean
mv  results*
#for i in  5 10 20 30
#do
#./calc-harmonic-h-var.sh  $i  $i 0.25 
##./calc-harmonic-std.sh  $i  $i
#done

#for i in 5 10 20 30 
#do
# cat >  plot_pics_nb$i.gp<< EOF
#set terminal pngcairo
##set terminal postscript enhanced
#set output 'Harmonic-pics-error-nb$i.png'
#set tics font "Times-news-romain,16"
# 
#set origin 0,0
##set size 1.5,1.75
##set multiplot
#set border linewidth 2.0
#
#set ylabel "Intensity (ua)" font "Times-news-romain,14"
#set xlabel "Time (ua)" font "Times-news-romain,14"
##set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
##set xtics 6
#set tics out nomirror
#set logscale y
#set boxwidth  3
#set key right bottom
#
#plot [:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb$i/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t "{/Symbol D}{C_1}(t)",\
#'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb$i/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "{R_C}(t)"
#
##pause -1
#   
#EOF

#done
#for i in  5 10 20 30 
#
#do
#gnuplot plot_pics_nb$i.gp
#done


cat >  plot_first-pics.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-first-pics-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel " {/Symbol D}{C_1}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
#set key right bottom
#set key outside

plot [4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t  "   nb=5 ",\
                 '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=10 ",\
                 '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=20",\
                 '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=30"

#pause -1
   
EOF




cat >  plot_r-pics.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-r-pics-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{R_C}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot[4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t " nb=5 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=10 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=20 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/file_norm_pics_hagedorn_taylor.txt' u 1:3 w p pt 8 ps 0.5  t "nb=30 "


#pause -1
   
EOF

cat >  plot_n3.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-N3-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{/Symbol D}{N_{3}}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot [4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/Norm_13_hagedorn_taylor.txt' u 1:(1-$2) w l lw 2 t " nb=5 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Norm_hagedorn_taylor.txt' u 1:(1-$2)  w l lw 2 t " nb=10",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/Norm_hagedorn_taylor.txt' u 1:(1-$2)  w l lw 2 t " nb=20",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/Norm_hagedorn_taylor.txt' u 1:(1-$2)  w l lw 2 t " nb=30 "

#pause -1
   
EOF


cat >  plot_n13.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-N13-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{/Symbol D}{N_{13}}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot [4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=10",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=20",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=30 "

#pause -1
   
EOF

gnuplot  plot_first-pics.gp
gnuplot  plot_r-pics.gp
gnuplot  plot_n3.gp
gnuplot  plot_n13.gp



for dt  in  0.4 0.3 0.2 0.1  0.01
do
./calc-harmonic-h-var.sh  5  5dt$dt $dt
#./calc-harmonic-std.sh  $i  $i
done




for i in 0.4 0.3 0.2 0.1  0.01
do
 cat >  plot_pics_nb5dt$i.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-pics-error-nb5dt$i.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "Intensity (ua)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom

plot [:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt$i/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=$i {/Symbol D}{C_1}(t)",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt$i/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t " nb=5 dt=$i {R_C}(t)"

#pause -1
   
EOF
done

for i in 0.4 0.3 0.2 0.1  0.01
do
gnuplot plot_pics_nb5dt$i.gp
done




 cat >  plot_first-picsnb5.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-first-pics-errornb5.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{/Symbol D}{C_1}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot[4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=0.4 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.3/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2  t " nb=5 dt=0.3  ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2  t " nb=5 dt=0.2  ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2  t " nb=5 dt=0.1  ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=0.01 "
#pause -1
   
EOF




 cat >  plot_r-picsnb5.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-r-pics-errornb5.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{R_C}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot[4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t " nb=5 dt=0.4 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.3/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=5 dt=0.3 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=5 dt=0.2 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=5 dt=0.1 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "nb=5 dt=0.01"


#pause -1
   
EOF





 cat >  plot_n13nb5.gp<< EOF
set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-N13nb5-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
#set multiplot
set border linewidth 2.0

set ylabel "{/Symbol D}{N_{13}}(t)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
#set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
#set xtics 6
set tics out nomirror
set logscale y
set boxwidth  3
set key right bottom
#set key outside

plot [4:][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt= 0.4 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.3/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=0.3 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=0.2 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=5 dt=0.1 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 2 t"nb=5 dt=0.01 "


#pause -1
   
EOF

gnuplot  plot_first-picsnb5.gp
gnuplot  plot_r-picsnb5.gp
gnuplot  plot_n3nb5.gp
gnuplot  plot_n13nb5.gp


 cat >  plot_ralphanb10nb0.25.gp<< EOF

set terminal pngcairo
set output 'ralpha-plotnb10dtO.25.png'
set tics font "Times-news-romain,12"
 
set origin 0,0
#set size 1.5,1.0
#set multiplot
set border linewidth 2.0

set ylabel "At" font "Times-news-romain,12"
set xlabel "Time (ua)" font "Times-news-romain,12"
#set xtics 6

set tics out nomirror


set boxwidth  3
set style fill solid 0.5 noborder
plot [-1 :][0.5:1.8] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' u 1:2 w l lw 2 t "nb=10 dt=0.25 " ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' u 1:3 w l lw 2  t " nb=10 dt=0.25  "



#pause -1

   
EOF



 cat >  plot_ialphanb10dt0.25.gp<< EOF

set terminal pngcairo
set output 'ialpha-plotnb10dt0.25.png'
set tics font "Times-news-romain,12"
 
set origin 0,0
#set size 1.5,1.0
set multiplot
set border linewidth 2.0

set ylabel "Bt" font "Times-news-romain,12"
set xlabel "Time (ua)" font "Times-news-romain,12"
#set xtics 6

set tics out nomirror


set boxwidth  3
set style fill solid 0.5 noborder
plot [-1:][-0.6:0.6] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' u 1:4 w l lw 2 t " nb=10 dt=0.25" ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' u 1:5 w  l lw 2 t " nb=10 dt=0.25  "



#pause -1

   
EOF



 cat >  plot_avqb10dt0.25.gp<< EOF

set terminal pngcairo
set output 'avqnb10dt0.25.png'
set tics font "Times-news-romain,12"
 
set origin 0,0
#set size 1.5,1.0
set multiplot
set border linewidth 2.0

set ylabel " {Qt}" font "Times-news-romain,12"
set xlabel "Time (ua)" font "Times-news-romain,12"
#set xtics 6

set tics out nomirror


set boxwidth  3
set style fill solid 0.5 noborder
plot [-0.5:][-2.5:3] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Qt_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=10 dt=0.25 " ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Qt_hagedorn_taylor.txt' u 1:3 w l lw 2 t " nb=10 dt=0.25  "



#pause -1

   
EOF



 cat >  plot_avpb10dt0.25.gp<< EOF

set terminal pngcairo
set output 'avpnb10dt0.25.png'
set tics font "Times-news-romain,12"
 
set origin 0,0
#set size 1.5,1.0
set multiplot
set border linewidth 2.0

set ylabel " {Pt}" font "Times-news-romain,12"
set xlabel "Time (ua)" font "Times-news-romain,12"
set xtics 6

set tics out nomirror


set boxwidth  3
set style fill solid 0.5 noborder
plot [-0.5:][-2.5:3] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/imp_k_hagedorn_taylor.txt' u 1:2 w l lw 2 t " nb=10 dt=0.25  " ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/imp_k_hagedorn_taylor.txt' u 1:3 w l lw 2 t " nb=10 dt=0.25 "

#pause -1

   
EOF


gnuplot plot_ralphanb10nb0.25.gp
gnuplot plot_ialphanb10dt0.25.gp

gnuplot plot_avqb10dt0.25.gp
gnuplot plot_avpb10dt0.25.gp
#gnuplot plot_energy-er.gp
#gnuplot plot_energy-ernb5.gp

mkdir -p IMG-A-VAR
mv  *.png ./IMG-A-VAR/.
#mv  results* ./IMG-A-VAR/.
