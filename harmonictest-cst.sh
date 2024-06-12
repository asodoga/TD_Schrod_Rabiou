#!/bin/bash
rm -rf  results*
rm   -rf ./IMG-A-CST
make clean
for i in  5 10 20 30
do
./calc-harmonic-h.sh  $i  $i 0.25 
#./calc-harmonic-std.sh  $i  $i
done


for dt  in  0.4 0.3 0.2 0.1  0.01
do
./calc-harmonic-h.sh  5  5dt$dt $dt 
#./calc-harmonic-std.sh  $i  $i
done


gnuplot plot_avqb10.gp
gnuplot plot_avpb10.gp
gnuplot plot_ralphanb10.gp
gnuplot plot_ialphanb10.gp
gnuplot plot_r-pics.gp
gnuplot plot_r-picsnb5.gp
gnuplot plot_n3.gp
gnuplot plot_n3nb5.gp
gnuplot plot_n13.gp
gnuplot plot_n13nb5.gp
mkdir -p IMG-A-CST
mv  *.pdf ./IMG-A-CST/.



rm -rf   results*
make clean



for i in  5 10 20 30
do
./calc-harmonic-h-var.sh  $i  $i 0.25 
#./calc-harmonic-std.sh  $i  $i
done


for dt  in  0.4 0.3 0.2 0.1  0.01
do
./calc-harmonic-h-var.sh  5  5dt$dt $dt 
#./calc-harmonic-std.sh  $i  $i
done


gnuplot plot_avqb10.gp
gnuplot plot_avpb10.gp
gnuplot plot_ralphanb10.gp
gnuplot plot_ialphanb10.gp
gnuplot plot_r-pics.gp
gnuplot plot_r-picsnb5.gp
gnuplot plot_n3.gp
gnuplot plot_n3nb5.gp
gnuplot plot_n13.gp
gnuplot plot_n13nb5.gp
mkdir -p IMG-A-VAR
mv  *.pdf ./IMG-A-VAR/.

#rm -rf   results*
make clean

