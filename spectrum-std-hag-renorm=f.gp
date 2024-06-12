reset
set term pdf size 40,30 font "Arial, 90"

set output 'spectrum-std-hag-renorm=f.pdf'
set border linewidth 18.0

set ylabel "spectrum " font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set  title "std-Hag dt=0.01 renorm=F"


plot [:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_stdnb4900.txt" every 1 u 1:3 w l lw 10  t"{ std nb_k=70}",\
"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb100.txt" every 1 u 1:(-$3) w lp lw 10 pt 7 ps 1  t"{Hag nb_k=10}" 

