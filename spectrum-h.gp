reset
set term pdf size 40,30 font "Arial, 90"

set output 'spectrum-h.pdf'
set border linewidth 18.0

set ylabel "spectrum " font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"


plot [:][:] "/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb25.txt" every 1 u 1:3 w l lw 10  t"{nb_k=5}",\
"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb100.txt" every 1 u 1:3 w lp lw 10 pt 7 ps 1  t"{nb_k=10}", \
"/Users/issa/Developpements/TD_Schrod_Rabiou/file_spectrum_Hnb400.txt" every 1 u 1:3 w lp lw 10 pt 12 ps 1  t"{nb_k=20}"

   
