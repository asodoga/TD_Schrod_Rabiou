gfortran -o spectrum_std.x spectre.f90

nb=$1
echo $nb
cd results_std_nb$nb

../spectrum_std.x << ** > spect_std.txt
&param Emin=0. Emax=5. file_auto='auto_cor_non_hagedorn_taylor.txt' option=2 dE=0.01 /

**

mv file_spectrum.txt ../file_spectrum_stdnb$nb.txt
cd ..
#gnuplot plot_spectrum_std.gp
