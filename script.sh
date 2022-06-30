#!/bin/bash
fpm build --flag  "-fopenmp /home/elprof/QuantumModelLib/libpot.a"
fpm run --flag    "-fopenmp /home/elprof/QuantumModelLib/libpot.a" --<DAT_files/dat_dia>resultat

