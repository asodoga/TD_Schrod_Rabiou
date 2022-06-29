#!/bin/bash
#fpm build --flag  "/home/elprof/QuantumModelLib/libpot.a"
fpm build
fpm run --<DAT_files/dat_dia>resultat

