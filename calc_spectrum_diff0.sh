#!/bin/bash

nb=$1
echo $nb

paste  file_spectrum_stdnb70.txt file_spectrum_Hnb$nb.txt |awk '{print $1,"    ", $3-$8}' > file_spectrum_stdnb70_Hnb$nb.txt