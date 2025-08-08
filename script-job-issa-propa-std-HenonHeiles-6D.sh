#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --job-name=issa-std
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=AMD3
### possible values AMD (includes both AMD3 and AMD4 nodes), AMD3, AMD4, Intel (includes all of the following non-AMD nodes), SKYLAKE, HASWELL, BROADWELL, IVYBRIDGE, SANDYBRIDGE
#SBATCH --hint=nomultithread
### nomultithread => no hyperthreading i.e. each thread receives its own core

#USEFUL DIAGNOSTICS
date
hostname
free -m
df -h


module purge
module load compilers
module load software-amd
module load gcc-13.1.0-amd
module load netlib-lapack-3.10.1-gcc-12.2.0-mxhc5e6
module load  openblas-0.3.20-gcc-12.2.0-u5nwqa4

#cd ./EXT_Lib/QDUtilLib/
#make clean
#make lib
#cd ../../
#cd ./EXT_Lib/QuantumModelLib/
#make clean
#make lib
#cd ../../
# Début du chronométrage
START_TIME=$(date +%s)

time ./calc_std-HenonHeiles-6D.sh 12
# Fin du chronométrage
END_TIME=$(date +%s)

# Calcul et affichage du temps écoulé
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Temps d'exécution : $ELAPSED_TIME secondes"
./calc_spectrum_std.sh 12
#./calc_spectrum_std.sh 5
# Début du chronométrage
#START_TIME=$(date +%s)
#time ./calc_H-HenonHeiles-6D.sh 5
## Fin du chronométrage
#END_TIME=$(date +%s)
#ELAPSED_TIME=$((END_TIME - START_TIME))
#echo "Temps d'exécution : $ELAPSED_TIME secondes"
#./calc_spectrum_H.sh 5
#time ./calc_std-HenonHeiles-6D.sh 10
#./calc_spectrum_std.sh 10
#make clean
#USEFUL DIAGNOSTICS
module purge
date
hostname
free -m
df -h



