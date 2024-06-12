gfortran -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -o diff.x diff.f90
nb1=$1
nb2=$2
 echo $nb1,$nb2
./diff.x << ** > diff.txt
 &res nbA=$nb1 nbB=$nb2 stA=t stB=f maxit=6000/
**