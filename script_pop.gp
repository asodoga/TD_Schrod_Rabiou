set yrange[0:1]
p 'pop_hagedorn_taylor.dat' u 1:2 w l lw 2 t'pop[1] ,H'
rep 'pop_hagedorn_taylor.dat' u 1:3 w l lw 2 t'pop[2],H'
rep 'pop_non_hagedorn_taylor.dat' u 1:2 w p pt 6 ps 0.25 t'pop[1],NH'
rep 'pop_non_hagedorn_taylor.dat' u 1:3 w p pt 6 ps 0.25 t'pop[2],NH'

