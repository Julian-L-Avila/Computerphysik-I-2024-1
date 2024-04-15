set term pdfcairo
set output 'velocity.pdf'
set grid
set auto xy
set tit 'Velocity vs Time'
set ylabel 'v(t) [ms^{-1}]'
set xlabel 't [s]' 
p 'dat-velocity.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Euler Method'
