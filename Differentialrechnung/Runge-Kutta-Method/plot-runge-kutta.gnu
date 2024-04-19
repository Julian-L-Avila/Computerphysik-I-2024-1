set term pdfcairo
set output 'runge-kutta.pdf'
set grid
set ylabel 'v [ms^{-1}]'
set xlabel 't [s]'
set auto xy
p 'dat-runge-kutta.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'
set output 'error-runge-kutta.pdf'
set ylabel 'Error %'
set logscale y
set auto xy
p 'dat-runge-kutta.dat' u 1:4 w l tit 'Error'
