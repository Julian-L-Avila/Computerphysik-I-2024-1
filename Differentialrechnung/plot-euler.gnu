set term pdfcairo
set output 'velocity-euler.pdf'
set grid
set ylabel 'v [ms^{-1}]'
set xlabel 't [s]'
set auto xy
p 'dat-velocity-euler.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'
set output 'error-velocity-euler.pdf'
set ylabel 'Error %'
set logscale y
set auto xy
p 'dat-velocity-euler.dat' u 1:4 w l tit 'Error'
