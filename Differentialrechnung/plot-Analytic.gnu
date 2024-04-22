set term pdfcairo
set output 'velocity-Analytic.pdf'
set grid
set ylabel 'v [ms^{-1}]'
set xlabel 't [s]'
set auto xy
p 'dat-velocity-Analytic.dat' u 1:2 w l tit 'Analytic'
