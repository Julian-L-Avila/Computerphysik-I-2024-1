set term pdfcairo
set output 'velocity-Taylor.pdf'
set grid
set ylabel 'v [ms^{-1}]'
set xlabel 't [s]'
set auto xy
set tit 'Taylor Method'
p 'dat-velocity-Taylor.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'
unset term
rep
pause -1
set term pdfcairo
set output 'error-velocity-Taylor.pdf'
set ylabel 'Error %'
set logscale y
set auto xy
set tit 'Error with Taylor Method'
p 'dat-velocity-Taylor.dat' u 1:4 w l tit 'Error'
unset term
rep
pause -1
