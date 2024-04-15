set term pdfcairo
set output 'position.pdf'
set grid
set auto xy
set tit 'Position vs Time'
set xlabel 'x(t) [m]'
set ylabel 't [s]'
p 'dat-position.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Simpson Method'
