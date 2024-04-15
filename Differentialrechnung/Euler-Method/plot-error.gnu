set term pdfcairo
set output 'error.pdf'
set grid
set logscale y
set tit 'Velocity and Position Error'
set xlabel 't [s]'
set ylabel 'Error %'
p 'dat-error.dat' u 1:2 w l tit 'Velocity Error', '' u 1:3 w l tit 'Position Error'
