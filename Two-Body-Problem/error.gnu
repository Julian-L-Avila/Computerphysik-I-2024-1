set term pdfcairo
set output 'error.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set auto xy
set logscale y
set tit 'Relative Percentage Error (Compared to Rectangular Rule)
p 'databserror.dat' u 1:2 w l tit 'Trapezoid', '' u 1:3 w l tit 'Simpson', '' u 1:4 w l tit 'Gauss'
