set term pdfcairo
set output 'Error.pdf'
set grid
set xlabel 'N'
set ylabel 'error (%)'
set logscale x
set logscale y
p 'results.dat' u 1:6 w l tit 'Rectangular', 'results.dat' u 1:7 w l tit 'Trapezoid', 'results.dat' u 1:8 w l tit 'Simpson'
