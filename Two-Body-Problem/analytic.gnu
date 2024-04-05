set term pdfcairo
set output 'analytic.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Analytic Solution'
p 'datanalytic.dat' u 3:4 w l tit 'Analytic solution'
