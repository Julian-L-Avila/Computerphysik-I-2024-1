set term pdfcairo
set output 'trapezoid.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Trapezoidal Rule'
p 'dattrapezoid.dat' u 4:5 w l tit 'Trapezoid Approx.'
