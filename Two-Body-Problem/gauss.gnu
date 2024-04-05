set term pdfcairo
set output 'gauss.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Gaussian Quadrature'
p 'datgauss.dat' u 4:5 w l tit 'Gaussian Approx.'
