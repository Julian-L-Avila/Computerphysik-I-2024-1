set term pdfcairo
set output 'rectangular.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Rectangular Rule'
p 'datrectangular.dat' u 4:5 w l tit 'Rectangular Approx.'
