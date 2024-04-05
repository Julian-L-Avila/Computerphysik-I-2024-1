set term pdfcairo
set output 'simpson.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Simpson Rule'
p 'datsimpson.dat' u 4:5 w l tit 'Simpson Approx.'
