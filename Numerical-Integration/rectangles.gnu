set term pdfcairo
set output 'Rectangles.pdf'
set grid
set xlabel 'x'
set ylabel 'f(x)'
f(x) = 100 * (x ** 2) * cos(20* x)
p f(x) tit 'Function', 'dataRectangular.dat' u 1:2 w boxes tit 'Rectangles'
