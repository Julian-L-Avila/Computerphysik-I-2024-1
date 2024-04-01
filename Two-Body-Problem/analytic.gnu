set term pdfcairo
set output 'analytic.pdf'
set grid
set parametric
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set tit 'Analytic Solution'
a = 5
e = 0.5
r(t) = a / (1 + e*cos(t))
p (r(t)*cos(t)),(r(t)*sin(t)) w l tit 'Analytic solution'
