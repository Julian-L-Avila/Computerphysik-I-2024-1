set term pdfcairo
set output 'analytic.pdf'
set grid
set xlabel 'x (AU)'
set ylabel 'y (AU)'
set parametric
a = 5
e = 0.5
r(t) = a * (1 - e**2) / (1 + e*cos(t))
set tit 'Analytic Solution'
p r(t)*cos(t), r(t)*sin(t) tit 'Analytic solution'
