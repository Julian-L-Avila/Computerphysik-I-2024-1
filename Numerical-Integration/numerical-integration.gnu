set term pdfcairo
set output 'Approximate-Integral.pdf'
set grid
set xlabel 'N'
set ylabel 'Numerical Approx'
set logscale x
p 'results.dat' u 1:2 w l tit 'Exact Integral', 'results.dat' u 1:3 w l tit 'Rectangular Approx.', 'results.dat' u 1:4 w l tit 'Trapezoid Approx.', 'results.dat' u 1:5 w l tit 'Simpson Approx.'
