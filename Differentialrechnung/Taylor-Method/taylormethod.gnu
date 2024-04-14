set term pdfcairo
set output 'Taylor-Method.pdf'
set multiplot layout 1,2 tit 'Taylor Method'
set grid 
set key l t
set xrange [0:1]
set auto y
set tit 'Taylor Method with Order 3'
f(x) = - 1 / ( x ** 2 /4 + x /2 - 1) 
set xlabel 'x'
set ylabel 'y(x)'
p f(x) tit 'Exact', 'results-h0.001.dat' u 1:2 w l dt 2 tit 'h = 0.001', 'results-h0.01.dat' u 1:2 w l dt 3 tit 'h = 0.01', 'results-h0.1.dat' u 1:2 w l dt 4 tit 'h = 0.1'
set grid
set key l t
set xrange [0:1]
set auto y
set tit 'Taylor Method with Order 1'
set xlabel 'x'
set ylabel 'y(x)'
p f(x) tit 'Exact', 'results-o1h0.001.dat' u 1:2 w l dt 2 tit 'h = 0.001', 'results-o1h0.01.dat' u 1:2 w l dt 3 tit 'h = 0.01', 'results-o1h0.1.dat' u 1:2 w l dt 4 tit 'h = 0.1'
unset multiplot
