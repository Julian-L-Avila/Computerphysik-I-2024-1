set grid
set xlabel 't [s]'
set auto xy
set term pdf
set output './Plot-Data/plot-Euler-270.pdf'
set tit 'Position'
set ylabel 'x [m]'
set auto xy
p './Approx-Data/Euler-270.dat' u 1:2 w l tit 'Euler','./Approx-Data/Analytic-270.dat' u 1:2 w l tit 'Analytic','./Experimental-Data/experimental-270.dat' u 1:2 w l tit 'Experimental'
set tit 'Velocity'
set ylabel 'v [ms^{-1}]'
set auto xy
p './Approx-Data/Euler-270.dat' u 1:3 w l tit 'Euler','./Approx-Data/Analytic-270.dat' u 1:3 w l tit 'Analytic','./Experimental-Data/experimental-270.dat' u 1:3 w l tit 'Experimental'
set tit 'Energy'
set ylabel 'E [J]'
set auto xy
E = 0.0310424
p './Approx-Data/Euler-270.dat' u 1:4 w l tit 'Euler','./Approx-Data/Analytic-270.dat' u 1:4 w l tit 'Analytic',E tit 'Experimental'
set tit 'Position Error'
set ylabel 'Absolute Error'
set key left
set log y
set auto xy
p './Error-Data/Error-Euler-270.dat' u 1:2 w l tit 'Exp. vs An.', '' u 1:3 w l tit 'Exp. vs Nu.', '' u 1:4 w l tit 'An vs Nu.'
exit