set grid
set xlabel 't [s]'
set auto xy
set term pdf
set output './Plot-Data/plot-Euler-250.pdf'
set tit 'Position'
set ylabel 'x [m]'
set auto xy
p './Approx-Data/Euler-250.dat' u 1:2 w l tit 'Euler','./Approx-Data/Analytic-250.dat' u 1:2 w l tit 'Analytic','./Experimental-Data/experimental-250.dat' u 1:2 w l tit 'Experimental'
set tit 'Velocity'
set ylabel 'v [ms^{-1}]'
set auto xy
p './Approx-Data/Euler-250.dat' u 1:3 w l tit 'Euler','./Approx-Data/Analytic-250.dat' u 1:3 w l tit 'Analytic','./Experimental-Data/experimental-250.dat' u 1:3 w l tit 'Experimental'
set tit 'Energy'
set ylabel 'E [J]'
set auto xy
E = 0.00137304
p './Approx-Data/Euler-250.dat' u 1:4 w l tit 'Euler','./Approx-Data/Analytic-250.dat' u 1:4 w l tit 'Analytic',E tit 'Experimental'
set tit 'Position Error'
set ylabel 'Absolute Error'
set key left
set log y
set auto xy
p './Error-Data/Error-Euler-250.dat' u 1:2 w l tit 'Exp. vs An.', '' u 1:3 w l tit 'Exp. vs Nu.', '' u 1:4 w l tit 'An vs Nu.'
set tit 'Velocity Error'
set auto xy
p './Error-Data/Error-Euler-250.dat' u 1:5 w l tit 'Exp. vs An.', '' u 1:6 w l tit 'Exp. vs Nu.', '' u 1:7 w l tit 'An vs Nu.'
set tit 'Energy Error'
set auto xy
p './Error-Data/Error-Euler-250.dat' u 1:8 w l tit 'An. vs Nu.'
exit