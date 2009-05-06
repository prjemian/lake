set terminal postscript landscape color
set output "plot.ps"
set logscale x 10
set logscale y 10
set data style linespoints
set xlabel "Q, 1/A" 0,0
set ylabel "intensity, a.u." 0,0
set autoscale xy
set zero 1e-08
plot   "test.smr" title "smeared data" with errorbars, \
   "test1.dsm" title "desmeared data (1x)" with errorbars, \
   "test10.dsm" title "desmeared data (10x)" with errorbars
