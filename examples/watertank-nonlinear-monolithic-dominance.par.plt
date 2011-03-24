example = "watertank-nl-mono-pr&watertank-nl-mono-hy[Kp,ref]-par"
set xrange [0.2:0.8]
set yrange [5.25:8.25]
set size 1.5, 1
set xlabel 'K_p'
set ylabel 'ref'
set terminal postscript portrait enhanced colour dashed lw 1 "Helvetica" 24
set output example.".eps"
set style fill solid 1.0 border lt -1
plot example.".indeterminate.dump" with filledcurves title '' lt rgb 'yellow', \
example.".true.dump" with filledcurves title '' ls 1 lt rgb 'green' lw 1, \
example.".false.dump" with filledcurves title '' ls 1 lt rgb 'red' lw 1
