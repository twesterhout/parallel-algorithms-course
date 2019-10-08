#!/usr/bin/env gnuplot

set style line 11 lt 1 lc rgb '#0072bd' # blue
set style line 12 lt 1 lc rgb '#d95319' # orange
set style line 13 lt 1 lc rgb '#edb120' # yellow
set style line 14 lt 1 lc rgb '#7e2f8e' # purple
set style line 15 lt 1 lc rgb '#77ac30' # green
set style line 16 lt 1 lc rgb '#4dbeee' # light-blue
set style line 17 lt 1 lc rgb '#a2142f' # red

set style line 101 lc rgb '#808080' lt 1
set border 3 back ls 101
set tics nomirror out scale 0.75

set lmargin 8
set xlabel "upper bound" offset 0,0.5
set xtics 0, 4e6
# set xtics ("0" 0, "1⋅10⁴" 10000, "2⋅10⁴" 20000, "3⋅10⁴" 30000, "4⋅10⁴" 40000, "5⋅10⁴" 50000)
set ylabel "time, seconds" offset 0,0
# set ytics ("0" 0, "2⋅10⁴" 20000, "4⋅10⁴" 40000, "6⋅10⁴" 60000, "8⋅10⁴" 80000, "10⋅10⁴" 100000)
set key top left

set terminal pngcairo size 480, 320 noenhanced font "Iosevka Light, 12"

set output "timings.png"
plot "timings_serial.dat"   u 1:2 w lp ls 15 pt 7 lw 2 title "serial", \
     "timings_parallel.dat" u 1:2 w lp ls 11 pt 7 lw 2 title "#proc = 1", \
     ""                     u 1:3 w lp ls 12 pt 7 lw 2 title "#proc = 2", \
     ""                     u 1:4 w lp ls 13 pt 7 lw 2 title "#proc = 3", \
     ""                     u 1:5 w lp ls 14 pt 7 lw 2 title "#proc = 4"
set output


set xtics 0, 2e8

set output "timings_large.png"
plot "timings_large_serial.dat"   u 1:2 w lp ls 15 pt 7 lw 2 title "serial", \
     "timings_large_parallel.dat" u 1:2 w lp ls 11 pt 7 lw 2 title "#proc = 1", \
     ""                           u 1:3 w lp ls 12 pt 7 lw 2 title "#proc = 2", \
     ""                           u 1:4 w lp ls 13 pt 7 lw 2 title "#proc = 3", \
     ""                           u 1:5 w lp ls 14 pt 7 lw 2 title "#proc = 4"
set output
