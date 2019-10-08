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

set lmargin 11
set xlabel "upper bound" offset 0,0.5
set xtics ("0" 0, "1⋅10⁴" 10000, "2⋅10⁴" 20000, "3⋅10⁴" 30000, "4⋅10⁴" 40000, "5⋅10⁴" 50000)
set ylabel "number operations" offset 3.5,0
set ytics ("0" 0, "2⋅10⁴" 20000, "4⋅10⁴" 40000, "6⋅10⁴" 60000, "8⋅10⁴" 80000, "10⋅10⁴" 100000)
set key top left

set terminal pngcairo size 480, 320 noenhanced font "Iosevka Light, 12"
set output "counts.png"

F(x) = x * (log(log(x)) + 0.261497212 - log(2.0) + 2.0 / (sqrt(x) * log(x)))
plot[][0:100000] "counts.txt" u 1:2 w l ls 11 lw 2 title "exact", "" u 1:(F($1)) w l ls 12 lw 2 title "approximation"

set output
