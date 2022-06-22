set terminal epslatex dashed  standalone  14 header "\\usepackage[T1]{fontenc}\n\\usepackage{mathptmx}\n\\usepackage{textcomp}\n\\usepackage{siunitx}"
set output "enerd.tex"


set lmargin  8.
set size 5/5, 6/5
set key left
#set key spacing 1.5
#set key at 3e10,5e6


set key spacing 1.5



set style line 2   lt 1 pt 7 ps 1.5 lc rgb  "red" lw 2
set style line 3 lt 1 pt 7 ps 1.0 lc rgb  "red" lw 3
set style line 4  lt 1 pt 7 ps 1.5 lc rgb  "blue" lw 2
set style line 5 lt 1 pt 7 ps 1.0 lc rgb  "black" lw 3
set style line 6   dt 4 pt 7 ps 1.0 lc rgb  "dark-green" lw 5
set style line 7 dt 5 pt 7 ps 1.0 lc rgb  "orange" lw 6
set style line 8 dt 3 pt 7 ps 1.0 lc rgb  "black" lw 8

# Make the bars semi-transparent so that the errorbars are easier to see.
#set style fill transparent solid 0.5 noborder
#set style data lines
set style fill  transparent solid 0.25 noborder


set xrange [0:35]
set yrange [0:500]

#set label 1 "\\large{$E_B=7.7$~MeV}" at 50,100

#set logscale y 10
#set logscale x 10
#set format y "$10^{\\scalebox{0.75}{%T}}$"
#set format y "\\large{$10^{\\scalebox{0.75}{%T}}$}"
#set format x "\\large{$10^{\\scalebox{0.75}{%T}}$}"
set ylabel "\\large{ $\\sigma$ [mb]}" offset 0,0 
set xlabel "\\large{ $E_B$ [MeV]}" 
#unset xlabel
#set xtics format ''
set mxtics 5
set mytics 5
set ytics 100 
set xtics 5

plot 'ener_17f_bu.dat' u 1:5:2:6 w xyerrorbars ls 2   t 'Inclusive', '' u 1:7:2:8 w xyerrorbars ls 4   t 'Exclusive'
