set terminal epslatex dashed  color  standalone  14 header "\\usepackage[T1]{fontenc}\n\\usepackage{mathptmx}\n\\usepackage{textcomp}\n\\usepackage{siunitx}"
set output "ang_80.tex"

set lmargin  8.
set size 5/5, 6/5
#set key at 40,49
set key spacing 1.5


set xrange [5:120]
set logscale y 10
#set format y "$10^{%L}$"
set yrange [0.1:2e1]
#set ylabel "ISGMR Strength [fm$^4$/MeV]" 
#set xlabel "$E_{x}$ [MeV]" 
#set label 1 "\\Large{${}^{24}$Mg}" at 30,30
#set label 2 "CDCC Calculation" at 14,0.035
#set label 1 "Grafico 3" at 0.1,1.8# rotate by 90 font "QuasiSwiss,22"
#set label 2 "^{152}Eu 244.7 keV" at 285,24000 center font "QuasiSwiss,26"
set mxtics 5
set mytics 5
#set nokey

#set xtics ("$10^o$" 10,"$20^o$" 20 , "$30^o$" 30, "$40^o$" 40, "$50^o$" 50, "$60^o$" 60)
#set ytics 300

#set arrow from 150, 1300 to 80,1250


a0 = 150
sigmaL = 3.2
x0 = 20.46

a1 = 15
sigmaL1 = 1.34
x1 = 10.5

a2 = 76
sigmaL2 = 3.5
x2 = 23.2

a3 = 24
sigmaL3 = 3.3
x3 = 26.0



y(x) = (a0*sigmaL/(2*3.1415))*(1/((x-x0)**2 + (0.5*sigmaL)**2))
y1(x) = (a1*sigmaL1/(2*3.1415))*(1/((x-x1)**2 + (0.5*sigmaL1)**2))
y2(x) = (a2*sigmaL2/(2*3.1415))*(1/((x-x2)**2 + (0.5*sigmaL2)**2))
y3(x) = (a3*sigmaL3/(2*3.1415))*(1/((x-x3)**2 + (0.5*sigmaL3)**2))

g(x)=y(x)+y1(x)+y2(x)+y3(x)

#fit [15:22] y(x) 'isgmr12C.dat' u 1:2 via a0, sigmaL, x0
#fit [8:30] g(x) 'isgmr12C.dat' u 1:2 via a0, sigmaL, x0, a2,x2, sigmaL2, a3, x3, sigmaL3, sigmaL1, x1, a1

binwidth=1
bin(x,width)=width*floor(x/width)

m1 = 4
z1 = 2
m2 = 17
z2 = 9
K1 = 3.95

ruth(x) = 10.*(z1*z2*197./137.)**2/( (4*K1*m2/(m1+m2))**2 *(sin(0.5*x*3.141592/180))**4 )

set boxwidth 0.5 relative
set bars small #no end error bars

set  style line 1 pt 7 lt 1 lc 7 ps 2.0
set style line 2   dt 1 pt 7 ps 1.5 lc rgb  "black" lw 3
set style line 3 dt 1 pt 7 ps 2.0 lc rgb  "gray" lw 2
set style line 4  dt 5 pt 7 ps 2.0 lc rgb  "web-green" lw 2
set  style line 5 dt 2 pt 7 ps 2.0 lc rgb  "slateblue1" lw 2
set  style line 6 dt 1 pt 7 ps 2.0 lc rgb  "red" lw 2


#plot 'alpha_17F_angdis_90.dat' u ($1):(0.4*$2) t 'in', 'alpha_17F_angdis_90.dat' u ($1):(0.24*$4) t 'out'#, 'spp_90.dat' u 1:($2) w l t 'SPP1'
plot 'alpha_17F_angdis_180.dat' u ($1):(0.44*$2/ruth($1)) t 'in', 'alpha_17F_angdis_180.dat' u ($1):(0.46*$4/ruth($1)) t 'out', 'spp_180.dat' u 1:($3) w l t 'SPP1'

#plot  'spp_120.dat' u 1:3 w l t 'SPP', 'spp_120_16O.dat' u 1:3 w l t 'SPP2', 'spp_120_sum.dat' u 1:3 w l t 'suma' 


#plot  'mono_24mg_texas.csv' u ($1):($3*0.69) w boxes ls 3 t '${}^{24}\text{Mg}(\alpha,\alpha^\prime)$ (scaled)', 'isgmr24mg.dat' u 1:2:6 w yerrorbars ls 2 t 'Exp. data' , 'qrpa_folded.dat' u ($1+2.5):($2*0.49)  smooth csplines ls 6  t 'RPA' #'amd_folded.dat' u ($1):($2*5.5)  smooth csplines ls 4  t 'AMD'#, g(x) w l ls 5 t 'fit', 'C12_gr_calctot_folded.dat' u 1:($2/3)  w l ls 4  t 'AMD', 'qrpa_c12_sigma3.3.dat' u ($1-1.0):($2*0.84)  smooth csplines ls 6  t 'QRPA'

#plot  'mono_24mg_texas.csv' u ($1):($3*0.69) w boxes ls 3 t '${}^{24}\text{Mg}(\alpha,\alpha^\prime)$ (scaled)', 'isgmr24mg.dat' u 1:2:5 w yerrorbars ls 2 t 'Exp. data' , 'amd_folded.dat' u ($1):($2*5.5)  smooth csplines ls 4  t 'AMD'#, g(x) w l ls 5 t 'fit', 'C12_gr_calctot_folded.dat' u 1:($2/3)  w l ls 4  t 'AMD', 'qrpa_c12_sigma3.3.dat' u ($1-1.0):($2*0.84)  smooth csplines ls 6  t 'QRPA'


replot
