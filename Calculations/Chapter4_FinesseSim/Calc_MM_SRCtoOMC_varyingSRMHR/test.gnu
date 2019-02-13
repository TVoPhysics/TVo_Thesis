reset
set xrange[5:5000]
set logscale x
set xlabel "f [Hz] (darm)"
set ylabel "Re "
set y2label "Im "
set y2tics nomirror
set my2tics 3
set ytics nomirror
set mxtics 2
set mytics 2
set zero 0.0
set title "test                Wed Sep  7 13:42:53 2016" offset 0, 2
set nolog y
set format y "%g"
set nolog y2
set format y2 "%g"
set term x11 
set size ratio .5 
set key below 
set grid xtics ytics 
set colors classic
set termoption noenhanced
plot\
'test.out' using ($1):($2) axes x1y1 title "NSR_with_RP nOMC_AROC_trans : Re  " w l lt 1 lw 2, \
'test.out' using ($1):($3) axes x1y2 title "NSR_with_RP nOMC_AROC_trans : Im  " w l lt 2 lw 2
