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
set title "VarySR3_95MM_wSqz_MMTTon_maxtem4                Tue Sep 13 22:09:04 2016" offset 0, 2
set nolog y
set nolog y2
set term x11 
set size ratio .5 
set key below 
set grid xtics ytics 
set colors classic
set termoption noenhanced
plot\
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($2) axes x1y1 title "NSR_with_RP nOMC_AROC_trans : Re  " w l lt 1 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($3) axes x1y2 title "NSR_with_RP nOMC_AROC_trans : Im  " w l lt 2 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($4) axes x1y1 title "nSRMHRaTEM00 nSRMHRa : Re  " w l lt 3 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($5) axes x1y2 title "nSRMHRaTEM00 nSRMHRa : Im  " w l lt 4 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($6) axes x1y1 title "nSRMHRaTEM01 nSRMHRa : Re  " w l lt 5 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($7) axes x1y2 title "nSRMHRaTEM01 nSRMHRa : Im  " w l lt 6 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($8) axes x1y1 title "nSRMHRaTEM02 nSRMHRa : Re  " w l lt 7 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($9) axes x1y2 title "nSRMHRaTEM02 nSRMHRa : Im  " w l lt 8 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($10) axes x1y1 title "SRCoutx nIBAin : Re  " w l lt 9 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($11) axes x1y2 title "SRCoutx nIBAin : Im  " w l lt 0 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($12) axes x1y1 title "SRCouty nIBAin : Re  " w l lt 11 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($13) axes x1y2 title "SRCouty nIBAin : Im  " w l lt 12 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($14) axes x1y1 title "SRMYqx nSRMHRa : Re  " w l lt 13 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($15) axes x1y2 title "SRMYqx nSRMHRa : Im  " w l lt 14 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($16) axes x1y1 title "SRMYqy nSRMHRa : Re  " w l lt 15 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($17) axes x1y2 title "SRMYqy nSRMHRa : Im  " w l lt 16 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($18) axes x1y1 title "ITMXqx nITMX2 : Re  " w l lt 17 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($19) axes x1y2 title "ITMXqx nITMX2 : Im  " w l lt 18 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($20) axes x1y1 title "ITMXqy nITMX2 : Re  " w l lt 19 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($21) axes x1y2 title "ITMXqy nITMX2 : Im  " w l lt 20 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($22) axes x1y1 title "ITMYqx nITMY2 : Re  " w l lt 21 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($23) axes x1y2 title "ITMYqx nITMY2 : Im  " w l lt 22 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($24) axes x1y1 title "ITMYqy nITMY2 : Re  " w l lt 23 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($25) axes x1y2 title "ITMYqy nITMY2 : Im  " w l lt 24 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($26) axes x1y1 title "OMCqx nOMC_HROC_refl : Re  " w l lt 25 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($27) axes x1y2 title "OMCqx nOMC_HROC_refl : Im  " w l lt 26 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($28) axes x1y1 title "OMCqy nOMC_HROC_refl : Re  " w l lt 27 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($29) axes x1y2 title "OMCqy nOMC_HROC_refl : Im  " w l lt 28 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($30) axes x1y1 title "OFIqx nIBAin : Re  " w l lt 29 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($31) axes x1y2 title "OFIqx nIBAin : Im  " w l lt 30 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($32) axes x1y1 title "OFIqy nIBAin : Re  " w l lt 31 lw 2, \
'VarySR3_95MM_wSqz_MMTTon_maxtem4.out' using ($1):($33) axes x1y2 title "OFIqy nIBAin : Im  " w l lt 32 lw 2
