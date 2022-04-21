#[ gnuplot file: gnu.flw ]
set style data lines
set title 'mod_ccoax.G; convecive air flow [ time unit: 1 TU = 0.2 seconds ]'
set xrange [-7.736271556184712e-02:1.919865229524329e-01]
set yrange [-1.346746192571400e-01:1.346746192571400e-01]
set xlabel 'x/m'
set ylabel 'y/m'
set grid
set size square
set style lines 1 lt 7 lw 1.3
set style arrow 1 head nofilled size screen 0.009, 30, 90 ls 1
set arrow from -7.06290e-02,-1.225539e-01 to -4.56457e-02,-1.22554e-01 as 1
set label "1.00000e-01 m/s " at -4.31474e-02,-1.22554e-01
#set label "max: " at -7.06290e-02,-1.29961e-01
#set label "5.89450e-02 m/s " at -4.31474e-02,-1.29961e-01
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'10_TU' dt 1 lt 6 lw 1.5
pause 5.
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'11_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'12_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'13_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'14_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'15_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'16_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'17_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'18_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'19_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'20_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'21_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'22_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'23_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'24_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'25_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'26_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'27_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'28_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'29_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'30_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'31_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'32_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'33_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'34_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'35_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'36_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'37_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'38_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'39_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'40_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'41_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'42_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'43_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'44_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'45_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'46_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'47_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'48_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'49_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'50_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'51_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'52_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'53_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'54_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'55_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'56_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'57_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'58_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'59_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'60_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'61_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'62_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'63_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'64_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'65_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'66_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'67_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'68_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'69_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'70_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'71_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'72_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'73_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'74_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'75_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'76_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'77_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'78_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'79_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'80_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'81_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'82_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'83_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'84_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'85_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'86_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'87_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'88_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'89_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'90_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'91_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'92_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'93_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'94_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'95_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'96_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'97_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'98_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'99_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'100_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'101_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'102_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'103_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'104_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'105_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'106_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'107_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'108_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'109_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'110_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'111_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'112_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'113_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'114_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'115_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'116_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'117_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'118_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'119_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'120_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'121_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'122_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'123_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'124_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'125_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'126_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'127_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'128_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'129_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'130_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'131_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'132_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'133_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'134_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'135_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'136_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'137_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'138_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'139_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'140_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'141_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'142_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'143_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'144_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'145_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'146_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'147_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'148_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'149_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'150_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'151_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'152_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'153_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'154_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'155_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'156_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'157_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'158_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'159_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'160_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'161_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'162_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'163_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'164_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'165_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'166_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'167_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'168_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'169_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'170_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'171_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'172_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'173_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'174_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'175_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'176_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'177_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'178_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'179_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'180_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'181_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'182_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'183_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'184_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'185_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'186_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'187_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'188_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'189_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'190_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'191_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'192_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'193_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'194_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'195_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'196_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'197_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'198_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'199_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'200_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'201_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'202_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'203_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'204_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'205_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'206_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'207_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'208_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'209_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'210_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'211_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'212_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'213_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'214_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'215_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'216_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'217_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'218_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'219_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'220_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'221_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'222_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'223_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'224_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'225_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'226_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'227_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'228_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'229_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'230_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'231_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'232_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'233_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'234_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'235_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'236_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'237_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'238_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'239_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'240_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'241_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'242_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'243_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'244_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'245_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'246_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'247_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'248_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'249_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'250_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'251_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'252_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'253_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'254_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'255_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'256_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'257_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'258_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'259_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'260_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'261_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'262_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'263_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'264_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'265_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'266_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'267_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'268_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'269_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'270_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'271_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'272_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'273_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'274_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'275_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'276_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'277_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'278_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'279_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'280_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'281_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'282_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'283_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'284_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'285_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'286_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'287_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'288_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'289_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'290_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'291_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'292_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'293_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'294_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'295_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'296_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'297_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'298_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'299_TU' dt 1 lt 6 lw 1.5
pause 0.2000000000
plot \
'mesh' dt 1 lt 8 lw 0.1,\
'no_slip' dt 1 lt 7 lw 2.0,\
'free_slip' dt 1 lt 2 lw 2.0,\
'300_TU' dt 1 lt 6 lw 1.5
pause -1 '[ hit return to quit ]'
