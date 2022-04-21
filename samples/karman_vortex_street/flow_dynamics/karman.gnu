#[ gnuplot file: gnu.flw ]
set style data lines
set title 'Karman vortex street [ time unit: 1 TU = 5.0E-4 seconds'
set xrange [-8.680000000000000e-02:2.468000000000000e-01]
set yrange [-1.668000000000000e-01:1.668000000000000e-01]
set xlabel 'x/m'
set ylabel 'y/m'
set grid
set size square
set style lines 1 lt 7 lw 1.3
set style arrow 1 head nofilled size screen 0.009, 30, 90 ls 1
set arrow from -7.84600e-02,-1.517880e-01 to -6.44239e-02,-1.51788e-01 as 1
set label "1.00000e+00 m/s " at -4.84360e-02,-1.51788e-01
#set label "max: " at -7.84600e-02,-1.60962e-01
#set label "1.48634e+00 m/s " at -4.84360e-02,-1.60962e-01
plot \
'0_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 10.
plot \
'1_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'2_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'3_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'4_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'5_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'6_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'7_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'8_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'9_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'10_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'11_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'12_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'13_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'14_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'15_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'16_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'17_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'18_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'19_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'20_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'21_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'22_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'23_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'24_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'25_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'26_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'27_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'28_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'29_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'30_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'31_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'32_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'33_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'34_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'35_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'36_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'37_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'38_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'39_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'40_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'41_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'42_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'43_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'44_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'45_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'46_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'47_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'48_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'49_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'50_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'51_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'52_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'53_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'54_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'55_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'56_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'57_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'58_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'59_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'60_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'61_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'62_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'63_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'64_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'65_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'66_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'67_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'68_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'69_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'70_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'71_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'72_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'73_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'74_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'75_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'76_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'77_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'78_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'79_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'80_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'81_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'82_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'83_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'84_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'85_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'86_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'87_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'88_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'89_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'90_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'91_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'92_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'93_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'94_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'95_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'96_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'97_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'98_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'99_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'100_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'101_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'102_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'103_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'104_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'105_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'106_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'107_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'108_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'109_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'110_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'111_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'112_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'113_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'114_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'115_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'116_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'117_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'118_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'119_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'120_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'121_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'122_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'123_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'124_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'125_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'126_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'127_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'128_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'129_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'130_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'131_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'132_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'133_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'134_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'135_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'136_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'137_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'138_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'139_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'140_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'141_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'142_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'143_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'144_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'145_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'146_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'147_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'148_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'149_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'150_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'151_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'152_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'153_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'154_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'155_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'156_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'157_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'158_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'159_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'160_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'161_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'162_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'163_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'164_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'165_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'166_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'167_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'168_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'169_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'170_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'171_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'172_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'173_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'174_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'175_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'176_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'177_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'178_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'179_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'180_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'181_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'182_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'183_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'184_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'185_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'186_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'187_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'188_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'189_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'190_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'191_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'192_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'193_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'194_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'195_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'196_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'197_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'198_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'199_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'200_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'201_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'202_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'203_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'204_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'205_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'206_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'207_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'208_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'209_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'210_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'211_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'212_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'213_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'214_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'215_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'216_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'217_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'218_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'219_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'220_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'221_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'222_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'223_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'224_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'225_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'226_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'227_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'228_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'229_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'230_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'231_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'232_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'233_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'234_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'235_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'236_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'237_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'238_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'239_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'240_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'241_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'242_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'243_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'244_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'245_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'246_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'247_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'248_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'249_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause 0.2000000000
plot \
'250_TU' with lines lt 6 lw 0.5,\
'no_slip' with lines lt 7 lw 2.0,\
'free_slip' with lines lt 2 lw 3.0
pause -1 '[ hit return to quit ]'
