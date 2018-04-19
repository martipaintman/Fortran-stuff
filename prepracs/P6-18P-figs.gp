
# apartat 1a)

set term png

set key outside
set key tmargin

set format x "%.e"
#set format y "%e"
#set key top right
set output 'P6-18P-fig1.png'

set title 'I1 en funció de N '
#set logscale y
#set logscale x

set xlabel 'N'
set ylabel 'Integral 1'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

plot 'P6-18Pres1.dat' using 1:2:3 with yerrorbars title 'I1 estimada', 122./27. -21.*pi**2/4


# apartat 1c)

set term png

set key outside
set key tmargin

set format x "%.e"
#set format y "%e"
#set key top right
set output 'P6-18P-fig2.png'

set title 'I2 en funció de N '
#set logscale y
#set logscale x

set xlabel 'N'
set ylabel 'Integral 2'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

plot 'P6-18Pres2.dat' using 1:2:3 with yerrorbars title 'I2 estimada'

# apartat 2)

set term png

set key outside
set key tmargin

set format x "%.e"
#set format y "%e"
#set key top right
set output 'P6-18P-fig3.png'

set title 'I3 en funció de N '
#set logscale y
#set logscale x

set xlabel 'N'
set ylabel 'Integral 3'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

plot 'P6-18Pres3.dat' using 1:2:3 with yerrorbars title 'I3 estimada'
