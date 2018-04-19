# apartat 1) histo

set term png

set key outside
set key tmargin
#set format y "%e"
#set key top right
set output 'P5-18P-fig1.png'

set title 'Apartat 1) histograma 120 caixes'
#set logscale y
#set logscale x

set xlabel 'x'
set ylabel 'P(x)'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

set style histogram rowstacked gap 0
set style fill solid 0.1

normal(x,mu,sd)=(1/(sd*sqrt(2*pi)))*exp(-(x-mu)**2/(2*sd**2))

plot 'P5-18P-res1.dat' with boxerrorbars ls 1 title 'Frequència dels numeros aleatoris',normal(x,0,1) t'Distribucio gaussiana estandaritzada'


# apartat 2) recorregut

set term png

set key outside
set key tmargin
#set format y "%e"
#set key top right
set output 'P5-18P-fig2.png'

set title 'Apartat 2a) recorregut 5 molecules '
#set logscale y
#set logscale x

set xlabel 'x(m)'
set ylabel 'y(m)'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

plot 'P5-18P-res3.dat' using 1:2 with l title 'Molec 1' ,\
'P5-18P-res3.dat' using 3:4 with l title 'Molec 2 ' ,\
'P5-18P-res3.dat' using 5:6 with l title 'Molec 3' ,\
'P5-18P-res3.dat' using 7:8 with l title 'Molec 4 ' ,\
'P5-18P-res3.dat' using 9:10 with l title 'Molec 5 ' ,\

# apartat 2c) recorregut

set term png

set key outside
set key tmargin
#set format y "%e"
#set key top right
set output 'P5-18P-fig3.png'

set title 'Apartat 2c) recorregut 5 molecules '
#set logscale y
#set logscale x

set xlabel 't'
set ylabel 'var(y(t))'
set zeroaxis

#set yrange [13:26]
#set xrange [0:3.5]

plot 'P5-18P-res4.dat' using 1:2 with l title 'Molec 1' ,\
