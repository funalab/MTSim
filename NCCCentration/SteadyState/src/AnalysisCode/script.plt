set tics font "Lucida-Grande,18"
set label font "Lucida-Grande,18"
set title font "Lucida-Grande,34"
unset key
set size square

set terminal pdf enhanced
plot '../../Result/Simulation/centering_vector_strain8_rad25_rads15.dat' using 1:2:(1.5*$3):(20*$4) with vector lw 1.5 lc rgb 'red'
set output "../../Result/Analysis/rads15.pdf"
replot
set output

set terminal pdf enhanced
plot '../../Result/Simulation/centering_vector_strain8_rad25_rads15.dat' using 1:2:(1.5*$3):(20*$4) with vector lw 1.5 lc rgb 'red'
replot '../../Result/Simulation/centering_vector_strain8_rad25_rads10.dat' using 1:2:(1.5*$3):(20*$4) with vector lw 1.5 lc rgb 'green'
set output "../../Result/Analysis/rads15_rads10.pdf"
replot
set output

set terminal pdf enhanced
plot '../../Result/Simulation/centering_vector_strain8_rad25_rads25.dat' using 1:2:(1.2*$3):(20*$4) with vector lw 1.5 lc rgb 'red'
set output "../../Result/Analysis/rads25.pdf"
replot
set output
