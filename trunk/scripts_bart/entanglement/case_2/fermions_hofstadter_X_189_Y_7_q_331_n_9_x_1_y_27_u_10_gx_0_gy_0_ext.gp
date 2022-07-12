set xrange [-1:28]
set yrange [-12.7067168175024:-12.7066126997043]
set xlabel "k_x*L_y + k_y" font "Helvetica,28"
set ylabel "E(k)" font "Helvetica,28"
#set size 1.0, 0.6
set terminal postscript eps size 3.5,3 enhanced colour "Helvetica" 16
set output "fermions_hofstadter_X_189_Y_7_q_331_n_9_x_1_y_27_u_10_gx_0_gy_0_ext.eps"
plot "fermions_hofstadter_X_189_Y_7_q_331_n_9_x_1_y_27_u_10_gx_0_gy_0_ext.dat" using ($1*27+$2):3 title "L_x=1, L_y=27, MUC: p=331, q=189 x 7" with points lc rgbcolor "red" lt 6 pt 1 ps 2
