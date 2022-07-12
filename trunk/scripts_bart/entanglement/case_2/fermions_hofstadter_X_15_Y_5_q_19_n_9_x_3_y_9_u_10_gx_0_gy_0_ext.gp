set xrange [-1:28]
set yrange [-12.3686457285112:-12.3513835977337]
set xlabel "k_x*L_y + k_y" font "Helvetica,28"
set ylabel "E(k)" font "Helvetica,28"
#set size 1.0, 0.6
set terminal postscript eps size 3.5,3 enhanced colour "Helvetica" 16
set output "fermions_hofstadter_X_15_Y_5_q_19_n_9_x_3_y_9_u_10_gx_0_gy_0_ext.eps"
plot "fermions_hofstadter_X_15_Y_5_q_19_n_9_x_3_y_9_u_10_gx_0_gy_0_ext.dat" using ($1*9+$2):3 title "L_x=3, L_y=9, MUC: p=19, q=15 x 5" with points lc rgbcolor "red" lt 6 pt 1 ps 2
