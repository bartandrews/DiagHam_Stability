set xrange [-1:28]
set yrange [-12.7108796794313:-12.7107570570802]
set xlabel "k_x*L_y + k_y" font "Helvetica,28"
set lmargin at screen 0.15
set ylabel "E(k)/10^{-5}" font "Helvetica,28"
#set size 1.0, 0.6
set terminal postscript eps size 3.5,3 enhanced colour "Helvetica" 16
set output "fermions_hofstadter_X_69_Y_23_q_397_n_9_x_3_y_9_u_10_gx_0_gy_0_ext.eps"
plot "fermions_hofstadter_X_69_Y_23_q_397_n_9_x_3_y_9_u_10_gx_0_gy_0_ext.dat" using ($1*9+$2):3 title "L_x=3, L_y=9, MUC: p=397, q=69 x 23" with points lc rgbcolor "red" lt 6 pt 1 ps 2
