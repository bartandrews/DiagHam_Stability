set xrange [-1:28]
set yrange [5.72195147485666:16.7808092607586]
unset key
set xlabel "k_x*L_y + k_y" font "Helvetica,28"
set ylabel "{/Symbol x}(k)" font "Helvetica,28"
set size 1.0, 0.6
set terminal postscript portrait enhanced "Helvetica" 14
set output "fermions_hofstadter_X_69_Y_23_q_397_n_9_x_3_y_9_u_10_gx_0_gy_0.na_3.parentspec.eps"
plot "fermions_hofstadter_X_69_Y_23_q_397_n_9_x_3_y_9_u_10_gx_0_gy_0.na_3.parentspec" using 4:6 title "N = 9 N_A = 3 3 x 9" with points lc rgbcolor "red" lt 6 pt 1 ps 2
