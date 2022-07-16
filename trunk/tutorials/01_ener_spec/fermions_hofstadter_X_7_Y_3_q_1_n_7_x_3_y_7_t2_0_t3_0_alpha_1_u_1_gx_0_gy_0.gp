set xrange [-1:22]
set yrange [-12.9916195931897:-12.9882773983036]
set xlabel "k_x*L_y + k_y"
set ylabel "E(k)"
set size 1.0, 0.6
set terminal postscript portrait enhanced "Helvetica" 14
set output "fermions_hofstadter_X_7_Y_3_q_1_n_7_x_3_y_7_t2_0_t3_0_alpha_1_u_1_gx_0_gy_0.ps"
plot "fermions_hofstadter_X_7_Y_3_q_1_n_7_x_3_y_7_t2_0_t3_0_alpha_1_u_1_gx_0_gy_0.dat" using ($1*7+$2):3 title "C=, {/Symbol n}=:  N = 7  L_x=3, L_y=7, MUC: p=1, q=7 x 3" with points lc rgbcolor "red" lt 6 pt 1 ps 2
