set xrange [0:46]
set yrange [0:46]
set xlabel "x"
set ylabel "y"
set zlabel "g(r)"
set key font ",10"
set terminal wxt
splot "fermions_hofstadter_X_15_Y_5_q_76_n_9_x_3_y_9_u_10_gx_0_gy_0_kx_0_ky_0.0.rhorho" title "|C|=1, {/Symbol n}=1/3:  N = 9  L_x=3, L_y=9, MUC: p=76, q=15 x 5, k=(0,0)" with points lt 6 pt 1 ps 2