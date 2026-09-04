# =============================================================
# Figure 6: Alpha dependence of relaxation and characteristic slip.
#
# Panel (a) shows R(u) on linear axes, with the same curves on logarithmic
# axes in the inset. Panel (b) shows D_c/(V2 tau) as a function of alpha.
# relaxation_alpha_scan.c generates the colored curves and points, while
# characteristic_slip.c generates the finely sampled gray reference curve.
#
# Run from this directory: mkdir -p fig && gnuplot load.plt
# Output: fig/relaxation_alpha.pdf
# =============================================================

reset

system "mkdir -p fig"

# ----- 出力設定 -----
set encoding utf8
set terminal pdfcairo enhanced color font "Times New Roman,22" size 18cm,24cm \
    background rgb "white" linewidth 1.5 rounded
set output "fig/relaxation_alpha.pdf"

# ----- データの場所 -----
dir = "data/"

# ----- モデルパラメータ -----
V2   = 1.0
TAU  = 1.0
XMIN = 0.01      # xi_min = Lmin/(V2 tau)
XMAX = 1000.0    # xi_max = Lmax/(V2 tau)

# ----- alpha のリスト -----
alphas = "1.20 1.40 1.60 1.80 2.00 2.20 2.40 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00"
na     = words(alphas)
A_MIN  = 1.2
A_MAX  = 4.0

fname(s) = dir.sprintf("R_u_alpha%s.txt", s)

# ----- 線形パネルの表示範囲 -----
U_LIN = 300.0    # (a) 主図の横軸上限

# =============================================================
#  カラーマップ (viridis)
# =============================================================
set palette defined (\
    0.0000 "#440154", \
    0.1429 "#46327e", \
    0.2857 "#365c8d", \
    0.4286 "#277f8e", \
    0.5714 "#1fa187", \
    0.7143 "#4ac16d", \
    0.8571 "#a0da39", \
    1.0000 "#fde725" )
set cbrange [A_MIN:A_MAX]
set cbtics 1.5,0.5,4.0 font "Times New Roman,24" offset -0.5,0

# ----- 共通の体裁 -----
set border lw 1.5
set tics in scale 1.2,0.6
set tics front
set tics font "Times New Roman,32"
unset key

# =============================================================
#  multiplot 開始 (パネル位置は screen 座標で直接指定)
# =============================================================
set multiplot

# =============================================================
#  (a) 主図 : 線形軸
# =============================================================
set lmargin at screen 0.155
set rmargin at screen 0.845
set bmargin at screen 0.575
set tmargin at screen 0.965

unset logscale
set xrange [0:U_LIN]
set yrange [0:1.05]
set format x "%.0f"
set format y "%.1f"
set xtics 0,100,300
set ytics 0,0.5,1.0
set mxtics 5
set mytics 5
set xlabel "u = t/{/Symbol t}" font "Times New Roman,32" offset 0,0.3
set ylabel "R(u)" font "Times New Roman,32" offset -1.5,0

set label 10 "(a)" at graph 0.04, graph 0.92 font "Times New Roman,32" front
set arrow 10 from first 150, first 0.27 to first 30, first 0.045 \
    head size screen 0.012,15 lw 2 lc rgb "gray30" front
set label 11 "increasing {/Symbol a}" at first 95, first 0.185 center \
    rotate by -22 tc rgb "gray30" font "Times New Roman,20" front

# カラーバー ((a) の右側)
set colorbox vertical user origin 0.875,0.575 size 0.025,0.39
set label 12 "{/Symbol a}" at screen 0.8875, screen 0.99 center font "Times New Roman,26"

plot \
    for [i=1:na] fname(word(alphas,i)) \
        u 1:2 w l lw 4 lc palette cb (word(alphas,i)+0) notitle

unset label 10
unset label 11
unset label 12
unset arrow 10
unset colorbox

# =============================================================
#  (a) 挿入図 : 両対数軸 (同じデータ)
# =============================================================
# 挿入図の下地 (主図の曲線を隠して目盛を読みやすくする)
set object 1 rectangle from screen 0.475,0.740 to screen 0.838,0.958 \
    fc rgb "white" fs solid 1.0 noborder behind

set lmargin at screen 0.545
set rmargin at screen 0.830
set bmargin at screen 0.775
set tmargin at screen 0.945

set logscale xy
set xrange [1e-2:2e3]
set yrange [1e-7:2]
set format x "10^{%T}"
set format y "10^{%T}"
set xtics 1e-2,100,1e2
set ytics 1e-7,100,1e0
set mxtics 10
set mytics 10
set tics font "Times New Roman,25"
unset xlabel
unset ylabel
set label 13 "u" at graph 0.5, graph -0.20 center font "Times New Roman,18"
set label 14 "R(u)" at graph -0.32, graph 0.5 center rotate by 90 font "Times New Roman,18"
set border lw 1.0

plot \
    for [i=1:na] fname(word(alphas,i)) \
        u 1:($2>0 ? $2 : 1/0) w l lw 2 lc palette cb (word(alphas,i)+0) notitle

unset label 13
unset label 14
unset object 1
set border lw 1.5
set tics font "Times New Roman,32"

# =============================================================
#  (b) alpha vs xi_rel = D_rel/(V2 tau)
# =============================================================
set lmargin at screen 0.155
set rmargin at screen 0.845
set bmargin at screen 0.075
set tmargin at screen 0.465

unset logscale x
set logscale y
set xrange [1.13:4.07]
set yrange [6e-3:2e3]
set format x "%.1f"
set format y "10^{%T}"
set xtics 1.5,0.5,4.0
set mxtics 5
set ytics 1e-2,10,1e3
set mytics 10
set xlabel "{/Symbol a}" font "Times New Roman,32" offset 0,0.3
set ylabel "{/Symbol x}_{rel} = D_{rel}/(V_2{/Symbol t})" font "Times New Roman,26" offset -1.5,0

set label 20 "(b)" at graph 0.04, graph 0.08 font "Times New Roman,26" front

# --- xi_min, xi_max : 点線 ---
set arrow 21 from graph 0, first XMIN to graph 1, first XMIN nohead lw 2 dt 3 lc rgb "gray50"
set arrow 22 from graph 0, first XMAX to graph 1, first XMAX nohead lw 2 dt 3 lc rgb "gray50"
set label 21 "{/Symbol x}_{min}" at graph 0.56, first XMIN*2.0 tc rgb "gray40" font "Times New Roman,20"
set label 22 "{/Symbol x}_{max}" at graph 0.80, first XMAX*0.42 tc rgb "gray40" font "Times New Roman,20"

# --- alpha = 2, 3 : スケーリング領域の境界 ---
set arrow 23 from first 2.0, graph 0 to first 2.0, graph 1 nohead lw 1.5 dt 2 lc rgb "gray80"
set arrow 24 from first 3.0, graph 0 to first 3.0, graph 1 nohead lw 1.5 dt 2 lc rgb "gray80"

# --- 漸近形のラベル ---
set label 23 "{/Symbol \265} {/Symbol x}_{max}"                at first 1.25, first 3.5e1  tc rgb "gray40" font "Times New Roman,18"
set label 24 "{/Symbol \265} {/Symbol x}_{max}^{3-{/Symbol a}}" at first 2.58, first 1.6e1 tc rgb "gray40" font "Times New Roman,18"
set label 25 "const."                                          at first 3.25, first 2.0e-1 tc rgb "gray40" font "Times New Roman,20"

plot \
    dir."Dc_ximax1e+03.txt" u 1:3   w l lw 4 lc rgb "gray30" notitle, \
    dir."Dc_vs_alpha.txt"   u 1:4   w p pt 7 ps 1.6 lc rgb "white" notitle, \
    dir."Dc_vs_alpha.txt"   u 1:4:1 w p pt 6 ps 1.6 lw 3 lc palette notitle

unset multiplot
unset output
