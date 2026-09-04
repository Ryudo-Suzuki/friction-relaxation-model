# =============================================================
# Figure 4: Dependence of the relaxation function on the upper cutoff.
#
# The three panels show alpha = 1.5, 2.0, and 3.0 for several values of
# xi_max. The input files are generated with relaxation_curve.c.
# Data columns: (1) u = t/tau, (2) R(u).
#
# Run from this directory: mkdir -p fig && gnuplot load.plt
# Output: fig/relaxation_all.pdf
# =============================================================

reset

system "mkdir -p fig"

# ----- 出力設定 -----
set terminal pdfcairo enhanced font "Times-New-Roman,25" size 46cm,13cm
set output "fig/relaxation_all.pdf"
set encoding utf8

# =============================================================
#  全パネル共通の設定
# =============================================================
set logscale xy
set xrange [1e-4:1e5]
set yrange [1e-5:2]          # 3図で y 軸を統一

set format x "10^{%T}"
set format y "10^{%T}"

set xtics 1e-4,10,1e5
set ytics 1e-5,10,1
set mxtics 10
set mytics 10

set tics font "Times New Roman,28"
set xlabel "{/=30 u = t/{/Symbol t}}"
set ylabel "{/=30 R(u)}"

set key left bottom font "Times-New-Roman,18"

set lmargin 10
set bmargin 4
set tmargin 2
set rmargin 3

# ----- 線色 (全図共通 : Lmax が大きいほど明るい緑側へ, viridis 系) -----
cL0 = "#482878"   # Lmax = 10^0
cL1 = "#2a788e"   # Lmax = 10^1
cL2 = "#22a884"   # Lmax = 10^2
cL3 = "#7ad151"   # Lmax = 10^3
cL4 = "#bddf26"   # Lmax = 10^4

# ----- ガイドライン共通スタイル -----
set style line 91 lc rgb "black" lw 3 dt 2

# alpha=2 用ガイドライン関数 : R_g(u) = 0.80 * log(Lref/u)/log(Lref)
Lref = 1e4
Rth(x, ximax) = 1.0 - log(x) / log(ximax)

# =============================================================
#  multiplot 開始 (1 行 3 列)
# =============================================================
set multiplot layout 1,3 rowsfirst

# =============================================================
#  (左) alpha = 1.5 :  R(u) ~ u^0  (const)
# =============================================================
set title "(a) {/Symbol a} = 1.5" font "Times-New-Roman,30"

unset label
set label 1 "{/=40 {/Symbol \265} u^{0} }" at graph 0.78,0.9 front

plot \
    "data/R_u_alpha1.500_Lmax1e+00_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL0 title "{/Symbol x}_{max}=10^{0}", \
    "data/R_u_alpha1.500_Lmax1e+01_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL1 title "{/Symbol x}_{max}=10^{1}", \
    "data/R_u_alpha1.500_Lmax1e+02_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL2 title "{/Symbol x}_{max}=10^{2}", \
    "data/R_u_alpha1.500_Lmax1e+03_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL3 title "{/Symbol x}_{max}=10^{3}", \
    "data/R_u_alpha1.500_Lmax1e+04_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL4 title "{/Symbol x}_{max}=10^{4}", \
    (x>=1e+0 && x<=1e3 ? 0.6 : 1/0) w l ls 91 notitle

# =============================================================
#  (中) alpha = 2.0 :  R(u) ~ log(xi_max/u)/log(xi_max)  (対数)
# =============================================================
set title "(b) {/Symbol a} = 2.0" font "Times-New-Roman,30"

unset label
#set label 1 "{/=28 {/Symbol \265} log\\({/Symbol x}_{max}/u\\)/log{/Symbol x}_{max}}" \
    at graph 0.55,0.93 tc rgb "black" front

plot \
    "data/R_u_alpha2.000_Lmax1e+00_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL0 title "{/Symbol x}_{max}=10^{0}", \
    "data/R_u_alpha2.000_Lmax1e+01_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL1 title "{/Symbol x}_{max}=10^{1}", \
    "data/R_u_alpha2.000_Lmax1e+02_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL2 title "{/Symbol x}_{max}=10^{2}", \
    "data/R_u_alpha2.000_Lmax1e+03_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL3 title "{/Symbol x}_{max}=10^{3}", \
    "data/R_u_alpha2.000_Lmax1e+04_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL4 title "{/Symbol x}_{max}=10^{4}", \
#    (x>=1e+0 && x<=1e1 ? Rth(x,1e1) : 1/0) w l lw 2 dt 3 lc rgb cL1 notitle, \
#    (x>=1e+0 && x<=1e2 ? Rth(x,1e2) : 1/0) w l lw 2 dt 3 lc rgb cL2 notitle, \
#    (x>=1e+0 && x<=1e3 ? Rth(x,1e3) : 1/0) w l lw 2 dt 3 lc rgb cL3 notitle, \
#    (x>=1e+0 && x<=1e4 ? Rth(x,1e4) : 1/0) w l lw 2 dt 3 lc rgb cL4 notitle

# =============================================================
#  (右) alpha = 3.0 :  R(u) ~ u^{-1}  (べき, gamma = alpha-2 = 1)
# =============================================================
set title "(c) {/Symbol a} = 3.0" font "Times-New-Roman,30"

unset label
set label 1 "{/=40 {/Symbol \265} u^{-1}}" at graph 0.6,0.65 front

plot \
    "data/R_u_alpha3.000_Lmax1e+00_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL0 title "{/Symbol x}_{max}=10^{0}", \
    "data/R_u_alpha3.000_Lmax1e+01_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL1 title "{/Symbol x}_{max}=10^{1}", \
    "data/R_u_alpha3.000_Lmax1e+02_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL2 title "{/Symbol x}_{max}=10^{2}", \
    "data/R_u_alpha3.000_Lmax1e+03_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL3 title "{/Symbol x}_{max}=10^{3}", \
    "data/R_u_alpha3.000_Lmax1e+04_tau1.txt" u 1:($2>0 ? $2 : 1/0) w l lw 5 lc rgb cL4 title "{/Symbol x}_{max}=10^{4}", \
    (x>=1e0 && x<=1e3 ? 0.2/x : 1/0) w l ls 91 notitle

unset multiplot
unset output
