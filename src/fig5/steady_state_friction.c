/*
 * Figure 5: Steady-state friction as a function of velocity.
 *
 * This program evaluates F_ss(V) in the N -> infinity model for alpha = 1.5,
 * 2.0, and 3.0. Each output file contains V and F_ss(V) in two columns.
 *
 * Build from the repository root:
 *   mkdir -p build
 *   cc -O2 -Wall -Wextra -std=c11 src/fig5/steady_state_friction.c -lm -o build/steady_state_friction
 * Run from the repository root:
 *   (cd figures/fig5 && ../../build/steady_state_friction)
 * Output: figures/fig5/data/Fss_alpha*.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ===================== Parameters ===================== */

static const double xi_min = 1.0e-2;
static const double xi_max = 1.0e3;

static const double c_age = 1.0;
static const double f0 = 1.0;
static const double tau0 = 1.0;

/* rho_g = rho_c, hence Phi = 1/2 */
static const double Phi = 0.5;

/* Velocity range */
static const double V_min = 1.0e-5;
static const double V_max = 1.0e6;

/* Number of output points per decade */
static const int n_per_decade = 100;

/* Adaptive Simpson tolerance */
static const double abs_tol = 1.0e-10;
static const double rel_tol = 1.0e-10;
static const int max_depth = 30;

/*
   Same numerical units as the reference calculation:
   tau = 1 and the reference velocity scale is 1,
   so that Lmin = xi_min and Lmax = xi_max.
*/
static const double Lmin = 1.0e-3;
static const double Lmax = 1.0e3;

static double alpha_now;
static double V_now;

/* ===================== Distribution ===================== */

static double norm_const(double alpha)
{
    if (fabs(alpha - 1.0) < 1.0e-12)
        return 1.0 / log(Lmax / Lmin);

    return (1.0 - alpha) /
           (pow(Lmax, 1.0 - alpha) - pow(Lmin, 1.0 - alpha));
}

static double mean_L(double alpha)
{
    const double A = norm_const(alpha);

    if (fabs(alpha - 2.0) < 1.0e-12)
        return A * log(Lmax / Lmin);

    return A *
           (pow(Lmax, 2.0 - alpha) - pow(Lmin, 2.0 - alpha)) /
           (2.0 - alpha);
}

static double survival_S(double L, double alpha)
{
    if (L <= Lmin)
        return 1.0;

    if (L >= Lmax)
        return 0.0;

    return (pow(Lmax, 1.0 - alpha) - pow(L, 1.0 - alpha)) /
           (pow(Lmax, 1.0 - alpha) - pow(Lmin, 1.0 - alpha));
}

/* ===================== Aging ===================== */

static double Z(double theta)
{
    return 1.0 + c_age * log1p(theta / tau0);
}

/*
   F_ss(V)
     = f0 Phi / <L>_c
       * int_0^Lmax dL Z(L/V) S(L).

   Split the integral into

      I1 = int_0^Lmin dL Z(L/V),

      I2 = int_Lmin^Lmax dL Z(L/V) S(L).

   For I2 use
      y = log(L/Lmin),
      L = Lmin exp(y).
*/

static double integrand_low(double L)
{
    return Z(L / V_now);
}

static double integrand_logL(double y)
{
    const double L = Lmin * exp(y);

    return Z(L / V_now) * survival_S(L, alpha_now) * L;
}

/* ===================== Adaptive Simpson ===================== */

typedef double (*func_t)(double);

static double simpson_basic(func_t f,
                            double a,
                            double b,
                            double fa,
                            double fm,
                            double fb)
{
    (void)f;

    return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
}

static double adaptive_simpson_rec(func_t f,
                                   double a,
                                   double b,
                                   double fa,
                                   double fm,
                                   double fb,
                                   double S,
                                   double tol,
                                   int depth)
{
    const double m = 0.5 * (a + b);

    const double lm = 0.5 * (a + m);
    const double rm = 0.5 * (m + b);

    const double flm = f(lm);
    const double frm = f(rm);

    const double Sleft =
        simpson_basic(f, a, m, fa, flm, fm);

    const double Sright =
        simpson_basic(f, m, b, fm, frm, fb);

    const double S2 = Sleft + Sright;

    if (depth <= 0 || fabs(S2 - S) <= 15.0 * tol)
        return S2 + (S2 - S) / 15.0;

    return adaptive_simpson_rec(
               f,
               a,
               m,
               fa,
               flm,
               fm,
               Sleft,
               0.5 * tol,
               depth - 1) +
           adaptive_simpson_rec(
               f,
               m,
               b,
               fm,
               frm,
               fb,
               Sright,
               0.5 * tol,
               depth - 1);
}

static double adaptive_simpson(func_t f,
                               double a,
                               double b,
                               double atol,
                               double rtol)
{
    if (b <= a)
        return 0.0;

    const double m = 0.5 * (a + b);

    const double fa = f(a);
    const double fm = f(m);
    const double fb = f(b);

    const double S =
        simpson_basic(f, a, b, fa, fm, fb);

    const double tol =
        fmax(atol, rtol * fabs(S));

    return adaptive_simpson_rec(
        f,
        a,
        b,
        fa,
        fm,
        fb,
        S,
        tol,
        max_depth);
}

/* ===================== F_ss(V) ===================== */

static double F_ss(double V,
                   double alpha,
                   double Lmean)
{
    V_now = V;
    alpha_now = alpha;

    /*
       0 <= L <= Lmin:
       S(L) = 1
    */
    const double I1 =
        adaptive_simpson(
            integrand_low,
            0.0,
            Lmin,
            abs_tol,
            rel_tol);

    /*
       Lmin <= L <= Lmax:
       y = log(L/Lmin)
    */
    const double ymax =
        log(Lmax / Lmin);

    const double I2 =
        adaptive_simpson(
            integrand_logL,
            0.0,
            ymax,
            abs_tol,
            rel_tol);

    return f0 * Phi * (I1 + I2) / Lmean;
}

/* ===================== Output ===================== */

static void output_alpha(double alpha)
{
    char filename[128];

    snprintf(
        filename,
        sizeof(filename),
        "data/Fss_alpha%.1f.txt",
        alpha);

    FILE *fp = fopen(filename, "w");

    if (fp == NULL)
    {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    const double Lmean = mean_L(alpha);

    const double logV_min = log10(V_min);
    const double logV_max = log10(V_max);

    const double decades =
        logV_max - logV_min;

    const int nout =
        (int)llround(decades * n_per_decade);

    fprintf(fp, "# alpha = %.16g\n", alpha);

    fprintf(fp,
            "# xi_min = %.16g, xi_max = %.16g\n",
            xi_min,
            xi_max);

    fprintf(fp,
            "# Lmin = %.16g, Lmax = %.16g\n",
            Lmin,
            Lmax);

    fprintf(fp,
            "# c = %.16g, f0 = %.16g, tau = %.16g, Phi = %.16g\n",
            c_age,
            f0,
            tau0,
            Phi);

    fprintf(fp,
            "# V range: %.16g -- %.16g\n",
            V_min,
            V_max);

    fprintf(fp,
            "# columns: V    F_ss(V)\n");

    for (int i = 0; i <= nout; ++i)
    {
        const double logV =
            logV_min +
            (logV_max - logV_min) * (double)i / (double)nout;

        const double V =
            pow(10.0, logV);

        const double F =
            F_ss(V, alpha, Lmean);

        fprintf(fp,
                "%.12e %.12e\n",
                V,
                F);
    }

    fclose(fp);

    printf("wrote %s\n", filename);
}

/* ===================== Main ===================== */

int main(void)
{
    const double alphas[] =
        {1.5, 2.0, 3.0};

    const int nalpha =
        (int)(sizeof(alphas) / sizeof(alphas[0]));

    for (int i = 0; i < nalpha; ++i)
        output_alpha(alphas[i]);

    return 0;
}
