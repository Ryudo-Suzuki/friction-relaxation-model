/*
 * Figure 2: Analytical friction relaxation in the N -> infinity limit.
 *
 * This program evaluates the continuum expression for the friction force
 * after an instantaneous velocity step from V1 to V2.
 *
 * Build from the repository root:
 *   mkdir -p build
 *   cc -O2 -Wall -Wextra -std=c11 src/fig2/analytic_relaxation.c -lm -o build/analytic_relaxation
 * Run from the repository root:
 *   (cd figures/fig2 && ../../build/analytic_relaxation)
 * Output: figures/fig2/data/F_Ninfty_*.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ===================== Parameters ===================== */

static const double alpha = 1.5;
static const double Lmin = 1.0e-2;
static const double Lmax = 1.0e3;

static const double c_age = 1.0;
static const double tau0 = 1.0;

static const double V1 = 0.1;
static const double V2 = 1.0;

static const double dt_out = 1.0;
static const double T_after = 2000.0;

/* Numerical integration step for theta */
static const double dtheta = 1.0e-1;

/* ===================== Distribution ===================== */

static double norm_const(void)
{
    if (fabs(alpha - 1.0) < 1.0e-12)
    {
        return 1.0 / log(Lmax / Lmin);
    }
    else
    {
        return (1.0 - alpha) /
               (pow(Lmax, 1.0 - alpha) - pow(Lmin, 1.0 - alpha));
    }
}

static double mean_L(void)
{
    const double A = norm_const();

    if (fabs(alpha - 2.0) < 1.0e-12)
    {
        return A * log(Lmax / Lmin);
    }
    else
    {
        return A * (pow(Lmax, 2.0 - alpha) - pow(Lmin, 2.0 - alpha)) / (2.0 - alpha);
    }
}

static double survival_S(double x)
{
    if (x <= Lmin)
    {
        return 1.0;
    }
    if (x >= Lmax)
    {
        return 0.0;
    }

    if (fabs(alpha - 1.0) < 1.0e-12)
    {
        return log(Lmax / x) / log(Lmax / Lmin);
    }
    else
    {
        return (pow(Lmax, 1.0 - alpha) - pow(x, 1.0 - alpha)) /
               (pow(Lmax, 1.0 - alpha) - pow(Lmin, 1.0 - alpha));
    }
}

/* ===================== Aging function ===================== */

static double Z(double theta)
{
    return 1.0 + c_age * log(1.0 + theta / tau0);
}

/* ===================== Integrands ===================== */

static double integrand_new(double theta)
{
    return Z(theta) * V2 * survival_S(V2 * theta);
}

static double integrand_old(double theta, double t)
{
    const double x = V2 * t + V1 * (theta - t);
    return Z(theta) * V1 * survival_S(x);
}

/* ===================== Simpson integration ===================== */

static double simpson_new(double a, double b)
{
    if (b <= a)
    {
        return 0.0;
    }

    int n = (int)ceil((b - a) / dtheta);
    if (n < 2)
    {
        n = 2;
    }
    if (n % 2 == 1)
    {
        n++;
    }

    const double h = (b - a) / (double)n;

    double sum = integrand_new(a) + integrand_new(b);

    for (int i = 1; i < n; i++)
    {
        const double x = a + h * (double)i;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * integrand_new(x);
    }

    return sum * h / 3.0;
}

static double simpson_old(double a, double b, double t)
{
    if (b <= a)
    {
        return 0.0;
    }

    int n = (int)ceil((b - a) / dtheta);
    if (n < 2)
    {
        n = 2;
    }
    if (n % 2 == 1)
    {
        n++;
    }

    const double h = (b - a) / (double)n;

    double sum = integrand_old(a, t) + integrand_old(b, t);

    for (int i = 1; i < n; i++)
    {
        const double x = a + h * (double)i;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * integrand_old(x, t);
    }

    return sum * h / 3.0;
}

/* ===================== F_Ninfty(t) ===================== */

static double F_Ninfty(double t, double Lmean)
{
    /*
        F(t) = 1/(2<L>) [
            int_0^t dtheta Z(theta) V2 S(V2 theta)
          + int_t^inf dtheta Z(theta) V1 S(V2 t + V1(theta - t))
        ]
    */

    double I1 = 0.0;
    double I2 = 0.0;

    /* S(V2 theta)=0 for V2 theta >= Lmax */
    const double upper1 = fmin(t, Lmax / V2);
    I1 = simpson_new(0.0, upper1);

    /*
        S(V2 t + V1(theta-t))=0 when
        V2 t + V1(theta-t) >= Lmax.
    */
    if (V2 * t < Lmax)
    {
        const double upper2 = t + (Lmax - V2 * t) / V1;
        I2 = simpson_old(t, upper2, t);
    }
    else
    {
        I2 = 0.0;
    }

    return (I1 + I2) / (2.0 * Lmean);
}

/* ===================== Main ===================== */

int main(void)
{
    const double Lmean = mean_L();

    char filename[256];

    snprintf(filename, sizeof(filename),
             "data/F_Ninfty_alpha%.3g_lmin%.3g_lmax%.3g"
             "_V1%.3g_V2%.3g_dt%.3g_dtheta%.3g.txt",
             alpha, Lmin, Lmax,
             V1, V2, dt_out, dtheta);

    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        perror(filename);
        return 1;
    }

    fprintf(fp, "# N -> infinity result\n");
    fprintf(fp, "# alpha %.12g\n", alpha);
    fprintf(fp, "# Lmin %.12g\n", Lmin);
    fprintf(fp, "# Lmax %.12g\n", Lmax);
    fprintf(fp, "# c %.12g\n", c_age);
    fprintf(fp, "# tau %.12g\n", tau0);
    fprintf(fp, "# V1 %.12g\n", V1);
    fprintf(fp, "# V2 %.12g\n", V2);
    fprintf(fp, "# <L> %.12g\n", Lmean);
    fprintf(fp, "# dt_out %.12g\n", dt_out);
    fprintf(fp, "# dtheta %.12g\n", dtheta);
    fprintf(fp, "# t F_Ninfty\n");

    const int nout = (int)floor(T_after / dt_out + 0.5);

    for (int i = 0; i <= nout; i++)
    {
        const double t = dt_out * (double)i;
        const double F = F_Ninfty(t, Lmean);
        fprintf(fp, "%.10f %.12f\n", t, F);
    }

    fclose(fp);

    return 0;
}
