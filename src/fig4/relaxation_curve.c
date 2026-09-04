/*
 * Figure 4: Relaxation function for a selected alpha and upper cutoff.
 *
 * The program evaluates R(u) after a velocity step. To reproduce every curve
 * in Figure 4, set PARAM_ALPHA and PARAM_LMAX to each combination used by
 * load.plt and run the program once for each combination.
 *
 * Build from the repository root:
 *   mkdir -p build
 *   cc -O2 -Wall -Wextra -std=c11 src/fig4/relaxation_curve.c -lm -o build/relaxation_curve
 * Run from the repository root:
 *   (cd figures/fig4 && ../../build/relaxation_curve)
 * Output: figures/fig4/data/R_u_alpha*_Lmax*_tau*.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ============================================================
   Parameter settings
   ============================================================ */

#define NGAUSS 64

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Model parameters */
#define PARAM_ALPHA 1.5
#define PARAM_LMIN 0.01
#define PARAM_LMAX 10.0
#define PARAM_TAU 1.0
#define PARAM_V1 0.1
#define PARAM_V2 1.0
#define PARAM_C 1.0

/* Output range for u = t/tau */
#define U_MIN 1.0e-4
#define U_MAX 1.0e7

/* Number of output points per decade in u */
#define U_POINTS_PER_DECADE 160

/* Numerical integration parameter:
   number of logarithmic subintervals per decade */
#define N_PER_DECADE 20

/* Output directory */
#define OUTPUT_DIR "data"

/* ============================================================
   Main code
   ============================================================ */

typedef struct
{
    double alpha;
    double Lmin;
    double Lmax;
    double tau;
    double V1;
    double V2;
    double c;
} Params;

double xg[NGAUSS], wg[NGAUSS];

double min2(double a, double b)
{
    return a < b ? a : b;
}

void gauss_legendre_init(void)
{
    int n = NGAUSS;
    int m = (n + 1) / 2;
    double eps = 1e-14;

    for (int i = 0; i < m; i++)
    {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1;

        do
        {
            double p1 = 1.0;
            double p2 = 0.0;

            for (int j = 1; j <= n; j++)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }

            double pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > eps);

        double p1 = 1.0;
        double p2 = 0.0;

        for (int j = 1; j <= n; j++)
        {
            double p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
        }

        double pp = n * (z * p1 - p2) / (z * z - 1.0);

        xg[i] = -z;
        xg[n - 1 - i] = z;
        wg[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        wg[n - 1 - i] = wg[i];
    }
}

double S_survival(double x, Params *p)
{
    double a = p->alpha;
    double Lmin = p->Lmin;
    double Lmax = p->Lmax;

    if (x <= Lmin)
        return 1.0;
    if (x >= Lmax)
        return 0.0;

    return (pow(x, 1.0 - a) - pow(Lmax, 1.0 - a)) / (pow(Lmin, 1.0 - a) - pow(Lmax, 1.0 - a));
}

double Z_age(double theta, Params *p)
{
    return 1.0 + p->c * log(1.0 + theta / p->tau);
}

typedef struct
{
    Params *p;
    double t;
    int mode;
} Context;

/*
mode = 0: delta F(t)
  integrand = S(y) [ Z(t + (y - V2 t)/V1) - Z(y/V2) ]

mode = 1: Delta F_ss
  integrand = S(y) [ Z(y/V1) - Z(y/V2) ]
*/
double integrand_y(double y, void *ctx_void)
{
    Context *ctx = (Context *)ctx_void;
    Params *p = ctx->p;

    double Sy = S_survival(y, p);

    if (ctx->mode == 0)
    {
        double t = ctx->t;

        double theta1 = t + (y - p->V2 * t) / p->V1;
        double theta2 = y / p->V2;

        return Sy * (Z_age(theta1, p) - Z_age(theta2, p));
    }
    else
    {
        double theta1 = y / p->V1;
        double theta2 = y / p->V2;

        return Sy * (Z_age(theta1, p) - Z_age(theta2, p));
    }
}

double integrate_interval(double (*f)(double, void *), void *ctx,
                          double a, double b)
{
    if (b <= a)
        return 0.0;

    double mid = 0.5 * (a + b);
    double half = 0.5 * (b - a);
    double sum = 0.0;

    for (int i = 0; i < NGAUSS; i++)
    {
        double y = mid + half * xg[i];
        sum += wg[i] * f(y, ctx);
    }

    return half * sum;
}

double integrate_logsplit(double (*f)(double, void *), void *ctx,
                          double a, double b)
{
    if (b <= a)
        return 0.0;

    double total = 0.0;

    if (a <= 0.0)
    {
        double eps = min2(1e-10, b * 1e-12);
        if (eps <= 0.0)
            eps = 1e-12;

        total += integrate_interval(f, ctx, 0.0, eps);
        a = eps;
    }

    double loga = log10(a);
    double logb = log10(b);

    int nseg = (int)ceil((logb - loga) * N_PER_DECADE);
    if (nseg < 1)
        nseg = 1;

    for (int i = 0; i < nseg; i++)
    {
        double q0 = (double)i / (double)nseg;
        double q1 = (double)(i + 1) / (double)nseg;

        double left = pow(10.0, loga + (logb - loga) * q0);
        double right = pow(10.0, loga + (logb - loga) * q1);

        total += integrate_interval(f, ctx, left, right);
    }

    return total;
}

double deltaF_y_integral(double t, Params *p)
{
    double y0 = p->V2 * t;
    double y1 = p->Lmax;

    if (y0 >= y1)
        return 0.0;

    Context ctx;
    ctx.p = p;
    ctx.t = t;
    ctx.mode = 0;

    double total = 0.0;

    if (y0 < p->Lmin && p->Lmin < y1)
    {
        total += integrate_logsplit(integrand_y, &ctx, y0, p->Lmin);
        total += integrate_logsplit(integrand_y, &ctx, p->Lmin, y1);
    }
    else
    {
        total += integrate_logsplit(integrand_y, &ctx, y0, y1);
    }

    return total;
}

double DeltaFss_y_integral(Params *p)
{
    Context ctx;
    ctx.p = p;
    ctx.t = 0.0;
    ctx.mode = 1;

    double total = 0.0;

    total += integrate_logsplit(integrand_y, &ctx, 0.0, p->Lmin);
    total += integrate_logsplit(integrand_y, &ctx, p->Lmin, p->Lmax);

    return total;
}

int main(void)
{
    gauss_legendre_init();

    Params p;

    p.alpha = PARAM_ALPHA;
    p.Lmin = PARAM_LMIN;
    p.Lmax = PARAM_LMAX;
    p.tau = PARAM_TAU;
    p.V1 = PARAM_V1;
    p.V2 = PARAM_V2;
    p.c = PARAM_C;

    double u_physical_max = p.Lmax / (p.V2 * p.tau);
    double u_min = U_MIN;
    double u_max = U_MAX;

    if (u_max > u_physical_max)
        u_max = u_physical_max;

    if (u_min <= 0.0 || u_max <= u_min)
    {
        fprintf(stderr, "Invalid u range: u_min=%g, u_max=%g\n", u_min, u_max);
        fprintf(stderr, "Physical upper bound is Lmax/(V2*tau) = %g\n", u_physical_max);
        return 1;
    }

    double DFss = DeltaFss_y_integral(&p);

    char filename[256];
    snprintf(filename, sizeof(filename),
             "%s/R_u_alpha%.3f_Lmax%.0e_tau%.3g.txt",
             OUTPUT_DIR, p.alpha, p.Lmax, p.tau);

    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Cannot open output file: %s\n", filename);
        fprintf(stderr, "If the output directory does not exist, run: mkdir -p %s\n", OUTPUT_DIR);
        return 1;
    }

    fprintf(fp, "# alpha = %.15g\n", p.alpha);
    fprintf(fp, "# Lmin  = %.15g\n", p.Lmin);
    fprintf(fp, "# Lmax  = %.15g\n", p.Lmax);
    fprintf(fp, "# tau   = %.15g\n", p.tau);
    fprintf(fp, "# V1    = %.15g\n", p.V1);
    fprintf(fp, "# V2    = %.15g\n", p.V2);
    fprintf(fp, "# c     = %.15g\n", p.c);
    fprintf(fp, "# u_min_requested = %.15g\n", U_MIN);
    fprintf(fp, "# u_max_requested = %.15g\n", U_MAX);
    fprintf(fp, "# u_max_used      = %.15g\n", u_max);
    fprintf(fp, "# u_physical_max  = %.15g\n", u_physical_max);
    fprintf(fp, "# points_per_decade = %d\n", U_POINTS_PER_DECADE);
    fprintf(fp, "# DeltaFss = %.15e\n", DFss);
    fprintf(fp, "# columns: u=t/tau, R(u)\n");

    int nsteps = (int)ceil(log10(u_max / u_min) * U_POINTS_PER_DECADE);
    if (nsteps < 1)
        nsteps = 1;

    for (int i = 0; i <= nsteps; i++)
    {
        double q = (double)i / (double)nsteps;
        double u = u_min * pow(u_max / u_min, q);
        double t = u * p.tau;

        double R;

        if (t >= p.Lmax / p.V2)
        {
            R = 0.0;
        }
        else
        {
            double dF = deltaF_y_integral(t, &p);
            R = dF / DFss;
        }

        if (R < 0.0 && R > -1e-12)
            R = 0.0;

        fprintf(fp, "%.15e %.15e\n", u, R);
    }

    fclose(fp);

    printf("Output written to %s\n", filename);
    printf("alpha = %g\n", p.alpha);
    printf("Lmax  = %g\n", p.Lmax);
    printf("u range used = [%g, %g]\n", u_min, u_max);
    printf("u physical max = %g\n", u_physical_max);
    printf("points/decade = %d\n", U_POINTS_PER_DECADE);
    printf("DeltaFss = %.15e\n", DFss);

    return 0;
}
