/* ============================================================
   Figure 6: Relaxation curves and characteristic slip over alpha.

   Relaxation function R(u) of the frictional force after a velocity
   step, scanned over the exponent alpha of the contact-length
   distribution, together with the characteristic slip

       D_c = int_0^infinity R(D/V2) dD = V2 tau int_0^infinity R(u) du ,

   which reduces to the decay length when R is exponential.

   R(u) = deltaF(t) / DeltaFss  with  u = t/tau ,

       deltaF(t)  = int_{V2 t}^{Lmax} dy S(y)
                      [ Z( t + (y - V2 t)/V1 ) - Z( y/V2 ) ]
       DeltaFss   = int_0^{Lmax}     dy S(y)
                      [ Z( y/V1 ) - Z( y/V2 ) ]

   D_c is obtained without a double integral: exchanging the order of
   integration and using  int_0^y Z = y + c A(y),
   A(x) = (tau + x) log(1 + x/tau) - x,  gives

       D_c = V2 * int_0^{Lmax} dy S(y) g(y) / DeltaFss ,
       g(y) = c { [A(y/V2) - A(y/V1)] / (1 - V2/V1)
                  - (y/V2) log(1 + y/(V2 tau)) } .

   The leading terms cancel analytically in g(y), so no catastrophic
   cancellation occurs at small y.

   Build from the repository root:
     mkdir -p build
     cc -O2 -Wall -Wextra -std=c11 src/fig6/relaxation_alpha_scan.c -lm -o build/relaxation_alpha_scan
   Run from the repository root:
     (cd figures/fig6 && ../../build/relaxation_alpha_scan)
   Output: figures/fig6/data/R_u_alpha*.txt and figures/fig6/data/Dc_vs_alpha.txt
   ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NGAUSS 64

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------- Parameter settings ---------------- */

/* Model parameters (same as Fig. 2: xi_min = 1e-2, xi_max = 1e3) */
#define PARAM_LMIN 0.01
#define PARAM_LMAX 1000.0
#define PARAM_TAU 1.0
#define PARAM_V1 0.1
#define PARAM_V2 1.0
#define PARAM_C 1.0

/* Scan over alpha */
#define ALPHA_MIN 1.2
#define ALPHA_MAX 4.0
#define ALPHA_STEP 0.2

/* Output range for u = t/tau (clipped to the physical maximum) */
#define U_MIN 1.0e-4
#define U_MAX 1.0e7

/* Output points per decade in u */
#define U_POINTS_PER_DECADE 60

/* Extra grid resolving the terminal region u -> u_max,
   where R(u) vanishes as (u_max - u)^3.
   Points are placed logarithmically in (u_max - u). */
#define TAIL_FRACTION 0.1     /* tail starts at u = (1-f) u_max */
#define TAIL_MIN 1.0e-6       /* smallest (u_max - u)/u_max     */
#define TAIL_POINTS_PER_DECADE 30

/* Numerical integration: logarithmic subintervals per decade */
#define N_PER_DECADE 20

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

    return (pow(x, 1.0 - a) - pow(Lmax, 1.0 - a)) /
           (pow(Lmin, 1.0 - a) - pow(Lmax, 1.0 - a));
}

/* A(x) = (tau + x) log(1 + x/tau) - x = tau * H(x/tau) */
double A_aging(double x, Params *p)
{
    return (p->tau + x) * log1p(x / p->tau) - x;
}

typedef struct
{
    Params *p;
    double t;
    int mode;
} Context;

/*
mode = 0: deltaF(t)
  integrand = S(y) * c * log[ (1 + theta1/tau) / (1 + theta2/tau) ]
              theta1 = t + (y - V2 t)/V1 ,  theta2 = y/V2
              The ratio is evaluated inside a single logarithm to
              avoid cancellation when the two ages are close.

mode = 1: DeltaFss
  same with theta1 = y/V1

mode = 2: numerator of D_c
  integrand = S(y) * g(y)
*/
double integrand_y(double y, void *ctx_void)
{
    Context *ctx = (Context *)ctx_void;
    Params *p = ctx->p;

    double Sy = S_survival(y, p);
    if (Sy == 0.0)
        return 0.0;

    if (ctx->mode == 0)
    {
        double t = ctx->t;
        double theta1 = t + (y - p->V2 * t) / p->V1;
        double theta2 = y / p->V2;

        return Sy * p->c * log((1.0 + theta1 / p->tau) / (1.0 + theta2 / p->tau));
    }
    else if (ctx->mode == 1)
    {
        double theta1 = y / p->V1;
        double theta2 = y / p->V2;

        return Sy * p->c * log((1.0 + theta1 / p->tau) / (1.0 + theta2 / p->tau));
    }
    else
    {
        double a1 = y / p->V1;
        double a2 = y / p->V2;
        double beta = 1.0 - p->V2 / p->V1;

        double g = p->c * ((A_aging(a2, p) - A_aging(a1, p)) / beta
                           - a2 * log1p(a2 / p->tau));
        return Sy * g;
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

/* integral over [y0, Lmax], split at Lmin */
double integrate_over_y(Context *ctx, double y0, Params *p)
{
    double y1 = p->Lmax;

    if (y0 >= y1)
        return 0.0;

    if (y0 < p->Lmin && p->Lmin < y1)
        return integrate_logsplit(integrand_y, ctx, y0, p->Lmin) +
               integrate_logsplit(integrand_y, ctx, p->Lmin, y1);

    return integrate_logsplit(integrand_y, ctx, y0, y1);
}

double deltaF_y_integral(double t, Params *p)
{
    Context ctx;
    ctx.p = p;
    ctx.t = t;
    ctx.mode = 0;

    return integrate_over_y(&ctx, p->V2 * t, p);
}

double DeltaFss_y_integral(Params *p)
{
    Context ctx;
    ctx.p = p;
    ctx.t = 0.0;
    ctx.mode = 1;

    return integrate_over_y(&ctx, 0.0, p);
}

/* D_c = V2 * int_0^{Lmax} dy S(y) g(y) / DeltaFss */
double Dc_integral(Params *p, double DFss)
{
    Context ctx;
    ctx.p = p;
    ctx.t = 0.0;
    ctx.mode = 2;

    return p->V2 * integrate_over_y(&ctx, 0.0, p) / DFss;
}

/* ---------------- Output grid in u ---------------- */

int build_u_grid(double u_min, double u_max, double **out)
{
    double u_split = (1.0 - TAIL_FRACTION) * u_max;

    int n_main = (int)ceil(log10(u_split / u_min) * U_POINTS_PER_DECADE);
    if (n_main < 1)
        n_main = 1;

    int n_tail = (int)ceil(log10(TAIL_FRACTION / TAIL_MIN) *
                           TAIL_POINTS_PER_DECADE);
    if (n_tail < 1)
        n_tail = 1;

    int n = n_main + 1 + n_tail;
    double *u = (double *)malloc(sizeof(double) * n);
    if (u == NULL)
        return 0;

    for (int i = 0; i <= n_main; i++)
    {
        double q = (double)i / (double)n_main;
        u[i] = u_min * pow(u_split / u_min, q);
    }

    /* tail: (u_max - u) decreasing, hence u increasing */
    for (int i = 1; i <= n_tail; i++)
    {
        double q = (double)i / (double)n_tail;
        double gap = TAIL_FRACTION * pow(TAIL_MIN / TAIL_FRACTION, q);
        u[n_main + i] = u_max * (1.0 - gap);
    }

    *out = u;
    return n;
}

/* ---------------- Main ---------------- */

int main(void)
{
    gauss_legendre_init();

    Params p;
    p.Lmin = PARAM_LMIN;
    p.Lmax = PARAM_LMAX;
    p.tau = PARAM_TAU;
    p.V1 = PARAM_V1;
    p.V2 = PARAM_V2;
    p.c = PARAM_C;

    if (p.V1 == p.V2)
    {
        fprintf(stderr, "V1 must differ from V2.\n");
        return 1;
    }

    double u_physical_max = p.Lmax / (p.V2 * p.tau);
    double u_min = U_MIN;
    double u_max = min2(U_MAX, u_physical_max);

    if (u_min <= 0.0 || u_max <= u_min)
    {
        fprintf(stderr, "Invalid u range: u_min=%g, u_max=%g\n", u_min, u_max);
        return 1;
    }

    double *ugrid = NULL;
    int nu = build_u_grid(u_min, u_max, &ugrid);
    if (nu == 0)
    {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    char fname[512];

    snprintf(fname, sizeof(fname), "%s/Dc_vs_alpha.txt", OUTPUT_DIR);
    FILE *fdc = fopen(fname, "w");
    if (fdc == NULL)
    {
        fprintf(stderr, "Cannot open %s\n", fname);
        fprintf(stderr, "If the output directory does not exist, run: mkdir -p %s\n",
                OUTPUT_DIR);
        return 1;
    }
    fprintf(fdc, "# characteristic slip  D_c = int_0^inf R(D/V2) dD\n");
    fprintf(fdc, "# Lmin = %g, Lmax = %g, tau = %g, V1 = %g, V2 = %g, c = %g\n",
            p.Lmin, p.Lmax, p.tau, p.V1, p.V2, p.c);
    fprintf(fdc, "# columns: alpha  D_c  D_c/Lmax  u_c=D_c/(V2 tau)  R(u_c)\n");

    int nalpha = (int)floor((ALPHA_MAX - ALPHA_MIN) / ALPHA_STEP + 0.5) + 1;

    for (int ia = 0; ia < nalpha; ia++)
    {
        p.alpha = ALPHA_MIN + ALPHA_STEP * ia;

        double DFss = DeltaFss_y_integral(&p);
        double Dc = Dc_integral(&p, DFss);
        double uc = Dc / (p.V2 * p.tau);
        double Rc = (uc < u_physical_max)
                        ? deltaF_y_integral(uc * p.tau, &p) / DFss
                        : 0.0;

        snprintf(fname, sizeof(fname),
                 "%s/R_u_alpha%.2f.txt", OUTPUT_DIR, p.alpha);

        FILE *fp = fopen(fname, "w");
        if (fp == NULL)
        {
            fprintf(stderr, "Cannot open output file: %s\n", fname);
            return 1;
        }

        fprintf(fp, "# alpha = %.15g\n", p.alpha);
        fprintf(fp, "# Lmin  = %.15g\n", p.Lmin);
        fprintf(fp, "# Lmax  = %.15g\n", p.Lmax);
        fprintf(fp, "# tau   = %.15g\n", p.tau);
        fprintf(fp, "# V1    = %.15g\n", p.V1);
        fprintf(fp, "# V2    = %.15g\n", p.V2);
        fprintf(fp, "# c     = %.15g\n", p.c);
        fprintf(fp, "# xi_min = %.15g\n", p.Lmin / (p.V2 * p.tau));
        fprintf(fp, "# xi_max = %.15g\n", u_physical_max);
        fprintf(fp, "# DeltaFss = %.15e\n", DFss);
        fprintf(fp, "# D_c   = %.15e\n", Dc);
        fprintf(fp, "# u_c = D_c/(V2 tau) = %.15e\n", uc);
        fprintf(fp, "# R(u_c) = %.15e\n", Rc);
        fprintf(fp, "# columns: u=t/tau, R(u)\n");

        for (int i = 0; i < nu; i++)
        {
            double u = ugrid[i];
            double t = u * p.tau;
            double R;

            if (t >= p.Lmax / p.V2)
                R = 0.0;
            else
                R = deltaF_y_integral(t, &p) / DFss;

            if (R < 0.0 && R > -1e-12)
                R = 0.0;

            fprintf(fp, "%.15e %.15e\n", u, R);
        }

        fclose(fp);

        fprintf(fdc, "%.3f %.15e %.15e %.15e %.15e\n",
                p.alpha, Dc, Dc / p.Lmax, uc, Rc);

        printf("alpha = %.2f  D_c = %.6e  u_c = %.6e  R(u_c) = %.5f\n",
               p.alpha, Dc, uc, Rc);
    }

    fclose(fdc);
    free(ugrid);

    printf("\nOutput written to %s/\n", OUTPUT_DIR);
    printf("u range used = [%g, %g]  (physical max %g)\n",
           u_min, u_max, u_physical_max);
    return 0;
}
