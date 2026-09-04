/*
    Figure 6: Characteristic slip as a function of alpha.

    Characteristic slip of the relaxation,

        D_c = int_0^infinity R(D/V2) dD ,

    evaluated from the closed-form expression

        d_c = D_c / (V2 tau) = N(r) / [ 2 D(r) ]

    where, with  r = V1/V2,  xi = L/(V2 tau),

        ell1 = log(1 + xi/r),   ell2 = log(1 + xi)

        D(r) = < (xi + 1) ell2 - (xi + r) ell1 >
        N(r) = <xi^2> + <xi>
             + 1/(1-r) < (xi+r)^2 (ell2 - ell1) - (1-r)^2 ell2 >

    and <...> denotes the average over the truncated power law
        rho(xi) ∝ xi^{-alpha},   xi in [xi_min, xi_max].

    The integrals are evaluated by Simpson's rule in the logarithmic
    variable eta = log(xi):
        int f(xi) dxi = int f(e^eta) e^eta deta .

    Verified against a direct numerical evaluation of
    int_0^infinity R(u) du to 10 significant digits.

    Build from the repository root:
      mkdir -p build
      cc -O2 -Wall -Wextra -std=c11 src/fig6/characteristic_slip.c -lm -o build/characteristic_slip
    Run from the repository root:
      (cd figures/fig6 && ../../build/characteristic_slip)
    Output: figures/fig6/data/Dc_ximax*.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(dir) _mkdir(dir)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(dir) mkdir(dir, 0777)
#endif

/* ---------- parameters ---------- */

/* Same as Fig. 2 and as the R(u) scan:
   xi_min = 1e-2, xi_max = 1e3, r = V1/V2 = 0.1 */
#define LMIN 0.01
#define V1 0.1
#define V2 1.0
#define TAU 1.0

/* upper cutoffs to scan (set N_LMAX = 1 for a single curve) */
static const double LMAX_LIST[] = {1000.0};
#define N_LMAX (int)(sizeof(LMAX_LIST) / sizeof(LMAX_LIST[0]))

/* alpha loop (alpha > 1 required for a finite mean contact length) */
#define ALPHA_MIN 1.2
#define ALPHA_MAX 4.0
#define ALPHA_STEP 0.05

/* Simpson integration (must be even; 2000 already gives 10 digits) */
#define N_INTEGRAL 20000
#define EPS 1.0e-12

#define OUTPUT_DIR "data"

/* ---------- global variables ---------- */

static double ALPHA = 0.0;
static double XI_MIN = 0.0;
static double XI_MAX = 0.0;

static const double r = V1 / V2;

/* ---------- distribution rho(xi) ---------- */

static double rho_norm(void)
{
    if (fabs(ALPHA - 1.0) < EPS)
        return 1.0 / log(XI_MAX / XI_MIN);

    return (1.0 - ALPHA) /
           (pow(XI_MAX, 1.0 - ALPHA) - pow(XI_MIN, 1.0 - ALPHA));
}

static double rho_xi(double xi)
{
    if (xi < XI_MIN || xi > XI_MAX)
        return 0.0;

    return rho_norm() * pow(xi, -ALPHA);
}

/*
    Simpson integration in the logarithmic variable:
        int_{XI_MIN}^{XI_MAX} f(xi) dxi
      = int_{log XI_MIN}^{log XI_MAX} f(e^eta) e^eta deta
*/
static double simpson_log(double (*f)(double))
{
    int n = N_INTEGRAL;
    if (n % 2 != 0)
        n += 1;

    double a = log(XI_MIN);
    double b = log(XI_MAX);
    double h = (b - a) / (double)n;

    double s = f(XI_MIN) * XI_MIN + f(XI_MAX) * XI_MAX;

    for (int i = 1; i < n; i++)
    {
        double eta = a + h * (double)i;
        double xi = exp(eta);
        s += (i % 2 == 0 ? 2.0 : 4.0) * f(xi) * xi; /* Jacobian */
    }

    return s * h / 3.0;
}

/* ---------- averaged quantities ---------- */

static double integrand_norm_check(double xi)
{
    return rho_xi(xi);
}

static double integrand_m1(double xi)
{
    return rho_xi(xi) * xi;
}

static double integrand_m2(double xi)
{
    return rho_xi(xi) * xi * xi;
}

static double integrand_D(double xi)
{
    double ell1 = log1p(xi / r);
    double ell2 = log1p(xi);

    return rho_xi(xi) * ((xi + 1.0) * ell2 - (xi + r) * ell1);
}

static double integrand_N_extra(double xi)
{
    double ell1 = log1p(xi / r);
    double ell2 = log1p(xi);

    double term = (xi + r) * (xi + r) * (ell2 - ell1) -
                  (1.0 - r) * (1.0 - r) * ell2;

    return rho_xi(xi) * term;
}

/* ---------- main ---------- */

int main(void)
{
    if (fabs(1.0 - r) < EPS)
    {
        fprintf(stderr, "Error: r = V1/V2 is too close to 1;"
                        " the formula contains 1/(1-r).\n");
        return 1;
    }

    MKDIR(OUTPUT_DIR);

    int nalpha = (int)floor((ALPHA_MAX - ALPHA_MIN) / ALPHA_STEP + 0.5) + 1;

    for (int k = 0; k < N_LMAX; k++)
    {
        double Lmax = LMAX_LIST[k];

        XI_MIN = LMIN / (V2 * TAU);
        XI_MAX = Lmax / (V2 * TAU);

        if (XI_MIN <= 0.0 || XI_MAX <= XI_MIN)
        {
            fprintf(stderr, "Error: invalid xi range "
                            "(XI_MIN = %.6e, XI_MAX = %.6e).\n",
                    XI_MIN, XI_MAX);
            return 1;
        }

        char filename[256];
        snprintf(filename, sizeof(filename),
                 "%s/Dc_ximax%.0e.txt", OUTPUT_DIR, XI_MAX);

        FILE *fp = fopen(filename, "w");
        if (fp == NULL)
        {
            fprintf(stderr, "Error: cannot open file %s\n", filename);
            return 1;
        }

        fprintf(fp, "# characteristic slip  D_c = int_0^inf R(D/V2) dD\n");
        fprintf(fp, "# LMIN   %.15e\n", LMIN);
        fprintf(fp, "# LMAX   %.15e\n", Lmax);
        fprintf(fp, "# V1     %.15e\n", (double)V1);
        fprintf(fp, "# V2     %.15e\n", (double)V2);
        fprintf(fp, "# TAU    %.15e\n", (double)TAU);
        fprintf(fp, "# XI_MIN %.15e\n", XI_MIN);
        fprintf(fp, "# XI_MAX %.15e\n", XI_MAX);
        fprintf(fp, "# r      %.15e\n", r);
        fprintf(fp, "# columns: alpha  Dc  dc=Dc/(V2*tau)  Dc/LMAX\n");

        printf("Writing %s   (XI_MIN = %.3e, XI_MAX = %.3e, r = %.3f)\n",
               filename, XI_MIN, XI_MAX, r);

        for (int i = 0; i < nalpha; i++)
        {
            ALPHA = ALPHA_MIN + ALPHA_STEP * (double)i;

            double norm = simpson_log(integrand_norm_check);
            if (fabs(norm - 1.0) > 1.0e-6)
                fprintf(stderr, "Warning: normalization off at alpha = %.3f"
                                " (norm = %.10f)\n",
                        ALPHA, norm);

            double m1 = simpson_log(integrand_m1);
            double m2 = simpson_log(integrand_m2);

            double D_val = simpson_log(integrand_D);
            double N_extra = simpson_log(integrand_N_extra);

            double N_val = m2 + m1 + N_extra / (1.0 - r);

            double dc = N_val / (2.0 * D_val);
            double Dc = dc * V2 * TAU;

            fprintf(fp, "%.3f %.10e %.10e %.10e\n",
                    ALPHA, Dc, dc, Dc / Lmax);

            if (fabs(ALPHA * 10.0 - floor(ALPHA * 10.0 + 0.5)) < 1e-9 &&
                fabs(fmod(ALPHA + 1e-9, 0.2)) < 1e-6)
                printf("  alpha = %.2f   Dc = %.6e   Dc/Lmax = %.3e\n",
                       ALPHA, Dc, Dc / Lmax);
        }

        fclose(fp);
    }

    printf("Done.\n");
    return 0;
}
