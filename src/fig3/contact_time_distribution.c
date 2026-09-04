/*
 * Figure 3: Contact-time distributions after a velocity step.
 *
 * Computes P_t(theta) after a velocity step V1 -> V2:
 *
 *   P_t(theta) = (V2 / <L>) S(V2 theta)                         0 <= theta < t
 *              = (V1 / <L>) S(V1 (theta - t) + V2 t)             theta >= t
 *
 * where rho(L) is proportional to L^{-ALPHA} on [LMIN, LMAX], and
 * S(x) = int_x^infty rho(L) dL is the survival probability.
 *
 * Output times:
 *   t/tau = 0, 1e-3, 1e-1, 10, 1000
 *
 * Build from the repository root:
 *   mkdir -p build
 *   cc -O2 -Wall -Wextra -std=c11 src/fig3/contact_time_distribution.c -lm -o build/contact_time_distribution
 * Run from the repository root:
 *   (cd figures/fig3 && ../../build/contact_time_distribution)
 *
 * Output:
 *   figures/fig3/data/ptheta_t0.txt
 *   figures/fig3/data/ptheta_t1em3.txt
 *   figures/fig3/data/ptheta_t1em1.txt
 *   figures/fig3/data/ptheta_t10.txt
 *   figures/fig3/data/ptheta_t1000.txt
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

static const double ALPHA = 1.5;
static const double LMIN = 1.0e-2;
static const double LMAX = 1.0e3;

/* Not used in P_t(theta) itself, but kept for consistency with the model parameters. */
static const double C_AGING = 1.0;
static const double TAU0 = 1.0;

static const double V1 = 0.1;
static const double V2 = 1.0;

/*
 * Output theta grid.
 *
 * The theta grid is split into two parts so that small theta is resolved on a
 * log scale while the tail is sampled uniformly:
 *
 *   - N_GRID_LOG points, logarithmically spaced on [THETA_MIN_LOG, theta_split]
 *   - N_GRID_LIN points, linearly spaced on (theta_split, theta_max]
 *
 * A single theta = 0 point is written first.
 */
static const int N_GRID_LOG = 2000;
static const int N_GRID_LIN = 3000;
static const double THETA_MIN_LOG = 1.0e-4;
static const double THETA_SPLIT_FRAC = 0.1;

/* Fixed output times. */
enum
{
    N_T_OUTPUT = 5
};
static const double T_OVER_TAU_LIST[N_T_OUTPUT] = {
    0.0,
    1.0e-3,
    1.0e-1,
    1.0e+1,
    1.0e3};

static const char *T_LABEL_LIST[N_T_OUTPUT] = {
    "t0",
    "t1em3",
    "t1em1",
    "t10",
    "t1000"};

static double rho_norm_const(void)
{
    const double a = ALPHA;

    if (fabs(a - 1.0) < 1.0e-12)
    {
        return 1.0 / log(LMAX / LMIN);
    }

    return (a - 1.0) /
           (pow(LMIN, 1.0 - a) - pow(LMAX, 1.0 - a));
}

static double mean_L(void)
{
    const double a = ALPHA;
    const double c = rho_norm_const();

    if (fabs(a - 2.0) < 1.0e-12)
    {
        return c * log(LMAX / LMIN);
    }

    return c * (pow(LMAX, 2.0 - a) - pow(LMIN, 2.0 - a)) / (2.0 - a);
}

static double survival_S(double x)
{
    const double a = ALPHA;
    double denom;

    if (x <= LMIN)
    {
        return 1.0;
    }
    if (x >= LMAX)
    {
        return 0.0;
    }

    if (fabs(a - 1.0) < 1.0e-12)
    {
        return log(LMAX / x) / log(LMAX / LMIN);
    }

    denom = pow(LMIN, 1.0 - a) - pow(LMAX, 1.0 - a);
    return (pow(x, 1.0 - a) - pow(LMAX, 1.0 - a)) / denom;
}

static double P_theta(double t, double theta)
{
    const double mL = mean_L();

    if (theta < 0.0)
    {
        return 0.0;
    }

    if (theta < t)
    {
        return V2 * survival_S(V2 * theta) / mL;
    }

    return V1 * survival_S(V1 * (theta - t) + V2 * t) / mL;
}

static double theta_max_support(double t)
{
    double theta_max_new = 0.0;
    double theta_max_old = 0.0;

    if (t > 0.0)
    {
        theta_max_new = fmin(t, LMAX / V2);
    }

    if (V2 * t < LMAX)
    {
        theta_max_old = t + (LMAX - V2 * t) / V1;
    }

    return fmax(theta_max_new, theta_max_old);
}

static double output_time_over_tau(int idx)
{
    if (idx < 0 || idx >= N_T_OUTPUT)
    {
        fprintf(stderr, "Invalid time index: %d\n", idx);
        exit(EXIT_FAILURE);
    }

    return T_OVER_TAU_LIST[idx];
}

static void make_data_dir(void)
{
    if (mkdir("data", 0777) != 0 && errno != EEXIST)
    {
        perror("data");
        exit(EXIT_FAILURE);
    }
}

static void write_distribution(int idx)
{
    char filename[256];
    FILE *fp;
    const double t_over_tau = output_time_over_tau(idx);
    const double t = TAU0 * t_over_tau;
    const double theta_max = theta_max_support(t);
    double theta_split = THETA_SPLIT_FRAC * theta_max;
    double log_lo;
    double log_hi;
    double dlog;
    double dtheta_lin;
    int i;

    if (theta_max <= 0.0)
    {
        fprintf(stderr, "theta_max <= 0 at index %d\n", idx);
        exit(EXIT_FAILURE);
    }

    /*
     * Guard: theta_split must lie strictly between THETA_MIN_LOG and theta_max
     * so both sub-grids are well defined.
     */
    if (theta_split <= THETA_MIN_LOG)
    {
        theta_split = fmin(10.0 * THETA_MIN_LOG, 0.5 * theta_max);
    }
    if (theta_split <= THETA_MIN_LOG || theta_split >= theta_max)
    {
        theta_split = 0.5 * theta_max;
    }

    snprintf(filename, sizeof(filename), "data/ptheta_%s.txt", T_LABEL_LIST[idx]);

    fp = fopen(filename, "w");
    if (fp == NULL)
    {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "# Parameters: ALPHA=%g LMIN=%g LMAX=%g C_AGING=%g TAU0=%g V1=%g V2=%g\n",
            ALPHA, LMIN, LMAX, C_AGING, TAU0, V1, V2);
    fprintf(fp, "# index = %d\n", idx);
    fprintf(fp, "# t/tau = %.16e\n", t_over_tau);
    fprintf(fp, "# t = %.16e\n", t);
    fprintf(fp, "# theta_max = %.16e\n", theta_max);
    fprintf(fp, "# columns: theta  P_t(theta)\n");

    /* theta = 0. */
    fprintf(fp, "%.16e %.16e\n", 0.0, P_theta(t, 0.0));

    /* Log-spaced part: THETA_MIN_LOG .. theta_split. */
    log_lo = log10(THETA_MIN_LOG);
    log_hi = log10(theta_split);
    dlog = (log_hi - log_lo) / (double)(N_GRID_LOG - 1);
    for (i = 0; i < N_GRID_LOG; i++)
    {
        const double theta = pow(10.0, log_lo + dlog * (double)i);
        fprintf(fp, "%.16e %.16e\n", theta, P_theta(t, theta));
    }

    /* Linear part: just above theta_split .. theta_max. */
    dtheta_lin = (theta_max - theta_split) / (double)N_GRID_LIN;
    for (i = 1; i <= N_GRID_LIN; i++)
    {
        const double theta = theta_split + dtheta_lin * (double)i;
        fprintf(fp, "%.16e %.16e\n", theta, P_theta(t, theta));
    }

    fclose(fp);
}

int main(void)
{
    int idx;

    make_data_dir();

    for (idx = 0; idx < N_T_OUTPUT; idx++)
    {
        write_distribution(idx);
    }

    printf("Wrote %d files to ./data/\n", N_T_OUTPUT);
    printf("Time points t/tau:\n");
    for (idx = 0; idx < N_T_OUTPUT; idx++)
    {
        printf("  idx=%02d  t/tau=%.8g  t=%.8g  file=data/ptheta_%s.txt\n",
               idx, output_time_over_tau(idx), TAU0 * output_time_over_tau(idx),
               T_LABEL_LIST[idx]);
    }

    return EXIT_SUCCESS;
}
