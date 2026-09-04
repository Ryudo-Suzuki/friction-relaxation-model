/*
 * Figure 2: Finite-N simulation of friction relaxation.
 *
 * The program simulates independent asperities, applies a velocity step from
 * V1 to V2, and records F(t) = (1/N) sum_i f_i(t) for three system sizes.
 * The fixed random seed makes the generated datasets reproducible.
 *
 * Build from the repository root:
 *   mkdir -p build
 *   cc -O2 -Wall -Wextra -std=c11 src/fig2/finite_n_relaxation.c -lm -o build/finite_n_relaxation
 * Run from the repository root:
 *   (cd figures/fig2 && ../../build/finite_n_relaxation)
 * Output: figures/fig2/data/F_alpha*.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/*=============================== パラメータ ===============================*/

static const double ALPHA = 1.5;
static const double LMIN = 1.0e-2;
static const double LMAX = 1.0e3;

static const double C_AGING = 1.0;
static const double TAU0 = 1.0;

static const double V1 = 0.1;
static const double V2 = 1.0;

/* 速度ステップ前に定常状態へ近づけるための滑走時間（この中の最後 T_BEFORE も含む） */
static const double T_BURN = 10000.0;

/* 速度ステップ前に出力する時間範囲 (t = -T_BEFORE から 0 まで速度 V1 で出力) */
static const double T_BEFORE = 500.0;

/* 速度ステップ後に出力する時間範囲 */
static const double T_AFTER = 2000.0;

/* 出力間隔 */
static const double OUTPUT_DT = 0.1;

/* 計算したい N */
static const int N_TARGETS[] = {1000, 10000, 100000};
enum
{
    NUM_TARGETS = sizeof(N_TARGETS) / sizeof(N_TARGETS[0])
};

static const uint64_t BASE_SEED = 42ULL;

/*=============================== 乱数 ===============================*/

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z;

    *x += 0x9e3779b97f4a7c15ULL;
    z = *x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static double urand_open(uint64_t *state)
{
    /* (0,1) の一様乱数 */
    uint64_t z = splitmix64(state);
    return ((double)(z >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/*=============================== 分布と接触力 ===============================*/

static double Finv_powerlaw(double u, double alpha, double lmin, double lmax)
{
    const double eps = 1.0e-12;

    if (fabs(alpha - 1.0) < eps)
    {
        return lmin * pow(lmax / lmin, u);
    }
    else
    {
        const double one_minus_a = 1.0 - alpha;
        const double A = pow(lmin, one_minus_a);
        const double B = pow(lmax, one_minus_a);
        const double T = A + (B - A) * u;
        return pow(T, 1.0 / one_minus_a);
    }
}

static double sample_powerlaw(uint64_t *state)
{
    double u = urand_open(state);
    return Finv_powerlaw(u, ALPHA, LMIN, LMAX);
}

static double contact_force(double theta)
{
    return 1.0 + C_AGING * log(1.0 + theta / TAU0);
}

/*=============================== 1本の矢印の状態 ===============================*/

typedef struct
{
    uint64_t rng_state;

    int occ;      /* 1: アスペリティ占有領域, 0: ギャップ領域 */
    double rem;   /* 現在の領域の右端までの残り長さ */
    double theta; /* 現在のアスペリティに入ってからの接触時間 */
    double f;     /* 現在の接触力 */
} Arrow;

static void init_arrow(Arrow *a, int seed_index)
{
    /* seed_index をハッシュして state を作る。
       state を 42 + GOLDEN*(s+1) と等差で与えると、隣り合う矢印の状態差が
       splitmix64 の内部増分と一致し、各矢印の「最初の1サンプル」に系統的な
       相関（モーメントの偏り）が出る。下のように一度かき混ぜてから使う。 */
    uint64_t z = BASE_SEED + 0x9e3779b97f4a7c15ULL * (uint64_t)(seed_index + 1);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    a->rng_state = z;
    (void)splitmix64(&a->rng_state); /* 念のため 1 ステップ進める */

    /* 最初はアスペリティ占有領域から始める */
    a->occ = 1;
    a->rem = sample_powerlaw(&a->rng_state);
    a->theta = 0.0;
    a->f = contact_force(0.0);
}

static void enter_next_segment(Arrow *a)
{
    a->occ = 1 - a->occ;
    a->rem = sample_powerlaw(&a->rng_state);
    a->theta = 0.0;

    if (a->occ)
    {
        a->f = contact_force(0.0);
    }
    else
    {
        a->f = 0.0;
    }
}

static void advance_arrow(Arrow *a, double v, double dt)
{
    double h = dt;

    while (h > 0.0)
    {
        double time_to_boundary = a->rem / v;

        if (h < time_to_boundary)
        {
            /* 現在の領域内で止まる */
            a->rem -= v * h;

            if (a->occ)
            {
                a->theta += h;
                a->f = contact_force(a->theta);
            }
            else
            {
                a->theta = 0.0;
                a->f = 0.0;
            }

            h = 0.0;
        }
        else
        {
            /* 領域境界を越える */
            if (a->occ)
            {
                a->theta += time_to_boundary;
            }

            h -= time_to_boundary;
            enter_next_segment(a);

            /* 丸め誤差対策 */
            if (h < 1.0e-14)
            {
                h = 0.0;
            }
        }
    }
}

/*=============================== 出力 ===============================*/

static void write_result(const char *filename, const double *sumF, int nout, int N)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        perror(filename);
        exit(1);
    }

    fprintf(fp, "# alpha %.12g\n", ALPHA);
    fprintf(fp, "# Lmin  %.12g\n", LMIN);
    fprintf(fp, "# Lmax  %.12g\n", LMAX);
    fprintf(fp, "# c     %.12g\n", C_AGING);
    fprintf(fp, "# tau   %.12g\n", TAU0);
    fprintf(fp, "# V1    %.12g\n", V1);
    fprintf(fp, "# V2    %.12g\n", V2);
    fprintf(fp, "# N     %d\n", N);
    fprintf(fp, "# T_BURN    %.12g\n", T_BURN);
    fprintf(fp, "# T_BEFORE  %.12g\n", T_BEFORE);
    fprintf(fp, "# OUTPUT_DT %.12g\n", OUTPUT_DT);
    fprintf(fp, "# t F\n");

    for (int i = 0; i < nout; i++)
    {
        double t = -T_BEFORE + OUTPUT_DT * (double)i;
        double F = sumF[i] / (double)N;
        fprintf(fp, "%.10f %.10f\n", t, F);
    }

    fclose(fp);
}

/*=============================== 本体 ===============================*/

int main(void)
{
    int Nmax = N_TARGETS[NUM_TARGETS - 1];

    /* 出力点数: t = -T_BEFORE, ..., 0, ..., T_AFTER */
    int nbefore = (int)floor(T_BEFORE / OUTPUT_DT + 0.5); /* t<0 の点数 (t=0 を除く) */
    int nafter = (int)floor(T_AFTER / OUTPUT_DT + 0.5);   /* t>0 の点数 (t=0 を除く) */
    int nout = nbefore + nafter + 1;
    int istep = nbefore; /* t=0 (速度ステップ) に対応するインデックス */

    double *sumF[NUM_TARGETS];

    for (int k = 0; k < NUM_TARGETS; k++)
    {
        sumF[k] = (double *)calloc((size_t)nout, sizeof(double));
        if (sumF[k] == NULL)
        {
            fprintf(stderr, "calloc failed\n");
            return 1;
        }
    }

    for (int s = 0; s < Nmax; s++)
    {
        Arrow a;
        init_arrow(&a, s);

        /* 速度 V1 で滑らせて定常状態に近づける。
           最後の T_BEFORE の間は出力するので、その手前まで進める。 */
        advance_arrow(&a, V1, T_BURN - T_BEFORE);

        /* t = -T_BEFORE: 出力開始時点の状態を記録 */
        for (int k = 0; k < NUM_TARGETS; k++)
        {
            if (s < N_TARGETS[k])
            {
                sumF[k][0] += a.f;
            }
        }

        /* -T_BEFORE < t <= 0: 速度 V1 で時間発展しつつ記録 */
        for (int i = 1; i <= istep; i++)
        {
            advance_arrow(&a, V1, OUTPUT_DT);

            for (int k = 0; k < NUM_TARGETS; k++)
            {
                if (s < N_TARGETS[k])
                {
                    sumF[k][i] += a.f;
                }
            }
        }

        /* t > 0: 速度 V2 で時間発展 */
        for (int i = istep + 1; i < nout; i++)
        {
            advance_arrow(&a, V2, OUTPUT_DT);

            for (int k = 0; k < NUM_TARGETS; k++)
            {
                if (s < N_TARGETS[k])
                {
                    sumF[k][i] += a.f;
                }
            }
        }

        if ((s + 1) % 10000 == 0)
        {
            fprintf(stderr, "finished %d / %d arrows\n", s + 1, Nmax);
        }
    }

    for (int k = 0; k < NUM_TARGETS; k++)
    {
        char filename[256];

        snprintf(filename, sizeof(filename),
                 "data/F_alpha%.3g_lmin%.3g_lmax%.3g_N%d_dt%.3g.txt",
                 ALPHA, LMIN, LMAX, N_TARGETS[k], OUTPUT_DT);

        write_result(filename, sumF[k], nout, N_TARGETS[k]);
    }

    for (int k = 0; k < NUM_TARGETS; k++)
    {
        free(sumF[k]);
    }

    return 0;
}
