#include <stdio.h>
#include <math.h>
#include <complex.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int ilog2(int N)
{
    int result;
    for (result = 0; N != 1; N >>= 1, result++);
    return result;
}

typedef struct {
    int N;
    int log2_N;
    int stage;
    int gcount;
    int i1;
    int i2;
    int k1;
    int k2;
    int w;
    int dw;
} index_t;

void index_init(index_t *m, int N)
{
    m->N = N;
    m->log2_N = ilog2(N);

    m->stage = 1;
    m->gcount = 0;
    m->i1 = 0;
    m->i2 = N/2;
    m->k1 = N;
    m->k2 = N/2;
    m->w = 0;
    m->dw = 1;
}

int index_update(index_t *m)
{
    // Return 0 if successful or 1 when finished
    // From Meyer-Baese 2014, page 446

    m->i1 = m->i1 + m->k1; // Next bufferfly in group
    m->i2 = m->i1 + m->k2;

    if (m->i1 >= m->N - 1) { // All butterfiles in group done?
        m->gcount = m->gcount + 1;
        m->i1 = m->gcount;
        m->i2 = m->i1 + m->k2;
        if (m->gcount >= m->k2) { // All groups done in stage?
            m->gcount = 0;
            m->i1 = 0;
            m->i2 = m->k2;
            m->dw = 2*m->dw;
            m->stage = m->stage + 1;
            if (m->stage > m->log2_N) { // All stages done?
                return 1;
            }
            else { // Start new stage
                m->k1 = m->k2;
                m->k2 = m->k2/2;
                m->i1 = 0;
                m->i2 = m->k2;
                m->w = 0;
            }
        }
        else { // Start new group
            m->i1 = m->gcount;
            m->i2 = m->i1 + m->k2;
            m->w = m->w + m->dw;
        }
    }
    return 0;
}

int bitreverse(int x, int width)
{
    int n;
    int y = 0;
    for (n = 0; n < width; n++) {
        if ((x >> n) & 1) {
            y |= 1 << ((width - 1) - n);
        }
    }
    return y;
}

void butterfly_update_inline(complex float *x1, complex float *x2, complex float w)
{
    complex float new_x1, new_x2;
    new_x1 = *x1 + *x2;
    new_x2 = (*x1 - *x2)*w;
    *x1 = new_x1;
    *x2 = new_x2;
}

void twiddle_factor_init(complex float *w, int N)
{
    int n;
    complex float w_base = cexp(-I*2*M_PI/N);
    w[0] = 1.0;
    for (n = 1; n < N; n++) {
        w[n] = w[n-1] * w_base;
    }
}

int main(int argc, char *argv[])
{
    #define N 8

    int n;
    index_t m;

    // Example input data
    complex float x[N] = {0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -2.0, -1.0};

    // Output data
    complex float y[N];

    // Twiddle factors
    complex float W[N];

    // Initialize twiddle factors
    twiddle_factor_init(W, N);

    // Iterate over all butterflies
    index_init(&m, N);
    do {
        butterfly_update_inline(&x[m.i1], &x[m.i2], W[m.w]);
    } while (index_update(&m) == 0);

    // Apply bit-reverse to create the output array and print the result
    for (n = 0; n < N; n++) {
        int n_rev = bitreverse(n, m.log2_N);
        y[n] = x[n_rev];
        printf("y[%d] = x[%d] = %+.2f%+.2fj\n", n, n_rev, creal(y[n]), cimag(y[n]));
    }

    return 0;
}