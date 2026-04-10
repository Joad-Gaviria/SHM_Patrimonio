/*
 * fft_analysis.c
 *
 *  Created on: 27/03/2026
 *      Author: joadj
 */


#include "fft_analysis.h"
#include <string.h>   /* memcpy */

/* ── Constantes ────────────────────────────────────────── */
#define PI       3.14159265358979323846f
#define DOS_PI   6.28318530717958647692f

/* ── Helpers internos ──────────────────────────────────── */

/**
 * Verifica si n es potencia de 2.
 * Truco: n & (n-1) == 0 solo cuando n es potencia de 2.
 */
static int es_potencia_de_2(uint32_t n)
{
    return (n >= 2) && ((n & (n - 1)) == 0);
}

/**
 * Calcula log2 entero de n (n debe ser potencia de 2).
 */
static uint32_t log2_entero(uint32_t n)
{
    uint32_t log = 0;
    while (n > 1) { n >>= 1; log++; }
    return log;
}

/**
 * Bit-reversal permutation.
 * Reordena el array complejo buf[N] según el índice con bits invertidos.
 * Esto es necesario para que el algoritmo DIT trabaje in-place.
 *
 * Ejemplo con N=8 (3 bits):
 *   índice 0 (000) → 0 (000)   sin cambio
 *   índice 1 (001) → 4 (100)   intercambio
 *   índice 2 (010) → 2 (010)   sin cambio
 *   índice 3 (011) → 6 (110)   intercambio
 *   etc.
 */
static void bit_reversal(FFT_Complex *buf, uint32_t N)
{
    uint32_t bits = log2_entero(N);

    for (uint32_t i = 0; i < N; i++) {
        /* Calcula el índice con bits invertidos */
        uint32_t j = 0;
        uint32_t tmp = i;
        for (uint32_t b = 0; b < bits; b++) {
            j = (j << 1) | (tmp & 1);
            tmp >>= 1;
        }
        /* Intercambia solo una vez (i < j evita doble swap) */
        if (i < j) {
            FFT_Complex swap = buf[i];
            buf[i] = buf[j];
            buf[j] = swap;
        }
    }
}

/**
 * Núcleo FFT Cooley-Tukey radix-2 DIT in-place.
 *
 * Funciona en 3 niveles conceptuales:
 *
 *  1. BIT-REVERSAL: reordena las muestras.
 *
 *  2. BUTTERFLY STAGES: itera log2(N) etapas.
 *     En cada etapa el array se divide en grupos de tamaño 'len'.
 *     Cada grupo ejecuta len/2 "mariposas".
 *
 *  3. BUTTERFLY (mariposa):
 *     Dados dos elementos A y B, y un factor de giro W = e^(-j*2π*k/len):
 *       A' = A + W*B
 *       B' = A - W*B
 *     Esto combina los resultados de sub-FFTs más pequeñas.
 *
 *  Complejidad: O(N log N)  en lugar de O(N²) de la DFT directa.
 */
static void fft_cooley_tukey(FFT_Complex *buf, uint32_t N)
{
    /* Paso 1: reordenar por bit-reversal */
    bit_reversal(buf, N);

    /* Paso 2: etapas de mariposas */
    /* 'len' es el tamaño del grupo actual: 2, 4, 8, ..., N */
    for (uint32_t len = 2; len <= N; len <<= 1)
    {
        /* Ángulo del factor de giro para esta etapa */
        float ang = -DOS_PI / (float)len;

        /* Factor de giro unitario W_len = e^(j*ang) */
        FFT_Complex W_len = { cosf(ang), sinf(ang) };

        /* Itera sobre todos los grupos de tamaño 'len' */
        for (uint32_t i = 0; i < N; i += len)
        {
            /* Factor de giro acumulado W = W_len^k */
            FFT_Complex W = { 1.0f, 0.0f };

            /* Mariposas dentro del grupo */
            for (uint32_t k = 0; k < len / 2; k++)
            {
                /* Elementos de la mariposa */
                FFT_Complex *A = &buf[i + k];
                FFT_Complex *B = &buf[i + k + len / 2];

                /* W * B  (multiplicación compleja) */
                FFT_Complex WB = {
                    W.real * B->real - W.imag * B->imag,
                    W.real * B->imag + W.imag * B->real
                };

                /* Mariposa: A' = A + WB,  B' = A - WB */
                B->real = A->real - WB.real;
                B->imag = A->imag - WB.imag;
                A->real = A->real + WB.real;
                A->imag = A->imag + WB.imag;

                /* Avanza el factor de giro: W *= W_len */
                float Wr_new = W.real * W_len.real - W.imag * W_len.imag;
                float Wi_new = W.real * W_len.imag + W.imag * W_len.real;
                W.real = Wr_new;
                W.imag = Wi_new;
            }
        }
    }
}

/* ── Buffer interno estático ───────────────────────────── */
/*
 * Se usa un buffer estático para no hacer malloc.
 * No es reentrant, pero en un micro sin RTOS es seguro.
 */
static FFT_Complex _buf[FFT_MAX_SIZE];

/* ── Implementación de la API pública ──────────────────── */

int FFT_Calcular(float *entrada, FFT_Bin *salida, uint32_t N,
                 float fs_hz, float *freq_bins)
{
    /* Validación de entrada */
    if (!es_potencia_de_2(N) || N > FFT_MAX_SIZE) {
        return -1;
    }

    /* ── Copiar señal real al buffer complejo ── */
    /*
     * Para señal real, parte imaginaria = 0.
     * La FFT de una señal real tiene simetría conjugada:
     * X[N-k] = X[k]*, por eso solo usamos los primeros N/2 bins.
     */
    for (uint32_t i = 0; i < N; i++) {
        _buf[i].real = entrada[i];
        _buf[i].imag = 0.0f;
    }

    /* ── Ejecutar FFT ── */
    fft_cooley_tukey(_buf, N);

    /* ── Calcular magnitud y fase para cada bin útil ── */
    /*
     * Solo calculamos N/2 bins porque:
     *  - El bin 0 es la componente DC (frecuencia 0 Hz)
     *  - Los bins 1 .. N/2-1 son las frecuencias positivas
     *  - Los bins N/2 .. N-1 son el espejo (frecuencias negativas)
     *    y no aportan información nueva para señal real.
     *
     * Frecuencia de cada bin k:  f_k = k * fs / N
     */
    uint32_t N_2 = N / 2;

    for (uint32_t k = 0; k < N_2; k++) {
        float re = _buf[k].real;
        float im = _buf[k].imag;

        salida[k].magnitud = sqrtf(re * re + im * im);
        salida[k].fase     = atan2f(im, re);   /* radianes en [-π, π] */

        if (freq_bins != NULL) {
            freq_bins[k] = (float)k * fs_hz / (float)N;
        }
    }

    return 0;
}

float FFT_FrecuenciaDominante(const FFT_Bin *salida, uint32_t N_2, float fs_hz)
{
    uint32_t idx_max = 0;
    float    mag_max = salida[0].magnitud;

    /*
     * Ignoramos el bin 0 (DC) porque suele ser el mayor
     * y no representa una frecuencia de vibración real.
     * Empezamos desde k=1.
     */
    for (uint32_t k = 1; k < N_2; k++) {
        if (salida[k].magnitud > mag_max) {
            mag_max = salida[k].magnitud;
            idx_max = k;
        }
    }

    /* f = k * fs / N = k * fs / (2 * N_2) */
    return (float)idx_max * fs_hz / (float)(2 * N_2);
}

void FFT_Normalizar(FFT_Bin *salida, uint32_t N_2)
{
    /*
     * La FFT sin normalizar acumula la energía de N muestras.
     * Dividir por N/2 escala la magnitud al rango de la señal original.
     * El factor 2 compensa el espejo descartado (bins N/2..N-1).
     */
    float escala = 1.0f / (float)N_2;
    for (uint32_t k = 0; k < N_2; k++) {
        salida[k].magnitud *= escala;
    }
}
