/*
 * fft_analysis.h
 *
 *  Created on: 27/03/2026
 *      Author: joadj
 */

#ifndef INC_FFT_ANALYSIS_H_
#define INC_FFT_ANALYSIS_H_


#include <stdint.h>
#include <math.h>

/* =========================================================
 *  fft.h  —  Librería FFT portable para STM32
 *  Algoritmo: Cooley-Tukey radix-2 DIT (in-place)
 *  Entrada:   señal real (float) de N puntos, N = potencia de 2
 *  Salida:    magnitud y fase por bin de frecuencia
 * ========================================================= */

/* Tamaños soportados: 64, 128, 256, 512, 1024, 2048 */
#define FFT_MAX_SIZE  2048

/* Estructura que representa un número complejo */
typedef struct {
    float real;
    float imag;
} FFT_Complex;

/* Resultado por bin de frecuencia */
typedef struct {
    float magnitud;   /* Amplitud del bin (unidades de la señal entrada) */
    float fase;       /* Fase en radianes [-π, π] */
} FFT_Bin;

/* ── API principal ─────────────────────────────────────── */

/**
 * @brief  Calcula la FFT de una señal real.
 *
 * @param  entrada    Buffer de N muestras float (señal del acelerómetro).
 *                    IMPORTANTE: el buffer se usa como scratchpad y se modifica.
 * @param  salida     Buffer de salida con N/2 bins (solo mitad útil del espectro).
 * @param  N          Número de puntos. DEBE ser potencia de 2 (64..2048).
 * @param  fs_hz      Frecuencia de muestreo en Hz (para calcular frecuencia de cada bin).
 * @param  freq_bins  Array de N/2 floats donde se escriben las frecuencias de cada bin (Hz).
 *                    Puede ser NULL si no se necesita.
 * @return 0 si OK, -1 si N no es potencia de 2 o supera FFT_MAX_SIZE.
 */
int FFT_Calcular(float *entrada, FFT_Bin *salida, uint32_t N,
                 float fs_hz, float *freq_bins);

/**
 * @brief  Encuentra el bin de mayor magnitud (frecuencia dominante).
 *
 * @param  salida   Resultado de FFT_Calcular().
 * @param  N_2      Número de bins = N/2.
 * @param  fs_hz    Frecuencia de muestreo.
 * @return Frecuencia dominante en Hz.
 */
float FFT_FrecuenciaDominante(const FFT_Bin *salida, uint32_t N_2, float fs_hz);

/**
 * @brief  Normaliza las magnitudes dividiendo por N/2.
 *         Llama esto si quieres magnitudes en la misma unidad que la entrada.
 *
 * @param  salida  Buffer de bins a normalizar (in-place).
 * @param  N_2     Número de bins = N/2.
 */
void FFT_Normalizar(FFT_Bin *salida, uint32_t N_2);

#endif /* INC_FFT_ANALYSIS_H_ */
