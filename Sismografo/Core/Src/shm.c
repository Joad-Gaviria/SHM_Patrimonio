/*
 * smh.c
 *
 *  Created on: 8/04/2026
 *      Author: joadj
 *
 * =========================================================
 *  shm.c  —  Structural Health Monitoring
 *
 *  Implementación de:
 *    1. Ventana Hanning
 *    2. Método de Welch (promediado de espectros)
 *    3. Detección de picos modales
 *    4. Línea base y Damage Index
 * ========================================================= */

#include <string.h>   /* memcpy, memset */
#include <math.h>     /* fabsf           */
#include <shm.h>
#include <stdio.h>
#include "stm32h7xx_hal.h"


/* ── Constantes internas ───────────────────────────────── */
#define PI      3.14159265358979323846f

/* ── Buffers estáticos internos ────────────────────────── */
/*
 * Se usan buffers estáticos para evitar malloc().
 * No son reentrantes, pero en un micro sin RTOS es seguro.
 */

/* Buffer de trabajo para una ventana individual */
static float _ventana_work[SHM_N];

/* Buffer acumulador del espectro promediado */
static float _acum_mag[SHM_N / 2];
static float _acum_fase[SHM_N / 2];
static FFT_Bin _fft_bin_tmp[SHM_N / 2];   /* ← este es el que falta */

/* ── Buffer público del espectro Welch ── */
FFT_Bin shm_espectro_publico[SHM_N / 2];
static float _freq_bins[SHM_N / 2];

/* ══════════════════════════════════════════════════════════
 *  1. VENTANA HANNING
 * ══════════════════════════════════════════════════════════
 *
 * La ventana Hanning w(n) = 0.5 * (1 - cos(2π*n / (N-1)))
 * tiene forma de campana: pesa 0 en los extremos y 1 en
 * el centro. Al multiplicar la señal por esta ventana antes
 * de la FFT, las discontinuidades en los bordes desaparecen,
 * eliminando el "spectral leakage" que ensancha los picos.
 *
 * Efecto visual:
 *   Sin ventana:  pico en 2 Hz contamina bins de 1 y 3 Hz
 *   Con Hanning:  pico en 2 Hz queda contenido en su bin
 */
void SHM_AplicarVentanaHanning(float *señal, uint32_t N)
{
    for (uint32_t i = 0; i < N; i++) {
        float w = 0.5f * (1.0f - cosf(2.0f * PI * (float)i / (float)(N - 1)));
        señal[i] *= w;
    }
}

/* ══════════════════════════════════════════════════════════
 *  2. MÉTODO DE WELCH
 * ══════════════════════════════════════════════════════════
 *
 * Problema con una sola FFT:
 *   El espectro de una ventana corta tiene mucha varianza
 *   estadística — los picos "bailan" de medición en medición
 *   aunque la señal no cambie. Esto haría falsas alarmas en SHM.
 *
 * Solución — Welch 1967:
 *   1. Dividir la señal en K ventanas solapadas al 50%
 *   2. Aplicar Hanning a cada ventana
 *   3. Calcular FFT de cada ventana
 *   4. Promediar las magnitudes al cuadrado (densidad espectral)
 *
 * El solapamiento al 50% recupera la información "perdida"
 * por la ventana Hanning en los bordes, sin duplicar el costo
 * de captura (cada nueva ventana solo necesita N/2 muestras nuevas).
 *
 * Reducción de varianza: factor 1/K en la varianza del estimador,
 * equivale a reducir ruido en sqrt(K).
 * Con K=8: ruido reducido al 35% del original.
 *
 * Longitud de buffer necesaria con solapamiento del 50%:
 *   n = N + (K-1) * (N/2)  =  N * (1 + (K-1)/2)
 *   Con N=512, K=8: n = 512 + 7*256 = 2304 muestras
 */

uint32_t SHM_LongitudBufferWelch(void)
{
    uint32_t paso = SHM_N * (100 - SHM_SOLAPAMIENTO) / 100;  /* 256 muestras */
    return SHM_N + (SHM_WELCH_VENTANAS - 1) * paso;
}

uint32_t SHM_Welch(const float *señal_larga, uint32_t n_muestras,
                   FFT_Bin *espectro_out, float *freq_bins)
{
    uint32_t paso = SHM_N * (100 - SHM_SOLAPAMIENTO) / 100;
    uint32_t n_min = SHM_LongitudBufferWelch();

    printf("    [W] n_muestras=%lu n_min=%lu paso=%lu\r\n",
               (unsigned long)n_muestras,
               (unsigned long)n_min,
               (unsigned long)paso);

    if (n_muestras < n_min){
    	printf("    [W] ERROR: buffer insuficiente\r\n");
    	return 0;
    }

    /* Inicializar acumuladores a 0 */
    memset(_acum_mag,  0, sizeof(_acum_mag));
    memset(_acum_fase, 0, sizeof(_acum_fase));

    uint32_t ventanas_ok = 0;

    for (uint32_t v = 0; v < SHM_WELCH_VENTANAS; v++)
    {
    	uint32_t t0 = HAL_GetTick();
        uint32_t offset = v * paso;

        /* ── Copiar ventana al buffer de trabajo ── */
        for (uint32_t i = 0; i < SHM_N; i++) {
            _ventana_work[i] = señal_larga[offset + i];
        }

        /* ── Aplicar ventana Hanning ── */
        SHM_AplicarVentanaHanning(_ventana_work, SHM_N);

        /* ── Calcular FFT ── */
        int ret = FFT_Calcular(_ventana_work, _fft_bin_tmp,
                               SHM_N, SHM_FS_HZ, NULL);

        uint32_t t1 = HAL_GetTick();
                printf("    [W] ventana %lu: FFT_ret=%d  %lu ms\r\n",
                       (unsigned long)v, ret, (unsigned long)(t1 - t0));


        if (ret != 0) continue;

        /* ── Acumular magnitud y fase ──
         *
         * Para la magnitud acumulamos el CUADRADO (potencia espectral).
         * Esto es lo correcto estadísticamente en Welch: promediar
         * potencias, no amplitudes.
         * Al final sacamos raíz del promedio para volver a magnitud.
         *
         * Para la fase acumulamos directamente y promediamos.
         * Nota: promediar fases es una simplificación; para señales
         * con fase muy variable entre ventanas, la fase promediada
         * puede no ser representativa. La magnitud sí es robusta.
         */
        for (uint32_t k = 0; k < SHM_N / 2; k++) {
            float m = _fft_bin_tmp[k].magnitud;
            _acum_mag[k]  += m * m;
            _acum_fase[k] += _fft_bin_tmp[k].fase;
        }

        ventanas_ok++;
    }

    if (ventanas_ok == 0) return 0;

    /* ── Calcular promedio final ── */
    float inv_k = 1.0f / (float)ventanas_ok;
    float escala = 1.0f / (float)(SHM_N / 2);   /* normalización de magnitud */

    for (uint32_t k = 0; k < SHM_N / 2; k++) {
        /* Raíz del promedio de potencias → magnitud promediada */
        espectro_out[k].magnitud = sqrtf(_acum_mag[k] * inv_k) * escala;
        espectro_out[k].fase     = _acum_fase[k] * inv_k;

        if (freq_bins != NULL) {
            freq_bins[k] = (float)k * SHM_FS_HZ / (float)SHM_N;
        }
    }

    return ventanas_ok;
}

/* ══════════════════════════════════════════════════════════
 *  3. DETECCIÓN DE PICOS MODALES
 * ══════════════════════════════════════════════════════════
 *
 * Un "pico modal" es una frecuencia natural de la estructura.
 * Aparece en el espectro como un máximo local cuya magnitud
 * supera el nivel de ruido de fondo.
 *
 * Criterios para declarar un bin como pico:
 *   a) magnitud[k] > umbral              (supera el ruido)
 *   b) magnitud[k] > magnitud[k-1]       (mayor que vecino izquierdo)
 *   c) magnitud[k] > magnitud[k+1]       (mayor que vecino derecho)
 *
 * Se ignora el bin 0 (componente DC, frecuencia 0 Hz) porque
 * representa el valor medio de la aceleración, no una vibración.
 *
 * Los picos se ordenan de mayor a menor magnitud para que
 * picos_out[0] sea siempre el modo dominante.
 */
void SHM_DetectarPicos(const FFT_Bin *espectro, uint32_t N_2,
                       float fs_hz, float umbral,
                       SHM_Pico *picos_out, uint32_t *n_picos)
{
    *n_picos = 0;

    for (uint32_t k = 1; k < N_2 - 1 && *n_picos < SHM_MAX_PICOS; k++)
    {
        float m = espectro[k].magnitud;

        if (m > umbral &&
            m > espectro[k - 1].magnitud &&
            m > espectro[k + 1].magnitud)
        {
            picos_out[*n_picos].frecuencia_hz = (float)k * fs_hz / (float)(2 * N_2);
            picos_out[*n_picos].magnitud      = m;
            picos_out[*n_picos].fase          = espectro[k].fase;
            (*n_picos)++;
        }
    }

    /* ── FIX: guardar en variable con signo antes de ordenar ── */
    if (*n_picos < 2) return;  /* 0 o 1 picos → nada que ordenar */

    uint32_t np = *n_picos;
    for (uint32_t i = 0; i < np - 1; i++) {
        for (uint32_t j = 0; j < np - i - 1; j++) {
            if (picos_out[j].magnitud < picos_out[j + 1].magnitud) {
                SHM_Pico tmp     = picos_out[j];
                picos_out[j]     = picos_out[j + 1];
                picos_out[j + 1] = tmp;
            }
        }
    }
    printf("  [P] n_picos=%lu, mag_max=%.6f\r\n",
           (unsigned long)*n_picos,
           espectro[1].magnitud);  /* para ver el orden de magnitud real */
}

/* ══════════════════════════════════════════════════════════
 *  4. LÍNEA BASE
 * ══════════════════════════════════════════════════════════
 *
 * La línea base es el "estado sano" de referencia.
 * Se graba una sola vez cuando se instala el sistema y el
 * edificio está en condiciones normales.
 *
 * En una implementación completa, esta estructura debería
 * guardarse en memoria no volátil (Flash interna del H723
 * o EEPROM externa) para que sobreviva reinicios del micro.
 */
void SHM_GrabarLineaBase(SHM_LineaBase *lb,
                         const SHM_Pico *picos, uint32_t n_picos)
{
    lb->n_picos = (n_picos > SHM_MAX_PICOS) ? SHM_MAX_PICOS : n_picos;

    for (uint32_t i = 0; i < lb->n_picos; i++) {
        lb->picos[i] = picos[i];
    }

    lb->valida = 1;
}

/* ══════════════════════════════════════════════════════════
 *  5. DAMAGE INDEX
 * ══════════════════════════════════════════════════════════
 *
 * El Damage Index (DI) cuantifica cuánto se desplazaron
 * las frecuencias modales respecto a la línea base.
 *
 * Fundamento físico:
 *   La frecuencia natural de un modo es f = (1/2π) * sqrt(K/M)
 *   donde K es la rigidez y M la masa.
 *   Si hay daño → K disminuye → f disminuye.
 *   Un desplazamiento de frecuencia del 5% indica daño significativo
 *   en la mayoría de los estándares de SHM (ISO 13822, FEMA P-58).
 *
 * Algoritmo:
 *   Para cada pico modal de la línea base (f_ref):
 *     1. Buscar el pico más cercano en la medición actual (f_act)
 *     2. Calcular DI_modo = |f_act - f_ref| / f_ref
 *   DI global = max(DI_modo) sobre todos los modos
 *
 * El DI global representa el modo más afectado.
 */
float SHM_CalcularDamageIndex(const SHM_LineaBase *lb,
                              const SHM_Pico *picos, uint32_t n_picos)
{
    if (!lb->valida || lb->n_picos == 0 || n_picos == 0) return 0.0f;

    float di_global = 0.0f;

    for (uint32_t r = 0; r < lb->n_picos; r++)
    {
        float f_ref = lb->picos[r].frecuencia_hz;

        /* Buscar el pico actual más cercano a f_ref */
        float dist_min = 1e10f;
        float f_cercano = f_ref;

        for (uint32_t a = 0; a < n_picos; a++) {
            float dist = fabsf(picos[a].frecuencia_hz - f_ref);
            if (dist < dist_min) {
                dist_min   = dist;
                f_cercano  = picos[a].frecuencia_hz;
            }
        }

        /* DI de este modo */
        float di_modo = fabsf(f_cercano - f_ref) / f_ref;

        /* Actualizar el máximo global */
        if (di_modo > di_global) {
            di_global = di_modo;
        }
    }

    return di_global;
}

/* ══════════════════════════════════════════════════════════
 *  6. PIPELINE COMPLETO
 * ══════════════════════════════════════════════════════════ */
static FFT_Bin  _espectro_welch[SHM_N / 2];
static float    _freq_bins[SHM_N / 2];

int SHM_Procesar(float *señal_larga, uint32_t n_muestras,
                 const SHM_LineaBase *lb, SHM_Resultado *resultado)
{
	uint32_t t0, t1;

	    /* ── Paso 1: Welch ── */
	    t0 = HAL_GetTick();
	    printf("  [W] Iniciando Welch...\r\n");

	    uint32_t ventanas = SHM_Welch(señal_larga, n_muestras,
	                                  shm_espectro_publico, _freq_bins);

	    t1 = HAL_GetTick();
	    printf("  [W] Welch termino: %lu ventanas en %lu ms\r\n",
	           (unsigned long)ventanas, (unsigned long)(t1 - t0));

	    if (ventanas == 0) {
	        printf("  [W] ERROR: ventanas=0\r\n");
	        return -1;
	    }

	    /* ── Paso 2: Detectar picos ── */
	    t0 = HAL_GetTick();
	    printf("  [P] Detectando picos...\r\n");

	    SHM_DetectarPicos(shm_espectro_publico, SHM_N / 2, SHM_FS_HZ,
	                      SHM_UMBRAL_PICO,
	                      resultado->picos, &resultado->n_picos);

	    t1 = HAL_GetTick();
	    printf("  [P] Picos: %lu encontrados en %lu ms\r\n",
	           (unsigned long)resultado->n_picos, (unsigned long)(t1 - t0));

	    /* ── Paso 3: Damage Index ── */
	    resultado->damage_index = 0.0f;

	    /* ── Paso 4: Estado ── */
	    resultado->estado = SHM_ESTADO_SANO;

	    return 0;
}
