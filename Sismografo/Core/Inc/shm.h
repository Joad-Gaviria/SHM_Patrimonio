/*
 * smh.h
 *
 *  Created on: 8/04/2026
 *      Author: joadj
 */

#ifndef INC_SHM_H_
#define INC_SHM_H_


/* =========================================================
 *  shm.h  —  Structural Health Monitoring
 *
 *  Módulos incluidos:
 *    1. Ventana Hanning
 *    2. Promediado de espectros (método de Welch)
 *    3. Detección de picos modales
 *    4. Línea base y Damage Index
 *
 *  Dependencia: fft.h / fft.c (debe estar en el proyecto)
 * ========================================================= */

#include "fft_analysis.h"
#include <stdint.h>


/* ── Configuración del sistema ─────────────────────────── */

/* Frecuencia de muestreo del acelerómetro en Hz.
 * Para edificios: las frecuencias modales suelen estar
 * entre 0.5 Hz y 20 Hz, así que 200 Hz es más que suficiente.
 * Ajusta según tu sensor. */
#define SHM_FS_HZ   100.0f

/* Tamaño de cada ventana FFT (potencia de 2).
 * Con fs=200 Hz y N=512: resolución = 200/512 ≈ 0.39 Hz/bin
 * Suficiente para distinguir modos estructurales. */
#define SHM_N               512

/* Número de ventanas a promediar (método Welch).
 * Más ventanas = menos ruido, pero más tiempo de captura.
 * 8 ventanas reduce el ruido a ~35% del original. */
#define SHM_WELCH_VENTANAS  8

/* Solapamiento entre ventanas consecutivas (porcentaje).
 * 50% es el estándar — equilibrio entre resolución y velocidad. */
#define SHM_SOLAPAMIENTO    50

/* Número máximo de picos modales a detectar */
#define SHM_MAX_PICOS       10

/* Umbral mínimo de magnitud para considerar un pico válido.
 * Ajustar según el nivel de ruido de tu instalación. */
#define SHM_UMBRAL_PICO     0.01f

/* Umbrales del Damage Index
 *  < DI_NORMAL  : estructura sana
 *  < DI_ALERTA  : cambio moderado, vigilar
 *  >= DI_ALERTA : posible daño estructural */
#define SHM_DI_NORMAL       0.02f   /* 2%  */
#define SHM_DI_ALERTA       0.05f   /* 5%  */

/* ── Estructuras de datos ──────────────────────────────── */

/* Estado de salud calculado en una medición */
typedef enum {
    SHM_ESTADO_SANO    = 0,
    SHM_ESTADO_VIGILAR = 1,
    SHM_ESTADO_ALERTA  = 2
} SHM_Estado;

/* Información de un pico modal detectado */
typedef struct {
    float    frecuencia_hz;   /* Frecuencia del pico                    */
    float    magnitud;        /* Amplitud del pico (normalizada)        */
    float    fase;            /* Fase en radianes                       */
} SHM_Pico;

/* Línea base — snapshot del edificio cuando está sano */
typedef struct {
    SHM_Pico picos[SHM_MAX_PICOS];   /* Picos modales de referencia    */
    uint32_t n_picos;                /* Cuántos picos se encontraron    */
    uint8_t  valida;                 /* 1 si la línea base fue grabada  */
} SHM_LineaBase;

/* Resultado completo de una medición SHM */
typedef struct {
    SHM_Pico   picos[SHM_MAX_PICOS]; /* Picos detectados en esta medición */
    uint32_t   n_picos;              /* Cantidad de picos detectados       */
    float      damage_index;         /* DI global (peor modo detectado)    */
    SHM_Estado estado;               /* Diagnóstico de salud               */
} SHM_Resultado;

/* ── API ───────────────────────────────────────────────── */

/**
 * @brief  Aplica ventana Hanning in-place sobre una señal.
 *
 *   La ventana reduce el spectral leakage multiplicando cada
 *   muestra por un peso que va de 0 en los extremos a 1 en
 *   el centro. Imprescindible antes de la FFT en SHM.
 *
 * @param  señal  Buffer de N muestras float (se modifica).
 * @param  N      Número de muestras.
 */
void SHM_AplicarVentanaHanning(float *señal, uint32_t N);

/**
 * @brief  Calcula el espectro promediado por el método de Welch.
 *
 *   Divide la señal en ventanas solapadas, calcula la FFT de
 *   cada una, y promedia las magnitudes. Esto reduce el ruido
 *   estadístico sin perder resolución en frecuencia.
 *
 * @param  señal_larga   Buffer con todas las muestras capturadas.
 *                       Longitud mínima: SHM_N * (1 + (SHM_WELCH_VENTANAS-1)
 *                       * (1 - SHM_SOLAPAMIENTO/100))
 *                       Usar SHM_LongitudBufferWelch() para calcularla.
 * @param  n_muestras    Longitud real del buffer señal_larga.
 * @param  espectro_out  Buffer de salida con SHM_N/2 bins promediados.
 * @param  freq_bins     Array de SHM_N/2 floats con la frecuencia de
 *                       cada bin en Hz. Puede ser NULL.
 * @return Número de ventanas efectivamente promediadas (≥1), o 0 si error.
 */
uint32_t SHM_Welch(const float *señal_larga, uint32_t n_muestras,
                   FFT_Bin *espectro_out, float *freq_bins);

/**
 * @brief  Calcula la longitud mínima del buffer para SHM_Welch().
 * @return Número de muestras necesarias.
 */
uint32_t SHM_LongitudBufferWelch(void);

/**
 * @brief  Detecta picos modales en un espectro promediado.
 *
 *   Un bin es pico si su magnitud supera el umbral Y es mayor
 *   que sus dos vecinos inmediatos (máximo local).
 *
 * @param  espectro   Espectro de entrada (salida de SHM_Welch).
 * @param  N_2        Número de bins = SHM_N / 2.
 * @param  fs_hz      Frecuencia de muestreo.
 * @param  umbral     Magnitud mínima para considerar un pico.
 * @param  picos_out  Array de salida con los picos encontrados.
 * @param  n_picos    Número de picos encontrados (salida).
 */
void SHM_DetectarPicos(const FFT_Bin *espectro, uint32_t N_2,
                       float fs_hz, float umbral,
                       SHM_Pico *picos_out, uint32_t *n_picos);

/**
 * @brief  Graba la línea base (estado sano del edificio).
 *
 *   Debe llamarse una vez cuando el edificio está en buen estado.
 *   Guarda los picos modales actuales como referencia.
 *
 * @param  lb       Puntero a la estructura de línea base a llenar.
 * @param  picos    Picos detectados en la medición actual.
 * @param  n_picos  Número de picos.
 */
void SHM_GrabarLineaBase(SHM_LineaBase *lb,
                         const SHM_Pico *picos, uint32_t n_picos);

/**
 * @brief  Calcula el Damage Index comparando medición vs línea base.
 *
 *   Para cada pico modal de la línea base, busca el pico más
 *   cercano en la medición actual y calcula el desplazamiento
 *   relativo de frecuencia: DI = |f_actual - f_ref| / f_ref.
 *   Retorna el DI máximo entre todos los modos.
 *
 * @param  lb       Línea base de referencia.
 * @param  picos    Picos de la medición actual.
 * @param  n_picos  Número de picos actuales.
 * @return DI en el rango [0, 1]. 0 = sin cambio. >0.05 = alerta.
 */
float SHM_CalcularDamageIndex(const SHM_LineaBase *lb,
                              const SHM_Pico *picos, uint32_t n_picos);

/**
 * @brief  Pipeline completo de SHM en una sola llamada.
 *
 *   Ejecuta en orden: Welch → DetectarPicos → DamageIndex → Estado.
 *
 * @param  señal_larga  Buffer de muestras (ver SHM_LongitudBufferWelch()).
 * @param  n_muestras   Longitud del buffer.
 * @param  lb           Línea base de referencia (puede ser NULL si aún
 *                      no se grabó; en ese caso DI = 0).
 * @param  resultado    Estructura de salida con el diagnóstico completo.
 * @return 0 si OK, -1 si error.
 */
int SHM_Procesar(float *señal_larga, uint32_t n_muestras,
                 const SHM_LineaBase *lb, SHM_Resultado *resultado);

/* Buffer público del espectro Welch — disponible después de SHM_Procesar() */
extern FFT_Bin shm_espectro_publico[SHM_N / 2];

#endif /* INC_SHM_H_ */
