# Proyecto: Sismógrafo SHM — Estructuras Patrimoniales de Popayán

> **Autor:** Joad Julian Gaviria Chaves, Maria de los Angeles Rubio Yasno, Daniel Alejandro Paz Mera  
> **Última actualización:** Mayo 2026  
> **Plataforma objetivo:** STM32H723VGT6 + LIS3DH (I2C)  
> **Propósito:** Monitoreo de salud estructural (Structural Health Monitoring) no invasivo para patrimonio arquitectónico

---

## 1. Descripción General

El sistema captura aceleraciones triaxiales de una estructura con un acelerómetro MEMS, construye un espectro de vibración promediado mediante el método de Welch, detecta las frecuencias modales dominantes y calcula un índice de daño (Damage Index) comparando el estado actual contra una línea base de referencia. El objetivo final es detectar deterioro estructural temprano sin intervenir físicamente el inmueble.

---

## 2. Hardware

### Microcontrolador: STM32H723VGT6

| Parámetro | Valor |
|---|---|
| Familia | STM32H7 (Cortex-M7) |
| Frecuencia | ~500 MHz (PLL desde HSI 64 MHz: M=4, N=34, P=1, FRACN=3072) |
| Bus AHB | HCLK / 2 |
| Buses APB1/2/3/4 | PCLK / 2 |
| Flash latency | 3 ciclos de espera |
| MPU | Configurada: región 0 = 4 GB sin acceso (protección general) |
| UART debug | USART1 @ 115200 bps, PA9/PA10 |
| I2C sensor | I2C1, timing `0x60404E72` |

### Acelerómetro: LIS3DH (en el código referido también como LIS3DSH)

| Parámetro | Valor |
|---|---|
| Interfaz | I2C, dirección `LIS3DSH_I2C_ADDR` |
| WHO_AM_I esperado | `0x3F` |
| ODR configurado | 100 Hz (`CTRL_REG4 = 0x6F`) |
| Rango de medición | ±2 g (`CTRL_REG5 = 0x00`) |
| Ancho de banda | 800 Hz |
| Sensibilidad | 0.06 mg/LSB (~0.00006 g/LSB) |
| BDU | Activado (Block Data Update) |
| Auto-incremento registros | Activado (bit 7 del registro base en lectura multibyte) |

---

## 3. Arquitectura de Software

### 3.1 Módulos

```
main.c
 ├─ lis3dh_driver.c / .h    → Inicialización y lectura del acelerómetro
 ├─ shm.c / .h              → Pipeline SHM (Welch + detección de picos + DI)
 └─ fft_analysis.c / .h     → FFT Cooley-Tukey radix-2 DIT in-place
```

### 3.2 Flujo principal (main.c)

```
Arranque
  └── Inicialización periféricos (GPIO, USART1, I2C1)
  └── LIS3DSH_Init() → verifica WHO_AM_I y configura ODR/rango
  └── Bucle infinito
        ├── Warm-up del filtro HP: descarta las primeras WARMUP_MUESTRAS (300)
        │     muestras hasta que el transitorio del filtro se estabilice
        │
        ├── Temporización precisa cada 10 ms (INTERVALO_MUESTREO_MS)
        │     └── LIS3DSH_ReadAccel() → ax_g, ay_g, az_g
        │     └── magnitud = sqrt(ax² + ay² + az²)  [aceleración escalar]
        │     └── filtro HP: muestra_hp = α*(prev_out + raw - prev_in)
        │     └── buffer_muestras[idx++]  ← solo si warmup_ok=1
        │     └── printf "RAW,x,y,z"
        │
        └── Cuando buffer_muestras llena (BUFFER_LEN = 2304)
              └── SHM_Procesar() → espectro Welch + picos
              └── Imprime espectro 1–20 Hz (ESPECTRO_START / ESPECTRO_END)
              └── Imprime picos detectados
              └── Si calibrado=0: acumular mediciones; al llegar a 5
                    → SHM_GrabarLineaBase() → calibrado=1 → continue
                    (el ciclo donde se graba la LB no calcula DI)
              └── Si calibrado=1: imprime DI y estado
```

---

## 4. Módulo `fft_analysis.c`

### Algoritmo: Cooley-Tukey Radix-2 DIT in-place

**Requisitos de entrada:**
- N debe ser potencia de 2
- N ≤ `FFT_MAX_SIZE`

**Etapas internas:**
1. **Bit-reversal permutation** — reordena las muestras para el algoritmo DIT
2. **Butterfly stages** — `log2(N)` etapas, cada una con grupos de longitud `len = 2, 4, …, N`
3. Cada mariposa: `A' = A + W·B`, `B' = A - W·B`, donde `W = e^(-j·2π·k/len)`

**Complejidad:** O(N log N)

**Salida:** Solo los `N/2` bins positivos (señal real → simetría conjugada).

| Función | Descripción |
|---|---|
| `FFT_Calcular(entrada, salida, N, fs_hz, freq_bins)` | Ejecuta la FFT y devuelve magnitud + fase por bin |
| `FFT_FrecuenciaDominante(salida, N_2, fs_hz)` | Retorna la frecuencia del bin de mayor magnitud (ignora DC) |
| `FFT_Normalizar(salida, N_2)` | Escala magnitudes dividiendo por `N/2` |

**Buffer estático interno:** `_buf[FFT_MAX_SIZE]` (no reentrante, seguro sin RTOS)

---

## 5. Módulo `shm.c`

### 5.1 Ventana Hanning

```
w(n) = 0.5 · (1 − cos(2π·n / (N−1)))
```

Elimina el *spectral leakage* causado por discontinuidades en los bordes del bloque de muestras. Pesa 0 en extremos, 1 en el centro.

---

### 5.2 Método de Welch

**Parámetros configurados en `shm.h`:**

| Constante | Valor | Descripción |
|---|---|---|
| `SHM_N` | 512 | Tamaño de ventana FFT |
| `SHM_FS_HZ` | 100.0 | Frecuencia de muestreo (Hz) |
| `SHM_WELCH_VENTANAS` (K) | 8 | Número de ventanas promediadas |
| `SHM_SOLAPAMIENTO` | 50 % | Overlap entre ventanas |
| `SHM_UMBRAL_PICO` | ver sección 8 | Umbral mínimo de magnitud para picos |
| `SHM_MAX_PICOS` | (definido en .h) | Máximo número de picos a registrar |

**Longitud de buffer requerida:**

```
n = N + (K−1) · (N/2) = 512 + 7·256 = 2304 muestras
```

Eso equivale a **23.04 segundos** de señal a 100 Hz.

**Importante:** el buffer cubre 23 s. Si ocurre un evento de excitación y el buffer se llena poco después, las primeras ventanas de Welch todavía contendrán muestras del evento, inflando la magnitud espectral y el DI del siguiente ciclo aunque la estructura ya esté en reposo. Ver sección 9.

**Acumulación de potencia:**
- Se acumula `magnitud²` (potencia espectral) para cada bin
- Al finalizar: `magnitud_promedio = sqrt(acum / K) · escala`
- Escala de normalización: `1 / (N/2)`

**Reducción de varianza:** factor `1/K` → ruido reducido al `1/√K ≈ 35%` con K=8.

---

### 5.3 Detección de Picos Modales

Criterio para declarar un bin `k` como pico:
- `magnitud[k] > umbral`
- `magnitud[k] > magnitud[k−1]` (máximo local)
- `magnitud[k] > magnitud[k+1]`
- Se ignora el bin 0 (componente DC)

Los picos detectados se ordenan por magnitud descendente (burbuja).

**Nota crítica sobre el umbral:** si `SHM_UMBRAL_PICO` es demasiado bajo, los picos detectados durante la calibración serán bins de ruido blanco, no frecuencias naturales de la estructura. Esto produce falsas alarmas en reposo porque el DI de magnitud compara bins de ruido (que cambian de ciclo en ciclo) contra la línea base. Ver sección 8 y 9.

---

### 5.4 Línea Base y Damage Index

**Línea base (`SHM_LineaBase`):**
- Se graba **automáticamente después de 5 mediciones** cuando el edificio está en condiciones normales
- El flujo de calibración es disparado desde `main.c` mediante un contador estático
- Contiene los picos modales de referencia (`frecuencia_hz`, `magnitud`, `fase`)
- Flag `valida` indica si ya fue grabada (se establece a 1 tras `SHM_GrabarLineaBase()`)
- El ciclo en que se graba la LB termina con `continue` — el DI no se calcula en ese mismo ciclo para evitar DI=0 espurio
- **Pendiente:** persistencia en Flash interna o EEPROM externa para sobrevivir reinicios

**Damage Index (DI):**

```
DI_freq = |f_actual − f_referencia| / f_referencia
DI_mag  = |m_actual − m_referencia| / m_referencia
DI_modo = max(DI_freq, DI_mag)
DI_global = max(DI_modo) sobre todos los modos de la línea base
```

Fundamento físico: frecuencia natural ∝ `√(K/M)`. Daño → rigidez K disminuye → frecuencia baja.

Umbrales de diagnóstico:
- **DI < 2%** → `SHM_ESTADO_SANO` (estructura íntegra)
- **2% ≤ DI < 5%** → `SHM_ESTADO_VIGILAR` (cambio moderado, monitorear)
- **DI ≥ 5%** → `SHM_ESTADO_ALERTA` (posible daño estructural, requiere inspección)

**Limitación conocida del DI de magnitud:** el término `DI_mag` es sensible a variaciones de amplitud entre ciclos (ruido, cambios en el acoplamiento del sensor). En estructuras con vibración ambiental baja, el `DI_freq` es más confiable que `DI_mag`. Considerar ponderar o desactivar `DI_mag` según el entorno de despliegue.

---

### 5.5 Pipeline `SHM_Procesar()`

```
SHM_Procesar(señal_larga, n_muestras, lb, resultado)
  ├── SHM_Welch()            → espectro promediado en shm_espectro_publico[]
  ├── SHM_DetectarPicos()    → resultado->picos[], resultado->n_picos
  ├── SHM_CalcularDamageIndex() → resultado->damage_index  [solo si lb != NULL]
  └── resultado->estado      → SHM_ESTADO_SANO / VIGILAR / ALERTA (basado en DI)
```

**Variable global pública:** `FFT_Bin shm_espectro_publico[SHM_N/2]` — contiene el espectro Welch y es leída directamente en `main.c` para la impresión por UART.

**Parámetro `lb` (línea base):**
- Si `lb == NULL` (antes de calibración): `damage_index = 0.0` (sin referencia)
- Si `lb != NULL` (después de calibración): `damage_index` se calcula comparando frecuencias modales contra la línea base

---

## 6. Filtro Pasa-Alto (HP) en `main.c`

Se aplica un filtro HP de primer orden a la magnitud escalar antes de almacenar en el buffer, para eliminar la componente de gravedad residual:

```c
muestra_hp = HP_ALPHA * (hp_prev_out + muestra_raw - hp_prev_in);
hp_prev_in  = muestra_raw;
hp_prev_out = muestra_hp;
```

| Parámetro | Valor | Descripción |
|---|---|---|
| `HP_ALPHA` | 0.9950 | α = fc/(fc+fs) con fc≈0.5 Hz, fs=100 Hz |
| `WARMUP_MUESTRAS` | 300 | Muestras descartadas al inicio para estabilizar el transitorio |
| `warmup_ok` | flag uint8_t | Se activa cuando `warmup_cnt >= WARMUP_MUESTRAS` |

**Comportamiento observado:** con 300 muestras de warm-up, el pico DC residual bajó de 0.035 g (sin warm-up) a 0.0078 g. Con 500 muestras debería desaparecer completamente. Ver sección 9.

---

## 7. Comunicación y Diagnóstico (UART)

Toda la salida es por `USART1` a 115200 bps vía `printf` redirigido a `HAL_UART_Transmit`.

| Mensaje | Significado |
|---|---|
| `RAW,ax,ay,az` | Muestra cruda en tiempo real (cada 10 ms) |
| `[HP] Filtro en regimen. Iniciando captura.` | Warm-up completado, inicio de captura válida |
| `T1: Iniciando SHM_Procesar...` | Inicio de procesamiento |
| `T2: SHM_Procesar termino en X ms` | Tiempo total del pipeline |
| `ESPECTRO_START` / `ESPECTRO_END` | Delimitadores del volcado de espectro |
| `freq,magnitud` | Un bin del espectro (rango 1–20 Hz) |
| `PICOS:N` | Cantidad de picos encontrados |
| `PICO i: X.XX Hz mag=Y` | Pico modal individual |
| `DI:X.XXXXX Estado:Y` | Damage Index y estado de salud (0=SANO, 1=VIGILAR, 2=ALERTA) |
| `[CALIB] Medicion N de 5...` | Contador de calibración automática |
| `[CALIB] Sin picos detectados, reintentando...` | Calibración rechazada por falta de picos |
| `[CALIB] Linea base grabada. Sistema CALIBRADO.` | Calibración completada |
| `*** ALERTA ESTRUCTURAL ***` | Estado ALERTA detectado (impreso junto al DI) |
| `ERROR: SHM_Procesar fallo` | Fallo del pipeline |

**Nota:** el mensaje `DI:X.XXXXX Estado:Y` se imprime en todos los ciclos post-calibración para mantener compatibilidad con el parser Python de la interfaz de visualización, incluyendo el ciclo inmediatamente siguiente a la grabación de la LB (donde el DI puede ser artificialmente bajo).

---

## 8. Memoria RAM estimada

| Buffer | Tamaño |
|---|---|
| `buffer_muestras[2304]` | ~9.0 KB |
| `_ventana_work[512]` | ~2.0 KB |
| `_acum_mag[256]` + `_acum_fase[256]` | ~2.0 KB |
| `_fft_bin_tmp[256]` | ~2.0 KB |
| `shm_espectro_publico[256]` | ~2.0 KB |
| `_buf[FFT_MAX_SIZE]` (fft_analysis.c) | variable |
| **Total aproximado (sin FFT_MAX_SIZE)** | **~17 KB** |

El STM32H723VGT6 tiene 564 KB de RAM → margen amplio.

---

## 9. Estado Actual del Proyecto

### Lo que funciona
- [x] Inicialización y lectura del LIS3DH por I2C
- [x] Muestreo temporizado a 100 Hz con `HAL_GetTick()`
- [x] Cálculo de magnitud escalar de la aceleración
- [x] Filtro HP de primer orden sobre la magnitud escalar
- [x] FFT Cooley-Tukey radix-2 completa con bit-reversal
- [x] Ventana Hanning implementada
- [x] Método de Welch con K=8 ventanas, 50% overlap
- [x] Detección de picos modales con ordenamiento
- [x] Línea base y cálculo de Damage Index (DI de frecuencia + magnitud)
- [x] Volcado de espectro y picos por UART (formato parseable por interfaz Python)
- [x] Instrumentación de tiempos con `HAL_GetTick()` en cada etapa
- [x] **[CORREGIDO 23/04/2026]** Eliminación de variable duplicada `_espectro_welch[]` en `shm.c`
- [x] **[CORREGIDO 23/04/2026]** Paso de `NULL` en lugar de `_freq_bins` a `SHM_Welch()`
- [x] **[CORREGIDO 23/04/2026]** Cálculo real del Damage Index con estado automático (SANO/VIGILAR/ALERTA)
- [x] **[CORREGIDO 23/04/2026]** Eliminación de doble llamada a `SHM_Procesar()` en `main.c`
- [x] **[CORREGIDO 23/04/2026]** Flujo de calibración automática: primeras 5 mediciones → `SHM_GrabarLineaBase()`
- [x] **[CORREGIDO 23/04/2026]** Reporte de DI y Estado por UART
- [x] **[CORREGIDO 14/05/2026]** Warm-up del filtro HP: se descartan 300 muestras al inicio para evitar el transitorio DC en la calibración (`WARMUP_MUESTRAS=300`, flag `warmup_ok`)
- [x] **[CORREGIDO 14/05/2026]** Separación de ciclos: el ciclo donde se graba la LB termina con `continue` para evitar DI=0 espurio por comparación de picos contra sí mismos
- [x] **[CORREGIDO 14/05/2026]** Guarda de calibración sin picos: si `n_picos==0` en la medición 5, el contador retrocede a 4 y fuerza una medición más
- [x] **[CORREGIDO 14/05/2026]** Eliminación de variable `resultado_shm` no utilizada en `main.c`
- [x] **[CORREGIDO 14/05/2026]** Impresión de `*** ALERTA ESTRUCTURAL ***` cuando `estado == SHM_ESTADO_ALERTA`

### Bugs conocidos activos

- [ ] **Contaminación del buffer Welch post-evento:** el buffer de 2304 muestras cubre 23 s. Cuando ocurre un evento de excitación y el siguiente ciclo se procesa antes de que el buffer se limpie completamente de muestras del evento, las primeras ventanas de Welch incluyen energía residual, produciendo un DI muy elevado (observado: DI=54 en el primer ciclo post-excitación con estructura en reposo). Requiere invalidar o marcar el buffer cuando se detecta una transición ALERTA→reposo.

- [ ] **Falsas alarmas por línea base de ruido:** si durante la calibración no hay modos estructurales claramente distinguibles sobre el noise floor (protoboard en escritorio, estructura sin vibración ambiental real), los "picos" detectados son bins de ruido blanco con frecuencias aleatorias. El DI en reposo post-calibración resulta en 0.15–0.44 (15%–44%), que supera el umbral de ALERTA (5%). La solución requiere o bien elevar `SHM_UMBRAL_PICO` por encima del noise floor medido, o bien validar la consistencia de los picos entre las 5 mediciones de calibración antes de grabar la LB. Ver sección 10.

- [ ] **Warm-up insuficiente:** con `WARMUP_MUESTRAS=300` el pico DC bajó de 0.035 g a 0.0078 g pero no desapareció. Con 500 muestras (5 s) debería estabilizarse completamente.

- [ ] **La línea base no persiste en memoria no volátil:** se pierde al reiniciar el microcontrolador.

### Variables y constantes clave en `main.c`

```c
#define BUFFER_LEN            2304U   // muestras para Welch
#define INTERVALO_MUESTREO_MS   10U   // 100 Hz ODR
#define HP_ALPHA             0.9950f  // filtro HP, fc≈0.5 Hz
#define WARMUP_MUESTRAS        300U   // muestras de calentamiento del filtro HP

static uint32_t warmup_cnt  = 0;     // contador de warm-up
static uint8_t  warmup_ok   = 0;     // flag: 1 cuando el filtro está en régimen
static uint8_t  calibrado   = 0;     // flag: 1 cuando la LB fue grabada
static uint32_t idx_muestra = 0;     // índice en buffer_muestras[]
static float    hp_prev_in  = 0.0f;  // estado del filtro HP
static float    hp_prev_out = 0.0f;  // estado del filtro HP
static SHM_LineaBase linea_base;     // línea base de referencia
```

---

## 10. Pendientes y Mejoras

### 10.1 Correcciones prioritarias (bugs activos)

1. **Elevar `WARMUP_MUESTRAS` a 500** en `main.c` para eliminar completamente el pico DC transitorio del filtro HP. Con 300 muestras el pico bajó ×4.5 pero sigue visible en el espectro 0 de calibración.

2. **Ajustar `SHM_UMBRAL_PICO` en `shm.h`** al doble del noise floor medido en calibración (~0.00020 g para protoboard en escritorio). El noise floor observado en sesiones de prueba es ~0.00009–0.00018 g. Un umbral demasiado bajo captura ruido blanco como picos modales, generando falsas alarmas en reposo.

3. **Validación de consistencia de picos en calibración:** antes de grabar la LB, verificar que los picos detectados aparezcan en al menos 3 de las 5 mediciones de calibración. Un pico de ruido aleatorio cambia de bin en cada ciclo; un modo estructural real permanece estable. Esto requiere almacenar los picos de cada medición durante el proceso de calibración y hacer votación antes de llamar a `SHM_GrabarLineaBase()`.

4. **Invalidación del primer buffer post-excitación:** cuando el estado pasa de `SHM_ESTADO_ALERTA` a cualquier otro estado, el primer ciclo posterior debería descartarse (o marcarse como no confiable) porque el buffer Welch aún contiene muestras del evento. Implementación sugerida: flag `post_evento` que se activa cuando `estado == ALERTA` y se consume en el ciclo siguiente.

5. **Persistencia de la línea base en Flash interna** del STM32H723 (sector dedicado, ~128 KB mínimo) para que sobreviva reinicios. Sin esto, el sistema siempre arranca sin calibración y requiere 5 ciclos (~115 s incluyendo warm-up) antes de poder reportar DI reales.

### 10.2 Mejoras funcionales

- **Aumentar `WARMUP_MUESTRAS` a 500** (corrección inmediata del pico DC)
- **Doble buffer ping-pong** para captura continua sin pérdida de muestras mientras se procesa el buffer anterior
- **Promediado temporal del DI:** filtro de media móvil sobre las últimas N mediciones para suavizar variaciones aleatorias y reducir falsas alarmas
- **Alarma LED/buzzer:** indicador físico cuando DI supera umbral de ALERTA
- **Timestamp real:** el STM32H7 tiene RTC integrado — datar cada medición
- **Modo bajo consumo:** STOP2 entre ciclos de captura para alargar vida de batería
- **Interpolación parabólica de picos:** refinar estimación de frecuencia modal más allá de Δf = fs/N = 0.195 Hz

### 10.3 Mejoras de robustez

- Verificar retorno de `HAL_I2C_Mem_Read/Write` en `lis3dh_driver.c` (actualmente se ignora)
- Watchdog (IWDG) para recuperación automática ante cuelgues
- CRC sobre el buffer antes de procesar para detectar corrupción de memoria
- Manejo de pérdida de muestras: si el procesamiento tarda más que `BUFFER_LEN × 10 ms = 23 s`, se pierden muestras → usar doble buffer (ping-pong)

---

## 11. Historial de Sesiones de Prueba

### Sesión 1 — `session_20260514_002106.csv`

**Condición:** sin warm-up (código original), 7 espectros, una sola excitación (esp 5).

| Espectro | Condición | DI | Estado | Observaciones |
|---|---|---|---|---|
| 0 | calibración | 0 | SANO | Pico DC enorme: 0.035 g @ 0.195 Hz — filtro HP sin estabilizar |
| 1–4 | calibración | 0 | SANO | Noise floor ~0.00009–0.00019 g, sin picos definidos |
| 5 | excitación | 0 | SANO | 0.00822 g @ 14.8 Hz (×54 sobre noise floor); DI=0 porque la LB se comparaba contra sí misma en el mismo ciclo — **bug principal** |
| 6 | post-excitación | 0 | SANO | DI=0 — bug del ciclo sin separar aún presente |

**Conclusión:** el DI era sistemáticamente 0 en todos los espectros. Bugs identificados: (a) LB grabada y comparada en el mismo ciclo → DI=0 por construcción, (b) filtro HP sin warm-up → pico DC de 0.035 g contaminando calibración.

---

### Sesión 2 — `session_20260514_005930.csv`

**Condición:** con warm-up de 300 muestras y separación de ciclos (fixes del 14/05/2026), 10 espectros, dos excitaciones (esp 5 y 6), luego quieto.

| Espectro | Condición | DI | Estado | Observaciones |
|---|---|---|---|---|
| 0 | calibración | 0 | SANO | Pico DC residual: 0.0078 g @ 0.195 Hz — mejoró ×4.5 pero no desaparece |
| 1–4 | calibración | 0 | SANO | Noise floor limpio: 0.000135–0.000175 g ✓ |
| 5 | 1ra excitación | 0 | SANO | LB=NULL → DI correcto en 0 ✓ (fix funcionó) |
| 6 | 2da excitación | 7.41 | ALERTA | 0.024 g @ 19.3 Hz, 6365× energía — detección correcta ✓ |
| 7 | quieto post-evento | 54.31 | ALERTA | **Falsa alarma** — buffer Welch contiene residuo de excitación |
| 8 | quieto | 0.44 | ALERTA | **Falsa alarma** — energía = noise floor, DI inflado por LB de ruido |
| 9 | quieto | 0.15 | ALERTA | **Falsa alarma** — misma causa |

**Conclusión:** los dos bugs principales de la sesión 1 quedaron resueltos. El sistema detectó la excitación correctamente (esp 6, DI=7.41). Quedan tres problemas nuevos: (a) buffer Welch contaminado en el ciclo inmediato post-evento, (b) falsas alarmas por LB grabada con ruido blanco en lugar de modos reales, (c) warm-up insuficiente para eliminar el pico DC completamente.

---

## 12. Diseño del Dispositivo Físico para Patrimonio

### 12.1 Restricciones del contexto patrimonial

El Consejo de Monumentos Nacionales (y equivalentes locales/UNESCO) exige que cualquier instrumentación sobre patrimonio sea:
- **Reversible:** no puede dejar huella permanente
- **No invasiva:** sin perforaciones, anclajes con tornillos, adhesivos permanentes
- **Estéticamente discreta:** no debe alterar la percepción visual del inmueble
- **Seguimiento documentado:** todo lo instalado debe poder retirarse sin daño

### 12.2 Propuesta de acoplamiento

**Principio:** adhesión temporal por vacío + peso propio + banda de sujeción textil

```
┌────────────────────────────────────────────┐
│           CARCASA EXTERIOR                 │
│   Material: PLA/PETG impreso o ABS         │
│   Color: tono piedra / neutro              │
│   Dimensiones objetivo: ~80×60×35 mm       │
│                                            │
│   ┌──────────────┐    ┌─────────────────┐  │
│   │  PCB + MCU   │    │   Batería LiPo  │  │
│   │  LIS3DH      │    │   3.7V 2000mAh  │  │
│   └──────────────┘    └─────────────────┘  │
│                                            │
│  ┌──────────────────────────────────────┐  │
│  │  BASE DE ACOPLAMIENTO                │  │
│  │  Goma EPDM o silicona 3mm            │  │
│  │  → amortigua vibraciones propias     │  │
│  │  → confiere agarre por fricción      │  │
└──┴──────────────────────────────────────┴──┘
```

**Método de fijación (sin intervenir el muro):**

| Opción | Descripción | Invasividad |
|---|---|---|
| **Ventosa industrial** | Copa de vacío 50–80 mm, ≥15 kg de carga, reversible 100% | Ninguna |
| **Cinta VHB removible 3M** | Doble cara espumada, removible con calor/solvente sin dañar la piedra | Mínima (deja residuo removible) |
| **Banda textil ajustable** | Correa que rodea un pilar/columna sin contacto adhesivo | Ninguna |
| **Imán de neodimio** | Solo si hay elementos metálicos embebidos (refuerzo, ventanas) | Ninguna |

**Recomendación para muros de tapia/mampostería (típicos en Popayán):** ventosa industrial como anclaje primario + banda textil de seguridad secundaria. La ventosa no deja marca en piedra, ladrillo ni estuco.

### 12.3 Orientación del sensor

El LIS3DH debe montarse con uno de sus ejes principales alineado con la dirección de análisis preferente. Para una pared:

- **Eje Z perpendicular al muro** → capta vibraciones de flexión (las más indicativas de daño)
- **Eje X vertical** → capta ondas sísmicas verticales
- **Eje Y horizontal** → capta vibraciones en el plano del muro

El uso de la magnitud escalar `sqrt(ax²+ay²+az²)` que ya implementa el código es robusto a la orientación, aunque sacrifica información direccional. Para una versión avanzada, procesar los tres ejes por separado daría más información modal.

### 12.4 Conectividad sugerida

| Opción | Rango | Consumo | Uso recomendado |
|---|---|---|---|
| **BLE 5.0** (interno H7 o módulo) | ~50 m | Muy bajo | Lectura de picos y DI en tiempo real desde smartphone |
| **LoRa (SX1276)** | ~2 km urbano | Bajo | Envío periódico de resumen a gateway central del edificio |
| **USB-Serial** | Local | Alto | Volcado de datos en inspecciones técnicas presenciales |

### 12.5 Autonomía estimada

Con batería LiPo 3.7V 2000 mAh y consumo promedio del H723 en medición activa (~50 mA):

```
Autonomía ≈ 2000 mAh / 50 mA ≈ 40 horas continuas
```

Con modo sleep entre ciclos Welch (cada 23 s activo, el resto dormido):

```
Duty cycle ≈ 23s / (23s + tiempo_sleep)
→ Con sleep de 5 min: duty ≈ 7% → ~570 horas (~24 días)
```

---

## 13. Hoja de Ruta

| Prioridad | Tarea | Esfuerzo estimado | Estado |
|---|---|---|---|
| 🔴 Alta | ✅ Corregir bugs (`_freq_bins` duplicado, doble llamada SHM, DI=0) | 1–2 h | **COMPLETADO 23/04** |
| 🔴 Alta | ✅ Flujo de calibración automática en `main.c` | 2–3 h | **COMPLETADO 23/04** |
| 🔴 Alta | ✅ Warm-up del filtro HP (WARMUP_MUESTRAS=300) | 1 h | **COMPLETADO 14/05** |
| 🔴 Alta | ✅ Separar ciclo de grabación de LB del ciclo de cálculo de DI | 1 h | **COMPLETADO 14/05** |
| 🔴 Alta | ✅ Guarda de calibración sin picos detectados | 0.5 h | **COMPLETADO 14/05** |
| 🔴 Alta | Aumentar `WARMUP_MUESTRAS` a 500 (pico DC residual aún presente) | 0.1 h | Pendiente |
| 🔴 Alta | Ajustar `SHM_UMBRAL_PICO` al doble del noise floor medido (~0.00020 g) | 0.5 h | Pendiente |
| 🔴 Alta | Validación de consistencia de picos entre mediciones de calibración | 3–4 h | Pendiente |
| 🟡 Media | Invalidar primer buffer post-excitación (flag `post_evento`) | 2 h | Pendiente |
| 🟡 Media | Persistencia de línea base en Flash interna del H723 | 4–6 h | Pendiente |
| 🟡 Media | Doble buffer ping-pong para captura continua | 3–4 h | Pendiente |
| 🟡 Media | Diseño e impresión 3D de la carcasa con base de ventosa | 4–8 h | Pendiente |
| 🟡 Media | Indicador LED/salida visual del estado (SANO/VIGILAR/ALERTA) | 2–3 h | Pendiente |
| 🟢 Baja | Agregar RTC y timestamp a cada medición | 2 h | Pendiente |
| 🟢 Baja | Módulo BLE/LoRa para telemetría inalámbrica | 8–12 h | Pendiente |
| 🟢 Baja | Procesamiento por eje (x, y, z independientes) | 4 h | Pendiente |
| 🟢 Baja | Herramienta Python/PC para visualización en tiempo real del espectro | 4–6 h | Pendiente |

---

## 14. Parámetros Clave de Referencia Rápida

```c
#define BUFFER_LEN              2304    // muestras necesarias para Welch
#define INTERVALO_MUESTREO_MS     10    // 100 Hz ODR
#define WARMUP_MUESTRAS          300    // warm-up del filtro HP (aumentar a 500)
#define HP_ALPHA              0.9950f   // filtro HP primer orden, fc≈0.5 Hz
#define SHM_N                    512    // tamaño ventana FFT
#define SHM_WELCH_VENTANAS         8    // K ventanas
#define SHM_SOLAPAMIENTO          50    // % overlap
#define SHM_UMBRAL_PICO      0.00020f   // ajustar al 2× noise floor medido
// Resolución espectral: Δf = 100/512 ≈ 0.195 Hz
// Rango útil: 0 – 50 Hz (Nyquist)
// Rango de interés sísmico/estructural: 1 – 20 Hz
// Noise floor observado en protoboard: ~0.00009–0.00018 g
// Tiempo de captura por ciclo: 2304 × 10 ms = 23.04 s
// Tiempo de warm-up: 300 × 10 ms = 3 s (aumentar a 500 → 5 s)
```
