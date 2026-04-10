import serial
import matplotlib.pyplot as plt
import matplotlib.animation as animation

PORT = "COM3"  # Windows: COM3
BAUD = 115200

ser = serial.Serial(PORT, BAUD, timeout=2)

freqs = []
mags = []

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
fig.suptitle("SHM - Sismógrafo Popayán", fontsize=14)


def leer_espectro():
    global freqs, mags
    linea = ser.readline().decode("utf-8", errors="ignore").strip()

    if linea == "ESPECTRO_START":
        f_tmp, m_tmp = [], []
        while True:
            l = ser.readline().decode("utf-8", errors="ignore").strip()
            if l == "ESPECTRO_END":
                break
            if "," in l:
                try:
                    f, m = l.split(",")
                    f_tmp.append(float(f))
                    m_tmp.append(float(m))
                except:
                    pass
        if f_tmp:
            freqs, mags = f_tmp, m_tmp


def actualizar(frame):
    leer_espectro()
    if not freqs:
        return

    # Gráfico 1: Espectro FFT
    ax1.clear()
    ax1.plot(freqs, mags, "b-", linewidth=1)
    ax1.fill_between(freqs, mags, alpha=0.3)
    ax1.set_xlabel("Frecuencia (Hz)")
    ax1.set_ylabel("Magnitud (g)")
    ax1.set_title("Espectro de Vibración")
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, 20])

    # Marcar picos automáticamente
    if len(mags) > 2:
        max_mag = max(mags)
        umbral = max_mag * 0.3  # picos > 30% del máximo
        for i in range(1, len(mags) - 1):
            if mags[i] > umbral and mags[i] > mags[i - 1] and mags[i] > mags[i + 1]:
                ax1.annotate(
                    f"{freqs[i]:.2f} Hz",
                    xy=(freqs[i], mags[i]),
                    xytext=(freqs[i] + 0.2, mags[i] * 1.05),
                    fontsize=8,
                    color="red",
                    arrowprops=dict(arrowstyle="->", color="red"),
                )

    # Gráfico 2: Historial de magnitud total (energía)
    energia = sum(m * m for m in mags) ** 0.5
    if not hasattr(actualizar, "historial"):
        actualizar.historial = []
    actualizar.historial.append(energia)
    if len(actualizar.historial) > 50:
        actualizar.historial.pop(0)

    ax2.clear()
    ax2.plot(actualizar.historial, "g-", linewidth=2)
    ax2.set_xlabel("Ciclos de medición")
    ax2.set_ylabel("Energía total")
    ax2.set_title("Tendencia de energía de vibración")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()


ani = animation.FuncAnimation(fig, actualizar, interval=500)
plt.show()
ser.close()
