import serial
import serial.tools.list_ports
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation
import tkinter as tk
from tkinter import ttk, messagebox


class SismografoApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SHM - Sismógrafo Popayán")
        self.root.geometry("1000x800")

        self.ser = None
        self.running = False

        # Datos del espectro (se actualizan cuando llega ESPECTRO_END)
        self.freqs = []
        self.mags = []
        self.historial_energia = []

        # Estado del parser de protocolo
        self.en_espectro = False  # True entre ESPECTRO_START y ESPECTRO_END
        self.freqs_tmp = []
        self.mags_tmp = []

        # Buffer de líneas incompletas (puede llegar medio línea en un read)
        self.line_buffer = ""

        self.setup_ui()

    def setup_ui(self):
        control_frame = ttk.LabelFrame(
            self.root, text="Configuración de Conexión", padding=10
        )
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        ttk.Label(control_frame, text="Puerto COM:").pack(side=tk.LEFT, padx=5)
        self.combo_puertos = ttk.Combobox(
            control_frame, width=12, postcommand=self.actualizar_puertos
        )
        self.combo_puertos.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            control_frame, text="Actualizar", command=self.actualizar_puertos
        ).pack(side=tk.LEFT, padx=2)

        ttk.Label(control_frame, text="Baudrate:").pack(side=tk.LEFT, padx=5)
        self.baudrate = ttk.Combobox(
            control_frame, values=[9600, 115200, 230400, 921600], width=8
        )
        self.baudrate.set(115200)
        self.baudrate.pack(side=tk.LEFT, padx=5)

        self.btn_conectar = ttk.Button(
            control_frame, text="Conectar", command=self.toggle_conexion
        )
        self.btn_conectar.pack(side=tk.LEFT, padx=10)

        self.lbl_status = ttk.Label(
            control_frame, text="Desconectado", foreground="red"
        )
        self.lbl_status.pack(side=tk.LEFT, padx=10)

        # Label de diagnóstico — muestra la última línea recibida
        self.lbl_debug = ttk.Label(
            control_frame, text="—", foreground="gray", font=("Courier", 9)
        )
        self.lbl_debug.pack(side=tk.RIGHT, padx=10)

        # Gráficos
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, figsize=(9, 6))
        self.fig.tight_layout(pad=4.0)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.ani = animation.FuncAnimation(
            self.fig, self.actualizar_grafico, interval=300, cache_frame_data=False
        )

    def actualizar_puertos(self):
        puertos = [p.device for p in serial.tools.list_ports.comports()]
        self.combo_puertos["values"] = puertos
        if puertos:
            self.combo_puertos.set(puertos[0])

    def toggle_conexion(self):
        if not self.running:
            puerto = self.combo_puertos.get()
            if not puerto:
                messagebox.showwarning("Advertencia", "Selecciona un puerto COM.")
                return
            try:
                self.ser = serial.Serial(puerto, int(self.baudrate.get()), timeout=0)
                self.running = True
                self.line_buffer = ""
                self.en_espectro = False
                self.btn_conectar.config(text="Desconectar")
                self.lbl_status.config(text=f"Conectado: {puerto}", foreground="green")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo abrir el puerto:\n{e}")
        else:
            self.running = False
            if self.ser and self.ser.is_open:
                self.ser.close()
            self.btn_conectar.config(text="Conectar")
            self.lbl_status.config(text="Desconectado", foreground="red")

    def leer_datos(self):
        """
        Lee TODOS los bytes disponibles en el buffer serial y los procesa
        línea por línea. Nunca bloquea — si no hay datos, retorna de inmediato.

        El truco clave es acumular en self.line_buffer los bytes que llegan
        parcialmente y solo procesar cuando encontramos un '\n' completo.
        """
        if not (self.ser and self.ser.is_open):
            return

        try:
            # Leer todo lo disponible de una vez (no bloqueante)
            n = self.ser.in_waiting
            if n == 0:
                return

            raw = self.ser.read(n)
            self.line_buffer += raw.decode("utf-8", errors="ignore")

        except Exception:
            return

        # Procesar todas las líneas completas que haya en el buffer
        while "\n" in self.line_buffer:
            linea, self.line_buffer = self.line_buffer.split("\n", 1)
            linea = linea.strip()
            if linea:
                self.procesar_linea(linea)

    def procesar_linea(self, linea):
        """Máquina de estados del protocolo serial."""

        # Mostrar última línea recibida en el label de debug
        self.lbl_debug.config(text=linea[:60])

        if linea == "ESPECTRO_START":
            self.en_espectro = True
            self.freqs_tmp = []
            self.mags_tmp = []
            return

        if linea == "ESPECTRO_END":
            self.en_espectro = False
            # Solo actualizar si llegaron datos válidos
            if len(self.freqs_tmp) > 2:
                self.freqs = list(self.freqs_tmp)
                self.mags = list(self.mags_tmp)
            return

        if self.en_espectro:
            # Línea de dato: "0.195,0.000123"
            if "," in linea:
                try:
                    partes = linea.split(",")
                    f = float(partes[0])
                    m = float(partes[1])
                    self.freqs_tmp.append(f)
                    self.mags_tmp.append(m)
                except (ValueError, IndexError):
                    pass
            return

        # Líneas fuera del espectro (RAW, PICOS, debug) — solo mostrar en label
        # No hacer nada más para no interferir con el plot

    def actualizar_grafico(self, frame):
        """Llamado por FuncAnimation cada 300ms."""
        self.leer_datos()

        if not self.freqs:
            return

        # ── Gráfico 1: Espectro FFT ──
        self.ax1.clear()
        self.ax1.plot(self.freqs, self.mags, color="#1f77b4", linewidth=1.5)
        self.ax1.fill_between(self.freqs, self.mags, alpha=0.25, color="#1f77b4")
        self.ax1.set_title(
            "Espectro de Frecuencia (FFT)", fontsize=10, fontweight="bold"
        )
        self.ax1.set_xlabel("Frecuencia (Hz)")
        self.ax1.set_ylabel("Magnitud (g)")
        self.ax1.set_xlim([0, 20])
        self.ax1.set_ylim(bottom=0)
        self.ax1.grid(True, linestyle="--", alpha=0.5)

        # Marcar picos automáticamente
        if len(self.mags) > 2:
            max_m = max(self.mags)
            if max_m > 0:
                umbral = max_m * 0.3
                for i in range(1, len(self.mags) - 1):
                    if (
                        self.mags[i] > umbral
                        and self.mags[i] > self.mags[i - 1]
                        and self.mags[i] > self.mags[i + 1]
                    ):
                        self.ax1.annotate(
                            f"{self.freqs[i]:.2f} Hz",
                            xy=(self.freqs[i], self.mags[i]),
                            xytext=(4, 6),
                            textcoords="offset points",
                            fontsize=8,
                            color="red",
                            arrowprops=dict(arrowstyle="->", color="red", lw=0.8),
                        )

        # ── Gráfico 2: Energía acumulada ──
        energia = sum(m * m for m in self.mags) ** 0.5
        self.historial_energia.append(energia)
        if len(self.historial_energia) > 60:
            self.historial_energia.pop(0)

        self.ax2.clear()
        self.ax2.plot(self.historial_energia, color="#2ca02c", linewidth=2)
        self.ax2.set_title("Tendencia de Energía", fontsize=10, fontweight="bold")
        self.ax2.set_xlabel("Ciclos de medición")
        self.ax2.set_ylabel("Energía (g·rms)")
        self.ax2.grid(True, linestyle="--", alpha=0.5)

        self.fig.tight_layout(pad=3.0)


if __name__ == "__main__":
    root = tk.Tk()
    app = SismografoApp(root)

    def on_closing():
        if app.ser and app.ser.is_open:
            app.ser.close()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
