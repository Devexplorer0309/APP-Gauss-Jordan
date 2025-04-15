import tkinter as tk
from tkinter import ttk, messagebox
import re
import numpy as np

# Función para hacer los pivotes (forma escalonada reducida)
def hacer_pivotes_completo(matriz):
    matriz = np.array(matriz, dtype=float)
    num_filas, num_cols = matriz.shape
    fila_actual = 0
    columna_actual = 0

    while fila_actual < num_filas and columna_actual < num_cols - 1:
        # Encuentra la primera fila con un valor no nulo en la columna actual (o debajo)
        fila_pivote = -1
        for i in range(fila_actual, num_filas):
            if matriz[i][columna_actual] != 0:
                fila_pivote = i
                break

        if fila_pivote != -1:
            # Intercambia la fila actual con la fila del pivote si es necesario
            if fila_pivote != fila_actual:
                matriz[[fila_actual, fila_pivote]] = matriz[[fila_pivote, fila_actual]]

            # Haz que el pivote sea 1
            pivote = matriz[fila_actual][columna_actual]
            matriz[fila_actual] = matriz[fila_actual] / pivote

            # Haz ceros arriba y abajo del pivote
            for i in range(num_filas):
                if i != fila_actual:
                    factor = matriz[i][columna_actual]
                    matriz[i] = matriz[i] - factor * matriz[fila_actual]

            fila_actual += 1
        columna_actual += 1

    return matriz

# Función para convertir una ecuación a matriz aumentada
def matriz_aum(ec):
    ec = ec.strip().replace(" ", "")
    izq, der = ec.split('=')
    izq = izq.replace('-', '+-')  # Solo aplicar a la parte izquierda

    terminos = re.findall(r'[+-]?\d*\.?\d*[a-zA-Z_]\w*', izq)
    coeficientes = []
    variables = {}

    for termino in terminos:
        numero = re.match(r'[+-]?\d*\.?\d*', termino).group()
        if numero in ('', '+'):
            coeficiente = 1.0
        elif numero == '-':
            coeficiente = -1.0
        else:
            coeficiente = float(numero)
        variable = termino.replace(numero, '')
        if variable:
            if variable not in variables:
                variables[variable] = len(variables)
            coeficientes.append((variables[variable], coeficiente))
        else:
            coeficientes.append((None, coeficiente))

    fila = [0.0] * len(variables) + [float(der)]  # Ya no fallará aquí
    for variable_index, coeficiente in coeficientes:
        if variable_index is not None:
            fila[variable_index] = coeficiente
        else:
            fila[-1] -= coeficiente

    return fila, variables

# Interfaz gráfica
class SistemaEcuacionesGUI:
    def __init__(self, master):
        self.master = master
        master.title("Sistema de Ecuaciones - Método Gauss-Jordan")
        master.configure(bg="#2E3440")

        style = ttk.Style()
        style.theme_use("clam")

        # Estilo para los botones
        style.configure("TButton",
                        background="#5E81AC",
                        foreground="white",
                        font=('Helvetica', 10, 'bold'),
                        padding=6)
        style.map("TButton",
                  background=[("active", "#81A1C1")])

        # Estilo para las etiquetas
        style.configure("TLabel",
                        background="#2E3440",
                        foreground="white",
                        font=('Helvetica', 11, 'bold'))

        # Estilo para las entradas
        style.configure("TEntry",
                        fieldbackground="#D8DEE9",
                        foreground="#2E3440",
                        padding=5)

        self.ecuaciones = []
        self.variables_encontradas = {}

        self.ecuacion_label = ttk.Label(master, text="Ingrese la ecuación:")
        self.ecuacion_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.ecuacion_entry = ttk.Entry(master, width=50)
        self.ecuacion_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        self.agregar_button = ttk.Button(master, text="Agregar Ecuación", command=self.agregar_ecuacion)
        self.agregar_button.grid(row=1, column=0, columnspan=2, pady=5)

        self.ecuaciones_listbox_label = ttk.Label(master, text="Ecuaciones Ingresadas:")
        self.ecuaciones_listbox_label.grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.ecuaciones_listbox = tk.Listbox(master, width=60, height=7, bg="#3B4252", fg="white",
                                            selectbackground="#81A1C1", selectforeground="black",
                                            highlightbackground="#4C566A", relief="flat", font=('Helvetica', 10))
        self.ecuaciones_listbox.grid(row=3, column=0, columnspan=2, padx=5, pady=5)

        self.mostrar_button = ttk.Button(master, text="Mostrar Sistema Resuelto", command=self.mostrar_sistema_resuelto)
        self.mostrar_button.grid(row=4, column=0, columnspan=2, pady=10)

        self.sistema_label = ttk.Label(master, text="Sistema Resuelto (Forma Escalonada Reducida):")
        self.sistema_label.grid(row=5, column=0, padx=5, pady=5, sticky="w")

        self.sistema_text = tk.Text(master, width=60, height=7, bg="#3B4252", fg="white",
                                    insertbackground="white", relief="flat", font=('Consolas', 10))
        self.sistema_text.grid(row=6, column=0, columnspan=2, padx=5, pady=5)
        self.sistema_text.config(state=tk.DISABLED)

        self.variables_label = ttk.Label(master, text="Soluciones:")
        self.variables_label.grid(row=7, column=0, padx=5, pady=5, sticky="w")

        self.variables_text = tk.Text(master, width=60, height=5, bg="#3B4252", fg="white",
                                     insertbackground="white", relief="flat", font=('Consolas', 10))
        self.variables_text.grid(row=8, column=0, columnspan=2, padx=5, pady=5)
        self.variables_text.config(state=tk.DISABLED)

        master.grid_columnconfigure(1, weight=1)

    def agregar_ecuacion(self):
        ecuacion = self.ecuacion_entry.get()
        if ecuacion:
            self.ecuaciones.append(ecuacion)
            self.ecuaciones_listbox.insert(tk.END, ecuacion)
            self.ecuacion_entry.delete(0, tk.END)
        else:
            messagebox.showerror("Error", "Por favor, ingrese una ecuación.")

    def mostrar_sistema_resuelto(self):
        if not self.ecuaciones:
            messagebox.showinfo("Información", "No se han ingresado ecuaciones.")
            return

        sistema_de_ecuaciones = []
        todas_las_variables = {}

        for ecuacion in self.ecuaciones:
            try:
                fila, variables_ecuacion = matriz_aum(ecuacion)
                sistema_de_ecuaciones.append(fila)
                todas_las_variables.update(variables_ecuacion)
            except Exception as e:
                messagebox.showerror("Error", f"Error al procesar la ecuación: '{ecuacion}'.\n\nDetalle: {e}")
                return

        if not sistema_de_ecuaciones:
            messagebox.showinfo("Información", "No se pudieron convertir las ecuaciones a una matriz.")
            return

        matriz_res = hacer_pivotes_completo(sistema_de_ecuaciones)

        if matriz_res is not None:
            self.mostrar_matriz_en_text(self.sistema_text, matriz_res)
            self.mostrar_soluciones(matriz_res, list(todas_las_variables.keys()))
        else:
            messagebox.showerror("Error", "Ocurrió un error al realizar las operaciones de pivote.")

    def mostrar_matriz_en_text(self, text_widget, matriz):
        text_widget.config(state=tk.NORMAL)
        text_widget.delete("1.0", tk.END)
        for fila in matriz:
            text_widget.insert(tk.END, "  ".join(f"{num:8.3f}" for num in fila) + "\n")
        text_widget.config(state=tk.DISABLED)

    def mostrar_soluciones(self, matriz_res, nombres_variables):
        num_filas, num_cols = matriz_res.shape
        num_variables = len(nombres_variables)
        soluciones = {}

        self.variables_text.config(state=tk.NORMAL)
        self.variables_text.delete("1.0", tk.END)

        if num_variables > 0:
            for i in range(num_filas):
                # Encuentra el primer coeficiente no nulo (pivote)
                pivot_col = -1
                for j in range(num_variables):
                    if abs(matriz_res[i][j]) > 1e-9:  # Usar una tolerancia para comparar con cero
                        pivot_col = j
                        break

                if pivot_col != -1:
                    variable_nombre = nombres_variables[pivot_col]
                    soluciones[variable_nombre] = matriz_res[i][num_cols - 1]
                    self.variables_text.insert(tk.END, f"{variable_nombre} = {soluciones[variable_nombre]:.3f}\n")
                elif all(abs(x) < 1e-9 for x in matriz_res[i][:num_variables]) and abs(matriz_res[i][num_cols - 1]) > 1e-9:
                    self.variables_text.insert(tk.END, "El sistema no tiene solución (inconsistente).\n")
                    self.variables_text.config(state=tk.DISABLED)
                    return
                elif all(abs(x) < 1e-9 for x in matriz_res[i]) and num_variables > 0:
                    self.variables_text.insert(tk.END, "El sistema tiene infinitas soluciones.\n")
                    self.variables_text.config(state=tk.DISABLED)
                    return
        elif num_cols > 1 and any(abs(matriz_res[i][0]) > 1e-9 for i in range(num_filas)):
            self.variables_text.insert(tk.END, "El sistema tiene una solución constante.\n")
            for i in range(num_filas):
                if all(abs(x) < 1e-9 for x in matriz_res[i][:num_variables]) and num_variables == 0:
                    self.variables_text.insert(tk.END, f"Constante {i+1} = {matriz_res[i][num_cols - 1]:.3f}\n")
        elif num_cols > 1 and all(abs(matriz_res[i][0]) < 1e-9 for i in range(num_filas)) and any(abs(matriz_res[i][num_cols - 1]) > 1e-9 for i in range(num_filas)):
            self.variables_text.insert(tk.END, "El sistema no tiene solución (inconsistente).\n")
        elif num_cols > 1 and all(abs(x) < 1e-9 for fila in matriz_res for x in fila):
            self.variables_text.insert(tk.END, "El sistema tiene infinitas soluciones (solo constantes).\n")
        elif num_variables == 0 and num_cols == 1:
            self.variables_text.insert(tk.END, "Sistema sin variables.\n")


        self.variables_text.config(state=tk.DISABLED)


if __name__ == "__main__":
    root = tk.Tk()
    app = SistemaEcuacionesGUI(root)
    root.mainloop()