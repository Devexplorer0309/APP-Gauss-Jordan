import tkinter as tk
from tkinter import ttk, messagebox
import re
import numpy as np


# Función para hacer los pivotes
def hacer_piv(matriz):
    matriz = np.array(matriz, dtype=float)
    n = len(matriz)


    for i in range(n):
        if matriz[i][i] != 0:
            matriz[i] = matriz[i] / matriz[i][i]
        else:
            print(f"No se puede hacer pivote en la fila {i} ya que el elemento es 0")
            return None
    return matriz


# Función para convertir una ecuación a matriz aumentada
def matriz_aum(ec):
    ec = ec.strip().replace(" ", "").replace('-', '+-')
    izq, der = ec.split('=')
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


    fila = [0.0] * len(variables) + [float(der)]
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
        master.configure(bg="#2E3440")  # Fondo elegante


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


        self.mostrar_button = ttk.Button(master, text="Mostrar Sistema y Variables", command=self.mostrar_sistema)
        self.mostrar_button.grid(row=4, column=0, columnspan=2, pady=10)


        self.sistema_label = ttk.Label(master, text="Sistema Guardado:")
        self.sistema_label.grid(row=5, column=0, padx=5, pady=5, sticky="w")


        self.sistema_text = tk.Text(master, width=60, height=7, bg="#3B4252", fg="white",
                                    insertbackground="white", relief="flat", font=('Consolas', 10))
        self.sistema_text.grid(row=6, column=0, columnspan=2, padx=5, pady=5)
        self.sistema_text.config(state=tk.DISABLED)


        master.grid_columnconfigure(1, weight=1)


    def agregar_ecuacion(self):
        ecuacion = self.ecuacion_entry.get()
        if ecuacion:
            self.ecuaciones.append(ecuacion)
            self.ecuaciones_listbox.insert(tk.END, ecuacion)
            self.ecuacion_entry.delete(0, tk.END)
        else:
            messagebox.showerror("Error", "Por favor, ingrese una ecuación.")


    def mostrar_sistema(self):
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


        matriz_pivoteada = hacer_piv(sistema_de_ecuaciones)


        if matriz_pivoteada is None:
            return


        self.sistema_text.config(state=tk.NORMAL)
        self.sistema_text.delete("1.0", tk.END)
        for fila in matriz_pivoteada:
            self.sistema_text.insert(tk.END, "  ".join(f"{num:8.3f}" for num in fila) + "\n")
        self.sistema_text.config(state=tk.DISABLED)


if __name__ == "__main__":
    root = tk.Tk()
    app = SistemaEcuacionesGUI(root)
    root.mainloop()
