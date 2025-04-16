# -*- coding: utf-8 -*- # Para asegurar compatibilidad con acentos
import tkinter as tk
from tkinter import ttk, messagebox
import re
import numpy as np

#--------------------------------------------------------------------------
# Función para realizar la Eliminación de Gauss-Jordan
#--------------------------------------------------------------------------
def eliminacion_gauss_jordan(matriz):
    """
    Convierte la matriz aumentada a su forma escalonada reducida por filas (RREF)
    utilizando el método de eliminación de Gauss-Jordan.

    Args:
        matriz (list or np.array): La matriz aumentada del sistema de ecuaciones.

    Returns:
        np.array or None: La matriz en forma escalonada reducida (RREF) o None si ocurre un error.
    """
    try:
        matriz = np.array(matriz, dtype=float)
        num_filas, num_cols = matriz.shape
        fila_actual = 0
        columna_actual = 0

        while fila_actual < num_filas and columna_actual < num_cols - 1:
            # Pivoteo Parcial
            fila_pivote = fila_actual
            for i in range(fila_actual + 1, num_filas):
                if abs(matriz[i][columna_actual]) > abs(matriz[fila_pivote][columna_actual]):
                    fila_pivote = i

            if abs(matriz[fila_pivote][columna_actual]) < 1e-9:
                columna_actual += 1
                continue

            if fila_pivote != fila_actual:
                matriz[[fila_actual, fila_pivote]] = matriz[[fila_pivote, fila_actual]]

            # Normalización
            pivote = matriz[fila_actual][columna_actual]
            if abs(pivote) > 1e-9:
                 matriz[fila_actual] = matriz[fila_actual] / pivote
            else:
                 columna_actual += 1
                 continue

            # Eliminación Gauss-Jordan (Ceros arriba y abajo)
            for i in range(num_filas):
                if i != fila_actual:
                    factor = matriz[i][columna_actual]
                    matriz[i] = matriz[i] - factor * matriz[fila_actual]

            fila_actual += 1
            columna_actual += 1

        matriz = np.round(matriz, 9)
        matriz[np.abs(matriz) < 1e-9] = 0.0
        return matriz

    except Exception as e:
        messagebox.showerror("Error en Cálculo", f"Ocurrió un error durante la eliminación Gauss-Jordan: {e}")
        return None

#--------------------------------------------------------------------------
# Función para convertir una ecuación a matriz aumentada
#--------------------------------------------------------------------------
def matriz_aum(ec, variables_dict):
    """
    Convierte una cadena de texto representando una ecuación lineal
    en una fila de la matriz aumentada y actualiza el diccionario de variables.
    """
    try:
        ec = ec.strip().replace(" ", "")
        if '=' not in ec: raise ValueError("La ecuación debe contener un signo '='.")
        partes = ec.split('=')
        if len(partes) != 2: raise ValueError("La ecuación debe tener solo un signo '='.")
        izq, der_str = partes

        try: der = float(der_str)
        except ValueError: raise ValueError(f"El lado derecho ('{der_str}') debe ser un número.")

        if izq.startswith('-'): izq = "-" + izq[1:].replace('-', '+-')
        else: izq = izq.replace('-', '+-')
        if izq.startswith('+'): izq = izq[1:]

        # Regex mejorada para capturar: [signo opcional][numero opcional][variable] O [signo obligatorio][numero]
        terminos_regex = r'([+-]?)(\d+\.?\d*|\.\d+)?([a-zA-Z_]\w*)' # Captura coef explícito o implícito 1
        constantes_regex = r'([+-])(\d+\.?\d*|\.\d+)(?![a-zA-Z_.\w])' # Captura constantes con signo
        constante_inicial_regex = r'^(\d+\.?\d*|\.\d+)(?![a-zA-Z_.\w])' # Captura constante inicial sin signo

        coeficientes_locales = {}
        constante_izq = 0.0
        procesado_hasta = 0

        # 1. Constante inicial sin signo
        match_inicio = re.match(constante_inicial_regex, izq)
        if match_inicio:
            constante_izq += float(match_inicio.group(1))
            procesado_hasta = match_inicio.end()

        # 2. Términos con variables y coeficientes
        for match in re.finditer(terminos_regex, izq[procesado_hasta:]):
            signo, num_str, var = match.groups()
            if not var: continue # No es una variable válida

            coeficiente = 1.0
            if num_str: coeficiente = float(num_str)
            if signo == '-': coeficiente *= -1.0

            if var not in variables_dict: variables_dict[var] = len(variables_dict)
            var_index = variables_dict[var]
            coeficientes_locales[var_index] = coeficientes_locales.get(var_index, 0.0) + coeficiente
            # Marcar como procesado para no contarlo como constante después
            izq = izq[:procesado_hasta + match.start()] + ' ' * (match.end() - match.start()) + izq[procesado_hasta + match.end():]


        # 3. Constantes restantes con signo explícito
        # Usar la versión original de izq para encontrar constantes que no eran parte de términos con variables
        izq_temp = izq[procesado_hasta:] # Analizar solo lo que no era constante inicial
        for match in re.finditer(constantes_regex, izq_temp):
             # Asegurarse de que este match no esté dentro de una parte ya procesada (marcada con espacios)
             if izq[procesado_hasta + match.start(): procesado_hasta + match.end()].strip():
                 signo, num_str = match.groups()
                 valor = float(num_str)
                 constante_izq += valor if signo == '+' else -valor
                 # Marcar como procesado
                 izq = izq[:procesado_hasta + match.start()] + ' ' * (match.end() - match.start()) + izq[procesado_hasta + match.end():]


        # Verificar si quedó algo no procesado que no sea solo +/- o espacios
        resto = izq.replace('+', '').replace('-', '').replace('.', '').strip()
        if resto:
             raise ValueError(f"Formato inválido. Parte no reconocida: '{resto}'")


        # Construir fila
        lado_derecho_final = der - constante_izq
        num_total_variables = len(variables_dict)
        fila = [0.0] * (num_total_variables + 1)

        for var_index, coef in coeficientes_locales.items():
            if var_index < num_total_variables: fila[var_index] = coef

        fila[-1] = lado_derecho_final
        return fila, variables_dict

    except ValueError as ve:
         messagebox.showerror("Error de Formato", f"Error en la ecuación: '{ec}'\nDetalle: {ve}")
         return None, variables_dict
    except Exception as e:
        messagebox.showerror("Error Inesperado", f"Error procesando ecuación: '{ec}'\nDetalle: {e}")
        return None, variables_dict

#--------------------------------------------------------------------------
# Interfaz Gráfica (GUI)
#--------------------------------------------------------------------------
class SistemaEcuacionesGUI:
    def __init__(self, master):
        self.master = master
        master.title("Sistema de Ecuaciones - Método Gauss-Jordan")
        master.configure(bg="#2E3440")

        # --- Estilos ttk ---
        style = ttk.Style()
        style.theme_use("clam")
        style.configure(".", background="#2E3440", foreground="white", font=('Helvetica', 10))
        style.map(".", foreground=[('disabled', '#6c757d')])
        style.configure("TButton", background="#5E81AC", foreground="#ECEFF4", font=('Helvetica', 10, 'bold'), padding=6, relief="flat", borderwidth=0)
        style.map("TButton", background=[("active", "#81A1C1"), ("disabled", "#4C566A")], foreground=[("disabled", "#D8DEE9")])
        style.configure("TLabel", background="#2E3440", foreground="#ECEFF4", font=('Helvetica', 11, 'bold'))
        style.configure("TEntry", fieldbackground="#4C566A", foreground="#ECEFF4", insertcolor="#ECEFF4", padding=5, relief="flat", borderwidth=1)
        style.map("TEntry", background=[("active", "#434C5E")], relief=[("focus", "solid")])

        # --- Variables de instancia ---
        self.ecuaciones_ingresadas = []
        self.variables_encontradas = {}

        # --- Widgets ---
        self.ecuacion_label = ttk.Label(master, text="Ingrese la ecuación (ej: 2x + 3y = 5):")
        self.ecuacion_label.grid(row=0, column=0, padx=10, pady=8, sticky="w")
        self.ecuacion_entry = ttk.Entry(master, width=60)
        self.ecuacion_entry.grid(row=0, column=1, padx=10, pady=8, sticky="ew")
        self.ecuacion_entry.bind("<Return>", self.agregar_ecuacion_event)

        self.agregar_button = ttk.Button(master, text="Agregar Ecuación", command=self.agregar_ecuacion)
        self.agregar_button.grid(row=1, column=0, columnspan=2, pady=8)

        self.ecuaciones_listbox_label = ttk.Label(master, text="Ecuaciones Ingresadas:")
        self.ecuaciones_listbox_label.grid(row=2, column=0, padx=10, pady=(10, 2), sticky="w")
        self.ecuaciones_listbox = tk.Listbox(master, width=70, height=7, bg="#3B4252", fg="#ECEFF4", selectbackground="#81A1C1", selectforeground="#2E3440", highlightthickness=1, highlightbackground="#4C566A", highlightcolor="#5E81AC", relief="flat", font=('Consolas', 10))
        self.ecuaciones_listbox.grid(row=3, column=0, columnspan=2, padx=10, pady=5, sticky="ew")

        self.resolver_button = ttk.Button(master, text="Resolver Sistema (Gauss-Jordan)", command=self.resolver_sistema)
        self.resolver_button.grid(row=4, column=0, columnspan=2, pady=10)

        self.matriz_resuelta_label = ttk.Label(master, text="Matriz en Forma Escalonada Reducida (Gauss-Jordan):")
        self.matriz_resuelta_label.grid(row=5, column=0, padx=10, pady=(10, 2), sticky="w")
        self.matriz_text = tk.Text(master, width=70, height=7, bg="#3B4252", fg="#ECEFF4", insertbackground="#ECEFF4", relief="flat", font=('Consolas', 10), highlightthickness=1, highlightbackground="#4C566A", highlightcolor="#5E81AC", wrap=tk.NONE)
        self.matriz_text.grid(row=6, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        self.matriz_text.config(state=tk.DISABLED)

        self.soluciones_label = ttk.Label(master, text="Interpretación de la Solución:")
        self.soluciones_label.grid(row=7, column=0, padx=10, pady=(10, 2), sticky="w")
        self.soluciones_text = tk.Text(master, width=70, height=7, bg="#3B4252", fg="#ECEFF4", insertbackground="#ECEFF4", relief="flat", font=('Consolas', 10), highlightthickness=1, highlightbackground="#4C566A", highlightcolor="#5E81AC", wrap=tk.WORD) # Height increased
        self.soluciones_text.grid(row=8, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        self.soluciones_text.config(state=tk.DISABLED)

        self.limpiar_button = ttk.Button(master, text="Limpiar Todo", command=self.limpiar_todo)
        self.limpiar_button.grid(row=9, column=0, columnspan=2, pady=15)

        master.grid_columnconfigure(1, weight=1)

    # --- Métodos de la GUI ---

    def agregar_ecuacion_event(self, event):
        self.agregar_ecuacion()

    def agregar_ecuacion(self):
        ecuacion = self.ecuacion_entry.get()
        if not ecuacion.strip():
            messagebox.showwarning("Entrada Vacía", "Por favor, ingrese una ecuación válida.")
            return
        if '=' not in ecuacion or ecuacion.count('=') > 1:
            messagebox.showwarning("Formato Inválido", "La ecuación debe contener exactamente un signo '='.")
            return

        # Intento preliminar de parseo para validación temprana
        temp_vars = self.variables_encontradas.copy() # No modificar el diccionario real aún
        fila_test, _ = matriz_aum(ecuacion, temp_vars)
        if fila_test is None:
             # El error ya se mostró en matriz_aum
             return

        self.ecuaciones_ingresadas.append(ecuacion)
        self.ecuaciones_listbox.insert(tk.END, ecuacion)
        self.ecuacion_entry.delete(0, tk.END)
        self.limpiar_resultados()

    def limpiar_resultados(self):
        for text_widget in [self.matriz_text, self.soluciones_text]:
            text_widget.config(state=tk.NORMAL)
            text_widget.delete("1.0", tk.END)
            text_widget.config(state=tk.DISABLED)

    def limpiar_todo(self):
        self.ecuacion_entry.delete(0, tk.END)
        self.ecuaciones_ingresadas.clear()
        self.variables_encontradas.clear()
        self.ecuaciones_listbox.delete(0, tk.END)
        self.limpiar_resultados()

    def resolver_sistema(self):
        if not self.ecuaciones_ingresadas:
            messagebox.showinfo("Sin Ecuaciones", "No se han ingresado ecuaciones para resolver.")
            return

        self.variables_encontradas = {}
        matriz_aumentada_filas = []
        error_parseo = False

        for ecuacion_str in self.ecuaciones_ingresadas:
            fila, self.variables_encontradas = matriz_aum(ecuacion_str, self.variables_encontradas)
            if fila is None:
                error_parseo = True
                break
            matriz_aumentada_filas.append(fila)

        if error_parseo:
             self.limpiar_resultados(); return

        num_vars_final = len(self.variables_encontradas)
        if num_vars_final == 0 and matriz_aumentada_filas:
            # Chequear si alguna fila es inconsistente (ej: 5=3 -> fila=[-2.0])
            if any(abs(f[0]) > 1e-9 for f in matriz_aumentada_filas):
                 messagebox.showerror("Error", "Sistema inconsistente detectado (constantes desiguales).")
            else: # Todas son identidades (ej: 5=5 -> fila=[0.0])
                 messagebox.showinfo("Info", "Las ecuaciones son identidades (ej: 5=5) o no contienen variables.")
            self.limpiar_resultados(); return
        elif not matriz_aumentada_filas:
             messagebox.showerror("Error", "No se pudo construir la matriz."); return

        matriz_completa = []
        for fila in matriz_aumentada_filas:
             fila_coef = fila[:-1]; fila_res = fila[-1]
             fila_ajustada = fila_coef + [0.0] * (num_vars_final - len(fila_coef)) + [fila_res]
             # Asegurar longitud correcta incluso si num_vars_final=0
             if num_vars_final == 0: fila_ajustada=[fila_res]

             # Corregir longitud si la fila original tenía más "variables" que las finales encontradas
             # (esto no debería pasar con la lógica actual de matriz_aum, pero por seguridad)
             matriz_completa.append(fila_ajustada[:num_vars_final+1])


        if not matriz_completa:
             messagebox.showerror("Error", "No se pudo construir la matriz final."); return

        matriz_np = np.array(matriz_completa, dtype=float)
        matriz_resuelta_rref = eliminacion_gauss_jordan(matriz_np)

        if matriz_resuelta_rref is not None:
            self.mostrar_matriz_en_text(self.matriz_text, matriz_resuelta_rref)
            nombres_variables_ordenados = sorted(self.variables_encontradas.keys(), key=lambda v: self.variables_encontradas[v])
            # Pasar la matriz RREF a la función de interpretación
            self.interpretar_y_mostrar_soluciones(matriz_resuelta_rref, nombres_variables_ordenados)
        else:
            self.limpiar_resultados()


    def mostrar_matriz_en_text(self, text_widget, matriz):
        text_widget.config(state=tk.NORMAL)
        text_widget.delete("1.0", tk.END)
        if matriz is None or matriz.size == 0:
            text_widget.insert(tk.END, "Matriz vacía o inválida.")
            text_widget.config(state=tk.DISABLED)
            return

        num_filas, num_cols = matriz.shape
        if num_cols == 0: # Manejar matriz sin columnas
             text_widget.insert(tk.END, "[Matriz sin columnas]\n" * num_filas)
             text_widget.config(state=tk.DISABLED)
             return

        max_len = 0
        for fila in matriz:
            for num in fila:
                 s = f"{num:.3f}"; max_len = max(max_len, len(s))
        col_width = max_len + 2

        matriz_str = ""
        for i in range(num_filas):
            fila_str_list = []
            for j in range(num_cols):
                num_str = f"{matriz[i, j]:.3f}"
                # Añadir barra antes de la última columna si hay más de una columna
                if j == num_cols - 1 and num_cols > 1:
                    fila_str_list.append("|")
                fila_str_list.append(num_str.rjust(col_width))
            matriz_str += " ".join(fila_str_list) + "\n"

        text_widget.insert(tk.END, matriz_str)
        text_widget.config(state=tk.DISABLED)


    # --- COMIENZO DE LA SECCIÓN CORREGIDA ---
    def interpretar_y_mostrar_soluciones(self, matriz_rref, nombres_variables):
        """Interpreta la matriz en RREF y muestra las soluciones."""
        self.soluciones_text.config(state=tk.NORMAL)
        self.soluciones_text.delete("1.0", tk.END)

        if matriz_rref is None or matriz_rref.size == 0:
             self.soluciones_text.insert(tk.END, "No se puede determinar la solución (matriz RREF inválida).") # CORREGIDO
             self.soluciones_text.config(state=tk.DISABLED)
             return

        num_filas, num_cols = matriz_rref.shape
        num_variables = len(nombres_variables)

        if num_cols == 0:
             self.soluciones_text.insert(tk.END, "Matriz RREF vacía.") # CORREGIDO
             self.soluciones_text.config(state=tk.DISABLED)
             return
        # Si hay variables, debe haber al menos num_variables + 1 columnas
        if num_variables > 0 and num_cols <= num_variables:
            self.soluciones_text.insert(tk.END, f"Error: Dimensiones inconsistentes en RREF (Variables: {num_variables}, Columnas: {num_cols}).") # CORREGIDO
            self.soluciones_text.config(state=tk.DISABLED)
            return
        # Si no hay variables (num_variables=0), debe haber 1 columna (resultado)
        if num_variables == 0 and num_cols != 1:
             self.soluciones_text.insert(tk.END, f"Error: Dimensiones inconsistentes en RREF (Sin variables, Columnas: {num_cols}).") # CORREGIDO
             self.soluciones_text.config(state=tk.DISABLED)
             return


        # --- Chequeo de Inconsistencia (Fila [0 0 ... 0 | k] con k != 0) ---
        inconsistente = False
        for i in range(num_filas):
            # Coeficientes son todas las columnas excepto la última (si existe)
            coeficientes_fila = matriz_rref[i, :num_variables]
            # Resultado es la última columna
            resultado_fila = matriz_rref[i, -1]
            # Comprobar si TODOS los coeficientes son ~0
            if np.all(np.abs(coeficientes_fila) < 1e-9):
                # Si además el resultado NO es ~0, es inconsistente
                if abs(resultado_fila) > 1e-9:
                    self.soluciones_text.insert(tk.END, "El sistema es INCONSISTENTE (no tiene solución).\n") # CORREGIDO
                    self.soluciones_text.insert(tk.END, f"(Se encontró una fila tipo [ 0 ... 0 | {resultado_fila:.3f} ≠ 0 ])\n") # CORREGIDO
                    inconsistente = True
                    break # No hay más que analizar

        if inconsistente:
            self.soluciones_text.config(state=tk.DISABLED)
            return

        # Si no hay variables y llegamos aquí, el sistema es consistente (e.g. 0=0)
        if num_variables == 0:
             self.soluciones_text.insert(tk.END, "El sistema es consistente (ej: 0 = 0) y no tiene variables.") # CORREGIDO
             self.soluciones_text.config(state=tk.DISABLED)
             return


        # --- Identificación de Pivotes y Rango ---
        pivote_cols = []
        for r in range(num_filas):
             primer_no_cero_col = -1
             # Buscar primer no-cero en la parte de coeficientes
             for c in range(num_variables):
                 if abs(matriz_rref[r, c]) > 1e-9:
                      primer_no_cero_col = c
                      break
             # Si se encontró un primer no-cero y es ~1 y es el único en su columna -> es pivote
             if primer_no_cero_col != -1 and abs(matriz_rref[r, primer_no_cero_col] - 1.0) < 1e-9:
                 if abs(np.sum(np.abs(matriz_rref[:, primer_no_cero_col])) - 1.0) < 1e-9:
                     # Asegurar que no añadamos la misma columna de pivote dos veces
                     if primer_no_cero_col not in pivote_cols:
                          pivote_cols.append(primer_no_cero_col)

        rango = len(pivote_cols) # Rango es el número de pivotes

        # --- Determinar Tipo de Solución ---
        if rango < num_variables:
            # --- Infinitas Soluciones ---
            self.soluciones_text.insert(tk.END, "El sistema tiene INFINITAS SOLUCIONES.\n") # CORREGIDO
            variables_libres_indices = [idx for idx in range(num_variables) if idx not in pivote_cols]
            variables_libres_nombres = [nombres_variables[idx] for idx in variables_libres_indices]

            if variables_libres_nombres:
                 self.soluciones_text.insert(tk.END, f"Variables libres: {', '.join(variables_libres_nombres)}\n\n") # CORREGIDO
            else:
                 self.soluciones_text.insert(tk.END, "(Sistema dependiente, sin variables libres explícitas)\n\n") # CORREGIDO

            # Expresar variables básicas (con pivote) en términos de las libres
            for var_idx_basica in pivote_cols:
                nombre_var_basica = nombres_variables[var_idx_basica]
                fila_pivote_actual = -1
                for r in range(num_filas):
                     if abs(matriz_rref[r, var_idx_basica] - 1.0) < 1e-9:
                          fila_pivote_actual = r; break

                if fila_pivote_actual != -1:
                    resultado = matriz_rref[fila_pivote_actual, -1]
                    expresion = f"{nombre_var_basica} = {resultado:.3f}"
                    for libre_idx in variables_libres_indices:
                        coef = matriz_rref[fila_pivote_actual, libre_idx]
                        if abs(coef) > 1e-9:
                            expresion += f" {-coef:+.3f}{nombres_variables[libre_idx]}"
                    self.soluciones_text.insert(tk.END, expresion + "\n") # CORREGIDO
                else:
                    self.soluciones_text.insert(tk.END, f"Error interno: No se encontró fila de pivote para {nombre_var_basica}\n") # CORREGIDO

        elif rango == num_variables:
             # --- Solución Única ---
             self.soluciones_text.insert(tk.END, "El sistema tiene SOLUCIÓN ÚNICA:\n") # CORREGIDO
             for var_idx in range(num_variables):
                  fila_pivote_actual = -1
                  for r in range(num_filas):
                       if abs(matriz_rref[r, var_idx] - 1.0) < 1e-9:
                            fila_pivote_actual = r; break
                  if fila_pivote_actual != -1:
                      nombre_variable = nombres_variables[var_idx]
                      valor_variable = matriz_rref[fila_pivote_actual, -1]
                      self.soluciones_text.insert(tk.END, f"{nombre_variable} = {valor_variable:.3f}\n") # CORREGIDO
                  else:
                       # Puede pasar si hay filas de ceros, pero debería haber un pivote para cada variable si rango=num_vars
                       # O si la RREF no es perfecta.
                       # Buscar si la variable tiene solución en alguna fila aunque no sea pivote perfecto
                       solucion_encontrada = False
                       for r in range(num_filas):
                           # Si esta fila define principalmente esta variable (coef es ~1) y otros coef son ~0
                           if abs(matriz_rref[r, var_idx]) > 1e-9: # Coef no cero para nuestra variable
                               # ¿Es el único coeficiente significativo en la fila?
                               otros_coef_significativos = False
                               for c in range(num_variables):
                                   if c != var_idx and abs(matriz_rref[r, c]) > 1e-9:
                                       otros_coef_significativos = True; break
                               if not otros_coef_significativos:
                                    # Podemos despejarla desde aquí
                                    valor_variable = matriz_rref[r, -1] / matriz_rref[r, var_idx]
                                    self.soluciones_text.insert(tk.END, f"{nombres_variables[var_idx]} ≈ {valor_variable:.3f} (calculado)\n") # CORREGIDO
                                    solucion_encontrada = True; break
                       if not solucion_encontrada:
                           self.soluciones_text.insert(tk.END, f"Advertencia: No se pudo determinar un valor claro para '{nombres_variables[var_idx]}'\n") # CORREGIDO


        self.soluciones_text.config(state=tk.DISABLED)
    # --- FIN DE LA SECCIÓN CORREGIDA ---


# --- Ejecución Principal ---
if __name__ == "__main__":
    root = tk.Tk()
    app = SistemaEcuacionesGUI(root)
    root.mainloop()