def newton_raphson(func, dfunc, x0, tol=1e-6, max_iter=100):
    """
    Método de Newton-Raphson para encontrar la raíz de una función.

    Args:
        func: La función de la cual encontrar la raíz.
        dfunc: La derivada de la función.
        x0: Valor inicial.
        tol: Tolerancia para la convergencia.
        max_iter: Número máximo de iteraciones.

    Returns:
        La raíz encontrada o un mensaje de error si no converge.
    """
    x = x0
    for i in range(max_iter):
        f_x = func(x)
        df_x = dfunc(x)
        if abs(df_x) < 1e-10:
            raise ValueError("La derivada es demasiado pequeña, posible división por cero.")
        
        x_new = x - f_x / df_x
        
        if abs(x_new - x) < tol:
            return x_new
        
        x = x_new

    raise ValueError("No se alcanzó la convergencia después de {} iteraciones.".format(max_iter))

# Datos del problema
z = [0.1,0.2,0.3,0.4]  # Fracción molar en la alimentación
K = [4.2,1.75,0.74,0.34]  # Constantes de equilibrio

# Definimos la función Rachford-Rice
def rachford_rice(V):
    return sum(z[i] * (K[i] - 1) / (1 + V * (K[i] - 1)) for i in range(len(z)))

# Derivada de la función Rachford-Rice
def d_rachford_rice(V):
    return sum(-z[i] * (K[i] - 1)**2 / (1 + V * (K[i] - 1))**2 for i in range(len(z)))

# Estimación inicial
V0 = 0

# Resolver utilizando Newton-Raphson
try:
    V_solution = newton_raphson(rachford_rice, d_rachford_rice, V0)
    print("La solución para V (fracción de fase vapor) es:", V_solution)
except ValueError as e:
    print("Error:", e)
