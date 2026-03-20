import math
import matplotlib.pyplot as plt

def f(x):
    return (x ** 3) - 3 * x + 1 

def df(x):
    return 3 * x ** 2 - 3

def exist_root(a, b):
    # Проверка знаков на концах отрезка
    if f(a) * f(b) >= 0:
        return False, "Функция не меняет знак на концах отрезка"
    
    # Проверка монотонности (производная не меняет знак)
    if df(a) * df(b) < 0:
        return False, "Производная меняет знак, возможны несколько корней"
    
    return True, "Корень существует и единственный на данном отрезке"

def bisection_method(a, b, eps=1e-6, max_iter=100):
    exists, message = exist_root(a, b)
    print(f"Проверка существования корня: {message}")
    if not exists:
        return None, None, []
    
    iter = 0
    errors = []
    
    while abs(b - a) > 2 * eps and iter < max_iter:
        x = (a + b) / 2
        fx = f(x)

        errors.append(abs(b - a))
        
        if abs(fx) < eps:
            return x, iter + 1, errors
        
        if f(a) * f(x) < 0:
            b = x
        else:
            a = x
        
        iter += 1
    
    return (a + b) / 2, iter + 1, errors

def plot_function(a, b, root=None, errors=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Первый график - функция
    x_points = [a + (b - a) * i / 100 for i in range(101)]
    y_points = [f(x) for x in x_points]
    
    ax1.plot(x_points, y_points, 'b-', linewidth=2, label='f(x) = x³ - 3x + 1')
    ax1.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax1.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax1.grid(True, alpha=0.3)
    
    if root is not None:
        ax1.plot(root, f(root), 'ro', markersize=8, label=f'Корень: x ≈ {root:.6f}')
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('f(x)')
    ax1.set_title('График функции f(x) = x³ - 3x + 1')
    ax1.legend()
    
    # Второй график - сходимость
    if errors is not None and len(errors) > 0:
        iterations = list(range(1, len(errors) + 1))
        ax2.semilogy(iterations, errors, 'go-', linewidth=2, markersize=6, label='Длина интервала')
        ax2.set_xlabel('Количество итераций')
        ax2.set_ylabel('Длина интервала (лог. шкала)')
        ax2.set_title('График сходимости метода')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
    
    plt.tight_layout()
    plt.show()

def main():
    a, b = 0, 1

    root, iter, errors = bisection_method(a, b)

    if root is not None:
        print(f"x =  {root:.6f}")
        print(f"Значение функции в корне: {f(root):.2e}")
        print(f"Количество итераций: {iter}")
        
        plot_function(a, b, root, errors)
    else:
        print("Корень не найден")

if __name__ == "__main__":
    main()