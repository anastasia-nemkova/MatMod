import math
import matplotlib.pyplot as plt
import numpy as np


def func(theta, k, beta, M1):
    a = math.sin(theta) ** 2
    b = ((k + 1) / 2) * math.sin(theta) * math.sin(beta) / math.cos(theta - beta)
    c = (1 / M1) ** 2
    return a - b - c

def dfunc(theta, k, beta, M1):
    da = 2 * math.cos(theta) * math.sin(theta)

    up = ((k + 1) / 2) * math.sin(theta) * math.sin(beta)
    down = math.cos(theta - beta)

    dup = ((k + 1) / 2) * math.sin(beta) * math.cos(theta)
    ddown = -math.sin(theta - beta)

    db = (dup* down - ddown * up) / (down ** 2)

    dc = 0

    return da - db - dc

def d2func(theta, k, beta):
    d2a = 2 * (math.cos(theta) ** 2 -  math.sin(theta) ** 2)

    A = ((k + 1) / 2) * math.sin(beta)

    # A * cos(theta) / cos(theta - beta)

    up1 = A * math.cos(theta)
    down1 = math.cos(theta - beta)

    dup1 = -A * math.sin(theta)
    ddown1 = -math.sin(theta - beta)
   
    dterm1 = (dup1 * down1 - up1 * ddown1) / (down1 ** 2)

    # A * sin(theta - beta) * sin(theta) / cos(theta - beta) ^ 2

    up2 = A * math.sin(theta - beta) * math.sin(theta)
    down2 = math.cos(theta - beta) ** 2

    dup2 = A * (math.cos(theta - beta) * math.sin(theta) + math.cos(theta) * math.sin(theta - beta))
    ddown2 = -2 * math.cos(theta - beta) * math.sin(theta - beta) 

    dterm2 = (dup2 * down2 - up2 * ddown2) / (down2 ** 2)

    return d2a - dterm1 - dterm2

def newton_convergence(theta, beta, M1, k):
    """
    Проверка условия сходимости метода Ньютона:
    |f(x)f''(x)| < (f'(x))²
    """
    f = func(theta, k, beta, M1)
    df = dfunc(theta, k, beta, M1)
    d2f = d2func(theta, k, beta)

    left = abs(f * d2f)
    right = df ** 2

    return left < right


def newton_method(beta, M1, k, max_iter=100, epsilon=1e-6):
    # Начальное приближение
    theta = math.radians(30)
    theta_history = [theta]
    f_history = [func(theta, k, beta, M1)]
    
    for i in range(max_iter):
        f = func(theta, k, beta, M1)
        df = dfunc(theta, k, beta, M1)

        converges = newton_convergence(theta, beta, M1, k)
        if not converges:
           return None, i + 1, theta_history, f_history
        
        theta_new = theta - f / df

        theta_history.append(theta_new)
        f_history.append(func(theta_new, k, beta, M1))

        if abs(theta_new - theta) < epsilon:
            return theta_new, i + 1, theta_history, f_history
        
        f_new = func(theta_new, k, beta, M1)
        if abs(f_new) < epsilon:
            return theta_new, i + 1, theta_history, f_history
        
        theta = theta_new
    return None, max_iter, theta_history, f_history

def calc_parameters(theta, beta, M1, p1, T1, R, k):
    rho1 = p1 / (R * T1)
    print(rho1)
    a1 = math.sqrt(k * R * T1)
    v1 = M1 * a1
    u1 = 1 - ((math.sin(theta - beta) * math.cos(theta)) / (math.sin(theta) * math.cos(theta - beta)))
    v1_n = v1 * math.sin(theta)

    rho2 = rho1 / (1 - u1)
    p2 = p1 + rho1 * u1 * (v1_n ** 2)
    T2 = p2 / (rho2 * R)

    a2 = math.sqrt(k * R * T2)
    v1_t = v1 * math.cos(theta)

    v2_t = v1_t
    v2_n = v1_n * math.tan(theta - beta) / math.tan(theta)

    v2 = math.sqrt(v2_n ** 2 + v2_t ** 2)
    M2 = v2 / a2

    if M2 > 1:
        flow_regime = "сверхзвуковое"
    elif M2 < 1:
        flow_regime = "дозвуковое"
    else:
        flow_regime = "звуковое"
    
    return p2, T2, rho2, M2, flow_regime

def plot_func(beta, M1, k):
    theta_values = [math.radians(i) for i in range(1, 90)]
    f_values = [func(theta, k, beta, M1) for theta in theta_values]
    
    plt.figure(figsize=(10, 6))
    plt.plot([math.degrees(theta) for theta in theta_values], f_values, 'b-', linewidth=2)
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.7)
    plt.xlabel('Угол θ, градусы')
    plt.ylabel('f(θ)')
    plt.title('График функции f(θ) для определения угла наклона КСУ')
    plt.grid(True, alpha=0.3)
    plt.show()

def plot_newton_convergence(theta_history, f_history):
    """
    График сходимости метода Ньютона
    """
    iterations = range(len(theta_history))
    theta_degrees = [math.degrees(theta) for theta in theta_history]
    
    plt.figure(figsize=(12, 5))
    
    # График угла θ
    plt.subplot(1, 2, 1)
    plt.plot(iterations, theta_degrees, 'bo-', linewidth=2, markersize=6)
    plt.xlabel('Номер итерации')
    plt.ylabel('Угол θ, градусы')
    plt.title('Сходимость метода Ньютона: угол θ')
    plt.grid(True, alpha=0.3)
    
    # График функции f(θ)
    plt.subplot(1, 2, 2)
    plt.plot(iterations, f_history, 'ro-', linewidth=2, markersize=6)
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.7)
    plt.xlabel('Номер итерации')
    plt.ylabel('f(θ)')
    plt.title('Сходимость метода Ньютона: значение функции')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def main():
    k = 1.4 # аддиабата
    p1 = 63950
    T1 = 320
    M1 = 2 # число Маха набегающего потока
    beta_grad = 18 # угол клина
    R = 287 # газовая постоянная воздуха

    beta_rad = math.radians(beta_grad)

    theta_rad, iterations, theta_history, f_history = newton_method(beta_rad, M1, k)

    if theta_rad is None:
        print("NOT SOSHLIS")
        return

    theta = math.degrees(theta_rad)

    print(f"Угол наклона КСУ: θ = {theta:.3f}°")
    print(f"Количество итераций: {iterations}")
    print()

    p2, T2, rho2, M2, flow_regime = calc_parameters(theta_rad, beta_rad, M1, p1, T1, R, k)

    print("\n" + "=" * 80)
    print("ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 80)
    print(f"{'Параметр':<25} {'Значение':<20} {'Размерность':<15}")
    print("-" * 80)
    print(f"{'Угол наклона КСУ (θ)':<25} {theta:<20.3f} {'°':<15}")
    print(f"{'Давление за скачком (p2)':<25} {p2:<20.3f} {'Па':<15}")
    print(f"{'Температура за скачком (T2)':<25} {T2:<20.3f} {'К':<15}")
    print(f"{'Плотность за скачком (ρ2)':<25} {rho2:<20.3f} {'кг/м³':<15}")
    print(f"{'Число Маха за скачком (M2)':<25} {M2:<20.3f} {'':<15}")
    print(f"{'Режим течения':<25} {flow_regime:<20} {'':<15}")
    print("=" * 80)

    # Вывод информации о сходимости
    print("\nСХОДИМОСТЬ МЕТОДА НЬЮТОНА:")
    print(f"{'Итерация':<10} {'θ, град':<15} {'f(θ)':<15}")
    print("-" * 40)
    for i, (theta_val, f_val) in enumerate(zip(theta_history, f_history)):
        print(f"{i:<10} {math.degrees(theta_val):<15.3f} {f_val:<15.3e}")

    plot_func(beta_rad, M1, k)
    plot_newton_convergence(theta_history, f_history)

if __name__ == "__main__":
    main()