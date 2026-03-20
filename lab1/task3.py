import math
import matplotlib.pyplot as plt
import time

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

def exist_root(a, b, k, beta, M1):
    # Проверка знаков на концах отрезка
    if func(a, k, beta, M1) * func(b, k, beta, M1) >= 0:
        return False, "Функция не меняет знак на концах отрезка"
    
    # Проверка монотонности (производная не меняет знак)
    if dfunc(a, k, beta, M1) * dfunc(b, k, beta, M1) < 0:
        return False, "Производная меняет знак, возможны несколько корней"
    
    return True, "Корень существует и единственный на данном отрезке"

def newton_method(beta, M1, k, max_iter=100, epsilon=1e-6):
    # Начальное приближение
    theta = beta + math.radians(10)
    
    for i in range(max_iter):
        f = func(theta, k, beta, M1)
        df = dfunc(theta, k, beta, M1)

        converges = newton_convergence(theta, beta, M1, k)
        if not converges:
           return None, i + 1
        
        theta_new = theta - f / df

        if abs(theta_new - theta) < epsilon:
            return theta_new, i + 1
        
        f_new = func(theta_new, k, beta, M1)
        if abs(f_new) < epsilon:
            return theta_new, i + 1
        
        theta = theta_new
    return None, max_iter

def bisection_method_for_theta(k, beta, M1, eps=1e-6, max_iter=100):
    """
    Метод половинного деления для уравнения косого скачка уплотнения
    """
    a = math.radians(30)
    b = math.radians(60)
    
    exists, message = exist_root(a, b, k, beta, M1)
    print(f"Проверка существования корня: {message}")
    if not exists:
        return None, None
    
    iter_count = 0
    
    while abs(b - a) > 2 * eps and iter_count < max_iter:
        mid = (a + b) / 2
        f_mid = func(mid, k, beta, M1)
        
        if abs(f_mid) < eps:
            return mid, iter_count + 1
        
        if func(a, k, beta, M1) * f_mid < 0:
            b = mid
        else:
            a = mid
        
        iter_count += 1
    
    return (a + b) / 2, iter_count + 1

def calc_parameters(theta, beta, M1, p1, T1, R, k):
    rho1 = p1 / (R * T1)
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

def main():
    k = 1.4 # аддиабата
    p1 = 63950
    T1 = 320
    M1 = 2 # число Маха набегающего потока
    beta_grad = 18 # угол клина
    R = 287 # газовая постоянная воздуха

    beta_rad = math.radians(beta_grad)
    
    print("=" * 80)
    print("СРАВНЕНИЕ МЕТОДОВ НЬЮТОНА И ДИХОТОМИИ")
    print("=" * 80)

    print("\n1. МЕТОД НЬЮТОНА:")
    start_time_newton = time.perf_counter() * 1000 
    theta_newton_rad, newton_iter = newton_method(beta_rad, M1, k)
    end_time_newton = time.perf_counter() * 1000 
    newton_time = end_time_newton - start_time_newton
    
    if theta_newton_rad is None:
        print("Метод Ньютона не сошелся")
        return
    
    theta_newton_grad = math.degrees(theta_newton_rad)
    p2_newton, T2_newton, rho2_newton, M2_newton, flow_regime_newton = calc_parameters(
        theta_newton_rad, beta_rad, M1, p1, T1, R, k)
    
    print(f"Угол наклона КСУ: θ = {theta_newton_grad:.6f}°")
    print(f"Количество итераций: {newton_iter}")
    print(f"Время выполнения: {newton_time:.3f} мс")
    print(f"Число Маха за скачком: M2 = {M2_newton:.6f}")
    print(f"Режим течения: {flow_regime_newton}")
    
    print("\n2. МЕТОД ДИХОТОМИИ:")
    
    start_time_bisection = time.perf_counter() * 1000 
    theta_bisection_rad, bisection_iter = bisection_method_for_theta(k, beta_rad, M1)
    end_time_bisection = time.perf_counter() * 1000 
    bisection_time = end_time_bisection - start_time_bisection
    
    if theta_bisection_rad is None:
        print("Метод дихотомии не нашел корень на заданном интервале")
        return
    
    theta_bisection_grad = math.degrees(theta_bisection_rad)
    p2_bisection, T2_bisection, rho2_bisection, M2_bisection, flow_regime_bisection = calc_parameters(
        theta_bisection_rad, beta_rad, M1, p1, T1, R, k)
    
    print(f"Угол наклона КСУ: θ = {theta_bisection_grad:.6f}°")
    print(f"Количество итераций: {bisection_iter}")
    print(f"Время выполнения: {bisection_time:.3f} мс")
    print(f"Число Маха за скачком: M2 = {M2_bisection:.6f}")
    print(f"Режим течения: {flow_regime_bisection}")
    
    print("\n" + "=" * 80)
    print("СРАВНИТЕЛЬНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 80)
    print(f"{'Параметр':<25} {'Метод Ньютона':<20} {'Метод дихотомии':<20} {'Разность':<15}")
    print("-" * 80)
    print(f"{'Угол θ, °':<25} {theta_newton_grad:<20.6f} {theta_bisection_grad:<20.6f} {abs(theta_newton_grad - theta_bisection_grad):<15.6e}")
    print(f"{'Число Маха M2':<25} {M2_newton:<20.6f} {M2_bisection:<20.6f} {abs(M2_newton - M2_bisection):<15.6e}")
    print(f"{'Давление p2, Па':<25} {p2_newton:<20.3f} {p2_bisection:<20.3f} {abs(p2_newton - p2_bisection):<15.3e}")
    print(f"{'Температура T2, K':<25} {T2_newton:<20.3f} {T2_bisection:<20.3f} {abs(T2_newton - T2_bisection):<15.3e}")
    print(f"{'Плотность ρ2, кг/м³':<25} {rho2_newton:<20.3f} {rho2_bisection:<20.3f} {abs(rho2_newton - rho2_bisection):<15.3e}")
    print(f"{'Количество итераций':<25} {newton_iter:<20} {bisection_iter:<20} {abs(newton_iter - bisection_iter):<15}")
    print(f"{'Время выполнения, мс':<25} {newton_time:<20.3f} {bisection_time:<20.3f} {abs(newton_time - bisection_time):<15.3f}")
    print("=" * 80)

if __name__ == "__main__":
    main()