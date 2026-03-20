import pandas as pd
import matplotlib.pyplot as plt

# Простое чтение данных
df = pd.read_csv(r".\grafic\residuals.csv")

# Минимальный код для построения
plt.figure(figsize=(10, 6))
plt.semilogy(df['iteration'], df['EPT'])
plt.xlabel('Итерация')
plt.ylabel('Невязка EPT')
plt.title('График сходимости')
plt.grid(True)
plt.show()