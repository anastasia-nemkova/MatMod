import pandas as pd
import matplotlib.pyplot as plt

# Чтение данных с указанием разделителя табуляции
df = pd.read_csv(r".\residuals.csv", sep=';')

# Проверьте, что данные загрузились правильно
print(df.head())
print(df.columns)

# Минимальный код для построения
plt.figure(figsize=(10, 6))
plt.semilogy(df['iteration'], df['EPT'])
plt.xlabel('Итерация')
plt.ylabel('Невязка EPT')
plt.title('График сходимости')
plt.grid(True)
plt.show()