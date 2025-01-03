# -*- coding: cp1251 -*-
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


errorFromItterNum1 = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromIterNum_gap_1.txt", 'r') as file:
    for line in file:
        errorFromItterNum1.append(float(line))

errorFromItterNum2 = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromIterNum_gap_2.txt", 'r') as file:
    for line in file:
        errorFromItterNum2.append(float(line))

errorFromItterNum3 = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromIterNum_gap_3.txt", 'r') as file:
    for line in file:
        errorFromItterNum3.append(float(line))

errorFromItterNum4 = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromIterNum_gap_4.txt", 'r') as file:
    for line in file:
        errorFromItterNum4.append(float(line))
errorFromItterNum5 = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromIterNum_gap_5.txt", 'r') as file:
    for line in file:
        errorFromItterNum5.append(float(line))
errorFromEps = []
with open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 4\Numercial Methods lab 6 matrix generation\Numercial Methods lab 6 matrix generation\matrices_output_modified\ErrorFromEps_gap_5.txt", 'r') as file:
    for line in file:
        errorFromEps.append(float(line))



plt.subplot(2, 1, 1) 

plt.plot([i for i in range(len(errorFromItterNum1))], errorFromItterNum1, marker = '.', lw = 2 ,color = 'red', label = "Отделимость = 0.9756")
plt.plot([i for i in range(len(errorFromItterNum2))], errorFromItterNum2, marker = '.', lw = 2 ,color = 'orange', label = "Отделимость = 0.9523")
plt.plot([i for i in range(len(errorFromItterNum3))], errorFromItterNum3, marker = '.', lw = 2 ,color = 'cyan', label = "Отделимость = 0.9302")
plt.plot([i for i in range(len(errorFromItterNum4))], errorFromItterNum4, marker = '.', lw = 2 ,color = 'blue', label = "Отделимость = 0.9090")
plt.plot([i for i in range(len(errorFromItterNum5))], errorFromItterNum5, marker = '.', lw = 2 ,color = 'green', label = "Отделимость = 0.8889")

plt.xlabel('Номер итерации')
plt.ylabel('Абсолютная погрешность')
plt.title(f'Зависимость абсолютной погрешности от номера итераций(n = 11, Eps = 10E-10) для разной отделимости 7 и 8 СЧ')
plt.grid(True)
plt.legend(fontsize=8)
plt.tight_layout()
plt.xscale('log')
plt.yscale('log')

plt.subplot(2, 1, 2) 
plt.plot([0.1**(i+1) for i in range(len(errorFromEps))], errorFromEps, marker = '.', lw = 2 ,color = 'blue')
x_values = [0.1**(i+1) for i in range(len(errorFromEps))]
plt.plot(x_values, x_values , color = 'red', label='Верхний предел погрешности')
plt.xlabel('Заданная точность')
plt.ylabel('Абсолютная погрешность')
plt.title(f'Достижимость заданной точности для 7 и 8 СЧ с отделимостью 0.8889 (n = 11, Eps = 10E-10)')
plt.grid(True)
plt.tight_layout()
plt.legend(fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.show()