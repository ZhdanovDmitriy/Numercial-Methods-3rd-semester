# -*- coding: cp1251 -*-
import matplotlib.pyplot as plt
import numpy as np

func = "x*exp(x) + 0.1"
def function(x):
    return x*np.exp(x) + 0.1;

func2 = "3x^3 + 5x^2 + 5x + 2"
def function2(x):
    return 3*x**3 + 5*x**2 + 5*x + 2;

a = -1.4;
b = 0.2;
ans = -0.111832559158963;
ans2 = -0.66666666666666667;
Eps = 0.00000001
sampling = 40;


ResultFile = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Results.txt")
ErrorRateFile = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\ErrorRate.txt")
CostFile = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Cost.txt")
achievabilityFile = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\achievability.txt")
choiceaprFile = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\choiceapr.txt")
ResultFile2 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Results2.txt")
ErrorRateFile2 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\ErrorRate2.txt")
CostFile2 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Cost2.txt")
achievabilityFile2 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\achievability2.txt")
choiceaprFile2 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\choiceapr2.txt")
ResultFile3 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Results3.txt")
ErrorRateFile3 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\ErrorRate3.txt")
CostFile3 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Cost3.txt")
achievabilityFile3 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\achievability3.txt")
choiceaprFile3 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\choiceapr3.txt")
ResultFile4 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Results4.txt")
ErrorRateFile4 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\ErrorRate4.txt")
CostFile4 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\Cost4.txt")
achievabilityFile4 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\achievability4.txt")
choiceaprFile4 = open(r"C:\Users\dzhda\OneDrive\Рабочий стол\политех\Численные методы\3 сем\lab 1\Numercial Methods lab 1 algorithm\choiceapr4.txt")
Result = [float(line.strip()) for line in ResultFile]
ErrorRate = [abs(float(line.strip())) for line in ErrorRateFile]
Cost = [int(line.strip()) for line in CostFile]
achievability = [abs(float(line.strip())) for line in achievabilityFile]
choiceapr = [int(line.strip()) for line in choiceaprFile]
Result2 = [float(line.strip()) for line in ResultFile2]
ErrorRate2 = [abs(float(line.strip())) for line in ErrorRateFile2]
Cost2 = [int(line.strip()) for line in CostFile2]
achievability2 = [abs(float(line.strip())) for line in achievabilityFile2]
choiceapr2 = [int(line.strip()) for line in choiceaprFile2]
Result3 = [float(line.strip()) for line in ResultFile3]
ErrorRate3 = [abs(float(line.strip())) for line in ErrorRateFile3]
Cost3 = [int(line.strip()) for line in CostFile3]
achievability3 = [abs(float(line.strip())) for line in achievabilityFile3]
choiceapr3 = [int(line.strip()) for line in choiceaprFile3]
Result4 = [float(line.strip()) for line in ResultFile4]
ErrorRate4 = [abs(float(line.strip())) for line in ErrorRateFile4]
Cost4 = [int(line.strip()) for line in CostFile4]
achievability4 = [abs(float(line.strip())) for line in achievabilityFile4]
choiceapr4 = [int(line.strip()) for line in choiceaprFile4]


print(Result[-1])
print(Result)
print(ErrorRate)
print(Cost)
print(achievability)
print(choiceapr)
print("=======================================")
print(Result2[-1])
print(Result2)
print(ErrorRate2)
print(achievability2)
print(choiceapr2)
print("=======================================")
print(Result3[-1])
print(Result3)
print(ErrorRate3)
print(achievability3)
print(choiceapr3)
print("=======================================")
print(Result4[-1])
print(Result4)
print(ErrorRate4)
print(achievability4)
print(choiceapr4)




plt.subplot(2, 2, 1) #строка, столбец, номер графика
x = np.arange(a, b, 0.01)
y = function(x);
plt.plot(x, y, label= func, color = 'green')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'График функции {func} на интервале [{a}, {b}]')
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.axvline(a, color='blue', linewidth=1, label='Границы отрезка')  
plt.axvline(b, color='blue', linewidth=1)
plt.grid(True)  # Включаем сетку
plt.scatter(ans, 0, marker = '.', linewidths=0.1, color='black', label='Истинное значение', zorder=5)
plt.scatter(Result[-1], function(Result[-1]), marker = '.', linewidths=0.1, color='red', label='Вычисленное значение МПД', zorder=5)
plt.scatter(Result2[-1], function(Result2[-1]), marker = '.', linewidths=0.1, color='blue', label='Вычисленное значение метод Эйткена', zorder=5)
plt.legend(fontsize=8)

plt.subplot(2, 2, 2) 
plt.plot(list(range(0, len(Result))), Result, marker = '.', color = 'green', label='МПД')
plt.plot(list(range(0, len(Result2))), Result2, marker = '.', color = 'orange', label='метод Эйткена')
plt.axhline(ans, color='red', label='Искомое значение', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=1)  # Линия x=0
plt.xlabel('Номер итерации')
plt.ylabel('Промежуточное значение')
plt.title('Сходимость')
plt.grid(True)
plt.legend(fontsize=8)

plt.subplot(2, 2, 3) 
plt.plot(list(range(0, len(Result))), ErrorRate, marker = '.', color = 'green', label='МПД')
plt.plot(list(range(0, len(Result2))), ErrorRate2, marker = '.', color = 'orange', label='метод Эйткена')
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.xlabel('Номер итерации')
plt.ylabel('Погрешность')
plt.yscale('log')
plt.title(f'Погрешность от номера итерации (E = {Eps})')
plt.grid(True)
plt.legend(fontsize=8)

plt.subplot(2, 2, 4)
x_values = [10**(-i) for i in range(len(Cost))]
plt.plot(x_values, Cost , marker = '.', color = 'green', label='МПД')
plt.plot(x_values, Cost2 , marker = '.', color = 'orange', label='метод Эйткена')
plt.xlabel('Требуемая точность')
plt.ylabel('Количество итераций')
plt.title('Вычислительные затраты')
plt.grid(True)
plt.xscale('log')
plt.legend(fontsize=8)
plt.tight_layout() #чтобы графики не накладывались




plt.figure()
plt.subplot(2, 2, 1) #строка, столбец, номер графика
x = np.arange(a, b, 0.01)
y = function2(x);
plt.plot(x, y, label= func2, color = 'green')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'График функции {func2} на интервале [{a}, {b}]')
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.axvline(a, color='blue', linewidth=1, label='Границы отрезка')  
plt.axvline(b, color='blue', linewidth=1)  
plt.axvline(-2.2667, color='orange', linewidth=1, label='Границы корня')  
plt.axvline(-0.286, color='orange', linewidth=1)  
plt.grid(True)  # Включаем сетку
plt.scatter(ans2, 0, marker = '.', linewidths=0.1, color='black', label='Истинное значение', zorder=5)
plt.scatter(Result3[-1], function2(Result3[-1]), marker = '.', linewidths=0.1, color='red', label='Вычисленное значение МПД', zorder=5)
plt.scatter(Result4[-1], function2(Result4[-1]), marker = '.', linewidths=0.1, color='blue', label='Вычисленное значение метод Эйткена', zorder=5)
plt.legend(fontsize=8)

plt.subplot(2, 2, 2) 
plt.plot(list(range(0, len(Result3))), Result3, marker = '.', color = 'green', label='МПД')
plt.plot(list(range(0, len(Result4))), Result4, marker = '.', color = 'orange', label='метод Эйткена')
plt.axhline(ans2, color='red', label='Искомое значение', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=1)  # Линия x=0
plt.xlabel('Номер итерации')
plt.ylabel('Промежуточное значение')
plt.title('Сходимость')
plt.grid(True)
plt.legend(fontsize=8)

plt.subplot(2, 2, 3) 
plt.plot(list(range(0, len(Result3))), ErrorRate3, marker = '.', color = 'green', label='МПД')
plt.plot(list(range(0, len(Result4))), ErrorRate4, marker = '.', color = 'orange', label='метод Эйткена')
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.xlabel('Номер итерации')
plt.ylabel('Погрешность')
plt.yscale('log')
plt.title(f'Погрешность от номера итерации (E = {Eps})')
plt.grid(True)
plt.legend(fontsize=8)

plt.subplot(2, 2, 4)
x_values = [10**(-i) for i in range(len(Cost3))]
plt.plot(x_values, Cost3 , marker = '.', color = 'green', label='МПД')
plt.plot(x_values, Cost4 , marker = '.', color = 'orange', label='метод Эйткена')
plt.xlabel('Требуемая точность')
plt.ylabel('Количество итераций')
plt.title('Вычислительные затраты')
plt.grid(True)
plt.xscale('log')
plt.legend(fontsize=8)
plt.tight_layout() #чтобы графики не накладывались


plt.figure()
plt.subplot(2, 2, 1)
x_values = [10**(-i) for i in range(len(achievability))]
plt.scatter(x_values, achievability ,  color = 'green', label='МПД')
plt.scatter(x_values, achievability2 , color = 'orange', label='метод Эйткена')
plt.plot(x_values, x_values , color = 'red', label='Верхний предел погрешности')
plt.xlabel('Требуемая точность')
plt.ylabel('Погрешность')
plt.title(f'Достижимость точности для {func}')
plt.grid(True)
plt.legend(fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout() 

plt.subplot(2, 2, 2)
x_values = [10**(-i) for i in range(len(achievability3))]
plt.scatter(x_values, achievability3 ,linewidths=0.1, color = 'green', label='МПД', zorder=2)
plt.scatter(x_values, achievability4 ,linewidths=0.1, color = 'orange', label='метод Эйткена', zorder=2)
plt.plot(x_values, x_values , color = 'red', label='Верхний предел погрешности')
plt.xlabel('Требуемая точность')
plt.ylabel('Погрешность')
plt.title(f'Достижимость точности для {func2}')
plt.grid(True)
plt.legend(fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout() 

plt.subplot(2, 2, 3)
x_values = []
x_values2 = []
for i in range(sampling + 1):
    add = a + i*(b-a)/sampling
    if(add <= ans):
        x_values.append(add)
for i in range(sampling + 1):
        x_values2.append(a + i*(b-a)/sampling)
plt.axvline(ans, color='red', linewidth=1, label = 'корень')
plt.axvline(a, color='blue', linewidth=1, label = 'границы отрезка')
plt.axvline(b, color='blue', linewidth=1,)
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.scatter(x_values, choiceapr ,  color = 'green', label='МПД')
plt.scatter(x_values2, choiceapr2 , color = 'orange', label='метод Эйткена')
plt.xlabel('Начальное приближение')
plt.ylabel('Количество итераций')
plt.title(f'Количество итераций от приближения {func}')
plt.grid(True)
plt.legend(fontsize=8)
plt.tight_layout() 

plt.subplot(2, 2, 4)
x_values = []
for i in range(sampling + 1):
    add = a + i*(b-a)/sampling
    if(add <= ans2):
        x_values.append(add)
plt.axvline(ans2, color='red', linewidth=1, label = 'корень')
plt.axvline(a, color='blue', linewidth=1, label = 'границы отрезка')
plt.axvline(b, color='blue', linewidth=1,)
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.scatter(x_values, choiceapr3 ,  color = 'green', label='МПД')
plt.scatter(x_values2, choiceapr4 , color = 'orange', label='метод Эйткена')
plt.xlabel('Начальное приближение')
plt.ylabel('Количество итераций')
plt.title(f'Количество итераций от приближения {func2}')
plt.grid(True)
plt.legend(fontsize=8)
plt.tight_layout() 

plt.show(block=True)




















"""
import sympy as sp
X = sp.Symbol('x')
equation = X ** (sp.sin(X) - 6) -2;
ans = sp.solveset(equation, X, domain=sp.Interval(a, b, left_open=False, right_open=False))
print("Roots:", ans)





plt.figure()  # Создаём второе окно
plt.subplot(1, 1, 1) #строка, столбец, номер графика
x_values = [10**i for i in range(len(ErrorRate))]
plt.plot(x_values, ErrorRate , color = 'green', label='МПД')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'График функции {func} на интервале [{a}, {b}]')
plt.axhline(0, color='black', linewidth=0.5)  # Линия y=0
plt.axvline(0, color='black', linewidth=0.5)  # Линия x=0
plt.grid(True)  # Включаем сетку
plt.scatter(ans, 0, marker = '.', linewidths=0.1, color='blue', label='Истинное значение', zorder=5)
plt.scatter(Result2[-1], function(Result2[-1]), marker = '.', linewidths=0.1, color='red', label='Вычисленное значение', zorder=5)
plt.legend(fontsize=8)
"""