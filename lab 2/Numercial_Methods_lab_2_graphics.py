
# -*- coding: cp1251 -*-
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

n = 10; #������ ������ � ������ ������
cnt = 15; #���������� ������ � ������ ������
cnt2 = 200; #���������� ������ �� ������ ������
cond1 = 10; #����� ���������������� ��� ������ �������
cond2 = 10**3;#����� ���������������� ��� ������ �������
cond3 = 10**6;#����� ���������������� ��� ������ �������
N = 100; #������ ������ � ������� ������

deltas = []
residual = []
for i in range(1, cnt+1):
    with open(f"C:/Users/dzhda/OneDrive/������� ����/�������/��������� ������/3 ���/lab 2/Numercial Methods lab 2 algorithm/1matrix_{i}.txt", 'r') as file:
        for line in file:
            first, second = line.split()
            deltas.append(float(first))
            residual.append(float(second)/10**(i-1))
print(deltas)
print(residual)


duration = []
for i in range(1, cnt2):
    with open(f"C:/Users/dzhda/OneDrive/������� ����/�������/��������� ������/3 ���/lab 2/Numercial Methods lab 2 algorithm/2matrix{i+1}.txt", 'r') as file:
        for line in file:
            duration.append(float(line))

dX1 = [];
dA1 = [];
sup1 = []
with open(f"C:/Users/dzhda/OneDrive/������� ����/�������/��������� ������/3 ���/lab 2/Numercial Methods lab 2 algorithm/3matrix{cond1}.txt", 'r') as file:
        for line in file:
            first, second = line.split()
            dX1.append(float(first));
            dA1.append(float(second));
            sup1.append(cond1 * float(second) / (1 - cond1 * float(second)))
dX2 = [];
dA2 = [];
sup2 = []
with open(f"C:/Users/dzhda/OneDrive/������� ����/�������/��������� ������/3 ���/lab 2/Numercial Methods lab 2 algorithm/3matrix{cond2}.txt", 'r') as file:
        for line in file:
            first, second = line.split()
            dX2.append(float(first));
            dA2.append(float(second));
            sup2.append(cond2 * float(second) / (1 - cond2 * float(second)))
dX3 = [];
dA3 = [];
sup3 = []
with open(f"C:/Users/dzhda/OneDrive/������� ����/�������/��������� ������/3 ���/lab 2/Numercial Methods lab 2 algorithm/3matrix{cond3}.txt", 'r') as file:
        for line in file:
            first, second = line.split()
            dX3.append(float(first));
            dA3.append(float(second));
            sup3.append((cond3 * float(second) / (1 - cond3 * float(second))) + 2)



plt.subplot(2, 1, 1) 
x = [10**(i) for i in range(cnt)]
plt.scatter(x, deltas ,linewidths=0.1, color = 'green', zorder=2, label = "�����������")
plt.scatter(x, residual ,linewidths=0.1, color = 'blue', zorder=2, label = "�������")
plt.xlabel('����� ����������������')
plt.ylabel('�����������')
plt.title(f'����������� ����������� � ������� �� ����� ���������������� (n = {n})')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout() 


plt.subplot(2, 1, 2) 
x = [i+1 for i in range(1, cnt2)]
y = [0.0000000140 * i**3 for i in x]
plt.scatter(x, duration ,linewidths=0.1, color = 'green', zorder=2)
plt.plot(x, y, lw = 2 ,color = 'red')
plt.xlabel('����������� �������')
plt.ylabel('����� ������� (���)')
plt.title('����������� ������� ������� �� ����������� �������')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout() 


plt.figure()
plt.subplot(2, 2, 1) 
plt.scatter(dA1, dX1, s=10, linewidths=0.05, color = 'green', label = cond1)
plt.plot(dA1, sup1, lw = 2 ,color = 'green')
plt.xlabel(r'$\frac{||\delta A||}{||A||}$')
plt.ylabel(r'$\frac{||\delta X||}{||X||}$')
plt.title('')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=8)
plt.tight_layout() 

plt.subplot(2, 2, 2) 
plt.scatter(dA2, dX2, s=10, linewidths=0.05, color = 'red',  label = cond2)
plt.plot(dA2, sup2, lw = 2 ,color = 'red')
plt.xlabel(r'$\frac{||\delta A||}{||A||}$')
plt.ylabel(r'$\frac{||\delta X||}{||X||}$')
plt.title('')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=8)
plt.tight_layout() 

plt.subplot(2, 2, 3) 
plt.scatter(dA3, dX3, s=10, linewidths=0.05, color = 'blue',  label = cond3)
plt.plot(dA3, sup3, lw = 2 ,color = 'blue')
plt.xlabel(r'$\frac{||\delta A||}{||A||}$')
plt.ylabel(r'$\frac{||\delta X||}{||X||}$')
plt.title('')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=8)
plt.tight_layout() 

plt.subplot(2, 2, 4) 
plt.scatter(dA1, dX1, s=10, linewidths=0.05, color = 'green', label = cond1)
plt.plot(dA1, sup1, lw = 2 ,color = 'green')
plt.scatter(dA2, dX2, s=10,linewidths=0.05, color = 'red',  label = cond2)
plt.plot(dA2, sup2, lw = 2 ,color = 'red')
plt.scatter(dA3, dX3, s=10, linewidths=0.05, color = 'blue',  label = cond3)
plt.plot(dA3, sup3, lw = 2 ,color = 'blue')
plt.xlabel(r'$\frac{||\delta A||}{||A||}$')
plt.ylabel(r'$\frac{||\delta X||}{||X||}$')
plt.title('')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=8)
plt.tight_layout() 

plt.show()