
# -*- coding: cp1251 -*-
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

n = 100;

errorFromNum = []
with open(r"C:\Users\dzhda\OneDrive\������� ����\�������\��������� ������\3 ���\lab 3\Numercial Methods lab 4 algorithm\ErrorFromIterNum.txt", 'r') as file:
    for line in file:
        errorFromNum.append(float(line))

norms = []
errorFromStartPosition = []
with open(r"C:\Users\dzhda\OneDrive\������� ����\�������\��������� ������\3 ���\lab 3\Numercial Methods lab 4 algorithm\ErrorFromStartPosition.txt", 'r') as file:
    for line in file:
            first, second = line.split()
            norms.append(float(first))
            errorFromStartPosition.append(float(second))



plt.subplot(2, 1, 1) 
plt.plot([i+1 for i in range(len(errorFromNum))], errorFromNum, marker = '.',lw = 2 ,color = 'green')
plt.xlabel('����� ��������')
plt.ylabel('�����������')
plt.title(f'����������� ���������� ����������� �� ������ ��������(n = 100, cond = 5, Eps = 10E-15)')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()

plt.subplot(2, 1, 2) 
plt.plot(norms, errorFromStartPosition, marker = '.', lw = 2 ,color = 'green')
plt.xlabel('����� ������� �� x0 �� ans')
plt.ylabel('����� �������� ��� �������')
plt.title(f'����������� ���������� �������� �� ����� ������� ���������� ����������� ���������� �����������(n = 100, cond = 5, Eps = 10E-15)')
plt.grid(True)
plt.tight_layout() 
plt.show()