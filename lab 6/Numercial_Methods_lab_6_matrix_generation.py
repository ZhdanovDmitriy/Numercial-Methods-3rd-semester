# -*- coding: utf-8 -*-
import numpy as np
import os

output_dir = "matrices_output_modified"
os.makedirs(output_dir, exist_ok=True)

def generate_eigenvalues(base_distance, gap_between_4th_and_5th, size):
    eigenvalues = np.arange(10, 10 + base_distance * size, base_distance, dtype=float)
    eigenvalues[4] = eigenvalues[3] + gap_between_4th_and_5th
    return eigenvalues[:size]

base_distance = 10
N = 11
Matrixes = []
gap_values = [5, 4, 3, 2, 1]

np.random.seed(N)
M = np.random.rand(N, N)

for idx, gap in enumerate(gap_values):
    eigenvalues = generate_eigenvalues(base_distance, gap, N)

    # Создание диагональной матрицы собственных значений
    D = np.diag(eigenvalues)
    # Генерация ортогональной матрицы Q с помощью разложения QR
    Q, R = np.linalg.qr(M)
    # Создание матрицы A
    A = (Q.dot(D)).dot(Q.T)
    Matrixes.append(A)

    # Сортировка собственных значений (на всякий случай)
    sorted_indices = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[sorted_indices]

    filename = os.path.join(output_dir, f"matrix_data_gap_{gap}.txt")
    with open(filename, "w") as f:
        f.write(f"Matrix Index: {idx + 1}\n")
        f.write(f"Otd_number: {gap}\n")
        f.write(f"Matrix:\n")
        for row in A:
            f.write(", ".join(map(lambda val: f"{val:.16e}", row)) + "\n")
        f.write(f"Eigenvalues:\n")
        f.write("\n".join(map(lambda ev: f"{ev:.16e}", eigenvalues)) + "\n")

print(f"All matrices and related data have been saved to the '{output_dir}' directory.")
