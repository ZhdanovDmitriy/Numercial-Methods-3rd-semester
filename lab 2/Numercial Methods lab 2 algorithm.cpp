#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <iomanip>  // для std::setprecision

using namespace std;
std::string basePath = "C:/Users/dzhda/OneDrive/Рабочий стол/политех/матлаб/";
int cnt = 15;
int cnt2 = 200;
int cnt3 = 100;
int N = 200;

void printVector(const vector<long double>& vec) {
    for (long double val : vec) {
        cout << val << "\t";
    }
    cout << endl;
}

void printMatrix(const vector<vector<long double>>& matrix) {
    for (const auto& row : matrix) {
        for (long double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

bool LU_Decomposition(const vector<vector<long double>>& A, vector<vector<long double>>& L, vector<vector<long double>>& U) {
    int n = A.size();
    L = vector<vector<long double>>(n, vector<long double>(n, 0));
    U = vector<vector<long double>>(n, vector<long double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            long double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k) {
                L[i][i] = 1;
            }
            else {
                long double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                if (U[i][i] == 0) {
                    cerr << "Ошибка: матрица не может быть разложена на LU (деление на ноль)." << endl;
                    return false;
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
    return true;
}

// Функция прямой подстановки для решения системы Ly = b
vector<long double> forwardSubstitution(const vector<vector<long double>>& L, const vector<long double>& b) {
    int n = L.size();
    vector<long double> y(n, 0);

    for (int i = 0; i < n; i++) {
        long double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    return y;
}

// Функция обратной подстановки для решения системы Ux = y
vector<long double> backSubstitution(const vector<vector<long double>>& U, const vector<long double>& y) {
    int n = U.size();
    vector<long double> x(n, 0);

    for (int i = n - 1; i >= 0; i--) {
        long double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

// Основная функция для решения СЛАУ Ax = b с использованием LU-разложения
vector<long double> solveLU(const vector<vector<long double>>& A, const vector<long double>& b) {
    int n = A.size();
    vector<vector<long double>> L, U;
    if (!LU_Decomposition(A, L, U)) {
        cerr << "Ошибка при выполнении LU-разложения." << endl;
        return {};
    }
    vector<long double> y = forwardSubstitution(L, b);
    vector<long double> x = backSubstitution(U, y);
    return x;
}

vector<long double> multiplyMatrixVector(const vector<vector<long double>>& A, const vector<long double>& x) {
    int n = A.size();
    vector<long double> result(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

vector<long double> computeResidual(const vector<vector<long double>>& A, const vector<long double>& x, const vector<long double>& b) {
    int n = A.size();

    // Вычисляем Ax
    vector<long double> Ax = multiplyMatrixVector(A, x);

    //Вычисляем r = Ax - b
    vector<long double> r(n, 0);
    for (int i = 0; i < n; i++) {
        r[i] = Ax[i] - b[i];
    }
    return r;
}
/*
long double norm(const std::vector<long double>& vec) {
    long double maxNorm = 0.;
    for (long double val : vec) {
        if (std::fabs(val) > maxNorm) {
            maxNorm = std::fabs(val);
        }
    }
    return maxNorm;
}
*/

long double norm(const std::vector<vector<long double>>& matrix) {
    long double norm = 0.;
    for (int i = 0; i != matrix.size(); ++i) {
        long double currentNorm = 0;
        for (int j = 0; j != matrix[0].size(); ++j) {
            currentNorm += std::fabs(matrix[i][j]);
        }
        if (currentNorm > norm) {
            norm = currentNorm;
        }
    }
    return norm;
}

long double norm(const std::vector<long double>& vec) {
    long double cur = 0.;
    for (long double val : vec) {
        cur += val * val;
    }
    cur = sqrt(cur);
    return cur;
}

std::vector<std::vector<long double>> delta(std::vector<std::vector<long double>>& matrix1, std::vector<std::vector<long double>>& matrix2) {
    std::vector<std::vector<long double>> delta;
    delta.resize(matrix1.size());
    for (int i = 0; i != matrix1.size(); ++i) {
        delta[i].resize(matrix1[i].size());
        for (int j = 0; j != matrix1[i].size(); ++j) {
            delta[i][j] = std::fabs(matrix1[i][j] - matrix2[i][j]);
        }
    }
    return delta;
}
std::vector<long double> delta(std::vector<long double>& vec1, std::vector<long double>& vec2) {
    std::vector<long double> delta;
    delta.resize(vec1.size());
    for (int i = 0; i != vec1.size(); ++i) {
        delta[i] = std::fabs(vec1[i] - vec2[i]);
    }
    return delta;
}

bool readDataFromFileSolveWriteToFile(const std::string& filename, const std::string& filename2) {
    ifstream file(filename);
    ofstream file2(filename2);
    vector<vector<long double>> A;
    vector<long double> b;
    vector<long double> ans;

    if (!file.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }

    int n;
    string line;
    file >> n;
    // Пропускаем строку заголовка "Matrix A:"
    getline(file, line); // для завершения строки после n
    getline(file, line);
    // Инициализируем матрицу A
    A.resize(n, std::vector<long double>(n));
    // Считываем матрицу A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    // Пропускаем строку заголовка "Vector b:"
    getline(file, line); // для завершения строки после матрицы
    getline(file, line);
    // Инициализируем вектор b
    b.resize(n);
    // Считываем вектор b
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    // Пропускаем строку заголовка "Solution x:"
    getline(file, line); // для завершения строки после вектора b
    getline(file, line);
    // Инициализируем вектор x
    ans.resize(n);
    // Считываем решение x
    for (int i = 0; i < n; ++i) {
        file >> ans[i];
    }

    vector<long double> x = solveLU(A, b);
    long double residual = norm(computeResidual(A, x, b)) / norm(b);
    long double error = norm(delta(x, ans)) / norm(ans);
    file2 << error << " " << residual << std::endl;

    file.close();
    file2.close();
    return true;
}

bool readDataFromFileSolveWriteTimeToFile(const std::string& filename, const std::string& filename2) {
    ifstream file(filename);
    ofstream file2(filename2);
    vector<vector<long double>> A;
    vector<long double> b;
    vector<long double> ans;

    if (!file.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }

    int n;
    string line;
    file >> n;
    // Пропускаем строку заголовка "Matrix A:"
    getline(file, line); // для завершения строки после n
    getline(file, line);
    // Инициализируем матрицу A
    A.resize(n, std::vector<long double>(n));
    // Считываем матрицу A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    // Пропускаем строку заголовка "Vector b:"
    getline(file, line); // для завершения строки после матрицы
    getline(file, line);
    // Инициализируем вектор b
    b.resize(n);
    // Считываем вектор b
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    // Пропускаем строку заголовка "Solution x:"
    getline(file, line); // для завершения строки после вектора b
    getline(file, line);
    // Инициализируем вектор x
    ans.resize(n);
    // Считываем решение x
    for (int i = 0; i < n; ++i) {
        file >> ans[i];
    }

    auto start = std::chrono::high_resolution_clock::now();
    vector<long double> x = solveLU(A, b);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<long double> duration = end - start;
    file2 << duration.count() << std::endl;

    file.close();
    file2.close();
    return true;
}

bool calculatingNoiseEffect(const std::string& filename, const std::string& filename2,long double noise) {
    ifstream file(filename);
    ofstream file2(filename2);
    if (!file.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        cerr << "Oppening error: " << filename << std::endl;
        return false;
    }

    vector<vector<long double>> A;
    vector<vector<long double>> A_;
    vector<long double> b;
    vector<long double> ans;
    int n;
    string line;
    file >> n;
    getline(file, line);
    getline(file, line);
    A.resize(n, std::vector<long double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    A_.resize(n, std::vector<long double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_[i][j] = A[i][j];
        }
    }
    getline(file, line);
    getline(file, line);
    b.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    getline(file, line);
    getline(file, line);
    ans.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> ans[i];
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1. - noise, 1. + noise);

    for (int iter = 0; iter != N; iter++) {
        for (int i = 0; i != A.size(); i++) {
            A[i][0] = A_[i][0] * dis(gen);
        }
        vector<long double> x = solveLU(A, b);
        long double aboutX = norm(delta(x, ans)) / norm(ans);
        long double aboutA = norm(delta(A, A_)) / norm(A_);
        file2 << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        file2 << aboutX << " " << aboutA << std::endl;
    }

    file.close();
    file2.close();
    return true;

}

int main() {
    for (int i = 1; i <= cnt; ++i) {
            std::string fileName = "1matrix_" + std::to_string(i) + ".txt";
            readDataFromFileSolveWriteToFile(basePath + fileName, fileName);
    }
    for (int i = 2; i <= cnt2; ++i) {
        std::string fileName = "2matrix" + std::to_string(i) + ".txt";
        readDataFromFileSolveWriteTimeToFile(basePath + fileName, fileName);
    }
    calculatingNoiseEffect(basePath + "3matrix10.txt", "3matrix10.txt", 0.001);
    calculatingNoiseEffect(basePath + "3matrix1000.txt", "3matrix1000.txt", 0.001);
    calculatingNoiseEffect(basePath + "3matrix1000000.txt", "3matrix1000000.txt", 0.001);
 
    return 0;
}






/*
bool calculatingNoiseEffect(const std::string& filename, const std::string& filename2, long double noise) {
    ifstream file(filename);
    ofstream file2(filename2);
    if (!file.is_open()) {
        cerr << "Opening error: " << filename << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        cerr << "Opening error: " << filename2 << std::endl;
        return false;
    }

    vector<vector<long double>> A;
    vector<vector<long double>> A_;
    vector<long double> b;
    vector<long double> ans;
    int n;
    string line;

    // Чтение размера матрицы
    file >> n;
    getline(file, line);
    getline(file, line);

    // Чтение матрицы A
    A.resize(n, vector<long double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }

    // Копия A для внесения шума
    A_ = A;

    // Чтение вектора b
    getline(file, line);
    getline(file, line);
    b.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }

    // Чтение правильного ответа
    getline(file, line);
    getline(file, line);
    ans.resize(n);
    for (int i = 0; i < n; ++i) {
        file >> ans[i];
    }

    // Печать исходного уравнения
    std::cout << "Matrix A:" << std::endl;
    for (const auto& row : A) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Vector b: ";
    for (const auto& val : b) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Correct solution (ans): ";
    for (const auto& val : ans) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Генерация случайного шума
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1 - noise, 1 + noise);

    for (int iter = 0; iter != N; iter++) {
        // Решение системы с текущей матрицей A
        vector<long double> x = solveLU(A, b);

        // Печать промежуточных решений
        std::cout << "Iteration " << iter + 1 << ": Solution vector x: ";
        for (const auto& val : x) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Рассчет относительных ошибок
        long double aboutX = (norm(delta(x, ans)) / norm(ans));
        long double aboutA = (norm(delta(A, A_)) / norm(A));
        file2 << aboutX << " " << aboutA << std::endl;

        // Печать относительных ошибок
        std::cout << "Relative error in x (||delta(x)|| / ||x||): " << aboutX << std::endl;
        std::cout << "Relative error in A (||delta(A)|| / ||A||): " << aboutA << std::endl;

        // Применение шума к матрице A
        std::cout << "Matrix A with noise:" << std::endl;
        for (int i = 0; i < A.size(); i++) {
            for (int j = 0; j < A[0].size(); j++) {
                A[i][j] = A_[i][j] * dis(gen);
                std::cout << A[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    file.close();
    file2.close();
    return true;
}
*/
