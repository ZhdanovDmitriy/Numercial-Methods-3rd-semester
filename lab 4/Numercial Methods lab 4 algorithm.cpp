#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <utility>
#include <iomanip>  // для std::setprecision

using namespace std;
long double Eps = 0.00000000000001;
int steps = 40;
long double condA = 5.0;
std::string basePath = "C:/Users/dzhda/OneDrive/Рабочий стол/политех/матлаб/lab4matrix5.txt";
//std::string basePath = "C:/Users/dzhda/OneDrive/Рабочий стол/политех/матлаб/lab4matrix_test.txt";
std::string fileName1 = "ErrorFromIterNum.txt";
std::string fileName2 = "ErrorFromStartPosition.txt";

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

std::vector<long double> multiplyVectorScalar(const std::vector<long double>& vec, long double scalar) {
    std::vector<long double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

std::vector<long double> addVectors(const std::vector<long double>& vec1, const std::vector<long double>& vec2) {
    if (vec1.size() != vec2.size()) {
        cout << "Vectors must have the same size";
    }
    std::vector<long double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

const std::vector<long double> subtractVectors(const std::vector<long double>& vec1,const std::vector<long double>& vec2) {
    std::vector<long double> delta;
    delta.resize(vec1.size());
    for (int i = 0; i != vec1.size(); ++i) {
        delta[i] = vec1[i] - vec2[i];
    }
    return delta;
}


long double dotVecs(const std::vector<long double>& vec1, const std::vector<long double>& vec2) {
    long double result = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}


std::vector<long double> multiplyMatrixVector(const vector<vector<long double>>& A, const vector<long double>& x) {
    int n = A.size();
    vector<long double> result(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

std::vector<long double> x_next(const std::vector<long double>& x, long double alpha, const std::vector<long double>& p) {
    return addVectors(x, multiplyVectorScalar(p, alpha));
}

long double alpha_k(const std::vector<long double>& r, const std::vector<long double>& p, const std::vector<vector<long double>>& A) {
    return dotVecs(r, p) / dotVecs(multiplyMatrixVector(A,p), p);
}

std::vector<long double> p_k(const std::vector<long double>& r, const std::vector<long double>& p_previous, long double betta_previous) {
    return subtractVectors(r, multiplyVectorScalar(p_previous, betta_previous));
}


long double betta_previous(const std::vector<long double>& r, const std::vector<long double>& p_previous, const std::vector<vector<long double>>& A) {
    return dotVecs(r, multiplyMatrixVector(A, p_previous)) / dotVecs(p_previous, multiplyMatrixVector(A, p_previous));
}

vector<long double> computeResidual(const vector<vector<long double>>& A, const vector<long double>& x, const vector<long double>& b) {
    int n = A.size();

    vector<long double> Ax = multiplyMatrixVector(A, x);

    vector<long double> r(n, 0);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - Ax[i];
    }
    return r;
}

vector<long double> copyVec(const vector<long double>& a) {
    vector<long double> r(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        r[i] = a[i];
    }
    return r;
}

long double norm(const std::vector<long double>& vec) {
    long double maxNorm = vec[0];
    for (long double val : vec) {
        if (std::fabs(val) > maxNorm) {
            maxNorm = std::fabs(val);
        }
    }
    return maxNorm;
}

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

std::vector<long double> delta(const std::vector<long double>& vec1,const std::vector<long double>& vec2) {
    std::vector<long double> delta;
    delta.resize(vec1.size());
    for (int i = 0; i != vec1.size(); ++i) {
        delta[i] = std::fabs(vec1[i] - vec2[i]);
    }
    return delta;
}


void solveConjugateGradientMethod(const vector<vector<long double>>& A, const vector<long double>& x0, const vector<long double>& b, const vector<long double>& ans, const std::string& fileName) {
    ofstream file2(fileName);

    if (!file2.is_open()) {
        cerr << "Oppening error: " << fileName << std::endl;
    }

    vector<long double> p0 = computeResidual(A, x0, b);
    vector<long double> r0 = copyVec(p0);
    long double a_k = alpha_k(r0, p0, A);
    vector<long double> x_new = x_next(x0, a_k, p0);
    file2 << norm(delta(x_new, ans)) << std::endl;
    vector<long double> r_new = subtractVectors(r0, multiplyVectorScalar(multiplyMatrixVector(A, p0), a_k));
    long double betta_prev = betta_previous(r_new, p0, A);
    vector<long double> p_new = p_k(r_new, p0, betta_prev);
    vector<long double> x_prev = copyVec(x_new);

    do {
        a_k = alpha_k(r_new, p_new, A);
        x_prev = copyVec(x_new);
        x_new = x_next(x_new, a_k, p_new);
        file2 << norm(delta(x_new, ans)) << std::endl;
        r_new = subtractVectors(r_new, multiplyVectorScalar(multiplyMatrixVector(A, p_new), a_k));
        betta_prev = betta_previous(r_new, p_new, A);
        p_new = p_k(r_new, p_new, betta_prev);
    } while ((condA * norm(computeResidual(A, x_new, b)) / norm(b)) >= Eps);

    file2.close();
}

pair<vector<long double>, int> solveConjugateGradientMethod(const vector<vector<long double>>& A, const vector<long double>& x0, const vector<long double>& b){
    vector<long double> p0 = computeResidual(A, x0, b);
    vector<long double> r0 = copyVec(p0);
    long double a_k = alpha_k(r0, p0, A);
    vector<long double> x_new = x_next(x0, a_k, p0);
    vector<long double> r_new = subtractVectors(r0, multiplyVectorScalar(multiplyMatrixVector(A, p0), a_k));
    long double betta_prev = betta_previous(r_new, p0, A);
    vector<long double> p_new = p_k(r_new, p0, betta_prev);
    int cnt = 1;
    vector<long double> x_prev = copyVec(x_new);

    do{
        a_k = alpha_k(r_new, p_new, A);
        x_prev = copyVec(x_new);
        x_new = x_next(x_new, a_k, p_new);
        r_new = subtractVectors(r_new, multiplyVectorScalar(multiplyMatrixVector(A, p_new), a_k));
        betta_prev = betta_previous(r_new, p_new, A);
        p_new = p_k(r_new, p_new, betta_prev);
        ++cnt;
    } while ((norm(computeResidual(A, x_new, b)) / norm(b)) >= Eps);

    return pair<vector<long double>, int> (x_new,cnt);
}

void ErrorFromStartPosition(const std::string& basePath, const std::string& fileName){
    ifstream file(basePath);
    ofstream file2(fileName);
    vector<vector<long double>> A;
    vector<long double> b;
    vector<long double> ans;

    if (!file.is_open()) {
        cerr << "Oppening error: " << basePath << std::endl;
    }
    if (!file2.is_open()) {
        cerr << "Oppening error: " << fileName << std::endl;
    }

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
    cout << "\n";
    vector<long double> x0(ans.size(), 0);
    for (int i = 0; i != steps; i++) {
        for (int j = 0; j != ans.size(); j++) {
            x0[j] = ans[j] / steps * (i+1);
        }
        pair<vector<long double>, int> solution = solveConjugateGradientMethod(A, x0, b);
        file2 << norm(delta(ans, x0)) << ' ' << solution.second << std::endl;
    }
    /*
    for (int i = 0; i != steps; i++) {
        for (int j = 0; j != ans.size(); j++) {
            x0[j] = ans[j] + ans[j] / steps * i;
        }
        pair<vector<long double>, int> solution = solveConjugateGradientMethod(A, x0, b);
        file2 << norm(delta(ans, x0)) << ' ' << solution.second << std::endl;
    }*/
    file.close();
    file2.close();
}

void ErrorFromIterNum(const std::string& basePath, const std::string& fileName) {
    ifstream file(basePath);
    ofstream file2(fileName);
    vector<vector<long double>> A;
    vector<long double> b;
    vector<long double> ans;

    if (!file.is_open()) {
        cerr << "Oppening error: " << basePath << std::endl;
    }
    if (!file2.is_open()) {
        cerr << "Oppening error: " << fileName << std::endl;
    }

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
    vector<long double> x0(b.size(), 0);
    solveConjugateGradientMethod(A, x0, b, ans, fileName);
    file.close();
}

int main() {
    ErrorFromIterNum(basePath, fileName1);
    ErrorFromStartPosition(basePath, fileName2);
    return 0;
}



































/*
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

    vector<long double> x = solveConjugateGradientMethod(A, b);
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
    vector<long double> x = solveConjugateGradientMethod(A, b);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<long double> duration = end - start;
    file2 << duration.count() << std::endl;

    file.close();
    file2.close();
    return true;
}

bool calculatingNoiseEffect(const std::string& filename, const std::string& filename2, long double noise) {
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
        vector<long double> x = solveConjugateGradientMethod(A, b);
        long double aboutX = norm(delta(x, ans)) / norm(ans);
        //long double aboutA = norm(delta(A, A_)) / norm(A_);
        long double aboutA = 0;
        file2 << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
        file2 << aboutX << " " << aboutA << std::endl;
    }

    file.close();
    file2.close();
    return true;

}
*/

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

/*
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
*/