#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip> // Для std::setprecision

using namespace std;
const std::string basePath = "C:/Users/dzhda/OneDrive/Рабочий стол/политех/Численные методы/3 сем/lab 4/Numercial Methods lab 6 matrix generation/Numercial Methods lab 6 matrix generation/matrices_output_modified/";
long double Eps = 0.0000000001;
int maxIterations = 100;

struct MatrixData {
    int matrixIndex;
    long double otdNumber;
    long double shift;
    long double k;
    std::vector<std::vector<long double>> matrix;
    std::vector<long double> eigenvalues;
};

MatrixData readMatrixFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    MatrixData data;
    std::string line;

    while (std::getline(file, line)) {
        if (line.find("Matrix Index:") != std::string::npos) {
            data.matrixIndex = std::stoi(line.substr(14));
        }
        else if (line.find("Otd_number:") != std::string::npos) {
            data.otdNumber = std::stold(line.substr(12));
        }
        else if (line.find("Shift:") != std::string::npos) {
            data.shift = std::stold(line.substr(7));
        }
        else if (line.find("k:") != std::string::npos) {
            data.k = std::stold(line.substr(3));
        }
        else if (line.find("Matrix:") != std::string::npos) {
            while (std::getline(file, line) && line.find("Eigenvalues:") == std::string::npos) {
                std::istringstream rowStream(line);
                std::vector<long double> row;
                std::string value;
                while (std::getline(rowStream, value, ',')) {
                    row.push_back(std::stold(value));
                }
                data.matrix.push_back(row);
            }
        }

        if (line.find("Eigenvalues:") != std::string::npos) {
            // Чтение собственных значений
            while (std::getline(file, line) && !line.empty()) {
                data.eigenvalues.push_back(std::stold(line));
            }
        }
    }

    return data;
}
void chekFile(const std::string& filename) {
    try {
        MatrixData data = readMatrixFile(filename);

        // Установка фиксированного формата вывода и точности
        std::cout << std::fixed << std::setprecision(6);

        // Вывод данных
        std::cout << "Matrix Index: " << data.matrixIndex << std::endl;
        std::cout << "Gap between 4-th and 5-th: " << data.otdNumber << std::endl;

        std::cout << "Matrix:" << std::endl;
        for (const auto& row : data.matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Eigenvalues:" << std::endl;
        for (const auto& eigenvalue : data.eigenvalues) {
            std::cout << eigenvalue << " ";
        }
        std::cout << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

std::vector<std::vector<long double>> subtractMatrices(const std::vector<std::vector<long double>>& matrixA, const std::vector<std::vector<long double>>& matrixB) {
    std::vector<std::vector<long double>> result(matrixA.size(), std::vector<long double>(matrixA[0].size(), 0.0));
    for (size_t i = 0; i < matrixA.size(); ++i) {
        for (size_t j = 0; j < matrixA[i].size(); ++j) {
            result[i][j] = matrixA[i][j] - matrixB[i][j];
        }
    }
    return result;
}
std::vector<std::vector<long double>> multiplyMatrixScalar(const std::vector<std::vector<long double>>& A, long double scalar) {
    int n = A.size();
    std::vector<std::vector<long double>> result(n, std::vector<long double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] * scalar;
        }
    }
    return result;
}
std::vector<std::vector<long double>> generateSingleMatrix(size_t N) {
    std::vector<std::vector<long double>> identityMatrix(N, std::vector<long double>(N, 0.0));
    for (size_t i = 0; i < N; ++i) {
        identityMatrix[i][i] = 1.0;
    }
    return identityMatrix;
}
long double findMu(const std::vector<long double>& currentX) {
    long double mu = fabs(currentX[0]);
    int sign;
    if (currentX[0] >= 0) {
        sign = 1;
    }
    else {
        sign = -1;
    }
    for (int i = 1; i < currentX.size(); i++) {
        if (fabs(currentX[i]) > mu) {
            mu = fabs(currentX[i]);
            if (currentX[i] >= 0) {
                sign = 1;
            }
            else {
                sign = -1;
            }
        }
    }
    return mu * sign;
}
std::vector<long double> normalization(const std::vector<long double>& currentX) {
    long double mu = findMu(currentX);
    std::vector<long double> newX(currentX.size());
    for (int i = 0; i < currentX.size(); i++) {
        newX[i] = currentX[i] / mu;
    }
    return newX;
}

std::vector<long double> multiplyVectorScalar(const std::vector<long double>& vec, long double scalar) {
    std::vector<long double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}
const std::vector<long double> subtractVectors(const std::vector<long double>& vec1, const std::vector<long double>& vec2) {
    std::vector<long double> delta;
    delta.resize(vec1.size());
    for (int i = 0; i != vec1.size(); ++i) {
        delta[i] = vec1[i] - vec2[i];
    }
    return delta;
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
vector<long double> solveLU(const vector<vector<long double>>& A, const vector<long double>& b, vector<vector<long double>> L, vector<vector<long double>> U) {
    int n = A.size();
    vector<long double> y = forwardSubstitution(L, b);
    vector<long double> x = backSubstitution(U, y);
    return x;
}


std::pair<long double, std::vector<long double>> reverseIterationsMethod(const vector<vector<long double>>& A, long double ApproxLambda){
    long double lambda_n;
    std::vector<long double> omega_n;
    vector<vector<long double>> L;
    vector<vector<long double>> U;
    vector<long double> currentY(A.size(), 1.0);
    vector<long double> newY;
    long double mu;

    vector<vector<long double>> A_new = subtractMatrices(A, multiplyMatrixScalar(generateSingleMatrix(A.size()), ApproxLambda));
    LU_Decomposition(A_new, L, U);
    currentY = normalization(currentY);
    newY = solveLU(A_new, currentY, L, U);
    int num = 1;
    mu = findMu(newY);
    lambda_n = ApproxLambda + 1.0 / mu;
    while ((norm(subtractVectors(multiplyMatrixVector(A, newY), multiplyVectorScalar(newY, lambda_n))) / norm(newY)) >= Eps && num < maxIterations) {
        num++;
        currentY = normalization(newY);
        newY = solveLU(A_new, currentY, L, U);
        mu = findMu(newY);
    }

    lambda_n = ApproxLambda +  1.0 / mu;
    omega_n = newY;

    return std::pair<long double, std::vector<long double>>(lambda_n, omega_n);
}
void errorFromIterNum(const std::string& input_filename, const std::string& output_filename, int finding_num) {
    MatrixData data = readMatrixFile(input_filename);
    ofstream output(output_filename);
    output << std::setprecision(std::numeric_limits<long double>::digits10);

    long double lambda_n;
    std::vector<long double> omega_n;
    vector<vector<long double>> L;
    vector<vector<long double>> U;
    vector<vector<long double>> A = data.matrix;
    vector<long double> currentY(A.size(), 1.0);
    vector<long double> newY;
    long double mu;
    long double ApproxLambda = data.eigenvalues[data.eigenvalues.size() - (finding_num)] * 1.1;

    vector<vector<long double>> A_new = subtractMatrices(A, multiplyMatrixScalar(generateSingleMatrix(A.size()), ApproxLambda));
    LU_Decomposition(A_new, L, U);
    currentY = normalization(currentY);
    newY = solveLU(A_new, currentY, L, U);
    int num = 1;
    mu = findMu(newY);
    lambda_n = ApproxLambda + 1.0 / mu;
    output << abs(data.eigenvalues[data.eigenvalues.size() - (finding_num)] - lambda_n) << std::endl;
    while (( norm(subtractVectors(multiplyMatrixVector(A, newY), multiplyVectorScalar(newY, lambda_n))) / norm(newY) ) >= Eps && num < maxIterations) {
        num++;
        currentY = normalization(newY);
        newY = solveLU(A_new, currentY, L, U);
        mu = findMu(newY);
        lambda_n = ApproxLambda + 1.0 / mu;
        output << abs(data.eigenvalues[data.eigenvalues.size() - (finding_num)] - lambda_n) << std::endl;
    }

    lambda_n = ApproxLambda + 1.0 / mu;
    output << abs(data.eigenvalues[data.eigenvalues.size() - (finding_num)] - lambda_n) << std::endl;

    output.close();
}
void errrorFromEps(const std::string& input_filename, const std::string& output_filename, int finding_num) {
    MatrixData data = readMatrixFile(input_filename);
    ofstream output(output_filename);
    output << std::setprecision(std::numeric_limits<long double>::digits10);

    long double lambda_n;
    std::vector<long double> omega_n;
    vector<vector<long double>> L;
    vector<vector<long double>> U;
    vector<vector<long double>> A = data.matrix;
    vector<long double> newY;
    long double mu;
    long double ApproxLambda = data.eigenvalues[data.eigenvalues.size() - (finding_num)] * 1.1;

    vector<vector<long double>> A_new = subtractMatrices(A, multiplyMatrixScalar(generateSingleMatrix(A.size()), ApproxLambda));
    LU_Decomposition(A_new, L, U);

    for (int i = 0; i != 10; i++) {
        long double CurrentEps = std::pow(0.1, i+1);

        vector<long double> currentY(A.size(), 1.0);
        currentY = normalization(currentY);
        newY = solveLU(A_new, currentY, L, U);
        int num = 1;
        mu = findMu(newY);
        lambda_n = ApproxLambda + 1.0 / mu;
        while ((norm(subtractVectors(multiplyMatrixVector(A, newY), multiplyVectorScalar(newY, lambda_n))) / norm(newY)) >= CurrentEps && num < maxIterations) {
            num++;
            currentY = normalization(newY);
            newY = solveLU(A_new, currentY, L, U);
            mu = findMu(newY);
            lambda_n = ApproxLambda + 1.0 / mu;
        }

        lambda_n = ApproxLambda + 1.0 / mu;
        output << abs(data.eigenvalues[data.eigenvalues.size() - (finding_num)] - lambda_n) << std::endl;
    }
    output.close();
}

void example() {
    const std::string path = basePath + "matrix_data_gap_1.txt";
    chekFile(path);
    MatrixData matrix1 = readMatrixFile(path);
    int finding_num = 7;
    std::pair<long double, std::vector<long double>> result = reverseIterationsMethod(matrix1.matrix, matrix1.eigenvalues[matrix1.eigenvalues.size() - (finding_num)] * 1.001);
    long double lambda_n = result.first;
    std::vector<long double> omega_n = result.second;
    std::cout << "***************************\n" << "Found eigenvalue: " << lambda_n << '\n';
    std::cout << "Found eigenvector: " << '\n';
    for (auto num : omega_n) {
        std::cout << num << ' ';
    }
    std::cout << '\n';
}

int main() {
    example();
    errorFromIterNum(basePath + "matrix_data_gap_1.txt", basePath + "ErrorFromIterNum_gap_1.txt", 7);
    errorFromIterNum(basePath + "matrix_data_gap_2.txt", basePath + "ErrorFromIterNum_gap_2.txt", 7);
    errorFromIterNum(basePath + "matrix_data_gap_3.txt", basePath + "ErrorFromIterNum_gap_3.txt", 7);
    errorFromIterNum(basePath + "matrix_data_gap_4.txt", basePath + "ErrorFromIterNum_gap_4.txt", 7);
    errorFromIterNum(basePath + "matrix_data_gap_5.txt", basePath + "ErrorFromIterNum_gap_5.txt", 7);
    errrorFromEps(basePath + "matrix_data_gap_5.txt", basePath + "ErrorFromEps_gap_5.txt", 7);
    return 0;
}
