#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <iomanip>  // äëÿ std::setprecision

using namespace std;
std::string basePath = "C:/Users/dzhda/OneDrive/Ðàáî÷èé ñòîë/ïîëèòåõ/ìàòëàá/";
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
                    cerr << "Îøèáêà: ìàòðèöà íå ìîæåò áûòü ðàçëîæåíà íà LU (äåëåíèå íà íîëü)." << endl;
                    return false;
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
    return true;
}

// Ôóíêöèÿ ïðÿìîé ïîäñòàíîâêè äëÿ ðåøåíèÿ ñèñòåìû Ly = b
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

// Ôóíêöèÿ îáðàòíîé ïîäñòàíîâêè äëÿ ðåøåíèÿ ñèñòåìû Ux = y
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

// Îñíîâíàÿ ôóíêöèÿ äëÿ ðåøåíèÿ ÑËÀÓ Ax = b ñ èñïîëüçîâàíèåì LU-ðàçëîæåíèÿ
vector<long double> solveLU(const vector<vector<long double>>& A, const vector<long double>& b) {
    int n = A.size();
    vector<vector<long double>> L, U;
    if (!LU_Decomposition(A, L, U)) {
        cerr << "Îøèáêà ïðè âûïîëíåíèè LU-ðàçëîæåíèÿ." << endl;
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

    // Âû÷èñëÿåì Ax
    vector<long double> Ax = multiplyMatrixVector(A, x);

    //Âû÷èñëÿåì r = Ax - b
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
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Matrix A:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå n
    getline(file, line);
    // Èíèöèàëèçèðóåì ìàòðèöó A
    A.resize(n, std::vector<long double>(n));
    // Ñ÷èòûâàåì ìàòðèöó A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Vector b:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå ìàòðèöû
    getline(file, line);
    // Èíèöèàëèçèðóåì âåêòîð b
    b.resize(n);
    // Ñ÷èòûâàåì âåêòîð b
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Solution x:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå âåêòîðà b
    getline(file, line);
    // Èíèöèàëèçèðóåì âåêòîð x
    ans.resize(n);
    // Ñ÷èòûâàåì ðåøåíèå x
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
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Matrix A:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå n
    getline(file, line);
    // Èíèöèàëèçèðóåì ìàòðèöó A
    A.resize(n, std::vector<long double>(n));
    // Ñ÷èòûâàåì ìàòðèöó A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Vector b:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå ìàòðèöû
    getline(file, line);
    // Èíèöèàëèçèðóåì âåêòîð b
    b.resize(n);
    // Ñ÷èòûâàåì âåêòîð b
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }
    // Ïðîïóñêàåì ñòðîêó çàãîëîâêà "Solution x:"
    getline(file, line); // äëÿ çàâåðøåíèÿ ñòðîêè ïîñëå âåêòîðà b
    getline(file, line);
    // Èíèöèàëèçèðóåì âåêòîð x
    ans.resize(n);
    // Ñ÷èòûâàåì ðåøåíèå x
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
