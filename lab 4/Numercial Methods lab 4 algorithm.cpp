#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <utility>
#include <iomanip>  // äëÿ std::setprecision

using namespace std;
long double Eps = 0.00000000000001;
int steps = 40;
long double condA = 5.0;
std::string basePath = "C:/Users/dzhda/OneDrive/Ðàáî÷èé ñòîë/ïîëèòåõ/ìàòëàá/lab4matrix5.txt";
//std::string basePath = "C:/Users/dzhda/OneDrive/Ðàáî÷èé ñòîë/ïîëèòåõ/ìàòëàá/lab4matrix_test.txt";
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
