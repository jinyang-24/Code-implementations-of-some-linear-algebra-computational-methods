#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;
using namespace std::chrono;

// ==================== Matrix Class ====================
class Matrix {
private:
    int rows_;
    int cols_;
    vector<vector<double>> data_;
    
public:
    Matrix(int rows = 0, int cols = 0) : rows_(rows), cols_(cols) {
        data_.resize(rows, vector<double>(cols, 0.0));
    }
    
    Matrix(const vector<vector<double>>& realData) : rows_(realData.size()), cols_(realData[0].size()), data_(realData) {}
    
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    
    double& operator()(int i, int j) { return data_[i][j]; }
    const double& operator()(int i, int j) const { return data_[i][j]; }
    
    // Matrix multiplication
    Matrix multiply(const Matrix& other) const {
        if (cols_ != other.rows_) {
            cerr << "Matrix dimensions don't match for multiplication!" << endl;
            return Matrix(0, 0);
        }
        
        Matrix result(rows_, other.cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < other.cols_; j++) {
                double sum = 0.0;
                for (int k = 0; k < cols_; k++) {
                    sum += data_[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }
    
    // Matrix transpose
    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(j, i) = data_[i][j];
            }
        }
        return result;
    }
    
    // Matrix subtraction
    Matrix subtract(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            cerr << "Matrix dimensions don't match for subtraction!" << endl;
            return Matrix(0, 0);
        }
        
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j] - other(i, j);
            }
        }
        return result;
    }
    
    // Matrix addition
    Matrix add(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            cerr << "Matrix dimensions don't match for addition!" << endl;
            return Matrix(0, 0);
        }
        
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j] + other(i, j);
            }
        }
        return result;
    }
    
    // Scalar multiplication
    Matrix scalarMultiply(double scalar) const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j] * scalar;
            }
        }
        return result;
    }
    
    // Get column vector
    Matrix getColumn(int j) const {
        Matrix result(rows_, 1);
        for (int i = 0; i < rows_; i++) {
            result(i, 0) = data_[i][j];
        }
        return result;
    }
    
    // Set column vector
    void setColumn(int j, const Matrix& col) {
        for (int i = 0; i < rows_; i++) {
            data_[i][j] = col(i, 0);
        }
    }
    
    // Frobenius norm
    double frobeniusNorm() const {
        double sum = 0.0;
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                sum += data_[i][j] * data_[i][j];
            }
        }
        return sqrt(sum);
    }
    
    // Check if matrix is square
    bool isSquare() const {
        return rows_ == cols_;
    }
    
    // Generate identity matrix
    static Matrix identity(int n) {
        Matrix result(n, n);
        for (int i = 0; i < n; i++) {
            result(i, i) = 1.0;
        }
        return result;
    }
    
    // Copy matrix
    Matrix copy() const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j];
            }
        }
        return result;
    }
    
    // Print matrix
    void print() const {
        cout << fixed << setprecision(6);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                cout << setw(12) << data_[i][j] << " ";
            }
            cout << endl;
        }
    }
    
    // Get submatrix
    Matrix submatrix(int r1, int r2, int c1, int c2) const {
        Matrix result(r2 - r1 + 1, c2 - c1 + 1);
        for (int i = r1; i <= r2; i++) {
            for (int j = c1; j <= c2; j++) {
                result(i - r1, j - c1) = data_[i][j];
            }
        }
        return result;
    }
    
    // Set submatrix
    void setSubmatrix(int r, int c, const Matrix& sub) {
        for (int i = 0; i < sub.rows(); i++) {
            for (int j = 0; j < sub.cols(); j++) {
                data_[r + i][c + j] = sub(i, j);
            }
        }
    }
    
    // Get row vector
    Matrix getRow(int i) const {
        Matrix result(1, cols_);
        for (int j = 0; j < cols_; j++) {
            result(0, j) = data_[i][j];
        }
        return result;
    }
    
    // Set row vector
    void setRow(int i, const Matrix& row) {
        for (int j = 0; j < cols_; j++) {
            data_[i][j] = row(0, j);
        }
    }
};

// ==================== File I/O Functions ====================
Matrix readMatrixFromFile(const string& filename, int& method) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return Matrix(0, 0);
    }
    
    string line;
    
    // Read first line: Method
    getline(file, line);
    size_t eqPos = line.find('=');
    if (eqPos == string::npos) {
        cerr << "Invalid file format: first line should contain Method=" << endl;
        return Matrix(0, 0);
    }
    
    // Extract Method value
    size_t semicolonPos = line.find(';');
    if (semicolonPos == string::npos) {
        cerr << "Invalid file format: Method line should end with semicolon" << endl;
        return Matrix(0, 0);
    }
    
    string methodStr = line.substr(eqPos + 1, semicolonPos - eqPos - 1);
    method = stoi(methodStr);
    
    // Read second line: A = ...
    getline(file, line);
    
    // Read matrix data
    vector<vector<double>> matrixData;
    while (getline(file, line)) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Check for end of matrix
        if (line[0] == ']') break;
        
        // Remove leading '[' or space
        if (line[0] == '[') {
            line = line.substr(1);
        }
        
        // Remove trailing ';'
        if (line.back() == ';') {
            line.pop_back();
        }
        
        // Parse numbers in a row
        vector<double> row;
        stringstream ss(line);
        double val;
        while (ss >> val) {
            row.push_back(val);
        }
        
        if (!row.empty()) {
            matrixData.push_back(row);
        }
    }
    
    file.close();
    
    if (matrixData.empty()) {
        return Matrix(0, 0);
    }
    
    // Check if all rows have same number of columns
    int cols = matrixData[0].size();
    for (size_t i = 1; i < matrixData.size(); i++) {
        if (matrixData[i].size() != cols) {
            cerr << "Matrix rows have inconsistent lengths!" << endl;
            return Matrix(0, 0);
        }
    }
    
    return Matrix(matrixData);
}

void writeMatrixToFile(const string& filename, const string& matrixName, const Matrix& mat) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot create file: " << filename << endl;
        return;
    }
    
    file << fixed << setprecision(6);
    file << matrixName << " = ..." << endl;
    file << "[";
    
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            file << mat(i, j);
            if (j < mat.cols() - 1) {
                file << " ";
            }
        }
        if (i < mat.rows() - 1) {
            file << ";" << endl;
        } else {
            file << "];" << endl;
        }
    }
    
    file.close();
}

// ==================== Helper Functions ====================

// Vector dot product
double dotProduct(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != 1 || b.cols() != 1) {
        cerr << "Vector dimensions don't match!" << endl;
        return 0.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < a.rows(); i++) {
        sum += a(i, 0) * b(i, 0);
    }
    return sum;
}

// Vector norm
double vectorNorm(const Matrix& v) {
    return sqrt(dotProduct(v, v));
}

// ==================== QR Decomposition Algorithms ====================

// Method 1: QR via Cholesky decomposition (Lecture 9 page 2)
bool qrCholesky(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    
    // Check if matrix is non-singular (via Cholesky decomposition)
    Matrix AT = A.transpose();
    Matrix ATA = AT.multiply(A);
    
    // Perform Cholesky decomposition: ATA = L * L^T, where L is lower triangular
    Matrix L(n, n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L(i, k) * L(j, k);
            }
            
            if (i == j) {
                double diag = ATA(i, i) - sum;
                if (diag <= 0) {
                    // Matrix not positive definite, cannot use Cholesky
                    return false;
                }
                L(i, i) = sqrt(diag);
            } else {
                L(i, j) = (ATA(i, j) - sum) / L(j, j);
            }
        }
    }
    
    // Now ATA = L * L^T, let R = L^T (upper triangular)
    R = L.transpose();
    
    // Solve R^T * Q^T = A^T to get Q
    // Equivalent to solving L * Q^T = A^T
    Matrix QT = AT.copy();  // A^T
    
    // Forward substitution to solve L * X = A^T
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L(i, j) * QT(j, k);
            }
            QT(i, k) = (QT(i, k) - sum) / L(i, i);
        }
    }
    
    // Q = (Q^T)^T
    Q = QT.transpose();
    
    return true;
}

// Method 2: Householder reflections (Lecture 9 pages 3-4)
void qrHouseholder(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    R = A.copy();
    Q = Matrix::identity(n);
    
    for (int k = 0; k < n - 1; k++) {
        // Compute x = R[k:n, k]
        int m = n - k;
        Matrix x(m, 1);
        for (int i = 0; i < m; i++) {
            x(i, 0) = R(k + i, k);
        }
        
        // Compute vector norm
        double norm_x = vectorNorm(x);
        
        if (norm_x < 1e-15) {
            continue; // Column already zero, skip
        }
        
        // Compute reflection vector u
        Matrix u = x.copy();
        
        // Choose sign to avoid numerical instability
        double sign = (x(0, 0) >= 0) ? 1.0 : -1.0;
        u(0, 0) += sign * norm_x;
        
        // Normalize u
        double norm_u = vectorNorm(u);
        
        if (norm_u < 1e-15) {
            continue;
        }
        
        for (int i = 0; i < m; i++) {
            u(i, 0) /= norm_u;
        }
        
        // Apply Householder transformation to the lower right part of R
        for (int j = k; j < n; j++) {
            // Compute u^T * R[k:n, j]
            double u_dot_r = 0.0;
            for (int i = 0; i < m; i++) {
                u_dot_r += u(i, 0) * R(k + i, j);
            }
            
            // Update R[k:n, j]
            for (int i = 0; i < m; i++) {
                R(k + i, j) -= 2.0 * u(i, 0) * u_dot_r;
            }
        }
        
        // Accumulate Q: Q = Q * H
        for (int i = 0; i < n; i++) {
            // Compute Q[i, k:n] * u
            double q_dot_u = 0.0;
            for (int j = 0; j < m; j++) {
                q_dot_u += Q(i, k + j) * u(j, 0);
            }
            
            // Update Q[i, k:n]
            for (int j = 0; j < m; j++) {
                Q(i, k + j) -= 2.0 * q_dot_u * u(j, 0);
            }
        }
    }
}

// Method 3: Givens rotations (Lecture 9 page 5) - fixed version
void qrGivens(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    R = A.copy();
    Q = Matrix::identity(n);
    
    // According to Lecture 9 page 5, Givens rotations should eliminate elements in specific order
    // Typically start from bottom-left corner, eliminate column by column
    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double a = R(j, j);
            double b = R(i, j);
            
            // If already zero, skip
            if (fabs(b) < 1e-15) {
                continue;
            }
            
            // Compute rotation parameters
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = -b / r;  // Note the sign, lecture uses negative
            
            // Apply rotation to rows j and i of R
            for (int k = j; k < n; k++) {
                double r_jk = R(j, k);
                double r_ik = R(i, k);
                R(j, k) = c * r_jk - s * r_ik;  // Note minus sign here
                R(i, k) = s * r_jk + c * r_ik;
            }
            
            // Set R(i, j) to zero (eliminate rounding errors)
            R(i, j) = 0.0;
            
            // Apply same rotation to Q
            // Q = Q * G^T, where G is rotation matrix
            for (int k = 0; k < n; k++) {
                double q_kj = Q(k, j);
                double q_ki = Q(k, i);
                Q(k, j) = c * q_kj - s * q_ki;
                Q(k, i) = s * q_kj + c * q_ki;
            }
        }
    }
    
    // Ensure R diagonal elements are non-negative (optional, but helpful for stability)
    for (int i = 0; i < n; i++) {
        if (R(i, i) < 0) {
            for (int j = i; j < n; j++) {
                R(i, j) = -R(i, j);
            }
            for (int j = 0; j < n; j++) {
                Q(j, i) = -Q(j, i);
            }
        }
    }
}

// Method 4: Gram-Schmidt orthogonalization (Lecture 9 pages 7-8)
void qrGramSchmidt(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    Q = Matrix(n, n);
    R = Matrix(n, n);
    
    for (int j = 0; j < n; j++) {
        // Get j-th column of A
        Matrix v = A.getColumn(j);
        
        // Orthogonalize
        for (int i = 0; i < j; i++) {
            R(i, j) = dotProduct(Q.getColumn(i), A.getColumn(j));
            
            // v = v - R(i,j) * Q[:,i]
            for (int k = 0; k < n; k++) {
                v(k, 0) -= R(i, j) * Q(k, i);
            }
        }
        
        // Compute R(j,j)
        R(j, j) = vectorNorm(v);
        
        if (R(j, j) < 1e-15) {
            // Linearly dependent, set q_j to zero vector
            for (int k = 0; k < n; k++) {
                Q(k, j) = 0.0;
            }
        } else {
            // Normalize to get q_j
            for (int k = 0; k < n; k++) {
                Q(k, j) = v(k, 0) / R(j, j);
            }
        }
    }
}

// ==================== QR Decomposition Verification Function ====================
void verifyQR(const Matrix& A, const Matrix& Q, const Matrix& R, const string& methodName) {
    cout << "\nVerifying " << methodName << " method:" << endl;
    
    // 1. Verify A = QR
    Matrix QR = Q.multiply(R);
    Matrix diff_A_QR = A.subtract(QR);
    double error_A_QR = diff_A_QR.frobeniusNorm();
    cout << "  ||A - QR||_E = " << scientific << setprecision(6) << error_A_QR;
    cout << fixed << " (relative error: " << error_A_QR / A.frobeniusNorm() << ")" << endl;
    
    // 2. Verify orthogonality of Q: Q^T * Q = I
    Matrix QT = Q.transpose();
    Matrix QTQ = QT.multiply(Q);
    Matrix I = Matrix::identity(A.rows());
    Matrix diff_QTQ_I = QTQ.subtract(I);
    double error_orth = diff_QTQ_I.frobeniusNorm();
    cout << "  ||Q^T * Q - I||_E = " << scientific << error_orth << endl;
    
    // 3. Verify R is upper triangular
    bool isUpperTriangular = true;
    for (int i = 1; i < R.rows(); i++) {
        for (int j = 0; j < i; j++) {
            if (fabs(R(i, j)) > 1e-10) {
                isUpperTriangular = false;
                cout << "  Found non-zero lower triangular element: R(" << i << "," << j << ") = " << R(i, j) << endl;
                break;
            }
        }
        if (!isUpperTriangular) break;
    }
    cout << "  R is upper triangular: " << (isUpperTriangular ? "Yes" : "No") << endl;
    
    // 4. Show diagonal elements of R
    cout << "  R diagonal elements: ";
    for (int i = 0; i < min(5, R.rows()); i++) {
        cout << R(i, i) << " ";
    }
    if (R.rows() > 5) cout << "...";
    cout << endl;
}

// ==================== Main Function ====================
int main() {
    // Get input filename
    cout << "Enter input filename (e.g., Amat8.num): ";
    string inputFilename;
    cin >> inputFilename;
    
    // Extract number part from input filename
    size_t start = inputFilename.find("Amat");
    size_t end = inputFilename.find(".num");
    if (start == string::npos || end == string::npos) {
        cerr << "Invalid filename format!" << endl;
        return 1;
    }
    string num = inputFilename.substr(start + 4, end - start - 4);
    
    // Read matrix and method
    int method;
    Matrix A = readMatrixFromFile(inputFilename, method);
    
    if (A.rows() == 0) {
        cerr << "Failed to read matrix!" << endl;
        return 1;
    }
    
    // Check if matrix is square
    if (!A.isSquare()) {
        cerr << "Error: Matrix is not square!" << endl;
        return 1;
    }
    
    // Check if method value is valid
    if (method < 0 || method > 4) {
        cerr << "Error: Invalid method value! Should be 0-4" << endl;
        return 1;
    }
    
    // Create result file
    ofstream resFile("Res" + num + ".num");
    if (!resFile.is_open()) {
        cerr << "Cannot create result file!" << endl;
        return 1;
    }
    
    // Execute QR decomposition according to method
    if (method == 0) {
        // Use all four methods
        resFile << fixed << setprecision(15);
        
        for (int m = 1; m <= 4; m++) {
            Matrix Q, R;
            bool success = true;
            
            auto startTime = high_resolution_clock::now();
            
            string methodName;
            switch (m) {
                case 1:
                    methodName = "Cholesky";
                    success = qrCholesky(A, Q, R);
                    break;
                case 2:
                    methodName = "Householder";
                    qrHouseholder(A, Q, R);
                    break;
                case 3:
                    methodName = "Givens";
                    qrGivens(A, Q, R);
                    break;
                case 4:
                    methodName = "Gram-Schmidt";
                    qrGramSchmidt(A, Q, R);
                    break;
            }
            
            auto endTime = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(endTime - startTime);
            double timeSeconds = duration.count() / 1000000.0;
            
            if (!success) {
                resFile << "Method = " << m << "." << endl;
                resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
                resFile << "Погрешность полученного QR-разложения: Cholesky decomposition failed (matrix not positive definite)." << endl;
                resFile << "Погрешность ортогональности матрицы Q: N/A." << endl;
                
                cout << "\nMethod " << m << " (Cholesky) failed: Matrix not positive definite" << endl;
                continue;
            }
            
            // Compute error: ||A - QR||_E
            Matrix QR = Q.multiply(R);
            Matrix diff = A.subtract(QR);
            double qrError = diff.frobeniusNorm();
            
            // Compute orthogonality error of Q: ||Q^T * Q - I||_E
            Matrix QT = Q.transpose();
            Matrix QTQ = QT.multiply(Q);
            Matrix I = Matrix::identity(A.rows());
            Matrix orthDiff = QTQ.subtract(I);
            double orthError = orthDiff.frobeniusNorm();
            
            // Write results
            resFile << "Method = " << m << "." << endl;
            resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
            resFile << "Погрешность полученного QR-разложения: " << qrError << "." << endl;
            resFile << "Погрешность ортогональности матрицы Q: " << orthError << "." << endl;
            
            // Write Q and R matrices to files
            writeMatrixToFile("Qmat" + num + "_" + to_string(m) + ".num", "Q", Q);
            writeMatrixToFile("Rmat" + num + "_" + to_string(m) + ".num", "R", R);
            
            // Verify QR decomposition
            verifyQR(A, Q, R, methodName);
        }
    } else {
        // Use specified method
        Matrix Q, R;
        bool success = true;
        string methodName;
        
        auto startTime = high_resolution_clock::now();
        
        switch (method) {
            case 1:
                methodName = "Cholesky";
                success = qrCholesky(A, Q, R);
                break;
            case 2:
                methodName = "Householder";
                qrHouseholder(A, Q, R);
                break;
            case 3:
                methodName = "Givens";
                qrGivens(A, Q, R);
                break;
            case 4:
                methodName = "Gram-Schmidt";
                qrGramSchmidt(A, Q, R);
                break;
        }
        
        auto endTime = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(endTime - startTime);
        double timeSeconds = duration.count() / 1000000.0;
        
        if (!success) {
            resFile << "Method = " << method << "." << endl;
            resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
            resFile << "Погрешность полученного QR-разложения: Cholesky decomposition failed (matrix not positive definite)." << endl;
            resFile << "Погрешность ортогональности матрицы Q: N/A." << endl;
            
            cout << "\nMethod " << method << " (Cholesky) failed: Matrix not positive definite" << endl;
        } else {
            // Compute errors
            Matrix QR = Q.multiply(R);
            Matrix diff = A.subtract(QR);
            double qrError = diff.frobeniusNorm();
            
            // Compute orthogonality error of Q
            Matrix QT = Q.transpose();
            Matrix QTQ = QT.multiply(Q);
            Matrix I = Matrix::identity(A.rows());
            Matrix orthDiff = QTQ.subtract(I);
            double orthError = orthDiff.frobeniusNorm();
            
            // Write results
            resFile << fixed << setprecision(15);
            resFile << "Method = " << method << "." << endl;
            resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
            resFile << "Погрешность полученного QR-разложения: " << qrError << "." << endl;
            resFile << "Погрешность ортогональности матрицы Q: " << orthError << "." << endl;
            
            // Write Q and R matrices to files
            writeMatrixToFile("Qmat" + num + ".num", "Q", Q);
            writeMatrixToFile("Rmat" + num + ".num", "R", R);
            
            // Verify QR decomposition
            verifyQR(A, Q, R, methodName);
        }
    }
    
    resFile.close();
    
    cout << "\nQR decomposition completed! Results saved to files." << endl;
    cout << "Note: Error information written to Res" << num << ".num file." << endl;
    
    return 0;
}