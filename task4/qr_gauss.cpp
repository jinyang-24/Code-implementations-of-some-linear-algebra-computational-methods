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
    
    Matrix(const vector<double>& vectorData) : rows_(vectorData.size()), cols_(1) {
        data_.resize(rows_, vector<double>(1));
        for (int i = 0; i < rows_; i++) {
            data_[i][0] = vectorData[i];
        }
    }
    
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
    
    // Vector 2-norm (for column vectors)
    double vectorNorm2() const {
        if (cols_ != 1) {
            cerr << "Not a column vector!" << endl;
            return 0.0;
        }
        
        double sum = 0.0;
        for (int i = 0; i < rows_; i++) {
            sum += data_[i][0] * data_[i][0];
        }
        return sqrt(sum);
    }
    
    // Check if matrix is square
    bool isSquare() const {
        return rows_ == cols_;
    }
    
    // Check if matrix is real (for this task, we only handle real matrices)
    bool isReal() const {
        return true;
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
        cerr << "Invalid file format: First line should contain Method=" << endl;
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
    
    // Read second line: A = ... or similar
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
        
        // Remove leading '[' or spaces
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
    
    // Check if all rows have the same number of columns
    int cols = matrixData[0].size();
    for (size_t i = 1; i < matrixData.size(); i++) {
        if (matrixData[i].size() != cols) {
            cerr << "Matrix rows have inconsistent lengths!" << endl;
            return Matrix(0, 0);
        }
    }
    
    // Create matrix
    Matrix result(matrixData.size(), matrixData[0].size());
    for (size_t i = 0; i < matrixData.size(); i++) {
        for (size_t j = 0; j < matrixData[i].size(); j++) {
            result(i, j) = matrixData[i][j];
        }
    }
    
    return result;
}

Matrix readMatrixFromFileSimple(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return Matrix(0, 0);
    }
    
    string line;
    
    // Read first line: MatrixName = ...
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
        
        // Remove leading '[' or spaces
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
    
    // Check if all rows have the same number of columns
    int cols = matrixData[0].size();
    for (size_t i = 1; i < matrixData.size(); i++) {
        if (matrixData[i].size() != cols) {
            cerr << "Matrix rows have inconsistent lengths!" << endl;
            return Matrix(0, 0);
        }
    }
    
    // Create matrix
    Matrix result(matrixData.size(), matrixData[0].size());
    for (size_t i = 0; i < matrixData.size(); i++) {
        for (size_t j = 0; j < matrixData[i].size(); j++) {
            result(i, j) = matrixData[i][j];
        }
    }
    
    return result;
}

// Simplified vector reading function - only handles single-line format: b = [num1;num2;num3;]
Matrix readVectorFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return Matrix(0, 0);
    }
    
    string line;
    
    // Read entire line: b = [2.00000000;-3.00000000;4.00000000];
    getline(file, line);
    
    file.close();
    
    // Check if format is correct
    if (line.empty()) {
        cerr << "Vector file format error: Empty line" << endl;
        return Matrix(0, 0);
    }
    
    // Find equals sign position
    size_t eqPos = line.find('=');
    if (eqPos == string::npos) {
        cerr << "Vector file format error: No equals sign found" << endl;
        return Matrix(0, 0);
    }
    
    // Find left bracket position
    size_t leftBracketPos = line.find('[', eqPos);
    if (leftBracketPos == string::npos) {
        cerr << "Vector file format error: No left bracket found" << endl;
        return Matrix(0, 0);
    }
    
    // Find right bracket position
    size_t rightBracketPos = line.find(']', leftBracketPos);
    if (rightBracketPos == string::npos) {
        cerr << "Vector file format error: No right bracket found" << endl;
        return Matrix(0, 0);
    }
    
    // Extract content inside brackets (including semicolons)
    string vectorContent = line.substr(leftBracketPos + 1, rightBracketPos - leftBracketPos - 1);
    
    // Split string by semicolons
    vector<double> vectorData;
    size_t start = 0;
    size_t end = 0;
    
    while ((end = vectorContent.find(';', start)) != string::npos) {
        string token = vectorContent.substr(start, end - start);
        
        // Remove possible whitespace
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        
        if (!token.empty()) {
            try {
                double val = stod(token);
                vectorData.push_back(val);
            } catch (const std::exception& e) {
                cerr << "Vector file format error: Cannot parse number '" << token << "'" << endl;
                return Matrix(0, 0);
            }
        }
        
        start = end + 1;
    }
    
    // Check content after last semicolon (if any)
    if (start < vectorContent.length()) {
        string token = vectorContent.substr(start);
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        
        if (!token.empty()) {
            try {
                double val = stod(token);
                vectorData.push_back(val);
            } catch (const std::exception& e) {
                cerr << "Vector file format error: Cannot parse number '" << token << "'" << endl;
                return Matrix(0, 0);
            }
        }
    }
    
    if (vectorData.empty()) {
        cerr << "Vector file format error: No valid data found" << endl;
        return Matrix(0, 0);
    }
    
    // Verify number of elements read
    cout << "Read " << vectorData.size() << " vector elements" << endl;
    
    // Print first few elements for verification
    int printCount = min(3, (int)vectorData.size());
    for (int i = 0; i < printCount; i++) {
        cout << "Element " << i << ": " << vectorData[i] << endl;
    }
    
    return Matrix(vectorData);
}

void writeVectorToFile(const string& filename, const string& vectorName, const Matrix& vec) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot create file: " << filename << endl;
        return;
    }
    
    file << fixed << setprecision(6);
    file << vectorName << " = ";
    file << "[";
    
    for (int i = 0; i < vec.rows(); i++) {
        file << vec(i, 0);
        if (i < vec.rows() - 1) {
            file << ";";
        }
    }
    
    file << "];" << endl;
    
    file.close();
}

// ==================== QR Solver for Linear Systems ====================

// Solve linear system Ax = b using QR decomposition
// Input: Q (orthogonal matrix), R (upper triangular matrix), b (vector)
// Output: Solution vector x
Matrix solveQR(const Matrix& Q, const Matrix& R, const Matrix& b) {
    int n = Q.rows();
    
    // Step 1: Compute Q^T * b
    Matrix QT = Q.transpose();  // For real matrices, transpose is conjugate transpose
    Matrix QTb = QT.multiply(b);
    
    // Step 2: Back substitution to solve Rx = Q^T * b
    Matrix x(n, 1);
    
    // Start from last row for back substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        
        // Compute contribution from known parts
        for (int j = i + 1; j < n; j++) {
            sum += R(i, j) * x(j, 0);
        }
        
        // Compute x_i
        x(i, 0) = (QTb(i, 0) - sum) / R(i, i);
    }
    
    return x;
}

// Check if file exists
bool fileExists(const string& filename) {
    ifstream file(filename);
    return file.good();
}

// ==================== Main Function ====================
int main() {
    // Get input file names
    cout << "Enter matrix A filename (e.g., Amat10.num): ";
    string matrixFilename;
    cin >> matrixFilename;
    
    cout << "Enter vector b filename (e.g., bvec10.num): ";
    string vectorFilename;
    cin >> vectorFilename;
    
    // Extract number part from matrix filename
    size_t start = matrixFilename.find("Amat");
    size_t end = matrixFilename.find(".num");
    if (start == string::npos || end == string::npos) {
        cerr << "Invalid filename format!" << endl;
        return 1;
    }
    string num = matrixFilename.substr(start + 4, end - start - 4);
    
    // Read matrix A and method
    int method;
    Matrix A = readMatrixFromFile(matrixFilename, method);
    
    if (A.rows() == 0) {
        cerr << "Failed to read matrix A!" << endl;
        return 1;
    }
    
    // Read vector b
    Matrix b = readVectorFromFile(vectorFilename);
    
    if (b.rows() == 0) {
        cerr << "Failed to read vector b!" << endl;
        return 1;
    }
    
    cout << "Successfully read vector b, dimensions: " << b.rows() << "×" << b.cols() << endl;
    
    // Check conditions
    if (!A.isSquare()) {
        cerr << "Error: Matrix A is not square!" << endl;
        return 1;
    }
    
    if (A.rows() != b.rows()) {
        cerr << "Error: Number of rows in A (" << A.rows() 
             << ") does not match length of b (" << b.rows() << ")!" << endl;
        return 1;
    }
    
    if (b.cols() != 1) {
        cerr << "Error: b is not a column vector!" << endl;
        return 1;
    }
    
    // Check if method value is valid
    if (method < 1 || method > 4) {
        cerr << "Error: Invalid method value! Should be 1-4" << endl;
        return 1;
    }
    
    // Check if QR decomposition files exist
    string QFilename = "Qmat" + num + ".num";
    string RFilename = "Rmat" + num + ".num";
    
    Matrix Q, R;
    
    // If QR files don't exist, need to call calc_qr (simplified handling, assuming files exist)
    if (!fileExists(QFilename) || !fileExists(RFilename)) {
        cerr << "Error: QR decomposition files not found!" << endl;
        cerr << "Please run calc_qr program first to generate QR decomposition files." << endl;
        return 1;
    }
    
    // Read Q and R matrices
    Q = readMatrixFromFileSimple(QFilename);
    R = readMatrixFromFileSimple(RFilename);
    
    if (Q.rows() == 0 || R.rows() == 0) {
        cerr << "Failed to read QR decomposition files!" << endl;
        return 1;
    }
    
    // Check dimension compatibility
    if (Q.rows() != A.rows() || R.rows() != A.rows()) {
        cerr << "Error: QR decomposition matrix dimensions do not match A!" << endl;
        return 1;
    }
    
    // Start timing (only SLAU solving time)
    auto startTime = high_resolution_clock::now();
    
    // Solve linear system using QR decomposition
    Matrix x = solveQR(Q, R, b);
    
    auto endTime = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(endTime - startTime);
    double timeSeconds = duration.count() / 1000000.0;
    
    // Compute error: ||Ax - b||_2
    Matrix Ax = A.multiply(x);
    Matrix residual = Ax.subtract(b);
    double residualNorm = residual.vectorNorm2();
    
    // Create result file
    ofstream resFile("Res" + num + ".num");
    if (!resFile.is_open()) {
        cerr << "Cannot create result file!" << endl;
        return 1;
    }
    
    // Write results (Russian content remains unchanged)
    resFile << fixed << setprecision(15);
    resFile << "Время решения СЛАУ: " << timeSeconds << " секунд." << endl;
    resFile << "Погрешность решения СЛАУ: " << residualNorm << "." << endl;
    
    resFile.close();
    
    // Write solution vector
    writeVectorToFile("xvec" + num + ".num", "x", x);
    
    cout << "\nLinear system solving completed!" << endl;
    cout << "Results saved to:" << endl;
    cout << "  Solution vector: xvec" << num << ".num" << endl;
    cout << "  Results: Res" << num << ".num" << endl;
    cout << "\nSolution time: " << timeSeconds << " seconds" << endl;
    cout << "Error ||Ax - b||_2 = " << residualNorm << endl;
    
    // Verify results
    cout << "\nVerification:" << endl;
    cout << "Solution vector dimensions: " << x.rows() << "×" << x.cols() << endl;
    
    // Display first few elements of solution vector
    int showElements = min(5, x.rows());
    cout << "First " << showElements << " elements of solution vector:" << endl;
    for (int i = 0; i < showElements; i++) {
        cout << "  x[" << i << "] = " << x(i, 0) << endl;
    }
    
    return 0;
}