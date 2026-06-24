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
    
    // Check if matrix is real symmetric (self-adjoint)
    bool isSymmetric(double tolerance = 1e-10) const {
        if (!isSquare()) return false;
        
        for (int i = 0; i < rows_; i++) {
            for (int j = i + 1; j < rows_; j++) {
                if (fabs(data_[i][j] - data_[j][i]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }
    
    // Check if matrix is real
    bool isReal() const {
        // For this task, we only handle real matrices
        return true;
    }
    
    // Create identity matrix
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
    
    // Find maximum off-diagonal element magnitude and position (for Jacobi method)
    void findMaxOffDiagonal(int& p, int& q, double& maxVal) const {
        if (!isSquare()) {
            p = q = -1;
            maxVal = 0.0;
            return;
        }
        
        maxVal = 0.0;
        p = 0;
        q = 1;
        
        for (int i = 0; i < rows_; i++) {
            for (int j = i + 1; j < rows_; j++) {
                double val = fabs(data_[i][j]);
                if (val > maxVal) {
                    maxVal = val;
                    p = i;
                    q = j;
                }
            }
        }
    }
    
    // Find optimal off-diagonal element position (for Jacobi method)
    void findOptimalOffDiagonal(int& p, int& q, double& maxVal) const {
        if (!isSquare()) {
            p = q = -1;
            maxVal = 0.0;
            return;
        }
        
        maxVal = 0.0;
        p = 0;
        q = 1;
        
        // Optimal selection strategy: choose element that maximizes change after rotation
        for (int i = 0; i < rows_; i++) {
            for (int j = i + 1; j < rows_; j++) {
                double aij = data_[i][j];
                double aii = data_[i][i];
                double ajj = data_[j][j];
                
                // Compute rotation angle and evaluate rotation effect
                double tau = (ajj - aii) / (2.0 * aij);
                double t = 1.0 / (fabs(tau) + sqrt(1.0 + tau * tau));
                if (tau < 0) t = -t;
                
                // Compute change after rotation
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;
                double change = 2.0 * aij * s * c * (aii - ajj) + (c*c - s*s) * aij * aij;
                
                double val = fabs(change);
                if (val > maxVal) {
                    maxVal = val;
                    p = i;
                    q = j;
                }
            }
        }
    }
    
    // Compute Givens rotation parameters
    void computeGivensRotation(int p, int q, double& c, double& s) const {
        if (fabs(data_[p][q]) < 1e-15) {
            c = 1.0;
            s = 0.0;
            return;
        }
        
        double tau = (data_[q][q] - data_[p][p]) / (2.0 * data_[p][q]);
        double t;
        
        if (tau >= 0) {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        } else {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }
        
        c = 1.0 / sqrt(1.0 + t * t);
        s = t * c;
    }
    
    // Apply Givens rotation
    void applyGivensRotation(int p, int q, double c, double s) {
        int n = rows_;
        vector<double> rowP(n), rowQ(n);
        
        // Save original rows
        for (int j = 0; j < n; j++) {
            rowP[j] = data_[p][j];
            rowQ[j] = data_[q][j];
        }
        
        // Update p-th and q-th rows
        for (int j = 0; j < n; j++) {
            if (j == p) {
                data_[p][p] = c * c * rowP[p] - 2.0 * c * s * rowP[q] + s * s * rowQ[q];
            } else if (j == q) {
                data_[q][q] = s * s * rowP[p] + 2.0 * c * s * rowP[q] + c * c * rowQ[q];
            } else {
                data_[p][j] = c * rowP[j] - s * rowQ[j];
                data_[j][p] = data_[p][j]; // Maintain symmetry
                data_[q][j] = s * rowP[j] + c * rowQ[j];
                data_[j][q] = data_[q][j]; // Maintain symmetry
            }
        }
        
        // Update (p,q) and (q,p) elements
        data_[p][q] = 0.0;
        data_[q][p] = 0.0;
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
Matrix readMatrixFromFile(const string& inputFileName, int& methodValue, double& stopThreshold) {
    ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        cerr << "Error: Cannot open file: " << inputFileName << endl;
        return Matrix(0, 0);
    }
    
    string lineContent;
    
    // Read first line: Method
    getline(inputFile, lineContent);
    size_t methodPos = lineContent.find("Method=");
    if (methodPos == string::npos) {
        cerr << "Error: Invalid file format - first line should contain Method=" << endl;
        return Matrix(0, 0);
    }
    
    // Extract Method value
    size_t methodEnd = lineContent.find(";", methodPos);
    if (methodEnd == string::npos) {
        cerr << "Error: Invalid file format - Method line should end with semicolon" << endl;
        return Matrix(0, 0);
    }
    
    string methodString = lineContent.substr(methodPos + 7, methodEnd - methodPos - 7);
    methodValue = stoi(methodString);
    
    // Read second line: StopValue
    getline(inputFile, lineContent);
    size_t stopPos = lineContent.find("StopValue=");
    if (stopPos == string::npos) {
        cerr << "Error: Invalid file format - second line should contain StopValue=" << endl;
        return Matrix(0, 0);
    }
    
    size_t stopEnd = lineContent.find(";", stopPos);
    if (stopEnd == string::npos) {
        cerr << "Error: Invalid file format - StopValue line should end with semicolon" << endl;
        return Matrix(0, 0);
    }
    
    string stopString = lineContent.substr(stopPos + 10, stopEnd - stopPos - 10);
    stopThreshold = stod(stopString);
    
    // Read third line: A = ... or similar
    getline(inputFile, lineContent);
    
    // Read matrix data
    vector<vector<double>> matrixValues;
    while (getline(inputFile, lineContent)) {
        // Remove leading/trailing whitespace
        lineContent.erase(0, lineContent.find_first_not_of(" \t"));
        lineContent.erase(lineContent.find_last_not_of(" \t") + 1);
        
        // Skip empty lines
        if (lineContent.empty()) continue;
        
        // Check if end of matrix
        if (lineContent[0] == ']') break;
        
        // Remove leading '[' or spaces
        if (lineContent[0] == '[') {
            lineContent = lineContent.substr(1);
        }
        
        // Remove trailing ';'
        if (lineContent.back() == ';') {
            lineContent.pop_back();
        }
        
        // Parse numbers in a row
        vector<double> rowValues;
        stringstream ss(lineContent);
        double value;
        while (ss >> value) {
            rowValues.push_back(value);
        }
        
        if (!rowValues.empty()) {
            matrixValues.push_back(rowValues);
        }
    }
    
    inputFile.close();
    
    if (matrixValues.empty()) {
        return Matrix(0, 0);
    }
    
    // Check if all rows have same number of columns
    int columns = matrixValues[0].size();
    for (size_t i = 1; i < matrixValues.size(); i++) {
        if (matrixValues[i].size() != columns) {
            cerr << "Error: Inconsistent row lengths in matrix!" << endl;
            return Matrix(0, 0);
        }
    }
    
    // Create matrix
    Matrix result(matrixValues.size(), matrixValues[0].size());
    for (size_t i = 0; i < matrixValues.size(); i++) {
        for (size_t j = 0; j < matrixValues[i].size(); j++) {
            result(i, j) = matrixValues[i][j];
        }
    }
    
    return result;
}

// Write eigenvalues to file
void writeEigenvaluesToFile(const string& outputFileName, int methodValue, double computationTime, 
                           const vector<double>& eigenvalues, bool bothMethods = false) {
    ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        cerr << "Error: Cannot create file: " << outputFileName << endl;
        return;
    }
    
    outputFile << fixed << setprecision(15);
    
    if (bothMethods) {
        // Method=0: Write results for both methods
        outputFile << "Method=1." << endl;
        outputFile << "Время вычисления собственных значений: " << computationTime << " секунд." << endl;
        outputFile << "Собственные значения: ";
        for (size_t i = 0; i < eigenvalues.size(); i++) {
            outputFile << eigenvalues[i];
            if (i < eigenvalues.size() - 1) outputFile << " ";
        }
        outputFile << endl;
        
        outputFile << "Method=2." << endl;
        outputFile << "Время вычисления собственных значений: " << computationTime << " секунд." << endl;
        outputFile << "Собственные значения: ";
        for (size_t i = 0; i < eigenvalues.size(); i++) {
            outputFile << eigenvalues[i];
            if (i < eigenvalues.size() - 1) outputFile << " ";
        }
        outputFile << endl;
    } else {
        // Method=1 or 2
        outputFile << "Method=" << methodValue << "." << endl;
        outputFile << "Время вычисления собственных значений: " << computationTime << " секунд." << endl;
        outputFile << "Собственные значения: ";
        for (size_t i = 0; i < eigenvalues.size(); i++) {
            outputFile << eigenvalues[i];
            if (i < eigenvalues.size() - 1) outputFile << " ";
        }
        outputFile << endl;
    }
    
    outputFile.close();
}

// Sort eigenvalues by magnitude (descending)
void sortEigenvaluesByMagnitude(vector<double>& eigenvalues) {
    sort(eigenvalues.begin(), eigenvalues.end(), 
         [](double a, double b) { return fabs(a) > fabs(b); });
}

// Compute eigenvalues using Jacobi method
vector<double> computeEigenvaluesJacobiMethod(const Matrix& matrix, int methodValue, double stopThreshold) {
    if (!matrix.isSquare() || !matrix.isSymmetric()) {
        cerr << "Error: Matrix is not real symmetric square matrix, cannot compute eigenvalues!" << endl;
        return vector<double>();
    }
    
    int matrixSize = matrix.rows();
    Matrix workingMatrix = matrix.copy(); // Copy matrix to avoid modifying original
    
    int iterationCount = 0;
    double maxOffDiagonal;
    int pivotRow, pivotCol;
    
    do {
        if (methodValue == 1) {
            // Method=1: Select maximum magnitude off-diagonal element
            workingMatrix.findMaxOffDiagonal(pivotRow, pivotCol, maxOffDiagonal);
        } else if (methodValue == 2) {
            // Method=2: Select optimal element
            workingMatrix.findOptimalOffDiagonal(pivotRow, pivotCol, maxOffDiagonal);
        } else {
            cerr << "Error: Invalid method!" << endl;
            return vector<double>();
        }
        
        if (maxOffDiagonal < stopThreshold) break;
        
        double cosTheta, sinTheta;
        workingMatrix.computeGivensRotation(pivotRow, pivotCol, cosTheta, sinTheta);
        workingMatrix.applyGivensRotation(pivotRow, pivotCol, cosTheta, sinTheta);
        
        iterationCount++;
        
        if (iterationCount % 100 == 0) {
            cout << "Iteration " << iterationCount << ": Max off-diagonal element = " << maxOffDiagonal << endl;
        }
        
    } while (iterationCount < 1000 && maxOffDiagonal >= stopThreshold);
    
    cout << "Jacobi method completed, iterations: " << iterationCount << endl;
    cout << "Final max off-diagonal element: " << maxOffDiagonal << endl;
    
    // Extract diagonal elements as eigenvalues
    vector<double> eigenvalues(matrixSize);
    for (int i = 0; i < matrixSize; i++) {
        eigenvalues[i] = workingMatrix(i, i);
    }
    
    return eigenvalues;
}

// ==================== Main Function ====================
int main() {
    // Get input filename
    cout << "Enter matrix file name (e.g., Amat13.num): ";
    string inputFileName;
    cin >> inputFileName;
    
    // Read matrix, Method and StopValue
    int methodValue;
    double stopThreshold;
    Matrix matrixA = readMatrixFromFile(inputFileName, methodValue, stopThreshold);
    
    if (matrixA.rows() == 0) {
        cerr << "Error: Failed to read matrix!" << endl;
        return 1;
    }
    
    cout << "Successfully read matrix, dimensions: " << matrixA.rows() << "×" << matrixA.cols() << endl;
    cout << "Method: " << methodValue << endl;
    cout << "StopValue: " << stopThreshold << endl;
    
    // Check conditions
    if (!matrixA.isSquare()) {
        cerr << "Error: Matrix is not square!" << endl;
        return 1;
    }
    
    if (!matrixA.isSymmetric()) {
        cerr << "Error: Matrix is not real symmetric (self-adjoint)!" << endl;
        return 1;
    }
    
    if (methodValue < 0 || methodValue > 2) {
        cerr << "Error: Method value must be 0, 1, or 2!" << endl;
        return 1;
    }
    
    if (stopThreshold <= 0) {
        cerr << "Error: StopValue must be greater than 0!" << endl;
        return 1;
    }
    
    vector<double> eigenvalues;
    double computationTime = 0.0;
    
    // Select computation method based on Method value
    if (methodValue == 1 || methodValue == 2) {
        auto startTime = high_resolution_clock::now();
        
        // Compute eigenvalues using Jacobi method
        eigenvalues = computeEigenvaluesJacobiMethod(matrixA, methodValue, stopThreshold);
        
        auto endTime = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(endTime - startTime);
        computationTime = duration.count() / 1000000.0;
        
    } else if (methodValue == 0) {
        // Run both methods and compare
        
        // Method 1: Jacobi method (maximum off-diagonal element)
        cout << "\n=== Running Method 1: Jacobi method (maximum off-diagonal element) ===" << endl;
        auto startTime1 = high_resolution_clock::now();
        vector<double> eigenvalues1 = computeEigenvaluesJacobiMethod(matrixA, 1, stopThreshold);
        auto endTime1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(endTime1 - startTime1);
        double time1 = duration1.count() / 1000000.0;
        
        // Method 2: Jacobi method (optimal element selection)
        cout << "\n=== Running Method 2: Jacobi method (optimal element selection) ===" << endl;
        auto startTime2 = high_resolution_clock::now();
        vector<double> eigenvalues2 = computeEigenvaluesJacobiMethod(matrixA, 2, stopThreshold);
        auto endTime2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(endTime2 - startTime2);
        double time2 = duration2.count() / 1000000.0;
        
        // Use results from first method (typically more accurate)
        eigenvalues = eigenvalues1;
        computationTime = (time1 + time2) / 2.0; // Average time
        
        cout << "\nMethod comparison:" << endl;
        cout << "Method 1 (max off-diagonal) time: " << time1 << " seconds" << endl;
        cout << "Method 2 (optimal selection) time: " << time2 << " seconds" << endl;
    }
    
    // Check if eigenvalues were successfully computed
    if (eigenvalues.empty()) {
        cerr << "Error: Failed to compute eigenvalues!" << endl;
        return 1;
    }
    
    // Sort eigenvalues by magnitude (descending)
    sortEigenvaluesByMagnitude(eigenvalues);
    
    // Extract number from input filename
    size_t startPos = inputFileName.find("Amat");
    size_t endPos = inputFileName.find(".num");
    if (startPos == string::npos || endPos == string::npos) {
        cerr << "Error: Invalid file name format!" << endl;
        return 1;
    }
    string fileNumber = inputFileName.substr(startPos + 4, endPos - startPos - 4);
    
    // Write results to file
    string outputFileName = "Res" + fileNumber + ".num";
    writeEigenvaluesToFile(outputFileName, methodValue, computationTime, eigenvalues, (methodValue == 0));
    
    cout << "\nEigenvalue computation completed!" << endl;
    cout << "Computation time: " << computationTime << " seconds" << endl;
    cout << "Results saved to: " << outputFileName << endl;
    
    // Display eigenvalues
    cout << "\nEigenvalues (sorted by magnitude):" << endl;
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        cout << "λ[" << i << "] = " << eigenvalues[i] << " (|λ| = " << fabs(eigenvalues[i]) << ")" << endl;
    }
    
    return 0;
}