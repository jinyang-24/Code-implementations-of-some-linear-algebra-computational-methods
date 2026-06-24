// calc_chol.cpp - 任务2：Cholesky分解（使用标准输入）
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <cmath>
#include <algorithm>
#include <string>
#include <cctype>
#include <iomanip>

using namespace std;
using namespace chrono;

// ==============================
// 类型定义和复数运算
// ==============================

typedef double Real;  // 使用双精度浮点数

struct Complex {
    Real re, im;
    
    Complex(Real r = 0, Real i = 0) : re(r), im(i) {}
    
    // 基本运算
    inline Complex operator+(const Complex& c) const { 
        return Complex(re + c.re, im + c.im); 
    }
    
    inline Complex operator-(const Complex& c) const { 
        return Complex(re - c.re, im - c.im); 
    }
    
    inline Complex operator*(const Complex& c) const {
        return Complex(re*c.re - im*c.im, re*c.im + im*c.re);
    }
    
    inline Complex operator/(const Complex& c) const {
        Real denom = c.re*c.re + c.im*c.im;
        return Complex(
            (re*c.re + im*c.im)/denom,
            (im*c.re - re*c.im)/denom
        );
    }
    
    inline Complex operator/(Real d) const { 
        return Complex(re/d, im/d); 
    }
    
    inline Complex operator*(Real d) const { 
        return Complex(re*d, im*d); 
    }
    
    // 共轭
    inline Complex conj() const {
        return Complex(re, -im);
    }
    
    // 模的平方
    inline Real norm2() const {
        return re*re + im*im;
    }
    
    // 模
    inline Real abs() const {
        return sqrt(norm2());
    }
    
    // 是否接近零（考虑数值误差）
    inline bool isZero(Real eps = 1e-12) const {
        return norm2() < eps*eps;
    }
    
    // 是否为正实数（实部为正，虚部接近0）
    inline bool isPositiveReal(Real eps = 1e-12) const {
        return re > 0 && im*im < eps*eps;
    }
};

// ==============================
// 矩阵类
// ==============================

class Matrix {
private:
    vector<Complex> data;
    int n;  // 方阵尺寸
    
public:
    Matrix() : n(0) {}
    
    Matrix(int size) : n(size), data(size*size) {}
    
    // 从文件读取矩阵
    bool loadFromFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open file " << filename << endl;
            return false;
        }
        
        // 读取文件内容
        stringstream buffer;
        buffer << file.rdbuf();
        file.close();
        
        string content = buffer.str();
        
        // 解析文件
        size_t start_pos = content.find("A = ...");
        if (start_pos == string::npos) {
            cerr << "Error: Cannot find 'A = ...' in file" << endl;
            return false;
        }
        
        start_pos += 7;  // 跳过"A = ..."
        
        // 检查是否为复数矩阵
        size_t complex_pos = content.find("complex", start_pos);
        
        if (complex_pos != string::npos && complex_pos < content.find('[', start_pos)) {
            // 复数矩阵
            return parseComplexMatrix(content, complex_pos);
        } else {
            // 实数矩阵
            return parseRealMatrix(content, start_pos);
        }
    }
    
    // 保存矩阵到文件 - 按照Cmat1.num格式
    bool saveToFile(const string& filename, const string& matrixName = "C") const {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open file " << filename << " for writing" << endl;
            return false;
        }
        
        file << matrixName << " = ... \n";
        
        // 设置输出格式：固定小数点，精度6位
        file << fixed << setprecision(6);
        
        // 检查是否为实数矩阵（虚部全为0）
        bool isReal = true;
        for (const auto& elem : data) {
            if (abs(elem.im) > 1e-12) {
                isReal = false;
                break;
            }
        }
        
        if (isReal) {
            // 保存为实数矩阵
            file << "[";
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    file << (*this)(i, j).re;
                    if (j < n-1) file << " ";
                }
                file << ";" << endl;
            }
            file << "]" << endl;
        } else {
            // 保存为复数矩阵 - 严格按照Cmat1.num格式
            file << "complex([";
            for (int i = 0; i < n; ++i) {
                // 实部行
                for (int j = 0; j < n; ++j) {
                    file << (*this)(i, j).re;
                    if (j < n-1) file << " ";
                }
                file << ";";
                if (i < n-1) file << endl;  // 除了最后一行，每行后换行
            }
            file << "]," << endl << " [";
            
            // 虚部矩阵
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    file << (*this)(i, j).im;
                    if (j < n-1) file << " ";
                }
                file << ";";
                if (i < n-1) file << endl;  // 除了最后一行，每行后换行
            }
            file << "]);" << endl;  // 注意：按照Cmat1.num格式，这里应该有分号
        }
        
        file.close();
        return true;
    }
    
    // 矩阵元素访问
    inline Complex& operator()(int i, int j) { 
        return data[i*n + j]; 
    }
    
    inline const Complex& operator()(int i, int j) const { 
        return data[i*n + j]; 
    }
    
    int size() const { return n; }
    
    // 检查矩阵是否为方阵
    bool isSquare() const {
        return n > 0;  // 我们只存储方阵
    }
    
    // 检查矩阵是否为厄米特矩阵（A = A^H）
    bool isHermitian(Real eps = 1e-10) const {
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                Complex diff = (*this)(i, j) - (*this)(j, i).conj();
                if (diff.norm2() > eps*eps) {
                    return false;
                }
            }
        }
        return true;
    }
    
    // 计算Frobenius范数
    Real frobeniusNorm() const {
        Real sum = 0;
        for (const auto& elem : data) {
            sum += elem.norm2();
        }
        return sqrt(sum);
    }
    
    // 矩阵乘法
    Matrix multiply(const Matrix& B) const {
        if (n != B.n) {
            cerr << "Error: Matrix size mismatch for multiplication" << endl;
            return Matrix();
        }
        
        Matrix result(n);
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Complex sum(0, 0);
                for (int k = 0; k < n; ++k) {
                    sum = sum + (*this)(i, k) * B(k, j);
                }
                result(i, j) = sum;
            }
        }
        
        return result;
    }
    
    // 矩阵减法
    Matrix subtract(const Matrix& B) const {
        if (n != B.n) {
            cerr << "Error: Matrix size mismatch for subtraction" << endl;
            return Matrix();
        }
        
        Matrix result(n);
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result(i, j) = (*this)(i, j) - B(i, j);
            }
        }
        
        return result;
    }
    
    // 共轭转置
    Matrix conjugateTranspose() const {
        Matrix result(n);
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result(j, i) = (*this)(i, j).conj();
            }
        }
        
        return result;
    }
    
    // Cholesky分解（返回下三角矩阵C，使得A = C * C^H）
    Matrix choleskyDecomposition(Real& error, Real& timeTaken, bool& success, Real eps = 1e-12) {
        success = false;
        Matrix C(n);
        
        auto start = high_resolution_clock::now();
        
        for (int j = 0; j < n; ++j) {
            // 计算C[j][j]
            Complex sum(0, 0);
            for (int k = 0; k < j; ++k) {
                sum = sum + C(j, k) * C(j, k).conj();
            }
            
            Complex diag = (*this)(j, j) - sum;
            
            // 检查对角线元素是否为正值实数
            if (!diag.isPositiveReal(eps) || diag.re <= eps) {
                cerr << "Error: Matrix is not positive definite at position (" 
                     << j << ", " << j << ")" << endl;
                return Matrix();
            }
            
            C(j, j) = Complex(sqrt(diag.re), 0);
            
            // 计算C[i][j] (i > j)
            for (int i = j + 1; i < n; ++i) {
                sum = Complex(0, 0);
                for (int k = 0; k < j; ++k) {
                    sum = sum + C(i, k) * C(j, k).conj();
                }
                
                C(i, j) = ((*this)(i, j) - sum) / C(j, j);
            }
        }
        
        auto end = high_resolution_clock::now();
        timeTaken = duration<double>(end - start).count();
        
        // 验证分解：计算误差 ||A - C * C^H||
        Matrix C_conj = C.conjugateTranspose();
        Matrix CC_H = C.multiply(C_conj);
        Matrix diff = this->subtract(CC_H);
        error = diff.frobeniusNorm();
        
        success = true;
        return C;
    }
    
private:
    // 解析实数矩阵
    bool parseRealMatrix(const string& content, size_t start_pos) {
        size_t bracket_start = content.find('[', start_pos);
        if (bracket_start == string::npos) {
            cerr << "Error: Cannot find '[' in matrix definition" << endl;
            return false;
        }
        
        size_t bracket_end = content.find(']', bracket_start);
        if (bracket_end == string::npos) {
            cerr << "Error: Cannot find ']' in matrix definition" << endl;
            return false;
        }
        
        string matrix_str = content.substr(bracket_start + 1, bracket_end - bracket_start - 1);
        
        // 按行分割
        vector<string> rows;
        size_t pos = 0;
        while (pos < matrix_str.length()) {
            size_t semicolon_pos = matrix_str.find(';', pos);
            if (semicolon_pos == string::npos) {
                rows.push_back(matrix_str.substr(pos));
                break;
            }
            rows.push_back(matrix_str.substr(pos, semicolon_pos - pos));
            pos = semicolon_pos + 1;
        }
        
        n = rows.size();
        data.resize(n * n);
        
        // 解析每一行
        for (int i = 0; i < n; ++i) {
            stringstream row_ss(rows[i]);
            for (int j = 0; j < n; ++j) {
                Real val;
                row_ss >> val;
                (*this)(i, j) = Complex(val, 0);
            }
        }
        
        return true;
    }
    
    // 解析复数矩阵
    bool parseComplexMatrix(const string& content, size_t complex_pos) {
        // 找到第一个'['（实部开始）
        size_t real_start = content.find('[', complex_pos);
        if (real_start == string::npos) {
            cerr << "Error: Cannot find '[' for real part" << endl;
            return false;
        }
        
        size_t real_end = content.find(']', real_start);
        if (real_end == string::npos) {
            cerr << "Error: Cannot find ']' for real part" << endl;
            return false;
        }
        
        // 找到第二个'['（虚部开始）
        size_t imag_start = content.find('[', real_end);
        if (imag_start == string::npos) {
            cerr << "Error: Cannot find '[' for imaginary part" << endl;
            return false;
        }
        
        size_t imag_end = content.find(']', imag_start);
        if (imag_end == string::npos) {
            cerr << "Error: Cannot find ']' for imaginary part" << endl;
            return false;
        }
        
        string real_str = content.substr(real_start + 1, real_end - real_start - 1);
        string imag_str = content.substr(imag_start + 1, imag_end - imag_start - 1);
        
        // 解析实部
        vector<vector<Real>> real_data;
        if (!parseMatrixString(real_str, real_data)) {
            return false;
        }
        
        // 解析虚部
        vector<vector<Real>> imag_data;
        if (!parseMatrixString(imag_str, imag_data)) {
            return false;
        }
        
        // 检查尺寸是否匹配
        if (real_data.size() != imag_data.size() || real_data.empty()) {
            cerr << "Error: Real and imaginary parts have different dimensions" << endl;
            return false;
        }
        
        n = real_data.size();
        data.resize(n * n);
        
        for (int i = 0; i < n; ++i) {
            if (real_data[i].size() != n || imag_data[i].size() != n) {
                cerr << "Error: Inconsistent row size in matrix" << endl;
                return false;
            }
            
            for (int j = 0; j < n; ++j) {
                (*this)(i, j) = Complex(real_data[i][j], imag_data[i][j]);
            }
        }
        
        return true;
    }
    
    // 解析矩阵字符串（辅助函数）
    bool parseMatrixString(const string& matrix_str, vector<vector<Real>>& result) {
        result.clear();
        
        vector<string> rows;
        size_t pos = 0;
        while (pos < matrix_str.length()) {
            size_t semicolon_pos = matrix_str.find(';', pos);
            if (semicolon_pos == string::npos) {
                rows.push_back(matrix_str.substr(pos));
                break;
            }
            rows.push_back(matrix_str.substr(pos, semicolon_pos - pos));
            pos = semicolon_pos + 1;
        }
        
        for (const auto& row_str : rows) {
            vector<Real> row;
            stringstream row_ss(row_str);
            Real val;
            
            while (row_ss >> val) {
                row.push_back(val);
            }
            
            if (!row.empty()) {
                result.push_back(row);
            }
        }
        
        return !result.empty();
    }
};

// ==============================
// 主程序（使用标准输入）
// ==============================

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cout << "========================================" << endl;
    cout << "Cholesky Decomposition Program" << endl;
    cout << "========================================" << endl;
    cout << "Enter the input matrix filename: ";
    
    string inputFile;
    getline(cin, inputFile);
    
    // 移除可能的引号或空格
    inputFile.erase(remove_if(inputFile.begin(), inputFile.end(), 
                             [](char c) { return c == '"' || c == '\''; }), 
                   inputFile.end());
    
    // 去除首尾空格
    size_t start = inputFile.find_first_not_of(" \t");
    size_t end = inputFile.find_last_not_of(" \t");
    if (start != string::npos && end != string::npos) {
        inputFile = inputFile.substr(start, end - start + 1);
    }
    
    // 如果用户直接按回车，使用默认文件名
    if (inputFile.empty()) {
        inputFile = "Amat2.num";
        cout << "No input provided, using default: " << inputFile << endl;
    }
    
    // 检查文件名是否包含预期的后缀
    if (inputFile.find(".num") == string::npos) {
        cout << "Note: Input file should have .num extension" << endl;
    }
    
    // 从文件名提取数字部分（用于输出文件命名）
    size_t start_pos = inputFile.find("Amat");
    string number = "";
    
    if (start_pos != string::npos) {
        size_t dot_pos = inputFile.find(".num", start_pos);
        if (dot_pos != string::npos) {
            number = inputFile.substr(start_pos + 4, dot_pos - start_pos - 4);
        }
    }
    
    if (number.empty()) {
        // 如果不能从文件名提取数字，提示用户输入
        cout << "Could not extract matrix number from filename." << endl;
        cout << "Please enter the matrix number (e.g., '2' for Amat2.num): ";
        getline(cin, number);
        
        if (number.empty()) {
            number = "2";  // 默认值
            cout << "Using default number: " << number << endl;
        }
    }
    
    cout << "\nInput file: " << inputFile << endl;
    cout << "Matrix number: " << number << endl;
    cout << "========================================\n" << endl;
    
    // 加载矩阵
    cout << "Loading matrix from file..." << endl;
    Matrix A;
    
    if (!A.loadFromFile(inputFile)) {
        cerr << "Error: Failed to load matrix from file" << endl;
        return 1;
    }
    
    cout << "Matrix loaded successfully." << endl;
    cout << "Matrix size: " << A.size() << "x" << A.size() << endl;
    
    // 检查矩阵是否为方阵
    if (!A.isSquare()) {
        cerr << "Error: Matrix is not square" << endl;
        return 1;
    }
    
    cout << "\nChecking if matrix is Hermitian..." << endl;
    
    // 检查矩阵是否为厄米特矩阵
    if (!A.isHermitian()) {
        cerr << "Error: Matrix is not Hermitian" << endl;
        return 1;
    }
    
    cout << "Matrix is Hermitian." << endl;
    
    // 执行Cholesky分解
    cout << "\nPerforming Cholesky decomposition..." << endl;
    
    Real decompositionError, decompositionTime;
    bool success;
    
    Matrix C = A.choleskyDecomposition(decompositionError, decompositionTime, success);
    
    if (!success) {
        cerr << "Error: Cholesky decomposition failed" << endl;
        cerr << "The matrix may not be positive definite" << endl;
        return 1;
    }
    
    cout << "Cholesky decomposition completed successfully." << endl;
    cout << "Decomposition time: " << fixed << setprecision(10) << decompositionTime << " seconds" << endl;
    cout << "Error ||A - C*C^H||_F: " << fixed << setprecision(10) << decompositionError << endl;
    
    // 保存结果矩阵C
    string outputCFile = "Cmat" + number + ".num";
    cout << "\nSaving result matrix C to " << outputCFile << "..." << endl;
    
    if (!C.saveToFile(outputCFile, "C")) {
        cerr << "Error: Failed to save result matrix C" << endl;
        return 1;
    }
    
    cout << "Result matrix C saved successfully." << endl;
    
    // 保存结果文件（时间和误差）- 设置为10位小数
    string outputResFile = "Res" + number + ".num";
    cout << "Saving result data to " << outputResFile << "..." << endl;
    
    ofstream resFile(outputResFile);
    if (!resFile.is_open()) {
        cerr << "Error: Cannot open result file " << outputResFile << endl;
        return 1;
    }
    
    // 设置输出精度为10位小数
    resFile << fixed << setprecision(10);
    resFile << "Время вычисления разложения Холецкого: " << decompositionTime << " секунд." << endl;
    resFile << "Погрешность полученного разложения Холецкого: " << decompositionError << "." << endl;
    resFile.close();
    
    cout << "Result data saved successfully." << endl;
    cout << "\n========================================" << endl;
    cout << "Cholesky decomposition completed." << endl;
    cout << "Generated files:" << endl;
    cout << "  - " << outputCFile << " (matrix C)" << endl;
    cout << "  - " << outputResFile << " (results)" << endl;
    cout << "========================================" << endl;
    
    return 0;
}