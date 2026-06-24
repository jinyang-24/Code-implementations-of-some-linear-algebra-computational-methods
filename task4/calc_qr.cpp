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

// ==================== 矩阵类 ====================
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
    
    // 矩阵乘法
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
    
    // 矩阵转置
    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(j, i) = data_[i][j];
            }
        }
        return result;
    }
    
    // 矩阵减法
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
    
    // 矩阵加法
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
    
    // 标量乘法
    Matrix scalarMultiply(double scalar) const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j] * scalar;
            }
        }
        return result;
    }
    
    // 获取列向量
    Matrix getColumn(int j) const {
        Matrix result(rows_, 1);
        for (int i = 0; i < rows_; i++) {
            result(i, 0) = data_[i][j];
        }
        return result;
    }
    
    // 设置列向量
    void setColumn(int j, const Matrix& col) {
        for (int i = 0; i < rows_; i++) {
            data_[i][j] = col(i, 0);
        }
    }
    
    // Frobenius范数
    double frobeniusNorm() const {
        double sum = 0.0;
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                sum += data_[i][j] * data_[i][j];
            }
        }
        return sqrt(sum);
    }
    
    // 检查是否为方阵
    bool isSquare() const {
        return rows_ == cols_;
    }
    
    // 生成单位矩阵
    static Matrix identity(int n) {
        Matrix result(n, n);
        for (int i = 0; i < n; i++) {
            result(i, i) = 1.0;
        }
        return result;
    }
    
    // 复制矩阵
    Matrix copy() const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                result(i, j) = data_[i][j];
            }
        }
        return result;
    }
    
    // 打印矩阵
    void print() const {
        cout << fixed << setprecision(6);
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                cout << setw(12) << data_[i][j] << " ";
            }
            cout << endl;
        }
    }
    
    // 获取子矩阵
    Matrix submatrix(int r1, int r2, int c1, int c2) const {
        Matrix result(r2 - r1 + 1, c2 - c1 + 1);
        for (int i = r1; i <= r2; i++) {
            for (int j = c1; j <= c2; j++) {
                result(i - r1, j - c1) = data_[i][j];
            }
        }
        return result;
    }
    
    // 设置子矩阵
    void setSubmatrix(int r, int c, const Matrix& sub) {
        for (int i = 0; i < sub.rows(); i++) {
            for (int j = 0; j < sub.cols(); j++) {
                data_[r + i][c + j] = sub(i, j);
            }
        }
    }
    
    // 获取行向量
    Matrix getRow(int i) const {
        Matrix result(1, cols_);
        for (int j = 0; j < cols_; j++) {
            result(0, j) = data_[i][j];
        }
        return result;
    }
    
    // 设置行向量
    void setRow(int i, const Matrix& row) {
        for (int j = 0; j < cols_; j++) {
            data_[i][j] = row(0, j);
        }
    }
};

// ==================== 文件读写函数 ====================
Matrix readMatrixFromFile(const string& filename, int& method) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        return Matrix(0, 0);
    }
    
    string line;
    
    // 读取第一行：Method
    getline(file, line);
    size_t eqPos = line.find('=');
    if (eqPos == string::npos) {
        cerr << "无效的文件格式: 第一行应包含Method=" << endl;
        return Matrix(0, 0);
    }
    
    // 提取Method值
    size_t semicolonPos = line.find(';');
    if (semicolonPos == string::npos) {
        cerr << "无效的文件格式: Method行应以分号结尾" << endl;
        return Matrix(0, 0);
    }
    
    string methodStr = line.substr(eqPos + 1, semicolonPos - eqPos - 1);
    method = stoi(methodStr);
    
    // 读取第二行：A = ...
    getline(file, line);
    
    // 读取矩阵数据
    vector<vector<double>> matrixData;
    while (getline(file, line)) {
        // 移除行首尾的空格
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        
        // 跳过空行
        if (line.empty()) continue;
        
        // 检查是否结束
        if (line[0] == ']') break;
        
        // 移除行首的'['或空格
        if (line[0] == '[') {
            line = line.substr(1);
        }
        
        // 移除行尾的';'
        if (line.back() == ';') {
            line.pop_back();
        }
        
        // 解析一行数字
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
    
    // 检查所有行是否有相同的列数
    int cols = matrixData[0].size();
    for (size_t i = 1; i < matrixData.size(); i++) {
        if (matrixData[i].size() != cols) {
            cerr << "矩阵行长度不一致!" << endl;
            return Matrix(0, 0);
        }
    }
    
    return Matrix(matrixData);
}

void writeMatrixToFile(const string& filename, const string& matrixName, const Matrix& mat) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "无法创建文件: " << filename << endl;
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

// ==================== 辅助函数 ====================

// 向量点积
double dotProduct(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != 1 || b.cols() != 1) {
        cerr << "向量维度不匹配!" << endl;
        return 0.0;
    }
    
    double sum = 0.0;
    for (int i = 0; i < a.rows(); i++) {
        sum += a(i, 0) * b(i, 0);
    }
    return sum;
}

// 向量范数
double vectorNorm(const Matrix& v) {
    return sqrt(dotProduct(v, v));
}

// ==================== QR分解算法 ====================

// 方法1：通过Cholesky分解计算QR（讲义9第2页）
bool qrCholesky(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    
    // 检查矩阵是否非奇异（通过Cholesky分解）
    Matrix AT = A.transpose();
    Matrix ATA = AT.multiply(A);
    
    // 进行Cholesky分解：ATA = L * L^T，其中L是下三角
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
                    // 矩阵不正定，无法使用Cholesky分解
                    return false;
                }
                L(i, i) = sqrt(diag);
            } else {
                L(i, j) = (ATA(i, j) - sum) / L(j, j);
            }
        }
    }
    
    // 现在ATA = L * L^T，令R = L^T（上三角）
    R = L.transpose();
    
    // 解方程 R^T * Q^T = A^T 得到Q
    // 等价于解 L * Q^T = A^T
    Matrix QT = AT.copy();  // A^T
    
    // 前向替换解 L * X = A^T
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

// 方法2：Householder反射法（讲义9第3-4页）
void qrHouseholder(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    R = A.copy();
    Q = Matrix::identity(n);
    
    for (int k = 0; k < n - 1; k++) {
        // 计算x = R[k:n, k]
        int m = n - k;
        Matrix x(m, 1);
        for (int i = 0; i < m; i++) {
            x(i, 0) = R(k + i, k);
        }
        
        // 计算向量的范数
        double norm_x = vectorNorm(x);
        
        if (norm_x < 1e-15) {
            continue; // 列已经为零，跳过
        }
        
        // 计算反射向量u
        Matrix u = x.copy();
        
        // 选择符号以避免数值不稳定
        double sign = (x(0, 0) >= 0) ? 1.0 : -1.0;
        u(0, 0) += sign * norm_x;
        
        // 归一化u
        double norm_u = vectorNorm(u);
        
        if (norm_u < 1e-15) {
            continue;
        }
        
        for (int i = 0; i < m; i++) {
            u(i, 0) /= norm_u;
        }
        
        // 应用Householder变换到R的右下部分
        for (int j = k; j < n; j++) {
            // 计算u^T * R[k:n, j]
            double u_dot_r = 0.0;
            for (int i = 0; i < m; i++) {
                u_dot_r += u(i, 0) * R(k + i, j);
            }
            
            // 更新R[k:n, j]
            for (int i = 0; i < m; i++) {
                R(k + i, j) -= 2.0 * u(i, 0) * u_dot_r;
            }
        }
        
        // 累积Q：Q = Q * H
        for (int i = 0; i < n; i++) {
            // 计算Q[i, k:n] * u
            double q_dot_u = 0.0;
            for (int j = 0; j < m; j++) {
                q_dot_u += Q(i, k + j) * u(j, 0);
            }
            
            // 更新Q[i, k:n]
            for (int j = 0; j < m; j++) {
                Q(i, k + j) -= 2.0 * q_dot_u * u(j, 0);
            }
        }
    }
}

// 方法3：Givens旋转法（讲义9第5页）- 修复版本
void qrGivens(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    R = A.copy();
    Q = Matrix::identity(n);
    
    // 根据讲义9第5页，Givens旋转应该按特定顺序消元
    // 通常是从左下角开始，按列消元
    for (int j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double a = R(j, j);
            double b = R(i, j);
            
            // 如果已经为零，跳过
            if (fabs(b) < 1e-15) {
                continue;
            }
            
            // 计算旋转参数
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = -b / r;  // 注意符号，讲义中使用负号
            
            // 应用旋转到R的第j行和第i行
            for (int k = j; k < n; k++) {
                double r_jk = R(j, k);
                double r_ik = R(i, k);
                R(j, k) = c * r_jk - s * r_ik;  // 注意这里是减号
                R(i, k) = s * r_jk + c * r_ik;
            }
            
            // 将R(i, j)设为零（消除舍入误差）
            R(i, j) = 0.0;
            
            // 应用相同的旋转到Q
            // Q = Q * G^T，其中G是旋转矩阵
            for (int k = 0; k < n; k++) {
                double q_kj = Q(k, j);
                double q_ki = Q(k, i);
                Q(k, j) = c * q_kj - s * q_ki;
                Q(k, i) = s * q_kj + c * q_ki;
            }
        }
    }
    
    // 确保R的对角线元素为非负（可选，但有助于稳定性）
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

// 方法4：Gram-Schmidt正交化（讲义9第7-8页）
void qrGramSchmidt(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.rows();
    Q = Matrix(n, n);
    R = Matrix(n, n);
    
    for (int j = 0; j < n; j++) {
        // 获取A的第j列
        Matrix v = A.getColumn(j);
        
        // 正交化
        for (int i = 0; i < j; i++) {
            R(i, j) = dotProduct(Q.getColumn(i), A.getColumn(j));
            
            // v = v - R(i,j) * Q[:,i]
            for (int k = 0; k < n; k++) {
                v(k, 0) -= R(i, j) * Q(k, i);
            }
        }
        
        // 计算R(j,j)
        R(j, j) = vectorNorm(v);
        
        if (R(j, j) < 1e-15) {
            // 线性相关，设置q_j为零向量
            for (int k = 0; k < n; k++) {
                Q(k, j) = 0.0;
            }
        } else {
            // 归一化得到q_j
            for (int k = 0; k < n; k++) {
                Q(k, j) = v(k, 0) / R(j, j);
            }
        }
    }
}

// ==================== 验证QR分解函数 ====================
void verifyQR(const Matrix& A, const Matrix& Q, const Matrix& R, const string& methodName) {
    cout << "\n验证 " << methodName << " 方法:" << endl;
    
    // 1. 验证A = QR
    Matrix QR = Q.multiply(R);
    Matrix diff_A_QR = A.subtract(QR);
    double error_A_QR = diff_A_QR.frobeniusNorm();
    cout << "  ||A - QR||_E = " << scientific << setprecision(6) << error_A_QR;
    cout << fixed << " (相对误差: " << error_A_QR / A.frobeniusNorm() << ")" << endl;
    
    // 2. 验证Q的正交性：Q^T * Q = I
    Matrix QT = Q.transpose();
    Matrix QTQ = QT.multiply(Q);
    Matrix I = Matrix::identity(A.rows());
    Matrix diff_QTQ_I = QTQ.subtract(I);
    double error_orth = diff_QTQ_I.frobeniusNorm();
    cout << "  ||Q^T * Q - I||_E = " << scientific << error_orth << endl;
    
    // 3. 验证R是上三角矩阵
    bool isUpperTriangular = true;
    for (int i = 1; i < R.rows(); i++) {
        for (int j = 0; j < i; j++) {
            if (fabs(R(i, j)) > 1e-10) {
                isUpperTriangular = false;
                cout << "  发现非零下三角元素: R(" << i << "," << j << ") = " << R(i, j) << endl;
                break;
            }
        }
        if (!isUpperTriangular) break;
    }
    cout << "  R是上三角矩阵: " << (isUpperTriangular ? "是" : "否") << endl;
    
    // 4. 显示R矩阵的对角线元素
    cout << "  R对角线元素: ";
    for (int i = 0; i < min(5, R.rows()); i++) {
        cout << R(i, i) << " ";
    }
    if (R.rows() > 5) cout << "...";
    cout << endl;
}

// ==================== 主函数 ====================
int main() {
    // 获取输入文件名
    cout << "请输入输入文件名 (例如: Amat8.num): ";
    string inputFilename;
    cin >> inputFilename;
    
    // 从输入文件名提取数字部分
    size_t start = inputFilename.find("Amat");
    size_t end = inputFilename.find(".num");
    if (start == string::npos || end == string::npos) {
        cerr << "无效的文件名格式!" << endl;
        return 1;
    }
    string num = inputFilename.substr(start + 4, end - start - 4);
    
    // 读取矩阵和方法
    int method;
    Matrix A = readMatrixFromFile(inputFilename, method);
    
    if (A.rows() == 0) {
        cerr << "读取矩阵失败!" << endl;
        return 1;
    }
    
    // 检查矩阵是否为方阵
    if (!A.isSquare()) {
        cerr << "错误: 矩阵不是方阵!" << endl;
        return 1;
    }
    
    // 检查方法值是否有效
    if (method < 0 || method > 4) {
        cerr << "错误: 无效的方法值! 应为0-4" << endl;
        return 1;
    }
    
    // 创建结果文件
    ofstream resFile("Res" + num + ".num");
    if (!resFile.is_open()) {
        cerr << "无法创建结果文件!" << endl;
        return 1;
    }
    
    // 根据方法执行QR分解
    if (method == 0) {
        // 使用所有四种方法
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
                resFile << "Погрешность полученного QR-разложения: Cholesky分解失败（矩阵不正定）." << endl;
                resFile << "Погрешность ортогональности матрицы Q: N/A." << endl;
                
                cout << "\n方法 " << m << " (Cholesky) 失败: 矩阵不正定" << endl;
                continue;
            }
            
            // 计算误差：||A - QR||_E
            Matrix QR = Q.multiply(R);
            Matrix diff = A.subtract(QR);
            double qrError = diff.frobeniusNorm();
            
            // 计算Q的正交性误差：||Q^T * Q - I||_E
            Matrix QT = Q.transpose();
            Matrix QTQ = QT.multiply(Q);
            Matrix I = Matrix::identity(A.rows());
            Matrix orthDiff = QTQ.subtract(I);
            double orthError = orthDiff.frobeniusNorm();
            
            // 写入结果
            resFile << "Method = " << m << "." << endl;
            resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
            resFile << "Погрешность полученного QR-разложения: " << qrError << "." << endl;
            resFile << "Погрешность ортогональности матрицы Q: " << orthError << "." << endl;
            
            // 写入Q和R矩阵到文件
            writeMatrixToFile("Qmat" + num + "_" + to_string(m) + ".num", "Q", Q);
            writeMatrixToFile("Rmat" + num + "_" + to_string(m) + ".num", "R", R);
            
            // 验证QR分解
            verifyQR(A, Q, R, methodName);
        }
    } else {
        // 使用指定方法
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
            resFile << "Погрешность полученного QR-разложения: Cholesky分解失败（矩阵不正定）." << endl;
            resFile << "Погрешность ортогональности матрицы Q: N/A." << endl;
            
            cout << "\n方法 " << method << " (Cholesky) 失败: 矩阵不正定" << endl;
        } else {
            // 计算误差
            Matrix QR = Q.multiply(R);
            Matrix diff = A.subtract(QR);
            double qrError = diff.frobeniusNorm();
            
            // 计算Q的正交性误差
            Matrix QT = Q.transpose();
            Matrix QTQ = QT.multiply(Q);
            Matrix I = Matrix::identity(A.rows());
            Matrix orthDiff = QTQ.subtract(I);
            double orthError = orthDiff.frobeniusNorm();
            
            // 写入结果
            resFile << fixed << setprecision(15);
            resFile << "Method = " << method << "." << endl;
            resFile << "Время вычисления QR-разложения: " << timeSeconds << " секунд." << endl;
            resFile << "Погрешность полученного QR-разложения: " << qrError << "." << endl;
            resFile << "Погрешность ортогональности матрицы Q: " << orthError << "." << endl;
            
            // 写入Q和R矩阵到文件
            writeMatrixToFile("Qmat" + num + ".num", "Q", Q);
            writeMatrixToFile("Rmat" + num + ".num", "R", R);
            
            // 验证QR分解
            verifyQR(A, Q, R, methodName);
        }
    }
    
    resFile.close();
    
    cout << "\nQR分解完成! 结果已保存到文件。" << endl;
    cout << "注意: 误差信息已写入Res" << num << ".num文件。" << endl;
    
    return 0;
}