#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;
using namespace chrono;

// 使用双精度，但保持性能优化
typedef double Real;
struct Complex {
    Real re, im;
    
    Complex(Real r = 0, Real i = 0) : re(r), im(i) {}
    
    inline bool isZero(Real eps = 1e-12) const {
        return re*re + im*im < eps*eps;
    }
    
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
    
    // 共轭
    inline Complex conj() const {
        return Complex(re, -im);
    }
    
    // 模平方
    inline Real norm2() const {
        return re*re + im*im;
    }
};

class Matrix {
private:
    vector<Complex> data;
    int n;
    
public:
    Matrix(int size) : n(size), data(size*size) {}
    
    inline Complex& operator()(int i, int j) { 
        return data[i*n + j]; 
    }
    
    inline const Complex& operator()(int i, int j) const { 
        return data[i*n + j]; 
    }
    
    int size() const { return n; }
    
    void fillRandom(mt19937& gen, uniform_real_distribution<Real>& dist) {
        for(auto& elem : data) {
            elem = Complex(dist(gen), dist(gen));
        }
    }
    
    // LU分解检查
    bool luDecompositionStrict(Real eps = 1e-12) {
        // 复制矩阵进行原地分解
        vector<Complex> A = data;
        Complex* matrix = A.data();
        
        for(int k = 0; k < n; ++k) {
            // 检查当前主子式是否非零
            Complex& pivot = matrix[k*n + k];
            if(pivot.isZero(eps)) return false;
            
            // 更新L列
            for(int i = k+1; i < n; ++i) {
                matrix[i*n + k] = matrix[i*n + k] / pivot;
            }
            
            // 更新子矩阵
            for(int i = k+1; i < n; ++i) {
                Complex factor = matrix[i*n + k];
                for(int j = k+1; j < n; ++j) {
                    matrix[i*n + j] = matrix[i*n + j] - factor * matrix[k*n + j];
                }
            }
        }
        return true;
    }
};

void run_lu_benchmark(int N, int K, const string& filename) {
    ofstream out(filename);
    
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<Real> dist(-10.0, 10.0);
    
    cout << "Starting LU decomposition benchmark" << endl;
    cout << "N = " << N << ", K = " << K << endl;
    cout << "Using double precision complex numbers" << endl;
    
    for(int n = 1; n <= N; ++n) {
        int success = 0;
        double total_time = 0;
        int attempts = 0;
        const int max_attempts = K * 20;
        
        while(success < K && attempts < max_attempts) {
            attempts++;
            
            // 创建并填充随机矩阵
            Matrix A(n);
            A.fillRandom(gen, dist);
            
            auto start = high_resolution_clock::now();
            bool ok = A.luDecompositionStrict();
            auto end = high_resolution_clock::now();
            
            if(ok) {
                success++;
                total_time += duration<double>(end - start).count();
            }
        }
        
        double avg_time = (success > 0) ? total_time / success : 0.0;
        out << avg_time << "\n";
        
        // 输出进度
        if(n % 100 == 0 || n <= 10) {
            cout << "n = " << n << ": " << avg_time << " s (success: " 
                 << success << "/" << K << ", attempts: " << attempts << ")" << endl;
        }
        
        // 如果连续失败，可能需要调整随机数范围
        if(success == 0 && n > 10) {
            cout << "Warning: no success for n = " << n << endl;
        }
    }
    
    out.close();
    cout << "Results saved to " << filename << endl;
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int N = 1000;
    int K = 10;
    
    // 从命令行读取参数
    if(argc >= 3) {
        try {
            N = stoi(argv[1]);
            K = stoi(argv[2]);
        } catch(...) {
            cout << "Using default parameters: N=1000, K=10" << endl;
        }
    }
    
    if(N <= 0 || K <= 0) {
        cerr << "Error: N and K must be positive integers" << endl;
        return 1;
    }
    
    auto start = high_resolution_clock::now();
    run_lu_benchmark(N, K, "results1.num");
    auto end = high_resolution_clock::now();
    
    cout << "Total time: " << duration<double>(end-start).count() << " s" << endl;
    return 0;
}