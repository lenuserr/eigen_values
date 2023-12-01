#include <cmath>
#include "inc.h"
#define EPS 1e-14

// 13:10. Приведение к трех.диаг виду готово.
// 14:10. Residual1 и Residual2 отличные, всё правильно работает.
// 16:50. QR пишу.

void solution(int n, double* matrix, double* x, double* y, double* r1, double* r2, double* r3) {
    double a_norm = matrix_norm(n, matrix);
    double trace1 = trace(n, matrix);
    double length1 = matrix_length(n, matrix);
    for (int j = 0; j < n - 2; ++j) { // n - 2 шага надо для приведения к трех.диаг виду.
        int size = n - j - 1;
        get_column(n, j, j, matrix, y);
        if (!build_x(size, x, y)) {
            continue;
        }

        // умножение слева на U(x)
        for (int k = j; k < n; ++k) {
            get_column(n, j, k, matrix, y);
            product(size, x, y, r1);
            put_column(n, j, k, matrix, r1);
        }

        // умножение справа на U(x).T = U(x)
        for (int k = j; k < n; ++k) {
           get_row(n, j, k, matrix, y);
           product(size, x, y, r1);
           put_row(n, j, k, matrix, r1);
        }
    }

    
    double trace2 = trace(n, matrix);
    [[maybe_unused]] double residual1 = std::fabs(trace1 - trace2) / a_norm;
    double length2 = matrix_length(n, matrix);
    [[maybe_unused]] double residual2 = std::fabs(length1 - length2) / a_norm;
    //printf("Residual1 = %e Residual2 = %e\n", residual1, residual2);

    // считаем матрицу R в QR разложении по методу вращений.
    // R храним в трех векторах r1, r2, r3 как в книжке.
    
    // инициализация при k = 0:
    r1[0] = matrix[0]; r2[0] = matrix[1]; r3[0] = matrix[2];
    r1[1] = matrix[n + 1]; r2[1] = matrix[n + 2];

    for (int k = 0; k < n - 1; ++k) {
        double u = r1[k];
        // это не баг, этот элемент действительно можно из matrix брать, т.к. он не меняется.
        double v = matrix[n*(k + 1) + k]; 
        double denominator = std::sqrt(u * u + v * v);
        double cos = u / denominator;
        double sin = -v / denominator;

        // храним матрицу Q в двух векторах как в книжке написано.
        x[k] = cos;
        y[k] = sin; 
        
        double a = cos * matrix[n*k + k] - sin * matrix[n*(k+1) + k]; // r_kk
        double b = cos * matrix[n*k + k + 1] - sin * matrix[n*(k + 1) + k + 1]; // r_k_k+1

        if (k + 2 < n) {
            matrix[n*k + k + 2] = -sin * matrix[n*(k + 1) + k + 2];
            matrix[n*(k + 1) + k + 2] = cos * matrix[n*(k + 1) + k + 2];
            r2[k + 1] = matrix[n*(k + 1) + k + 2];
            r3[k] = matrix[n*k + k + 2];
        }
        matrix[n*(k + 1) + k] = sin * matrix[n*k + k] + cos * matrix[n*(k + 1) + k];
        matrix[n*(k + 1) + k + 1] = sin * matrix[n*k + k + 1] + cos * matrix[n*(k + 1) + k + 1];

        matrix[n*k + k] = a;
        matrix[n*k + k + 1] = b;

        r1[k] = matrix[n*k + k];
        r1[k + 1] = matrix[n*(k + 1) + k + 1];

        r2[k] = matrix[n*k + k + 1];
    } 
    // Еще R лежит в matrix... Я это пока не фиксил, но можно.

    print_y(n, r1);
    print_y(n - 1, r2);
    print_y(n - 2, r3);
}

void get_column(int n, int j, int k, double* matrix, double* y) {
    int ind = 0;
    for (int i = j + 1; i < n; ++i) {
        y[ind] = matrix[n*i + k];
        ind++;
    }
}

void put_column(int n, int j, int k, double* matrix, double* c) {
    int ind = 0;
    for (int i = j + 1; i < n; ++i) {
        matrix[n*i + k] = c[ind];
        ind++;
    }
}

void get_row(int n, int j, int k, double* matrix, double* y) {
    int ind = 0; 
    for (int i = j + 1; i < n; ++i) {
        y[ind] = matrix[n*k + i];
        ind++;  
    }   
}

void put_row(int n, int j, int k, double* matrix, double* c) {
    int ind = 0; 
    for (int i = j + 1; i < n; ++i) {
        matrix[n*k + i] = c[ind];
        ind++;  
    }   
}

bool is_zero(int size, double* y) {
    for (int i = 0; i < size; ++i) {
        if (std::fabs(y[i]) > EPS) {
            return false;
        }
    } 

    return true;
}

bool build_x(int size, double* x, double* y) {
    bool flag = true;
    double norm_y = vector_norm(size, y);
    y[0] -= norm_y;
    
    if (!is_zero(size, y)) {
        double norm_denominator = vector_norm(size, y);
        for (int i = 0; i < size; ++i) {
            x[i] = y[i] / norm_denominator;
        }
    } else {
        flag = false;
    }

    y[0] += norm_y; // чтобы y стал таким каким был.
    // это супер важно, иначе дальше неправильно считаться всё будет.
    return flag;
}

void product(int size, double* x, double* y, double* c) {
    double scalar_prod = scalar_product(size, x, y);
    for (int i = 0; i < size; ++i) {
        c[i] = y[i] - 2*scalar_prod*x[i];        
    }
}

double vector_norm(int size, double* vector) {
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += vector[i] * vector[i];        
    }

    return std::sqrt(sum);
}

double scalar_product(int size, double* a, double* b) {
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += a[i] * b[i];    
    }

    return sum;
}

double trace(int n, double* matrix) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += matrix[n*i + i];
    }

    return sum;
}

double matrix_norm(int n, double* matrix) {
    double norm = -1;
    for (int j = 0; j < n; ++j) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += std::fabs(matrix[n*i + j]);           
        }

        norm = std::max(norm, sum);
    }

    return norm;
}

double matrix_length(int n, double* matrix) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += matrix[n*i + j] * matrix[n*i + j];
        }
    }

    return std::sqrt(sum);
}

void print_y(int size, double* y) {
    std::cout << "\ny:\n";
    for (int q = 0; q < size; ++q) {
        printf("%e ", y[q]);
    }
    std::cout << "\n";
}


/*
WRONG:
y:
1.962142e+01 3.640604e+01 4.271541e+00 1.733508e+00 9.741234e-01 6.461675e-01 4.756767e-01 3.764488e-01 3.141168e-01 2.709010e-01 

y:
3.974972e+01 4.285055e+00 1.463439e+00 7.171818e-01 4.052346e-01 2.423755e-01 1.450517e-01 8.149180e-02 3.731916e-02 

y:
3.327805e+00 1.145048e-01 1.202875e-01 6.908339e-02 3.675512e-02 1.806739e-02 7.686232e-03 2.371138e-03 

RIGHT:
y:
1.962142e+01 5.509991e+00 1.896175e+00 1.049475e+00 7.064626e-01 5.290485e-01 4.233567e-01 3.542917e-01 3.061182e-01 2.708488e-01 

y:
3.974972e+01 4.305380e+00 1.578885e+00 7.769641e-01 4.332691e-01 2.545264e-01 1.497108e-01 8.293891e-02 3.761167e-02 

y:
3.327805e+00 7.565650e-01 2.709734e-01 1.141109e-01 5.068070e-02 2.206709e-02 8.636125e-03 2.519427e-03

*/