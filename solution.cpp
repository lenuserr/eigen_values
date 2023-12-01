#include <cmath>
#include "inc.h"
#define EPS 1e-14

// 13:10. Приведение к трех.диаг виду готово.

void solution(int n, double* matrix, double* x, double* y, double* c) {
    for (int j = 0; j < n - 2; ++j) { // n - 2 шага надо для приведения к трех.диаг виду.
        int size = n - j - 1;
        get_column(n, j, j, matrix, y);
        if (!build_x(size, x, y)) {
            continue;
        }

        // умножение слева на U(x)
        for (int k = j; k < n; ++k) {
            get_column(n, j, k, matrix, y);
            product(size, x, y, c);
            put_column(n, j, k, matrix, c);
        }

        // умножение справа на U(x).T = U(x)
        for (int k = j; k < n; ++k) {
           get_row(n, j, k, matrix, y);
           product(size, x, y, c);
           put_row(n, j, k, matrix, c);
        }
    }
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

void print_y(int size, double* y) {
    std::cout << "\ny:\n";
    for (int q = 0; q < size; ++q) {
        std::cout << y[q] << " ";
    }
    std::cout << "\n";
}
