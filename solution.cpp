#include <cmath>
#include "inc.h"
// 21:25. Приводим симметричную матрицу к трехдиагональному виду епта. Алгоритм прозрачный. 
// e = (1, 0, 0, ...., 0) в алгоритме по теореме из книжечки великолепной.

void solution(int n, double* matrix, double* x, double* y, double* c) {
    for (int j = 0; j < n - 2; ++j) { // n - 2 шага надо для приведения к трех.диаг виду.
        int size = n - j - 1;
        get_vector(n, j, matrix, y);
        build_x(size, x, y);
        product(size, x, y, c);
    }
}

void get_vector(int n, int j, double* matrix, double* y) {
    int ind = 0;
    for (int i = j + 1; i < n; ++i) {
        y[ind] = matrix[n*i + j];
        ind++;
    }
}

void build_x(int size, double* x, double* y) {
    double norm_y = vector_norm(size, y);
    y[0] -= norm_y;
    double norm_denominator = vector_norm(size, y);
    for (int i = 0; i < size; ++i) {
        x[i] = y[i] / norm_denominator;
    }

    y[0] += norm_y; // чтобы y стал таким каким был.
    // это супер важно, иначе дальше неправильно считаться всё будет.
}

void product(int size, double* x, double* y, double* c) {
    // крч по максимуму ща памяти заебашу. ну больше трех точно делать не буду. потом соптимизирую.
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
