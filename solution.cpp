#include <sys/resource.h>
#include "inc.h"

int solution(int n, double a_norm, double* matrix, double* x, double* y,
 double* r1, double* r2, double* r3, double* lambda, double EPS, double* t1, double* t2) {
    *t1 = get_cpu_time();
    for (int j = 0; j < n - 2; ++j) { 
        int size = n - j - 1;
        get_column(n, j, j, matrix, y);
        if (!build_x(size, x, y, a_norm, EPS)) {
            continue;
        }

        for (int k = j; k < n; ++k) {
            get_column(n, j, k, matrix, y);
            product(size, x, y, r1);
            put_column(n, j, k, matrix, r1);
        }

        for (int k = j; k < n; ++k) {
           get_row(n, j, k, matrix, y);
           product(size, x, y, r1);
           put_row(n, j, k, matrix, r1);
        }
    }

    *t1 = get_cpu_time() - *t1;

    *t2 = get_cpu_time();
    for (int k = 0; k < n; ++k) {
        r1[k] = matrix[n*k + k]; 
        if (k < n - 1) {
            r2[k] = matrix[n*k + k + 1]; 
            r3[k] = matrix[n*(k + 1) + k]; 
        }
    }

    double s = 0;
    int k = n;
    int its = 0;
    while (k > 2) { 
        while (std::fabs(r3[k - 2]) >= EPS * a_norm) {
            if ((r1[k - 1] > EPS * a_norm && r3[k - 2] > EPS * a_norm) || (r1[k - 1] < EPS * a_norm && r3[k - 2] < EPS * a_norm)) {
                s = r1[k - 1] + 0.5 * r3[k - 2];
            } else {
                s = r1[k - 1] - 0.5 * r3[k - 2];
            }

            
            for (int i = 0; i < k; ++i) {
                r1[i] -= s;
            }
            QR(k, x, y, r1, r2, r3);
            RQ_product(k, x, y, r1, r2, r3);
            for (int i = 0; i < k; ++i) {
                r1[i] += s;
            }

            its++;
        }
        
        while (k > 2 && std::fabs(r3[k - 2]) < EPS * a_norm) {
            lambda[k - 1] = r1[k - 1];
            --k;
        }
    }

    *t2 = get_cpu_time() - *t2;

    if (n == 1) {
        lambda[0] = matrix[0];
        return 0;
    }

    double D = r1[0]*r1[0] + r1[1]*r1[1] + 4*r3[0]*r3[0] -2*r1[0]*r1[1];
    if (r1[0] + r1[1] > EPS * a_norm) {
        lambda[0] = (r1[0] + r1[1] + std::sqrt(D)) / 2;
    } else {
        lambda[0] = (r1[0] + r1[1] - std::sqrt(D)) / 2;
    }

    lambda[1] = (r1[0]*r1[1] - r3[0]*r3[0]) / lambda[0];

    return its;
}

void QR(int n, double* x, double* y, double* r1, double* r2, double* r3) {
    for (int k = 0; k < n - 1; ++k) {
        double u = r1[k];
        double v = r3[k];
        double denominator = std::sqrt(u * u + v * v);
        double cos = u / denominator;
        double sin = -v / denominator; 

        x[k] = cos;
        y[k] = sin; 
        
        double a = cos * r1[k] - sin * r3[k]; 
        double b = cos * r2[k] - sin * r1[k + 1];

        if (k + 2 < n) {
            r2[k + 1] *= cos;
        }

        r3[k] = sin * r1[k] + cos * r3[k];
        r1[k + 1] = sin * r2[k] + cos * r1[k + 1];
        r1[k] = a;
        r2[k] = b;
    }
}

void RQ_product(int n, double* x, double* y, double* r1, double* r2, double* r3) {
    for (int k = 0; k < n - 1; ++k) {
        r1[k] = r1[k] * x[k] - r2[k] * y[k];
        r3[k] = -r1[k + 1] * y[k];
        r1[k + 1] *= x[k];
    }

    for (int i = 0; i < n - 1; ++i) {
        r2[i] = r3[i];
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

bool is_zero(int size, double* y, double a_norm, double EPS) {
    for (int i = 0; i < size; ++i) {
        if (std::fabs(y[i]) > EPS * a_norm) {
            return false;
        }
    } 

    return true;
}

bool build_x(int size, double* x, double* y, double a_norm, double EPS) {
    bool flag = true;
    double norm_y = vector_norm(size, y);
    y[0] -= norm_y;
    
    if (!is_zero(size, y, a_norm, EPS)) {
        double norm_denominator = vector_norm(size, y);
        for (int i = 0; i < size; ++i) {
            x[i] = y[i] / norm_denominator;
        }
    } else {
        flag = false;
    }

    y[0] += norm_y; 
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

double residual1(int n, double trace_a, double a_norm, double* lambda) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += lambda[i];
    }     

    return std::fabs(trace_a - sum) / a_norm;
}

double residual2(int n, double length_a, double a_norm, double* lambda) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += lambda[i] * lambda[i];
    }   

    sum = std::sqrt(sum);
    return std::fabs(length_a - sum) / a_norm;
}

double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}
