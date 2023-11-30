#include <fstream>
#include "inc.h"

double f(int s, int n, int i, int j) { 
    switch(s) {
        case 1:
            return n - std::max(i, j);
        case 2:
            if (i == j) {
                return 2.;
            } else if (std::abs(i - j) == 1) {
                return -1.;
            } else {
                return 0.;
            }
        case 3:
            if (i == j && j < n - 1) {
                return 1;
            } else if (j == n - 1) {
                return i;
            } else if (i == n - 1) {
                return j;
            } else {
                return 0;
            }
        case 4:
            return 1. / (i + j + 1);
        default:
            return 1;
    }
}

void input_matrix(int s, int n, double* matrix) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[n * i + j] = f(s, n, i, j);
        }
    }       
}

bool read_file(const std::string& name_file, int n, double* matrix) {
    std::ifstream fin(name_file);
    
    if (!fin.is_open()) {
        return false;
    }
    
    double x;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(fin >> x)) {
                return false;  
            }
            
            matrix[n * i + j] = x;
        }
    }
    return true;
}
