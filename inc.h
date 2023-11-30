#include <iostream>
#include <string>

double f(int s, int n, int i, int j);
void input_matrix(int s, int n, double* matrix);
bool read_file(const std::string& name_file, int n, double* matrix);
void output(int n, int r, int l, double* vec);
void solution(int n, double* matrix, double* x, double* y, double* c);
void get_vector(int n, int j, double* matrix, double* y);
void build_x(int size, double* x, double* y);
double vector_norm(int size, double* vector);
double scalar_product(int size, double* a, double* b);
void product(int size, double* x, double* y, double* c);
