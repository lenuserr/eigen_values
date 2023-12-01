#include <iostream>
#include <string>

double f(int s, int n, int i, int j);
void input_matrix(int s, int n, double* matrix);
bool read_file(const std::string& name_file, int n, double* matrix);
void output(int n, int r, int l, double* vec);
void print_y(int size, double* y);
void solution(int n, double* matrix, double* x, double* y, double* r1, double* r2, double* r3);
void get_column(int n, int j, int k, double* matrix, double* y);
void put_column(int n, int j, int k, double* matrix, double* c);
void get_row(int n, int j, int k, double* matrix, double* y);
void put_row(int n, int j, int k, double* matrix, double* c);
bool build_x(int size, double* x, double* y);
bool is_zero(int size, double* y);
double trace(int n, double* matrix);
double matrix_norm(int n, double* matrix);
double matrix_length(int n, double* matrix);
double vector_norm(int size, double* vector);
double scalar_product(int size, double* a, double* b);
void product(int size, double* x, double* y, double* c);
