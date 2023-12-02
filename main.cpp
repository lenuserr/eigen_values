#include "inc.h"

int main(int argc, char* argv[]) {
    // ./a.out n m eps k filename
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    double eps = std::stod(argv[3]);
    int k = std::stoi(argv[4]);

    double* matrix = new double[n*n];

    if (k) {
        input_matrix(k, n, matrix);
    } else {
        std::string name_file = argv[argc - 1];
        if (!read_file(name_file, n, matrix)) {
            std::cout << "Проблемы с чтением файла" << "\n";

            delete[] matrix;
            return -1;
        }
    }

    std::cout << "A:\n";
    output(n, m, n, matrix);
    std::cout << "\n";

    double a_norm = matrix_norm(n, matrix);
    double trace_a = trace(n, matrix);
    double length_a = matrix_length(n, matrix);

    if (std::fabs(matrix[n] - matrix[1]) > eps * a_norm) {
        std::cout << "Метод работает только для симметричных матриц" << "\n";

        printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d \
        Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
        argv[0], -1., -1., 0, 0, 0., 0.);
        delete[] matrix;
        return -1;
    } 

    double* x = new double[n];
    double* y = new double[n];
    double* r1 = new double[n];
    double* r2 = new double[n];
    double* r3 = new double[n];
    double* lambda = new double[n];

    double t1;
    double t2;
    int its = solution(n, a_norm, matrix, x, y, r1, r2, r3, lambda, eps, &t1, &t2);

    std::cout << "\nlambda:\n";
    output(n, m, 1, lambda);
    std::cout << "\n";

    double res1 = residual1(n, trace_a, a_norm, lambda);
    double res2 = residual2(n, length_a, a_norm, lambda);
    
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d \
    Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
    argv[0], res1, res2, its, its / n, t1, t2);

    delete[] matrix;
    delete[] x;
    delete[] y;
    delete[] r1;
    delete[] r2;
    delete[] r3;
    delete[] lambda;
    return 0;
}
