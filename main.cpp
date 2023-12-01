#include "inc.h"
// 20:45. 30.11.23. Ебашим собственные значения. Приведение симм. м-цы к 3-хдиаг. виду. Let's go.

int main(int argc, char* argv[]) {
    // ./a.out n m eps k filename
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    [[maybe_unused]] double eps = std::stoi(argv[3]);
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

    std::cout << "\nA:\n";
    output(n, m, n, matrix);
    std::cout << "\n";

    double* x = new double[n];
    double* y = new double[n];
    double* c = new double[n];
    solution(n, matrix, x, y, c);

    std::cout << "\nA:\n";
    output(n, m, n, matrix);
    std::cout << "\n";

    /*
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d \
    Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
    argv[0], res1, res2, its, its / n, t1, t2);
    */

    delete[] matrix;
    delete[] x;
    delete[] y;
    delete[] c;
    return 0;
}
