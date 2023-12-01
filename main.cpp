#include "inc.h"

int main(int argc, char* argv[]) {
    // ./a.out n m eps k filename
    int n = std::stoi(argv[1]);
    [[maybe_unused]] int m = std::stoi(argv[2]);
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

    double* x = new double[n];
    double* y = new double[n];
    double* r1 = new double[n];
    double* r2 = new double[n];
    double* r3 = new double[n];

    solution(n, matrix, x, y, r1, r2, r3);

    std::cout << "\nR:\n";
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
    delete[] r1;
    delete[] r2;
    delete[] r3;
    return 0;
}
