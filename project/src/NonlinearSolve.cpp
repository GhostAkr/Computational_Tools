//
// Created by ighos on 08.01.2019.
//

#include "../include/NonlinearSolve.h"

double f1(double x) {
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double f2(double x) {
    return sqrt(x + 1) - 1;
}

double f3(double x) {
    double x2 = x * x;
    double x3 = x2 * x;
    return 35 * x3 - 67 * x2 - 3 * x + 3;
}