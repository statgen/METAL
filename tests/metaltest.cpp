#include <iostream>
#include <cmath>
#include "MathStats.h"

using namespace std;

inline int fcmp(double x, double y, double epsilon) {
    int max_exponent = 0;
    double delta = 0.0;
    double diff = 0.0;

    frexp(fabs(x) > fabs(y) ? x : y, &max_exponent);
    delta = ldexp(epsilon, max_exponent);

    diff = x - y;

    if (diff > delta) {
        return 1;
    } else if (diff < -delta) {
        return -1;
    } else {
        return 0;
    }
}

int test_normp() {
    double z[] = {-5, -1, 0, 2, 5};
    double p[] = {0.00000028665157187919, 0.15865525393145704647, 0.5, 0.97724986805182079141, 0.99999971334842807646};

    for (int i = 0; i < 5; ++i) {
        if (fcmp(p[i], normp(z[i]), 0.000000001) != 0) {
            return 1;
        }
    }
    return 0;
}

int test_binormp() {
    double x[] = {0.0, 1.0, 0.0, 0.0, 10.0, 0.0, 1.0};
    double y[] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    double rho[] = {0.5, 0.5, 0.5, 0.0, 0.0, -0.5, 1.0};
    double p[] = {0.18377630, 0.0943538977, 0.0943538977, 0.1591549431, 3.06970072291198343e-23, 0.18377630, 0.0};

    for (int i = 0; i < 6; ++i) {
        if (fcmp(p[i], binormp(x[i], y[i], rho[i]), 0.00000001) != 0) {
            return 1;
        }
    }

    if (numeric_limits<double>::infinity() != binormp(0.0, 0.0, 1.0)) {
        return 1;
    }

    if (numeric_limits<double>::infinity() != binormp(0.0, 0.0, -1.0)) {
        return 1;
    }

    return 0;
}

int test_binomq() {
    double lx[] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    double ly[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    double rho[] = {0.0, 0.5, 0.5, 0.5, -0.5, 0.98};
    double q[] = {0.25, 0.3333333333, 0.1273982066, 0.1273982066, 0.1666666667, 0.4681157196};

    for (int i = 0; i < 6; ++i) {
        if (fcmp(q[i], binormq(lx[i], ly[i], rho[i]), 0.00000001) != 0) {
            return 1;
        }
    }
    return 0;
}

int main(int argc, char ** argv) {

    if (argc != 2) {
        return 1;
    }

    if (strcmp(argv[1], "test_normp") == 0) {
        return test_normp();
    } else if (strcmp(argv[1], "test_binormp") == 0) {
        return test_binormp();
    } else if (strcmp(argv[1], "test_binormq") == 0) {
        return test_binomq();
    }

    return 1;
}