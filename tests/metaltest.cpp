#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <map>
#include <vector>
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

int read_column_values(const char* file_name, const char* column_name, map<string, double>& entries) {
    ifstream ifile_stream;
    string line("");
    int columns[2] = { numeric_limits<int>::min(), numeric_limits<int>::min() };
    int bounds[2][2] = { {0, 0}, {0, 0} };

    ifile_stream.open(file_name, ios::binary);
    if (ifile_stream.fail()) {
        return 1;
    }

    getline(ifile_stream, line);
    for (string::size_type s = 0, e = 0, c = 0; s < line.length(); ++e) {
        if ((e == line.length()) || (line[e] == '\t') || (line[e] == '\0')) {
            if (line.compare(s, e - s, "MarkerName") == 0) {
                columns[0] = c;
            } else if (line.compare(s, e - s, column_name) == 0) {
                columns[1] = c;
            }
            s = e + 1;
            ++c;
        }
    }

    if ((columns[0] < 0) || (columns[1] < 0)) {
        return 1;
    }

    if (ifile_stream.eof()) {
        return 1;
    }

    while (!getline(ifile_stream, line).eof()) {
        for (string::size_type s = 0, e = 0, c = 0; s < line.length(); ++e) {
            if ((e == line.length()) || (line[e] == '\t') || (line[e] == '\0')) {
                if (c == columns[0]) {
                    bounds[0][0] = s;
                    bounds[0][1] = e - s;
                } else if (c == columns[1]) {
                    bounds[1][0] = s;
                    bounds[1][1] = e - s;
                }
                s = e + 1;
                ++c;
            }
        }
        entries.insert(pair<string, double>(line.substr(bounds[0][0], bounds[0][1]), atof(line.substr(bounds[1][0], bounds[1][1]).c_str())));
    }

    if (!ifile_stream.eof() && ifile_stream.fail()) {
        return 1;
    }


    ifile_stream.close();

    return 0;
}


int test_file_correlation(const char* file_name1, const char* file_name2, const char* column_name, double corr_threshold, double tolerance) {
    map<string, double> entries1;
    map<string, double> entries2;
    vector<double> values1;
    vector<double> values2;
    int n = 0;
    double mean1 = 0.0, mean2 = 0.0;
    double std1 = 0.0, std2 = 0.0, cov = 0.0, corr = 0.0;
    if (read_column_values(file_name1, column_name, entries1) != 0) {
        return 1;
    }
    if (read_column_values(file_name2, column_name, entries2) != 0) {
        return 1;
    }
    if ((n = entries1.size()) != entries1.size()) {
        return 1;
    }
    for (map<string, double>::iterator it = entries1.begin(); it != entries1.end(); ++it) {
        if (entries2.find(it->first) == entries2.end()) {
            return 1;
        }
    }
    for (map<string, double>::iterator it = entries1.begin(); it != entries1.end(); ++it) {
        mean1 += it->second;
        values1.push_back(it->second);
    }
    mean1 /= values1.size();
    entries1.clear();
    for (map<string, double>::iterator it = entries2.begin(); it != entries2.end(); ++it) {
        mean2 += it->second;
        values2.push_back(it->second);
    }
    mean2 /= values2.size();
    entries2.clear();
    for (int i = 0; i < n; ++i) {
        if (abs(values1[i] - values2[i]) > tolerance) {
            return 1;
        }
        cov += (values1[i] - mean1) * (values2[i] - mean2);
        std1 += (values1[i] - mean1) * (values1[i] - mean1);
        std2 += (values2[i] - mean2) * (values2[i] - mean2);
    }
    corr = cov / (sqrt(std1) * sqrt(std2));
    if (signbit(corr) != signbit(corr_threshold)) {
        return 1;
    }
    if (abs(corr) < abs(corr_threshold)) {
        return 1;
    }
    return 0;
}


int main(int argc, char ** argv) {
    if (argc < 2) {
        return 1;
    }
    if (strcmp(argv[1], "test_normp") == 0) {
        return test_normp();
    } else if (strcmp(argv[1], "test_binormp") == 0) {
        return test_binormp();
    } else if (strcmp(argv[1], "test_binormq") == 0) {
        return test_binomq();
    } else if (strcmp(argv[1], "test_file_correlation") == 0) {
        if (argc != 7) {
            return 1;
        }
        return test_file_correlation(argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]));
    }
    return 1;
}
