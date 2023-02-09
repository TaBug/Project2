#pragma once
#include <iostream>
#include <vector>

std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        std::cout << "Error: vectors have different sizes" << std::endl;
        return std::vector<double>();
    }

    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}
