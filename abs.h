#pragma once
#include <iostream>
#include <vector>
using namespace std;

vector<double> abs(const vector<double>& v) {
	vector<double> output(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < 0) {
			output[i] = -v[i];
		} else {
			output[i] = v[i];
		}
	}
	return output;
}