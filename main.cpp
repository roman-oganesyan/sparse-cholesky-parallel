#include "CSR_Matrix.h"
#include <iostream>
#include <string>

int main() {
	std::string filename = "filename.mtx";
	Matrix test = load_mtx(filename);
	return 0;
}