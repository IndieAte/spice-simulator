#pragma once

#include <complex>
#include <Eigen/Dense>

namespace Eigen {
	typedef Matrix<std::complex<double>, Dynamic, Dynamic> ComplexMatrix;
	typedef Matrix<std::complex<double>, Dynamic, 1> ComplexVector;
}