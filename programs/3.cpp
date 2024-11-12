#include "RandomizedSVD.hpp"
#include <Eigen/Dense>
#include <array>
#include <format>
#include <iostream>
#include <memory>

namespace eg = Eigen;

using u64 = size_t;

class MatrixCompressor {
public:

	MatrixCompressor(
		const eg::MatrixXd& A,
		const u64 t_min,
		const u64 t_max,
		const u64 s_min,
		const u64 s_max,
		const u64 r,
		const double epsilon
	) noexcept(!MATRIX_DEBUG) {
		eg::MatrixXd A_block = A.block(t_min, s_min, t_max - t_min, s_max - s_min);
		/*if (A_block.rows() == 0 || A_block.cols() == 0) {
			sons = nullptr;
			singularvalues = eg::VectorXd::Zero(0);
			return;
		}*/
		// RandomizedSvd SVD;

		eg::MatrixXd D;

		eg::VectorXd singVal;

		if (A_block.size() != 0) {
			eg::JacobiSVD SVD = eg::JacobiSVD<eg::MatrixXd>(A_block, eg::ComputeThinU | eg::ComputeThinV);
			singVal = SVD.singularValues();
			D = SVD.singularValues().head(r).asDiagonal();
			U = SVD.matrixU().leftCols(r);
			V = SVD.matrixV().leftCols(r).transpose();
		}

		std::cout << "U:\n" << U << '\n';
		std::cout << "D:\n" << D << '\n';
		std::cout << "V:\n" << V << '\n';
		std::cout << "row_max = " << t_max << '\n';
		std::cout << "row_min = " << t_min << '\n';
		std::cout << '\n';

		std::cout << U * D * V << '\n';

		if (t_max == t_min + r or singVal.size() < r or singVal[r - 1] < epsilon) {
			sons = nullptr;
			singularvalues = D.diagonal();
			size = { t_min, t_max, s_min, s_max };

			/*if (A_block.isZero()) {
				rank = 0;
				size = { t_min, t_max, s_min, s_max };
			} else {
				auto sigma = D.diagonal();
				auto rank = r;

				this->rank = rank;
				singularvalues = sigma.head(rank);
				U = U.leftCols(rank);
				V = D.block(0, 0, rank, rank) * V.topRows(rank);
				sons = nullptr;
				size = { t_min, t_max, s_min, s_max };
			}*/
		} else {
			auto t_newmax = (t_min + t_max) / 2;
			auto s_newmax = (s_min + s_max) / 2;

			size = { t_min, t_max, s_min, s_max };

			sons = std::unique_ptr<MatrixCompressor[]>(new MatrixCompressor[4]{
				{ A,	 t_min, t_newmax,	  s_min, s_newmax, r, epsilon },
				{ A,	 t_min, t_newmax, s_newmax,	s_max, r, epsilon },
				{ A, t_newmax,	   t_max,	  s_min, s_newmax, r, epsilon },
				{ A, t_newmax,	   t_max, s_newmax,	s_max, r, epsilon },
			});
		}
	}

	eg::MatrixXd& decompress(eg::MatrixXd& dest) const noexcept {
		if (not sons) {
			eg::MatrixXd sigma = eg::MatrixXd::Zero(singularvalues.size(), singularvalues.size());
			sigma.diagonal() = singularvalues;

			auto Usigma = U * sigma;

			std::cout << std::format("{}, {} | {}, {}\n", Usigma.rows(), Usigma.cols(), V.rows(), V.cols());

			if (size[0] != size[1] and size[2] != size[3]) {
				dest.block(size[0], size[2], size[1] - size[0], size[3] - size[2]) = Usigma * V;
			}
		} else {
			for (u64 i = 0; i != 4; ++i) {
				sons[i].decompress(dest);
			}
		}
		return dest;
	}

	std::unique_ptr<MatrixCompressor[]> sons = nullptr;
	u64 rank = 0;
	std::array<u64, 4> size;
	eg::MatrixXd U;
	eg::VectorXd singularvalues;
	eg::MatrixXd V;
};

Eigen::MatrixXd generateRandomMatrix(int N, double min, double max) {
	// Generate NxN matrix with random values in the range [-1, 1]
	Eigen::MatrixXd random_matrix = Eigen::MatrixXd::Random(N, N);

	// Scale and shift values to the specified range [min, max]
	random_matrix = 0.5 * (random_matrix + Eigen::MatrixXd::Ones(N, N)); // Scale to [0, 1]
	random_matrix = (max - min) * random_matrix + Eigen::MatrixXd::Constant(N, N, min); // Scale to [min, max]

	return random_matrix;
}

int main() {
	u64 N = 8;

	eg::MatrixXd A = /*eg::MatrixXd{
		{ -0.997497,	 0.617481, -0.299417 },
		{  0.127171,	0.170019,  0.791925 },
		{ -0.613392, -0.0402539,	 0.64568 }
	};*/
		eg::MatrixXd::Random(N, N);
	// generateRandomMatrix(N, -1, 1);

	// Eigen::MatrixXd A(8, 8);

	//// Initialize the matrix with the given values
	// A << -0.997497, 0.64568, -0.817194, -0.982177, -0.0984222, 0.751946, 0.724479, -0.467574, 0.127171, 0.49321,
	//	-0.271096, -0.24424, -0.295755, 0.453352, -0.580798, -0.405438, -0.613392, -0.651784, -0.705374, 0.0633259,
	//	-0.885922, 0.911802, 0.559313, 0.680288, 0.617481, 0.717887, -0.668203, 0.142369, 0.215369, 0.851436, 0.687307,
	//	-0.952513, 0.170019, 0.421003, 0.97705, 0.203528, 0.566637, 0.0787072, 0.993591, -0.248268, -0.0402539,
	//	0.0270699, -0.108615, 0.214331, 0.605213, -0.715323, 0.99939, -0.814753, -0.299417, -0.39201, -0.761834,
	//	-0.667531, 0.0397656, -0.0758385, 0.222999, 0.354411, 0.791925, -0.970031, -0.990661, 0.32609, -0.3961,
	//	-0.529344, -0.215125, -0.88757;

	std::cout << A << '\n';

	auto compressed_A = MatrixCompressor(A, 0, N, 0, N, 1, 0.1);

	eg::MatrixXd decompressed_A = eg::MatrixXd::Zero(N, N);
	compressed_A.decompress(decompressed_A);

	std::cout << (A - decompressed_A) << '\n';
}
