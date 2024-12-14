#include "RandomizedSVD.hpp"
#include <Eigen/Dense>
#include <array>
#include <format>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <opencv2/opencv.hpp>
#include <ranges>

namespace eg = Eigen;

using u64 = size_t;

struct MatrixCompressor {
	MatrixCompressor(
		const eg::MatrixXd& A,
		const u64 t_min,
		const u64 t_max,
		const u64 s_min,
		const u64 s_max,
		const u64 r,
		const double epsilon
	) {
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
			auto SVD = eg::BDCSVD<eg::MatrixXd>(A_block, eg::ComputeThinU | eg::ComputeThinV);
			singVal = SVD.singularValues();
			D = SVD.singularValues().head(r).asDiagonal();
			U = SVD.matrixU().leftCols(r);
			V = SVD.matrixV().leftCols(r).transpose();
		}

		/*std::cout << "U:\n" << U << '\n';
		std::cout << "D:\n" << D << '\n';
		std::cout << "V:\n" << V << '\n';
		std::cout << "row_max = " << t_max << '\n';
		std::cout << "row_min = " << t_min << '\n';
		std::cout << '\n';

		std::cout << U * D * V << '\n';*/

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

			// auto Usigma = U * sigma;

			// std::cout << std::format("{}, {} | {}, {}\n", Usigma.rows(), Usigma.cols(), V.rows(), V.cols());

			if (size[0] != size[1] and size[2] != size[3]) {
				dest.block(size[0], size[2], size[1] - size[0], size[3] - size[2]) = (U * sigma * V);
			}
		} else {
			for (u64 i = 0; i != 4; ++i) {
				sons[i].decompress(dest);
			}
		}
		return dest;
	}

	eg::MatrixXd grid() const noexcept {
		eg::MatrixXd result = eg::MatrixXd::Ones(size[1], size[3]);
		grid(result);
		return result;
	}

	void grid(eg::MatrixXd& result) const noexcept {
		for (int x = size[2]; x != size[3]; ++x) {
			result(size[0], x) = 0;
		}
		for (int y = size[0]; y != size[1]; ++y) {
			result(y, size[2]) = 0;
		}
		if (sons) {
			for (int i = 0; i != 4; ++i) {
				sons[i].grid(result);
			}
		}
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

struct RGB {
	uint8_t red;
	uint8_t green;
	uint8_t blue;
};

// Function to read a BMP file and return image dimensions and pixel data
void saveToBMP(const MatrixXd& R, const MatrixXd& G, const MatrixXd& B, const std::string& filename) {
	int height = R.rows();
	int width = R.cols();

	// Ensure the input matrices are of the correct size
	if (G.rows() != height || G.cols() != width || B.rows() != height || B.cols() != width) {
		throw std::invalid_argument("Input matrices must have the same dimensions.");
	}

	// BMP file header (14 bytes)
	uint8_t fileHeader[14] = {
		0x42, 0x4D, // Signature "BM"
		0,	  0,	0, 0, // File size in bytes (will be filled later)
		0,	  0, // Reserved
		0,	  0, // Reserved
		54,	  0,	0, 0 // Pixel data offset (54 bytes header)
	};

	// BMP information header (40 bytes)
	uint8_t infoHeader[40] = {
		40, 0, 0, 0, // Header size (40 bytes)
		0,	0, 0, 0, // Image width (will be filled later)
		0,	0, 0, 0, // Image height (will be filled later)
		1,	0, // Planes (1)
		24, 0, // Bits per pixel (24 bpp)
		0,	0, 0, 0, // Compression (0 - none)
		0,	0, 0, 0, // Image size (can be 0 for no compression)
		0,	0, 0, 0, // X pixels per meter (not specified)
		0,	0, 0, 0, // Y pixels per meter (not specified)
		0,	0, 0, 0, // Total colors (0 - default)
		0,	0, 0, 0 // Important colors (0 - all)
	};

	// Fill in image dimensions in the info header
	int fileSize = 54 + 3 * width * height;
	fileHeader[2] = fileSize;
	fileHeader[3] = fileSize >> 8;
	fileHeader[4] = fileSize >> 16;
	fileHeader[5] = fileSize >> 24;

	infoHeader[4] = width;
	infoHeader[5] = width >> 8;
	infoHeader[6] = width >> 16;
	infoHeader[7] = width >> 24;

	infoHeader[8] = height;
	infoHeader[9] = height >> 8;
	infoHeader[10] = height >> 16;
	infoHeader[11] = height >> 24;

	// Create the file
	std::ofstream outFile(filename, std::ios::binary);
	if (!outFile) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	// Write headers
	outFile.write(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));
	outFile.write(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	// BMP pixel data (bottom-to-top row order, padded to multiples of 4 bytes)
	int rowPadding = (4 - (width * 3) % 4) % 4; // Each row must be a multiple of 4 bytes
	std::vector<uint8_t> padding(rowPadding, 0);

	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x < width; ++x) {
			outFile.put(static_cast<uint8_t>(B(y, x))); // Blue
			outFile.put(static_cast<uint8_t>(G(y, x))); // Green
			outFile.put(static_cast<uint8_t>(R(y, x))); // Red
		}
		outFile.write(reinterpret_cast<char*>(padding.data()), rowPadding);
	}

	outFile.close();
	std::cout << "Image saved to " << filename << "\n";
}

void saveGridToBMP(const MatrixXd& M, const std::string& filename) {
	int height = M.rows();
	int width = M.cols();

	// BMP file header (14 bytes)
	uint8_t fileHeader[14] = {
		0x42, 0x4D, // Signature "BM"
		0,	  0,	0, 0, // File size in bytes (will be filled later)
		0,	  0, // Reserved
		0,	  0, // Reserved
		54,	  0,	0, 0 // Pixel data offset (54 bytes header)
	};

	// BMP information header (40 bytes)
	uint8_t infoHeader[40] = {
		40, 0, 0, 0, // Header size (40 bytes)
		0,	0, 0, 0, // Image width (will be filled later)
		0,	0, 0, 0, // Image height (will be filled later)
		1,	0, // Planes (1)
		24, 0, // Bits per pixel (24 bpp)
		0,	0, 0, 0, // Compression (0 - none)
		0,	0, 0, 0, // Image size (can be 0 for no compression)
		0,	0, 0, 0, // X pixels per meter (not specified)
		0,	0, 0, 0, // Y pixels per meter (not specified)
		0,	0, 0, 0, // Total colors (0 - default)
		0,	0, 0, 0 // Important colors (0 - all)
	};

	// Fill in image dimensions in the info header
	int fileSize = 54 + 3 * width * height;
	fileHeader[2] = fileSize;
	fileHeader[3] = fileSize >> 8;
	fileHeader[4] = fileSize >> 16;
	fileHeader[5] = fileSize >> 24;

	infoHeader[4] = width;
	infoHeader[5] = width >> 8;
	infoHeader[6] = width >> 16;
	infoHeader[7] = width >> 24;

	infoHeader[8] = height;
	infoHeader[9] = height >> 8;
	infoHeader[10] = height >> 16;
	infoHeader[11] = height >> 24;

	// Create the file
	std::ofstream outFile(filename, std::ios::binary);
	if (!outFile) {
		throw std::ios_base::failure("Failed to open file for writing.");
	}

	// Write headers
	outFile.write(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));
	outFile.write(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	// BMP pixel data (bottom-to-top row order, padded to multiples of 4 bytes)
	int rowPadding = (4 - (width * 3) % 4) % 4; // Each row must be a multiple of 4 bytes
	std::vector<uint8_t> padding(rowPadding, 0);

	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x < width; ++x) {
			if (M(y, x) == 0) {
				outFile.put(0); // Blue
				outFile.put(0); // Green
				outFile.put(255); // Red
			} else {
				outFile.put(255); // Blue
				outFile.put(255); // Green
				outFile.put(255); // Red
			}
			// outFile.put(static_cast<uint8_t>(M(y, x))); // Blue
			// outFile.put(static_cast<uint8_t>(M(y, x))); // Green
			// outFile.put(static_cast<uint8_t>(M(y, x))); // Red
		}
		outFile.write(reinterpret_cast<char*>(padding.data()), rowPadding);
	}

	outFile.close();
	std::cout << "Image saved to " << filename << "\n";
}

// Function to read a BMP file and return Eigen Matrices as a tuple
std::tuple<MatrixXd, MatrixXd, MatrixXd> readFromBMP(const std::string& filename) {
	std::ifstream inFile(filename, std::ios::binary);
	if (!inFile) {
		throw std::ios_base::failure("Failed to open file for reading.");
	}

	// Read BMP file header
	uint8_t fileHeader[14];
	inFile.read(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));

	// Read BMP info header
	uint8_t infoHeader[40];
	inFile.read(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	// Extract width and height
	int width = infoHeader[4] | (infoHeader[5] << 8) | (infoHeader[6] << 16) | (infoHeader[7] << 24);
	int height = infoHeader[8] | (infoHeader[9] << 8) | (infoHeader[10] << 16) | (infoHeader[11] << 24);

	// Ensure 24 bpp (3 bytes per pixel)
	if (infoHeader[14] != 24) {
		throw std::runtime_error("Only 24 bpp BMP files are supported.");
	}

	// Initialize matrices
	MatrixXd R(height, width);
	MatrixXd G(height, width);
	MatrixXd B(height, width);

	// Pixel data offset
	int dataOffset = fileHeader[10] | (fileHeader[11] << 8) | (fileHeader[12] << 16) | (fileHeader[13] << 24);
	inFile.seekg(dataOffset, std::ios::beg);

	// Read pixel data (bottom-to-top row order, padded to multiples of 4 bytes)
	int rowPadding = (4 - (width * 3) % 4) % 4;
	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x < width; ++x) {
			B(y, x) = inFile.get();
			G(y, x) = inFile.get();
			R(y, x) = inFile.get();
		}
		inFile.ignore(rowPadding);
	}

	inFile.close();
	std::cout << "Image read from " << filename << "\n";
	return std::make_tuple(R, G, B);
}

// std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> loadImageAsMatrices(const std::string& imagePath) {
//	cv::Mat img = cv::imread(imagePath);
//	if (img.empty()) {
//		throw std::runtime_error("Image not found!");
//	}
//
//	cv::Mat imgFloat;
//	img.convertTo(imgFloat, CV_64F);
//
//	int rows = img.rows, cols = img.cols;
//	Eigen::MatrixXd R(rows, cols), G(rows, cols), B(rows, cols);
//
//	for (int i = 0; i < rows; ++i) {
//		for (int j = 0; j < cols; ++j) {
//			cv::Vec3b intensity = img.at<cv::Vec3b>(i, j);
//			R(i, j) = intensity[2]; // OpenCV uses BGR
//			G(i, j) = intensity[1];
//			B(i, j) = intensity[0];
//		}
//	}
//
//	return { R, G, B };
// }

// void drawMatrix(const Eigen::MatrixXd& M, const std::string& outputPath) {
//	Eigen::MatrixXd normalized = M;
//	normalized = (normalized - normalized.minCoeff()) / (normalized.maxCoeff() - normalized.minCoeff()) * 255.0;
//
//	cv::Mat img(M.rows(), M.cols(), CV_8UC1);
//	for (int i = 0; i < M.rows(); ++i) {
//		for (int j = 0; j < M.cols(); ++j) {
//			img.at<uchar>(i, j) = static_cast<uchar>(normalized(i, j));
//		}
//	}
//
//	cv::imwrite(outputPath, img);
// }

// void drawBitmap(
//	const Eigen::MatrixXd& R,
//	const Eigen::MatrixXd& G,
//	const Eigen::MatrixXd& B,
//	const std::string& outputPath
//) {
//	int rows = R.rows(), cols = R.cols();
//	cv::Mat img(rows, cols, CV_8UC3);
//
//	for (int i = 0; i < rows; ++i) {
//		for (int j = 0; j < cols; ++j) {
//			img.at<cv::Vec3b>(i, j)[2] = static_cast<uchar>(R(i, j));
//			img.at<cv::Vec3b>(i, j)[1] = static_cast<uchar>(G(i, j));
//			img.at<cv::Vec3b>(i, j)[0] = static_cast<uchar>(B(i, j));
//		}
//	}
//
//	cv::imwrite(outputPath, img);
// }

int main() {
	try {
		const u64 N = 3;

		eg::MatrixXd A = /*eg::MatrixXd{
			{ -0.997497,	 0.617481, -0.299417 },
			{  0.127171,	0.170019,  0.791925 },
			{ -0.613392, -0.0402539,	 0.64568 }
		};*/
			eg::MatrixXd::Random(N, N);

		A *= 255;

		eg::MatrixXd decompressed_A = eg::MatrixXd::Zero(N, N);

		auto compressed_a = MatrixCompressor(A, 0, N, 0, N, 1, 0.1);
		compressed_a.decompress(decompressed_A);

		std::cout << A << "\n\n";
		std::cout << decompressed_A << '\n';

		// generateRandomMatrix(N, -1, 1);

		// Eigen::MatrixXd A(8, 8);

		//// Initialize the matrix with the given values
		// A << -0.997497, 0.64568, -0.817194, -0.982177, -0.0984222, 0.751946, 0.724479, -0.467574,
		// 0.127171, 0.49321,
		//	-0.271096, -0.24424, -0.295755, 0.453352, -0.580798, -0.405438, -0.613392, -0.651784, -0.705374,
		// 0.0633259, 	-0.885922, 0.911802, 0.559313, 0.680288, 0.617481, 0.717887, -0.668203, 0.142369,
		// 0.215369, 0.851436,
		// 0.687307, 	-0.952513, 0.170019, 0.421003, 0.97705, 0.203528, 0.566637, 0.0787072, 0.993591,
		// -0.248268,
		//-0.0402539, 	0.0270699, -0.108615, 0.214331, 0.605213, -0.715323, 0.99939, -0.814753, -0.299417,
		//-0.39201, -0.761834, 	-0.667531, 0.0397656, -0.0758385, 0.222999, 0.354411, 0.791925,
		//-0.970031, -0.990661, 0.32609, -0.3961, 	-0.529344, -0.215125, -0.88757;

		// std::cout << A << '\n';

		int height, width;
		std::string path;
		int rank;
		double epsilon;
		std::cin >> path >> rank >> epsilon;
		auto [R, G, B] = readFromBMP(path);
		height = R.rows();
		width = R.cols();

		// auto saveTestFile = std::ofstream("test.bmp");
		saveToBMP(R, G, B, "test.bmp");

		/*auto Rv_ = std::vector<double>();
		std::ranges::transform(Rv, std::back_inserter(Rv_), [](const uint8_t v) { return (double)v; });
		auto Gv_ = std::vector<double>();
		std::ranges::transform(Gv, std::back_inserter(Gv_), [](const uint8_t v) { return (double)v; });
		auto Bv_ = std::vector<double>();
		std::ranges::transform(Bv, std::back_inserter(Bv_), [](const uint8_t v) { return (double)v; });

		eg::MatrixXd R = eg::Map<Eigen::MatrixXd>(Rv_.data(), height, width);
		eg::MatrixXd G = eg::Map<Eigen::MatrixXd>(Gv_.data(), height, width);
		eg::MatrixXd B = eg::Map<Eigen::MatrixXd>(Bv_.data(), height, width);*/

		auto compressed_R_future =
			std::async(std::launch::async, [&]() { return MatrixCompressor(R, 0, height, 0, width, rank, epsilon); });
		auto compressed_G_future =
			std::async(std::launch::async, [&]() { return MatrixCompressor(G, 0, height, 0, width, rank, epsilon); });
		auto compressed_B_future =
			std::async(std::launch::async, [&]() { return MatrixCompressor(B, 0, height, 0, width, rank, epsilon); });

		auto compressed_R = compressed_R_future.get();
		auto compressed_G = compressed_G_future.get();
		auto compressed_B = compressed_B_future.get();

		eg::MatrixXd decompressed_R = eg::MatrixXd::Zero(height, width);
		eg::MatrixXd decompressed_G = eg::MatrixXd::Zero(height, width);
		eg::MatrixXd decompressed_B = eg::MatrixXd::Zero(height, width);

		compressed_R.decompress(decompressed_R);
		compressed_G.decompress(decompressed_G);
		compressed_B.decompress(decompressed_B);

		// std::vector<uint8_t> decRV;
		// std::vector<uint8_t> decGV;
		// std::vector<uint8_t> decBV;

		// for (int i = 0; i < decompressed_R.rows(); ++i) {
		//	for (int j = 0; j < decompressed_R.cols(); ++j) {
		//		decRV.push_back(decompressed_R(i, j)); // Add element at (i, j) to the vector
		//		decGV.push_back(decompressed_G(i, j)); // Add element at (i, j) to the vector
		//		decBV.push_back(decompressed_B(i, j)); // Add element at (i, j) to the vector
		//	}
		//	std::cout << '\n';
		// }

		// std::cout << (decRV == Rv) << '\n';
		// std::cout << (decGV == Gv) << '\n';
		// std::cout << (decBV == Bv) << '\n';
		/*std::cout << decRV.size() << '\n';
		std::cout << decGV.size() << '\n';
		std::cout << decBV.size() << '\n';*/

		// auto ofile = std::ofstream("xd2.bmp");

		auto grid_R = compressed_R.grid();
		auto grid_G = compressed_G.grid();
		auto grid_B = compressed_B.grid();

		saveGridToBMP(grid_R, "grid_R.bmp");
		saveGridToBMP(grid_G, "grid_G.bmp");
		saveGridToBMP(grid_B, "grid_B.bmp");

		saveToBMP(decompressed_R, decompressed_G, decompressed_B, "xd2.bmp");
	} catch (std::runtime_error& e) {
		std::cout << e.what() << '\n';
	}
}
