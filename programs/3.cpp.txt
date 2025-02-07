#include <Eigen/Dense>
#include <array>
#include <filesystem>
#include <format>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
// #include <opencv2/opencv.hpp>
#include <ranges>

namespace eg = Eigen;
namespace fs = std::filesystem;
using namespace std::string_literals;
using namespace std::string_view_literals;

using u64 = size_t;
using u32 = uint32_t;

struct MatrixCompressor {
	MatrixCompressor(
		const eg::MatrixXd& A,
		const u32 row_min,
		const u32 row_max,
		const u32 col_min,
		const u32 col_max,
		const u32 rank,
		const double epsilon
	) {
		eg::MatrixXd A_block = A.block(row_min, col_min, row_max - row_min, col_max - col_min);

		eg::MatrixXd D;

		eg::VectorXd singVal;

		if (A_block.size() != 0) {
			auto SVD = eg::BDCSVD(A_block, eg::ComputeThinU | eg::ComputeThinV);
			singVal = SVD.singularValues();
			D = SVD.singularValues().head(std::min(rank, (u32)SVD.singularValues().size())).asDiagonal();
			U = SVD.matrixU().leftCols(std::min(rank, (u32)SVD.matrixU().cols()));
			V = SVD.matrixV().leftCols(std::min(rank, (u32)SVD.matrixV().cols())).transpose();
		}

		if (row_max == row_min + rank or singVal.size() < rank or singVal[rank - 1] < epsilon) {
			sons = nullptr;
			singularvalues = D.diagonal();
			size = { row_min, row_max, col_min, col_max };
		} else {
			auto row_newmax = (row_min + row_max) / 2;
			auto col_newmax = (col_min + col_max) / 2;

			size = { row_min, row_max, col_min, col_max };

			sons = std::unique_ptr<MatrixCompressor[]>(new MatrixCompressor[4]{
				{ A,	 row_min, row_newmax,	  col_min, col_newmax, rank, epsilon },
				{ A,	 row_min, row_newmax, col_newmax,	  col_max, rank, epsilon },
				{ A, row_newmax,	 row_max,	  col_min, col_newmax, rank, epsilon },
				{ A, row_newmax,	 row_max, col_newmax,	  col_max, rank, epsilon },
			});
		}
	}

	u64 bytes() const noexcept {
		return sizeof(*this) + U.size() * sizeof(double) + V.size() * sizeof(double) +
			singularvalues.size() * sizeof(double) + [this]() {
				u64 sons_sum = 0;
				if (sons) {
					for (u32 i = 0; i != 4; ++i) {
						sons_sum += sons[i].bytes();
					}
				}
				return sons_sum;
			}();
	}

	eg::MatrixXd& decompress(eg::MatrixXd& dest) const noexcept {
		if (not sons) {
			eg::MatrixXd sigma = eg::MatrixXd::Zero(singularvalues.size(), singularvalues.size());
			sigma.diagonal() = singularvalues;

			if (size[0] != size[1] and size[2] != size[3]) {
				dest.block(size[0], size[2], size[1] - size[0], size[3] - size[2]) = U * sigma * V;
			}
		} else {
			for (u32 i = 0; i != 4; ++i) {
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
	u32 rank = 0;
	std::array<u32, 4> size;
	eg::MatrixXd U;
	eg::VectorXd singularvalues;
	eg::MatrixXd V;
};

void saveToBMP(const eg::MatrixXd& R, const eg::MatrixXd& G, const eg::MatrixXd& B, const std::string& filename) {
	u32 height = R.rows();
	u32 width = R.cols();

	// BMP file header
	uint8_t fileHeader[14] = { 0x42, 0x4D, 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };

	// BMP information header
	uint8_t infoHeader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, 0, 0, 0, 0,
							   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0 };

	// dimensions
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

	auto file = std::ofstream(filename, std::ios::binary);

	file.write(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));
	file.write(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	u32 rowPadding = (4 - (width * 3) % 4) % 4;
	auto padding = std::vector<uint8_t>(rowPadding, 0);

	for (u32 y = height - 1; y != (u32)-1; --y) {
		for (u32 x = 0; x < width; ++x) {
			file.put(static_cast<uint8_t>(B(y, x)));
			file.put(static_cast<uint8_t>(G(y, x)));
			file.put(static_cast<uint8_t>(R(y, x)));
		}
		file.write(reinterpret_cast<char*>(padding.data()), rowPadding);
	}

	std::cout << std::format("RGB {}x{} saved to '{}'\n", width, height, filename);
}

void saveMToBMP(const eg::MatrixXd& M, const std::string& mode, const std::string& filename) {
	u32 height = M.rows();
	u32 width = M.cols();

	// BMP file header
	uint8_t fileHeader[14] = { 0x42, 0x4D, 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };

	// BMP information header
	uint8_t infoHeader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, 0, 0, 0, 0,
							   0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0 };

	// dimensions
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

	auto file = std::ofstream(filename, std::ios::binary);

	file.write(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));
	file.write(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	u32 rowPadding = (4 - (width * 3) % 4) % 4;
	auto padding = std::vector<uint8_t>(rowPadding, 0);

	for (u32 y = height - 1; y != (u32)-1; --y) {
		for (u32 x = 0; x < width; ++x) {
			auto val = M(y, x);
			if (mode == "R") {
				file.put(0);
				file.put(0);
				file.put(val);
			} else if (mode == "G") {
				file.put(0);
				file.put(val);
				file.put(0);
			} else if (mode == "B") {
				file.put(val);
				file.put(0);
				file.put(0);
			} else if (mode == "grid") {
				if (val == 0) {
					file.put(0);
					file.put(0);
					file.put(0);
				} else {
					file.put(255);
					file.put(255);
					file.put(255);
				}
			}
		}
		file.write(reinterpret_cast<char*>(padding.data()), rowPadding);
	}

	std::cout << std::format("Grid {}x{} saved to '{}'\n", width, height, filename);
}

std::tuple<eg::MatrixXd, eg::MatrixXd, eg::MatrixXd> readFromBMP(const std::string& filename) {
	auto file = std::ifstream(filename, std::ios::binary);

	uint8_t fileHeader[14]{};
	file.read(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));

	uint8_t infoHeader[40]{};
	file.read(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

	u32 width = infoHeader[4] | (infoHeader[5] << 8) | (infoHeader[6] << 16) | (infoHeader[7] << 24);
	u32 height = infoHeader[8] | (infoHeader[9] << 8) | (infoHeader[10] << 16) | (infoHeader[11] << 24);

	auto R = eg::MatrixXd(height, width);
	auto G = eg::MatrixXd(height, width);
	auto B = eg::MatrixXd(height, width);

	int dataOffset = fileHeader[10] | (fileHeader[11] << 8) | (fileHeader[12] << 16) | (fileHeader[13] << 24);
	file.seekg(dataOffset, std::ios::beg);

	// Read pixel data (bottom-to-top row order, padded to multiples of 4 bytes)
	int rowPadding = (4 - (width * 3) % 4) % 4;
	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x < width; ++x) {
			B(y, x) = file.get();
			G(y, x) = file.get();
			R(y, x) = file.get();
		}
		file.ignore(rowPadding);
	}

	std::cout << std::format("Read {}x{} image from '{}'\n", width, height, filename);
	return std::make_tuple(std::move(R), std::move(G), std::move(B));
}

u64 bytes(const eg::MatrixXd& M) noexcept {
	return sizeof(M) + M.size() * sizeof(double);
}

u64 bytes(const MatrixCompressor& C) noexcept {
	return C.bytes();
}

int main() {
	std::string mode;
	std::cout << "mode (manual or auto) = ";
	std::cin >> mode;
	if (mode != "manual" and mode != "auto") {
		return 1;
	}

	std::string path;
	std::cout << "File: ";
	std::cin >> path;

	auto [R, G, B] = readFromBMP(path);
	u32 height = R.rows();
	u32 width = R.cols();

	fs::create_directory(fs::current_path() / "results");

	if (mode == "manual") {
		u32 rank;
		double epsilon;

		std::cout << "rank = ";
		std::cin >> rank;
		std::cout << "epsilon = ";
		std::cin >> epsilon;

		// paralellize compression
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

		std::cout << std::format(
			"RGB = {}\ncompressed_RGB = {}\n",
			3 * bytes(R),
			compressed_R.bytes() + compressed_G.bytes() + compressed_B.bytes()
		);
		std::cout << std::format(
			"errors R:\n\tmean: {}\n\tMSE: {}\n",
			(R - decompressed_R).cwiseAbs().mean(),
			(R - decompressed_R).cwiseAbs2().mean()
		);
		std::cout << std::format(
			"errors G:\n\tmean: {}\n\tMSE: {}\n",
			(G - decompressed_G).cwiseAbs().mean(),
			(G - decompressed_G).cwiseAbs2().mean()
		);
		std::cout << std::format(
			"errors B:\n\tmean: {}\n\tMSE: {}\n",
			(B - decompressed_B).cwiseAbs().mean(),
			(B - decompressed_B).cwiseAbs2().mean()
		);

		saveMToBMP(decompressed_R, "R", std::format("results/manual_R.bmp"));
		saveMToBMP(decompressed_G, "G", std::format("results/manual_G.bmp"));
		saveMToBMP(decompressed_B, "B", std::format("results/manual_B.bmp"));
		saveMToBMP(compressed_R.grid(), "grid", std::format("results/manual_R_grid.bmp"));
		saveMToBMP(compressed_G.grid(), "grid", std::format("results/manual_G_grid.bmp"));
		saveMToBMP(compressed_B.grid(), "grid", std::format("results/manual_B_grid.bmp"));
		saveToBMP(decompressed_R, decompressed_G, decompressed_B, std::format("results/manual_decompressed_RGB.bmp"));

		return 0;
	}

	auto compression_size_file = std::ofstream("results/sizes.txt");
	auto errors_file = std::ofstream("results/errors.txt");
	errors_file << "M\tmean\tMSE\n";
	compression_size_file << std::format("original_RGB\t{}\n", 3 * bytes(R));
	compression_size_file << std::format("original_R\t{}\n", bytes(R));
	compression_size_file << std::format("original_G\t{}\n", bytes(G));
	compression_size_file << std::format("original_B\t{}\n", bytes(B));

	saveToBMP(R, G, B, "results/RBG.bmp");

	eg::VectorXd svs_R = eg::BDCSVD(R, eg::ComputeThinU | eg::ComputeThinV).singularValues();
	eg::VectorXd svs_G = eg::BDCSVD(G, eg::ComputeThinU | eg::ComputeThinV).singularValues();
	eg::VectorXd svs_B = eg::BDCSVD(B, eg::ComputeThinU | eg::ComputeThinV).singularValues();

	{
		auto file = std::ofstream("results/svs_R.txt");
		for (auto&& sv : svs_R) {
			file << std::format("{}\n", sv);
		}
	}
	{
		auto file = std::ofstream("results/svs_G.txt");
		for (auto&& sv : svs_G) {
			file << std::format("{}\n", sv);
		}
	}
	{
		auto file = std::ofstream("results/svs_B.txt");
		for (auto&& sv : svs_B) {
			file << std::format("{}\n", sv);
		}
	}

	for (u32 rank : { 1, 4 }) {
		for (auto epsilon_of : { 'R', 'G', 'B' }) {
			for (auto epsilon_mode : { "first"sv, "middle"sv, "last"sv }) {
				std::cout << std::format(
					"work = (rank: {}, epsilon_of: {}, epsilon_mode: {})\n",
					rank,
					epsilon_of,
					epsilon_mode
				);

				double epsilon = [&]() {
					auto&& svs = epsilon_of == 'R' ? svs_R : (epsilon_of == 'G' ? svs_G : svs_B);

					return svs
						[epsilon_mode == "first" ? 0 : (epsilon_mode == "last" ? svs.size() - 1 : svs.size() / 2)];
				}();

				// paralellize compression
				auto compressed_R_future = std::async(std::launch::async, [&]() {
					return MatrixCompressor(R, 0, height, 0, width, rank, epsilon);
				});
				auto compressed_G_future = std::async(std::launch::async, [&]() {
					return MatrixCompressor(G, 0, height, 0, width, rank, epsilon);
				});
				auto compressed_B_future = std::async(std::launch::async, [&]() {
					return MatrixCompressor(B, 0, height, 0, width, rank, epsilon);
				});

				auto compressed_R = compressed_R_future.get();
				auto compressed_G = compressed_G_future.get();
				auto compressed_B = compressed_B_future.get();

				eg::MatrixXd decompressed_R = eg::MatrixXd::Zero(height, width);
				eg::MatrixXd decompressed_G = eg::MatrixXd::Zero(height, width);
				eg::MatrixXd decompressed_B = eg::MatrixXd::Zero(height, width);

				compressed_R.decompress(decompressed_R);
				compressed_G.decompress(decompressed_G);
				compressed_B.decompress(decompressed_B);

				saveMToBMP(
					decompressed_R,
					"R",
					std::format("results/R_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveMToBMP(
					decompressed_G,
					"G",
					std::format("results/G_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveMToBMP(
					decompressed_B,
					"B",
					std::format("results/B_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveMToBMP(
					compressed_R.grid(),
					"grid",
					std::format("results/R_grid_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveMToBMP(
					compressed_G.grid(),
					"grid",
					std::format("results/G_grid_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveMToBMP(
					compressed_B.grid(),
					"grid",
					std::format("results/B_grid_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				saveToBMP(
					decompressed_R,
					decompressed_G,
					decompressed_B,
					std::format("results/decompressed_RGB_r{}_eps_{}_{}.bmp", rank, epsilon_of, epsilon_mode)
				);
				compression_size_file << std::format(
					"RGB_r{}_eps_{}_{}\t{}\n",
					rank,
					epsilon_of,
					epsilon_mode,
					compressed_R.bytes() + compressed_G.bytes() + compressed_B.bytes()
				);
				errors_file << std::format(
					"R_r{}_eps_{}_{}\t{}\t{}\n",
					rank,
					epsilon_of,
					epsilon_mode,
					(R - decompressed_R).cwiseAbs().mean(),
					(R - decompressed_R).cwiseAbs2().mean()
				);
				errors_file << std::format(
					"G_r{}_eps_{}_{}\t{}\t{}\n",
					rank,
					epsilon_of,
					epsilon_mode,
					(G - decompressed_G).cwiseAbs().mean(),
					(G - decompressed_G).cwiseAbs2().mean()
				);
				errors_file << std::format(
					"B_r{}_eps_{}_{}\t{}\t{}\n",
					rank,
					epsilon_of,
					epsilon_mode,
					(B - decompressed_B).cwiseAbs().mean(),
					(B - decompressed_B).cwiseAbs2().mean()
				);
			}
		}
	}
}
