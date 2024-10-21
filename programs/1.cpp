#pragma region operation_counting

#include <array>
#include <string_view>

namespace operation_counting {

enum class operation {
	add,
	sub,
	mul,
	div,
	unknown,
};

size_t counter[(size_t)operation::unknown]{};

void reset() noexcept {
	for (auto&& v : counter) {
		v = 0;
	}
}

consteval std::array<size_t, 4> parse_operation(std::string_view str) noexcept {
	std::array<size_t, 4> result{};
	for (auto&& c : str) {
		if (c == '+') {
			++result[(size_t)operation::add];
		} else if (c == '-') {
			++result[(size_t)operation::sub];
		} else if (c == '*') {
			++result[(size_t)operation::mul];
		} else if (c == '/') {
			++result[(size_t)operation::div];
		}
	}

	return result;
}

template<const std::array<size_t, 4> Add, auto Counter = counter>
void increment_counter() noexcept {
	if constexpr (Add[(size_t)operation::add]) {
		Counter[(size_t)operation::add] += Add[(size_t)operation::add];
	}
	if constexpr (Add[(size_t)operation::sub]) {
		Counter[(size_t)operation::sub] += Add[(size_t)operation::sub];
	}
	if constexpr (Add[(size_t)operation::mul]) {
		Counter[(size_t)operation::mul] += Add[(size_t)operation::mul];
	}
	if constexpr (Add[(size_t)operation::div]) {
		Counter[(size_t)operation::div] += Add[(size_t)operation::div];
	}
}

#define op_impl(...) \
	(operation_counting::increment_counter<operation_counting::parse_operation(#__VA_ARGS__)>(), __VA_ARGS__)

#define op(...) op_impl(__VA_ARGS__)
// #define op(...) __VA_ARGS__

} // namespace operation_counting

#pragma endregion

#include <chrono>
#include <format>
#include <fstream>
#include <iostream>
#include <random>
#include <stacktrace>

namespace chr = std::chrono;
using clk = chr::high_resolution_clock;

#ifndef MATRIX_DEBUG
#define MATRIX_DEBUG 0
#endif

#if MATRIX_DEBUG
#define ASSERT(...)                                                               \
	if (not(__VA_ARGS__)) {                                                       \
		std::cout << "line " << __LINE__ << ": " << #__VA_ARGS__ << " failed\n "; \
		throw ::std::logic_error("bruh");                                         \
	}
#else
#define ASSERT(...) ((void)0)
#endif

/// @brief Klasa macierzy, dane przechowywane w jednym ciągłym wektorze (cache-friendly)
class Matrix {
public:
	Matrix(size_t rows, size_t cols) noexcept: _rows{ rows }, _cols{ cols }, _storage(_rows * _cols) {}

	auto begin() noexcept { return _storage.begin(); }

	auto begin() const noexcept { return _storage.begin(); }

	auto cbegin() const noexcept { return _storage.cbegin(); }

	auto rbegin() noexcept { return _storage.rbegin(); }

	auto rbegin() const noexcept { return _storage.rbegin(); }

	auto crbegin() const noexcept { return _storage.crbegin(); }

	auto end() noexcept { return _storage.end(); }

	auto end() const noexcept { return _storage.end(); }

	auto cend() const noexcept { return _storage.cend(); }

	auto rend() noexcept { return _storage.rend(); }

	auto rend() const noexcept { return _storage.rend(); }

	auto crend() const noexcept { return _storage.crend(); }

	Matrix(std::initializer_list<std::initializer_list<double>> data):
		_rows{ data.size() },
		_cols{ data.begin()->size() },
		_storage(_rows * _cols) {
		auto i = begin();
		for (auto&& row_list : data) {
			for (auto&& v : row_list) {
				*(i++) = v;
			}
		}
	}

	double& operator[](const std::array<size_t, 2> idx) noexcept { return _storage[_cols * idx[0] + idx[1]]; }

	const double& operator[](const std::array<size_t, 2> idx) const noexcept {
		return _storage[_cols * idx[0] + idx[1]];
	}

	Matrix& operator+=(const Matrix& other) noexcept(!MATRIX_DEBUG) {
		ASSERT(_rows >= other._rows and _cols >= other._cols);
		const auto ending = other.end();
		if (_rows == other._rows and _cols == other._cols) {
			auto i2 = other.begin();
			for (auto i1 = begin(); i2 != ending; ++i1, ++i2) {
				auto& v1 = *i1;
				auto& v2 = *i2;
				op(v1 += v2);
			}
		} else {
			const auto coldiff = _cols - other._cols;
			auto i1 = begin();
			auto i2 = other.begin();
			for (size_t col2 = 0; i2 != ending; ++i1, ++i2) {
				auto& v1 = *i1;
				auto& v2 = *i2;
				op(v1 += v2);

				col2 = (col2 + 1) % other.cols();
				if (col2 == 0) {
					i1 += coldiff;
				}
			}
		}

		return *this;
	}

	Matrix operator+(const Matrix& other) const noexcept {
		if (_rows >= other._rows and _cols >= other._cols) {
			auto temp = *this;
			temp += other;
			return temp;
		} else {
			auto temp = other;
			temp += *this;
			return temp;
		}
	}

	Matrix& operator-=(const Matrix& other) noexcept(!MATRIX_DEBUG) {
		ASSERT(_rows >= other._rows and _cols >= other._cols);
		const auto ending = other.end();
		if (_rows == other._rows and _cols == other._cols) {
			auto i2 = other.begin();
			for (auto i1 = begin(); i2 != ending; ++i1, ++i2) {
				auto& v1 = *i1;
				auto& v2 = *i2;
				op(v1 -= v2);
			}
		} else {
			const auto coldiff = _cols - other._cols;
			auto i1 = begin();
			auto i2 = other.begin();
			for (size_t col2 = 0; i2 != ending; ++i1, ++i2) {
				auto& v1 = *i1;
				auto& v2 = *i2;
				op(v1 -= v2);

				col2 = (col2 + 1) % other.cols();
				if (col2 == 0) {
					i1 += coldiff;
				}
			}
		}
		return *this;
	}

	Matrix operator-(const Matrix& other) const noexcept {
		auto temp = *this;
		temp -= other;
		return temp;
		/*if (not(_rows < other._rows) or not(_cols < other._cols)) {
			auto temp = other;
			temp -= *this;
			return temp;
		} else {
		}*/
	}

	Matrix operator*(const Matrix& other) const noexcept(!MATRIX_DEBUG) {
		ASSERT(_cols == other._rows);

		auto result = Matrix(_rows, other._cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < other._cols; ++j) {
				for (int k = 0; k < _cols; ++k) {
					auto& in_result = result._storage[i * other._cols + j];
					auto& in_this = this->_storage[i * _cols + k];
					auto& in_other = other._storage[k * other._cols + j];

					op(in_result += in_this * in_other);
				}
			}
		}

		return result;
	}

	/// @brief Kopiuje dane z innej macierzy
	/// @param idx - miejsce gdzie skopiować
	/// @param other - macierz z której kopiować
	/// @param subidx - miejsce skąd kopiować
	Matrix& set_at(
		const std::array<size_t, 2> idx,
		const Matrix& other,
		const std::array<size_t, 2> subidx = { 0, 0 }
	) noexcept {
		for (size_t row = idx[0]; row < idx[0] + other._rows - subidx[0] && row < _rows; ++row) {
			for (size_t column = idx[1]; column < idx[1] + other._cols - subidx[1] && column < _cols; ++column) {
				(*this)[{ row, column }] = other[{ row - idx[0] + subidx[0], column - idx[1] + subidx[1] }];
			}
		}
		return *this;
	}

	size_t rows() const noexcept { return _rows; }

	size_t cols() const noexcept { return _cols; }

	/// @brief Generuje losową macierz
	/// @param rows - liczba wierszy
	/// @param cols - liczba kolumn
	/// @param seed - seed do generatora liczb pseudolosowych
	/// @param a - początek
	/// @param b
	/// @return
	static Matrix random(
		size_t rows,
		size_t cols,
		std::random_device::result_type seed = 0,
		double a = 0.00000001,
		double b = 1.0
	) noexcept {
		auto mt = std::mt19937(seed);
		auto random_double = std::uniform_real_distribution<double>(a, b);

		auto result = Matrix(rows, cols);
		for (auto&& v : result) {
			do {
				v = random_double(mt);
			} while (v == 1.0);
		}

		return result;
	}

	void print() const noexcept {
		for (size_t row = 0; row < _rows; ++row) {
			for (size_t col = 0; col < _cols; ++col) {
				std::cout << (*this)[{ row, col }] << '\t';
			}
			std::cout << '\n';
		}
	}

	std::array<Matrix, 4> split() const noexcept {
		const size_t half_rows = _rows / 2;
		const size_t half_cols = _cols / 2;

		std::array<Matrix, 4> result{
			Matrix(half_rows, half_cols),
			Matrix(half_rows, _cols - half_cols),
			Matrix(_rows - half_rows, half_cols),
			Matrix(_rows - half_rows, _cols - half_cols),
		};
		result[0].set_at({ 0, 0 }, *this, { 0, 0 });
		result[1].set_at({ 0, 0 }, *this, { 0, half_cols });
		result[2].set_at({ 0, 0 }, *this, { half_rows, 0 });
		result[3].set_at({ 0, 0 }, *this, { half_rows, half_cols });

		return result;
	}

	std::array<Matrix, 4> split(const bool A) const noexcept {
		const size_t rows = _rows % 2 == 0 ? _rows : _rows + 1;
		const size_t cols = _cols % 2 == 0 ? _cols : _cols + 1;
		const bool padding_rows_A = A and _rows % 2;
		const bool padding_cols_B = not A and _cols % 2;
		const size_t half_rows = rows / 2;
		const size_t half_cols = cols / 2;

		std::array<Matrix, 4> result{
			Matrix(half_rows, half_cols),
			Matrix(half_rows, cols - half_cols - padding_cols_B),
			Matrix(rows - half_rows, half_cols),
			Matrix(rows - half_rows - padding_rows_A, cols - half_cols - padding_cols_B),
		};
		result[0].set_at({ 0, 0 }, *this, { 0, 0 });
		result[1].set_at({ 0, 0 }, *this, { 0, half_cols });
		result[2].set_at({ 0, 0 }, *this, { half_rows, 0 });
		result[3].set_at({ 0, 0 }, *this, { half_rows, half_cols });

		return result;
	}

	bool operator==(const Matrix& other) const noexcept {
		return _rows == other._rows and _cols == other._cols and [&]() {
			const auto iend = this->end();
			for (auto i1 = this->begin(), i2 = other.begin(); i1 != iend; ++i1, ++i2) {
				if (std::fabs(*i1 - *i2) > 1e-12) {
					return false;
				}
			}
			return true;
		}
		();
	}

	double diffsum(const Matrix& other) const noexcept {
		double result = 0;

		const auto ending = end();
		for (auto i1 = begin(), i2 = other.begin(); i1 != ending; ++i1, ++i2) {
			result += std::fabs(*i1 - *i2);
		}

		return result;
	}

private:

	const size_t _rows;
	const size_t _cols;
	std::vector<double> _storage;
};

/// @brief Rozdziela zadaną macierz na 4 mniejsze
/// @brief ich rozmiar jest potęgą 2
// std::array<Matrix, 4> split_matrix(const Matrix& A) noexcept {
//	const auto subrows = std::bit_ceil(A.rows()) >> 1;
//	const auto subcols = std::bit_ceil(A.cols()) >> 1;
//
//	std::array<Matrix, 4> result{
//		Matrix(subrows, subcols),
//		Matrix(subrows, A.cols() - subcols),
//		Matrix(A.rows() - subrows, subcols),
//		Matrix(A.rows() - subrows, A.cols() - subcols),
//	};
//
//	// kopiowanie wartości z A
//	result[0].set_at({ 0, 0 }, A, { 0, 0 });
//	result[1].set_at({ 0, 0 }, A, { 0, subcols });
//	result[2].set_at({ 0, 0 }, A, { subrows, 0 });
//	result[3].set_at({ 0, 0 }, A, { subrows, subcols });
//
//	return result;
// }

Matrix binet_recursive(const Matrix& A, const Matrix& B) noexcept(!MATRIX_DEBUG) {
	if (A.rows() == 1 or A.cols() == 1 or B.rows() == 1 or B.cols() == 1) {
		return A * B;
	}

	auto&& [A_11, A_12, A_21, A_22] = A.split();
	auto&& [B_11, B_12, B_21, B_22] = B.split();

	auto C = Matrix(A.rows(), B.cols());
	{
		auto temp = binet_recursive(A_11, B_11);
		temp += binet_recursive(A_12, B_21);
		C.set_at({ 0, 0 }, temp);
	}
	{
		auto temp = binet_recursive(A_11, B_12);
		temp += binet_recursive(A_12, B_22);
		C.set_at({ 0, B_11.cols() }, temp);
	}
	{
		auto temp = binet_recursive(A_21, B_11);
		temp += binet_recursive(A_22, B_21);
		C.set_at({ A_11.rows(), 0 }, temp);
	}
	{
		auto temp = binet_recursive(A_21, B_12);
		temp += binet_recursive(A_22, B_22);
		C.set_at({ A_11.rows(), B_11.cols() }, temp);
	}

	return C;
}

Matrix strassen_recursive(const Matrix& A, const Matrix& B) noexcept(!MATRIX_DEBUG) {
	if (A.rows() == 1 or A.cols() == 1 or B.rows() == 1 or B.cols() == 1) {
		return A * B;
	}

	auto&& [A_11, A_12, A_21, A_22] = A.split(true);
	auto&& [B_11, B_12, B_21, B_22] = B.split(false);

	auto M_1 = strassen_recursive(A_11 + A_22, B_11 + B_22);
	auto M_4 = strassen_recursive(A_22, B_21 - B_11);
	auto M_5 = strassen_recursive(A_11 + A_12, B_22);
	auto M_7 = strassen_recursive(A_12 - A_22, B_21 + B_22);

	auto C = Matrix(A.rows(), B.cols());
	{
		auto C_11 = M_1 + M_4 - M_5 + M_7;
		C.set_at({ 0, 0 }, C_11);
	}
	auto M_3 = strassen_recursive(A_11, B_12 - B_22);
	{
		auto C_12 = M_3 + M_5;
		C.set_at({ 0, B_11.cols() }, C_12);
	}
	auto M_2 = strassen_recursive(A_21 + A_22, B_11);
	{
		auto C_21 = M_2 + M_4;
		C.set_at({ A_11.rows(), 0 }, C_21);
	}
	auto M_6 = strassen_recursive(A_21 - A_11, B_11 + B_12);
	{
		auto C_22 = M_1 - M_2 + M_3 + M_6;
		C.set_at({ A_11.rows(), B_11.cols() }, C_22);
	}

	return C;
}

// Matrix binet_recursive2(
//	const Matrix& A,
//	const Matrix& B,
//	const size_t size,
//	const size_t A_curr_y = 0,
//	const size_t A_curr_x = 0,
//	const size_t B_curr_y = 0,
//	const size_t B_curr_x = 0
//) noexcept {
//	ASSERT(A.size() == B.size());
//	const auto subsize = std::bit_ceil(A.size()) >> 1;
//	if (A_curr_x >= size or A_curr_y >= size or B_curr_x >= size or B_curr_y >= size) {
//		return Matrix(A.size());
//	} else if (A.size() == 1) {
//		return op(Matrix({ { (A[{ 0, 0 }]) * (B[{ 0, 0 }]) } }));
//	}
//
//	auto&& [A_11, A_12, A_21, A_22] = split_matrix(A);
//	auto&& [B_11, B_12, B_21, B_22] = split_matrix(B);
//
//	auto C = Matrix(A.size());
//	// const auto subsize = A_11.size();
//
//	{
//		auto C_11 = binet_recursive(A_11, B_11, size, A_curr_y, A_curr_x, B_curr_y, B_curr_x);
//		C_11 += binet_recursive(A_12, B_21, size, A_curr_y, A_curr_x + subsize, B_curr_y + subsize, B_curr_x);
//		C.set_at({ 0, 0 }, C_11);
//	}
//	{
//		auto C_12 = binet_recursive(A_11, B_12, size, A_curr_y, A_curr_x, B_curr_y, B_curr_x + subsize);
//		C_12 += binet_recursive(A_12, B_22, size, A_curr_y, A_curr_x + subsize, B_curr_y + subsize, B_curr_x + subsize);
//		C.set_at({ 0, subsize }, C_12);
//	}
//	{
//		auto C_21 = binet_recursive(A_21, B_11, size, A_curr_y + subsize, A_curr_x, B_curr_y, B_curr_x);
//		C_21 += binet_recursive(A_22, B_21, size, A_curr_y + subsize, A_curr_x + subsize, B_curr_y + subsize, B_curr_x);
//		C.set_at({ subsize, 0 }, C_21);
//	}
//	{
//		auto C_22 = binet_recursive(A_21, B_12, size, A_curr_y + subsize, A_curr_x, B_curr_y, B_curr_x + subsize);
//		C_22 += binet_recursive(
//			A_22,
//			B_22,
//			size,
//			A_curr_y + subsize,
//			A_curr_x + subsize,
//			B_curr_y + subsize,
//			B_curr_x + subsize
//		);
//		C.set_at({ subsize, subsize }, C_22);
//	}
//
//	return C;
// }
//

int main() {
	// auto A11 = Matrix(2, 2);
	// auto A12 = Matrix(2, 2); // 2, 1 -> 2, 2
	// auto A21 = Matrix(2, 2); // 1, 2 -> 2, 2
	// auto A22 = Matrix(1, 2); // 1, 1 -> 1, 2
	// auto B11 = Matrix(2, 2);
	// auto B12 = Matrix(2, 1);
	// auto B21 = Matrix(2, 2); // 1, 2 -> 2, 2
	// auto B22 = Matrix(2, 1); // 1, 1 -> 2, 1

	// auto M1 = (A11 + A22) * (B11 + B22);
	// auto M2 = (A21 + A22) * (B11);
	// auto M3 = (A11) * (B12 - B22);
	// auto M4 = (A22) * (B21 - B11);
	// auto M5 = (A11 + A12) * (B22);
	// auto M6 = (A21 - A11) * (B11 + B12);
	// auto M7 = (A12 - A22) * (B21 + B22);

	// std::cout << std::format("M1: {}, {}\n", M1.rows(), M1.cols());
	// std::cout << std::format("M2: {}, {}\n", M2.rows(), M2.cols());
	// std::cout << std::format("M3: {}, {}\n", M3.rows(), M3.cols());
	// std::cout << std::format("M4: {}, {}\n", M4.rows(), M4.cols());
	// std::cout << std::format("M5: {}, {}\n", M5.rows(), M5.cols());
	// std::cout << std::format("M6: {}, {}\n", M6.rows(), M6.cols());
	// std::cout << std::format("M7: {}, {}\n", M7.rows(), M7.cols());

	// auto C11 = M1 + M4 - M5 + M7;
	// auto C12 = M3 + M5;
	// auto C21 = M2 + M4;
	// auto C22 = M1 - M2 + M3 + M6;

	// return 0;

	/*auto mat1 = Matrix{
		{ 1, 2, 3 },
		{ 1, 2, 3 },
		{ 1, 2, 3 }
	};
	auto mat2 = Matrix{
		{ 4, 5, 6 },
		{ 4, 5, 6 },
		{ 4, 5, 6 }
	};*/

	/*auto mat1 = Matrix{
		{ 1, 2 },
		{ 1, 2 }
	};
	auto mat2 = Matrix{
		{ 4, 5 },
		{ 4, 5 }
	};

	auto mat_normal = mat1 * mat2;
	auto mat_binet = binet_recursive(mat1, mat2);
	auto mat_strassen = strassen_recursive(mat1, mat2);

	ASSERT(mat_normal == mat_binet);
	ASSERT(mat_normal == mat_strassen);*/

	/*for (size_t i = 2; i != 100 + 1; ++i) {
		auto mat1 = Matrix::random(i, i, 0);
		auto mat2 = Matrix::random(i, i, 1);

		auto mat_normal = mat1 * mat2;
		auto mat_binet = binet_recursive(mat1, mat2);
		auto mat_strassen = strassen_recursive(mat1, mat2);

		std::cout
			<< std::format("{}:\n\t{}\n\t{}\n", i, mat_normal.diffsum(mat_binet), mat_normal.diffsum(mat_strassen));
	}

	return 0;*/

	// const size_t N = 256;

	auto mat1 = Matrix::random(10, 10, 0);
	auto mat2 = Matrix::random(10, 10, 1);

	auto mat_binet = binet_recursive(mat1, mat2);
	auto mat_strassen = strassen_recursive(mat1, mat2);

	std::cout << "[\n";
	mat1.print();
	std::cout << "]\n";
	std::cout << "[\n";
	mat2.print();
	std::cout << "]\n";
	std::cout << "[\n";
	mat_binet.print();
	std::cout << "]\n";
	std::cout << "[\n";
	mat_strassen.print();
	std::cout << "]\n";

	return 0;

	auto times_binet = std::ofstream("times_binet4.txt");
	auto times_strassen = std::ofstream("times_strassen4.txt");
	auto ops_binet = std::ofstream("ops_binet4.txt");
	auto ops_strassen = std::ofstream("ops_strassen4.txt");
	times_binet << "N\ttime(ms)\n";
	times_strassen << "N\ttime(ms)\n";
	ops_binet << "N\t+\t-\t*\t/\n";
	ops_strassen << "N\t+\t-\t*\t/\n";

	for (size_t i = 1; i <= 1'000; i += [i]() {
			 if (i < 150) {
				 return 1;
			 } else if (i < 250) {
				 return 5;
			 } else if (i < 500) {
				 return 10;
			 } else {
				 return 100;
			 }
		 }()) {
		const size_t N = i;
		//(size_t)1 << i;

		std::cout << i << '\n';

		auto mat1 = Matrix::random(N, N, 0);
		std::cout << i << '\n';
		auto mat2 = Matrix::random(N, N, 1);
		std::cout << i << '\n';

		operation_counting::reset();
		{
			auto start = clk::now();
			auto mat3 = binet_recursive(mat1, mat2);
			auto end = clk::now();

			/*std::cout << std::format(
				"BINET:\n+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
				operation_counting::counter[0],
				operation_counting::counter[1],
				operation_counting::counter[2],
				operation_counting::counter[3],
				chr::duration_cast<chr::milliseconds>(end - start)
			);*/
			times_binet << std::format("{}\t{}\n", N, chr::duration_cast<chr::milliseconds>(end - start).count());
			ops_binet << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(size_t)operation_counting::operation::add],
				operation_counting::counter[(size_t)operation_counting::operation::sub],
				operation_counting::counter[(size_t)operation_counting::operation::mul],
				operation_counting::counter[(size_t)operation_counting::operation::div]
			);
			operation_counting::reset();
		}
		std::cout << i << '\n';
		{
			auto start = clk::now();
			auto mat3 = strassen_recursive(mat1, mat2);
			auto end = clk::now();

			/*std::cout << std::format(
				"STRASSEN:\n+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
				operation_counting::counter[0],
				operation_counting::counter[1],
				operation_counting::counter[2],
				operation_counting::counter[3],
				chr::duration_cast<chr::milliseconds>(end - start)
			);*/
			times_strassen << std::format("{}\t{}\n", N, chr::duration_cast<chr::milliseconds>(end - start).count());
			ops_strassen << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(size_t)operation_counting::operation::add],
				operation_counting::counter[(size_t)operation_counting::operation::sub],
				operation_counting::counter[(size_t)operation_counting::operation::mul],
				operation_counting::counter[(size_t)operation_counting::operation::div]
			);
			operation_counting::reset();
		}
		// std::cout << '\n';
	}
}
