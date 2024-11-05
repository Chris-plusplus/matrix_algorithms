#include <algorithm>
#include <chrono>
#include <format>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <ratio>

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

using u32 = uint32_t;

class Matrix {
public:
	Matrix(u32 rows, u32 cols) noexcept: _rows{ rows }, _cols{ cols }, _storage(_rows * _cols) {}

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
		_rows{ (u32)data.size() },
		_cols{ (u32)data.begin()->size() },
		_storage(_rows * _cols) {
		auto i = begin();
		for (auto&& row_list : data) {
			for (auto&& v : row_list) {
				*(i++) = v;
			}
		}
	}

	double& operator[](const std::array<u32, 2> idx) noexcept { return _storage[_cols * idx[0] + idx[1]]; }

	double operator[](const std::array<u32, 2> idx) const noexcept {
		return idx[0] < _rows ? (idx[1] < _cols ? _storage[_cols * idx[0] + idx[1]] : 0) : 0;
	}

	double at_nocheck(const std::array<u32, 2> idx) const noexcept { return _storage[_cols * idx[0] + idx[1]]; }

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
			for (u32 col2 = 0; i2 != ending; ++i1, ++i2) {
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
			for (u32 col2 = 0; i2 != ending; ++i1, ++i2) {
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

	Matrix operator-(const Matrix& other) const noexcept(!MATRIX_DEBUG) {
		if (_rows >= other._rows and _cols >= other._cols) {
			auto temp = *this;
			temp -= other;
			return temp;
		} else {
			auto temp = Matrix(other.rows(), other.cols());
			for (u32 row = 0; row != temp._rows; ++row) {
				for (u32 col = 0; col != temp._cols; ++col) {
					temp[{ row, col }] = -other.at_nocheck({ row, col });
				}
			}
			temp += *this;
			return temp;
		}
		/*if (not(_rows < other._rows) or not(_cols < other._cols)) {
			auto temp = other;
			temp -= *this;
			return temp;
		} else {
		}*/
	}

	Matrix operator*(const Matrix& other) const noexcept(!MATRIX_DEBUG) {
		// ASSERT(_cols == other._rows);

		auto result = Matrix(_rows, other._cols);
		for (u32 i = 0; i < _rows; ++i) {
			for (u32 k = 0; k < _cols; ++k) {
				auto& in_this = this->_storage[i * _cols + k];
				for (u32 j = 0; j < other._cols; ++j) {
					auto& in_result = result._storage[i * other._cols + j];
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
		const std::array<u32, 2> idx,
		const Matrix& other,
		const std::array<u32, 2> subidx = { 0, 0 }
	) noexcept {
		for (u32 row = idx[0]; row < idx[0] + other._rows - subidx[0] && row < _rows; ++row) {
			for (u32 column = idx[1]; column < idx[1] + other._cols - subidx[1] && column < _cols; ++column) {
				(*this)[{ row, column }] = other[{ row - idx[0] + subidx[0], column - idx[1] + subidx[1] }];
			}
		}
		return *this;
	}

	u32 rows() const noexcept { return _rows; }

	u32 cols() const noexcept { return _cols; }

	/// @brief Generuje losową macierz
	/// @param rows - liczba wierszy
	/// @param cols - liczba kolumn
	/// @param seed - seed do generatora liczb pseudolosowych
	/// @param a - początek
	/// @param b
	/// @return
	static Matrix random(
		u32 rows,
		u32 cols,
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

	static Matrix eye(const u32 size) noexcept {
		auto result = Matrix(size, size);
		for (u32 row_col = 0; row_col != size; ++row_col) {
			result[{ row_col, row_col }] = 1;
		}
		return result;
	}

	void print(std::string_view msg = "") const noexcept {
		std::cout << msg << (msg.length() ? " " : "") << "[\n";
		for (u32 row = 0; row < _rows; ++row) {
			for (u32 col = 0; col < _cols; ++col) {
				std::cout << (*this)[{ row, col }];
				if (col < _cols - 1) {
					std::cout << ",\t";
				}
			}
			std::cout << ";\n";
		}
		std::cout << "]\n";
	}

	std::array<Matrix, 4> split() const noexcept {
		const u32 half_rows = _rows / 2;
		const u32 half_cols = _cols / 2;

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
		const u32 rows = _rows % 2 == 0 ? _rows : _rows + 1;
		const u32 cols = _cols % 2 == 0 ? _cols : _cols + 1;
		const bool padding_rows_A = A and _rows % 2;
		const bool padding_cols_B = not A and _cols % 2;
		const u32 half_rows = rows / 2;
		const u32 half_cols = cols / 2;

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

	Matrix operator-() const noexcept {
		auto temp = Matrix(_rows, _cols);
		temp -= *this;
		return temp;
	}

private:

	const u32 _rows;
	const u32 _cols;
	std::vector<double> _storage;
};

Matrix smul(const Matrix& A, const Matrix& B) noexcept(!MATRIX_DEBUG) {
	if (A.rows() == 1 or A.cols() == 1 or B.rows() == 1 or B.cols() == 1) {
		return A * B;
	}

	auto&& [A_11, A_12, A_21, A_22] = A.split(true);
	auto&& [B_11, B_12, B_21, B_22] = B.split(false);

	auto M_1 = smul(A_11 + A_22, B_11 + B_22);
	auto M_4 = smul(A_22, B_21 - B_11);
	auto M_5 = smul(A_11 + A_12, B_22);
	auto M_7 = smul(A_12 - A_22, B_21 + B_22);

	auto C = Matrix(A.rows(), B.cols());
	{
		auto C_11 = M_1 + M_4 - M_5 + M_7;
		C.set_at({ 0, 0 }, C_11);
	}
	auto M_3 = smul(A_11, B_12 - B_22);
	{
		auto C_12 = M_3 + M_5;
		C.set_at({ 0, B_11.cols() }, C_12);
	}
	auto M_2 = smul(A_21 + A_22, B_11);
	{
		auto C_21 = M_2 + M_4;
		C.set_at({ A_11.rows(), 0 }, C_21);
	}
	auto M_6 = smul(A_21 - A_11, B_11 + B_12);
	{
		auto C_22 = M_1 - M_2 + M_3 + M_6;
		C.set_at({ A_11.rows(), B_11.cols() }, C_22);
	}

	return C;
}

Matrix inverse(const Matrix& A) noexcept(!MATRIX_DEBUG) {
	if (A.rows() == A.cols() and A.rows() == 1) {
		return op(Matrix({ { 1.0 / A.at_nocheck({ 0, 0 }) } }));
	}

	auto&& B = Matrix(A.rows(), A.cols());

	auto&& [A11, A12, A21, A22] = A.split();

	auto&& A11_inv = inverse(A11);

	auto&& I11 = Matrix::eye(A11.rows());
	auto&& I22 = Matrix::eye(A22.rows());

	// auto&& M1 = A21 * A11_inv;
	auto&& M1 = smul(A21, A11_inv);

	// auto&& S22 = A22 - M1 * A12;
	auto&& S22 = A22 - smul(M1, A12);
	auto&& S22_inv = inverse(S22);

	// auto&& M2 = S22_inv * M1;
	auto&& M2 = smul(S22_inv, M1);

	// auto&& B11 = A11_inv * (I11 + A12 * M2);
	auto&& B11 = smul(A11_inv, I11 + smul(A12, M2));
	// auto&& B12 = -A11_inv * A12 * S22_inv;
	auto&& B12 = smul(smul(-A11_inv, A12), S22_inv);
	auto&& B21 = -M2;
	auto&& B22 = S22_inv;

	B.set_at({ 0, 0 }, B11);
	B.set_at({ 0, B11.cols() }, B12);
	B.set_at({ B11.rows(), 0 }, B21);
	B.set_at({ B11.rows(), B11.cols() }, B22);

	return B;
}

std::array<Matrix, 2> LU(const Matrix& A) noexcept(!MATRIX_DEBUG) {
	if (A.rows() == A.cols() and A.rows() == 1) {
		return { Matrix({ { 1.0 } }), Matrix({ { A.at_nocheck({ 0, 0 }) } }) };
	}

	auto&& [A11, A12, A21, A22] = A.split();

	auto&& [L11, U11] = LU(A11);

	auto&& L11_inv = inverse(L11);
	auto&& U11_inv = inverse(U11);

	auto&& L21 = smul(A21, U11_inv);
	auto&& U12 = smul(L11_inv, A12);
	auto&& S = A22 - smul(L21, U12);

	auto&& [Ls, Us] = LU(S);
	auto&& L22 = Ls;
	auto&& U22 = Us;

	auto&& result = std::array{ Matrix(A.rows(), A.cols()), Matrix(A.rows(), A.cols()) };

	// L
	result[0].set_at({ 0, 0 }, L11);
	// result[0].set_at({ 0, A.cols() }, L12);
	result[0].set_at({ A11.rows(), 0 }, L21);
	result[0].set_at({ A11.rows(), A11.cols() }, L22);

	// U
	result[1].set_at({ 0, 0 }, U11);
	result[1].set_at({ 0, A11.cols() }, U12);
	// result[1].set_at({ A.rows(), 0 }, U21);
	result[1].set_at({ A11.rows(), A11.cols() }, U22);

	return result;
}

Matrix gauss_elimination_recursive(const Matrix& A, const Matrix& b) noexcept(!MATRIX_DEBUG) {
	const u32 N = A.rows();

	/*if (N != A.cols()) {
		throw std::invalid_argument("Matrix A must be square.");
	}

	if (b.rows() != N || b.cols() != 1) {
		throw std::invalid_argument("Matrix b must have the same number of rows as A and exactly one column.");
	}*/

	if (N == 1) {
		Matrix x(1, 1);
		op(x[{ 0, 0 }] = b[{ 0, 0 }] / A[{ 0, 0 }]);
		return x;
	}

	auto&& blocks = A.split();
	const Matrix& A11 = blocks[0];
	const Matrix& A12 = blocks[1];
	const Matrix& A21 = blocks[2];
	const Matrix& A22 = blocks[3];

	Matrix b1(A11.rows(), 1);
	Matrix b2(A22.rows(), 1);
	for (u32 i = 0; i < A11.rows(); ++i) {
		b1[{ i, 0 }] = b[{ i, 0 }];
	}
	for (u32 i = 0; i < A22.rows(); ++i) {
		b2[{ i, 0 }] = b[{ i + A11.rows(), 0 }];
	}

	auto&& [L11, U11] = LU(A11);

	Matrix L11_inv = inverse(L11);
	Matrix U11_inv = inverse(U11);

	auto temp1 = smul(L11_inv, A12);

	Matrix S = A22 - smul(A21, smul(U11_inv, temp1));

	auto&& [Ls, Us] = LU(S);
	auto&& Ls_inv = inverse(Ls);

	auto&& C11 = U11;
	auto& C12 = temp1;
	// auto&& C21 = Matrix(A21.rows(), A21.cols());
	auto&& C22 = Us;

	Matrix RHS1 = smul(L11_inv, b1);
	Matrix RHS2 = smul(Ls_inv, b2) - smul(Ls_inv, smul(A21, smul(U11_inv, RHS1)));

	auto&& result = Matrix(A.rows(), A.cols() + 1);

	result.set_at({ 0, 0 }, C11);
	result.set_at({ 0, A11.cols() }, C12);
	// result.set_at({ A11.rows(), 0 }, C21);
	result.set_at({ A11.rows(), A11.cols() }, C22);
	result.set_at({ 0, A11.cols() + A12.cols() }, RHS1);
	result.set_at({ A11.rows(), A11.cols() + A12.cols() }, RHS2);

	return result;
}

double determinant_recursive(const Matrix& A) noexcept(!MATRIX_DEBUG) {
	const u32 N = A.rows();

	/*if (N != A.cols()) {
		throw std::runtime_error("Only square matrices have a determinant");
	}*/

	auto [L, U] = LU(A);

	double det = 1.0;
	for (u32 i = 0; i < N; ++i) {
		op(det *= L[{ i, i }] * U[{ i, i }]);
	}
	return det;
}

constexpr double to_ms =
	(double)std::ratio_divide<std::nano, std::milli>::num / (double)std::ratio_divide<std::nano, std::milli>::den;

int main() noexcept(!MATRIX_DEBUG) {
	/*auto A = Matrix({ { 4 } });
	auto B = Matrix({
		{ 1, 2 },
		{ 2, 1 }
	 });
	auto C = Matrix({
		{ 1, 2, 3 },
		{ 2, 1, 2 },
		{ 3, 2, 1 }
	});

	auto&& [L, U] = LU(B);
	L.print();
	U.print();*/

	auto&& A = Matrix::random(10, 10, 0);
	auto&& b = Matrix::random(10, 1, 0);

	auto gauss_result = gauss_elimination_recursive(A, b);
	gauss_result.print();

	return 0;

	/*for (u32 i = 1; i <= 100; ++i) {
		auto&& mat = Matrix::random(i, i, 0);

		auto&& [L, U] = LU(mat);
		auto&& inv = inverse(mat);

		std::cout
			<< std::format("{}:\n\tLU: {}\n\tinv: {}\n", i, mat.diffsum(L * U), Matrix::eye(1).diffsum(mat * inv));
	}
	return 0;*/

	getchar();

	auto times = std::ofstream("times.txt");
	auto ops_inverse = std::ofstream("ops_inverse.txt");
	auto ops_LU = std::ofstream("ops_LU.txt");
	auto ops_gauss = std::ofstream("ops_gauss.txt");
	auto ops_det = std::ofstream("ops_det.txt");
	times << "N\tinverse\tLU\tgauss\tdet\n";
	ops_inverse << "N\t+\t-\t*\t/\tsum\n";
	ops_LU << "N\t+\t-\t*\t/\tsum\n";
	ops_gauss << "N\t+\t-\t*\t/\tsum\n";
	ops_det << "N\t+\t-\t*\t/\tsum\n";

	for (u32 i = 5; i <= 500; i += 5/* [i]() {
			 if (i < 150) {
				 return 1;
			 } else if (i < 250) {
				 return 5;
			 } else if (i < 500) {
				 return 10;
			 } else {
				 return 100;
			 }
		 }()*/) {
		const u32 N = i;

		auto&& A = Matrix::random(N, N, 0);
		auto&& b = Matrix::random(N, 1, 0);

		times << std::format("{}\t", N);

		operation_counting::reset();
		{
			auto start = clk::now();
			auto&& A_inv = inverse(A);
			auto end = clk::now();

			times << std::format("{}\t", (end - start).count() * to_ms);
			ops_inverse << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(u32)operation_counting::operation::add],
				operation_counting::counter[(u32)operation_counting::operation::sub],
				operation_counting::counter[(u32)operation_counting::operation::mul],
				operation_counting::counter[(u32)operation_counting::operation::div],
				std::reduce(operation_counting::counter, operation_counting::counter + 4)
			);
			operation_counting::reset();
		}
		std::cout << i << '\n';
		{
			auto start = clk::now();
			auto&& [L, U] = LU(A);
			auto end = clk::now();

			times << std::format("{}\t", (end - start).count() * to_ms);
			ops_LU << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(u32)operation_counting::operation::add],
				operation_counting::counter[(u32)operation_counting::operation::sub],
				operation_counting::counter[(u32)operation_counting::operation::mul],
				operation_counting::counter[(u32)operation_counting::operation::div],
				std::reduce(operation_counting::counter, operation_counting::counter + 4)
			);
			operation_counting::reset();
		}
		std::cout << i << '\n';
		{
			auto start = clk::now();
			auto&& gauss = gauss_elimination_recursive(A, b);
			auto end = clk::now();

			times << std::format("{}\t", (end - start).count() * to_ms);
			ops_gauss << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(u32)operation_counting::operation::add],
				operation_counting::counter[(u32)operation_counting::operation::sub],
				operation_counting::counter[(u32)operation_counting::operation::mul],
				operation_counting::counter[(u32)operation_counting::operation::div],
				std::reduce(operation_counting::counter, operation_counting::counter + 4)
			);
			operation_counting::reset();
		}
		std::cout << i << '\n';
		{
			auto start = clk::now();
			auto&& det = determinant_recursive(A);
			auto end = clk::now();

			times << std::format("{}\n", (end - start).count() * to_ms);
			ops_det << std::format(
				"{}\t{}\t{}\t{}\t{}\n",
				N,
				operation_counting::counter[(u32)operation_counting::operation::add],
				operation_counting::counter[(u32)operation_counting::operation::sub],
				operation_counting::counter[(u32)operation_counting::operation::mul],
				operation_counting::counter[(u32)operation_counting::operation::div],
				std::reduce(operation_counting::counter, operation_counting::counter + 4)
			);
			operation_counting::reset();
		}
		std::cout << i << '\n';
	}
}
