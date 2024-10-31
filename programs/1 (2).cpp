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

namespace chr = std::chrono;
using clk = chr::high_resolution_clock;

#if MATRIX_DEBUG
#define ASSERT(...)                                                               \
	if (not(__VA_ARGS__)) {                                                       \
		std::cout << "line " << __LINE__ << ": " << #__VA_ARGS__ << " failed\n "; \
		::std::abort();                                                           \
	}
#else
#define ASSERT(...) ((void)0)
#endif

class MatrixView;

/// @brief Klasa macierzy, dane przechowywane w jednym ciągłym wektorze (cache-friendly)
class Matrix {
public:

	Matrix(size_t size) noexcept: _size{ size }, _storage(_size * _size) {}

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
		_size{ data.size() },
		_storage{ decltype(_storage)(_size * _size) } {
		auto i = begin();
		for (auto&& list : data) {
			for (auto&& v : list) {
				*(i++) = v;
			}
		}
	}

	double& operator[](const std::array<size_t, 2> idx) noexcept { return _storage[_size * idx[0] + idx[1]]; }

	const double& operator[](const std::array<size_t, 2> idx) const noexcept {
		return _storage[_size * idx[0] + idx[1]];
	}

	Matrix& operator+=(const Matrix& other) noexcept {
		ASSERT(size() == other.size());
		const auto ending = end();
		auto i2 = other.begin();
		for (auto i1 = begin(); i1 != ending; ++i1, ++i2) {
			auto&& v1 = *i1;
			auto&& v2 = *i2;
			op(v1 += v2);
		}
		return *this;
	}

	Matrix operator+(const Matrix& other) const noexcept {
		auto temp = *this;
		temp += other;
		return temp;
	}

	Matrix& operator-=(const Matrix& other) noexcept {
		ASSERT(size() == other.size());
		const auto ending = end();
		auto i2 = other.begin();
		for (auto i1 = begin(); i1 != ending; ++i1, ++i2) {
			auto&& v1 = *i1;
			auto&& v2 = *i2;
			op(v1 -= v2);
		}
		return *this;
	}

	Matrix operator-(const Matrix& other) const noexcept {
		auto temp = *this;
		temp -= other;
		return temp;
	}

	// operator MatrixView() const noexcept;
	MatrixView as_view() const noexcept;

	/// @brief Kopiuje dane z innej macierzy
	/// @param idx - miejsce gdzie skopiować
	/// @param other - macierz z której kopiować
	/// @param subidx - miejsce skąd kopiować
	Matrix& set_at(
		const std::array<size_t, 2> idx,
		const MatrixView& other
		//, const std::array<size_t, 2> subidx = { 0, 0 }
	) noexcept;

	Matrix& add_at(
		const std::array<size_t, 2> idx,
		const MatrixView& other
		//, const std::array<size_t, 2> subidx = { 0, 0 }
	) noexcept;
	Matrix& sub_at(
		const std::array<size_t, 2> idx,
		const MatrixView& other
		//, const std::array<size_t, 2> subidx = { 0, 0 }
	) noexcept;

	size_t size() const noexcept { return _size; }

	size_t rows() const noexcept { return _size; }

	size_t cols() const noexcept { return _size; }

	/// @brief Generuje losową macierz
	/// @param size - rozmiar macierzy (wiersze i kolumny)
	/// @param seed - seed do generatora liczb pseudolosowych
	/// @param a - początek
	/// @param b
	/// @return
	static Matrix random(
		size_t size,
		std::random_device::result_type seed = 0,
		double a = 0.00000001,
		double b = 1.0
	) noexcept {
		auto mt = std::mt19937(seed);
		auto random_double = std::uniform_real_distribution<double>(a, b);

		auto result = Matrix(size);
		for (auto&& v : result) {
			do {
				v = random_double(mt);
			} while (v == 1.0);
		}

		return result;
	}

	void print() const noexcept {
		for (size_t row = 0; row != size(); ++row) {
			for (size_t col = 0; col != size(); ++col) {
				std::cout << (*this)[{ row, col }] << '\t';
			}
			std::cout << '\n';
		}
	}

private:

	friend class MatrixView;

	const size_t _size;
	std::vector<double> _storage{};
};

class MatrixViewIterator {
public:

	MatrixViewIterator(const double* ptr, size_t pos, const size_t orig_cols, const size_t sub_cols) noexcept:
		_ptr{ ptr },
		_pos{ pos },
		_orig_cols{ orig_cols },
		_sub_cols{ sub_cols } {}

	bool operator++(int) noexcept {
		_pos = (_pos + 1) % _sub_cols;
		if (_pos == 0) {
			_ptr += _orig_cols - _sub_cols + 1;
			return true;
		} else {
			++_ptr;
			return false;
		}
	}

	const double& operator*() const noexcept { return *_ptr; }

	bool operator==(const MatrixViewIterator& other) const noexcept { return _ptr == other._ptr; }

private:

	const size_t _orig_cols;
	const size_t _sub_cols;
	const double* _ptr;
	size_t _pos;
};

class MatrixView {
public:

	MatrixView(const Matrix& matrix, const std::array<size_t, 2> pos, const std::array<size_t, 2> size) noexcept:
		_storage{ &matrix[pos] },
		_orig_size{ matrix.size() },
		_pos{ pos },
		_size{ size } {}

	MatrixView(const MatrixView& matrix, const std::array<size_t, 2> pos, const std::array<size_t, 2> size) noexcept:
		_storage{ matrix._storage + pos[1] },
		_orig_size{ matrix._orig_size },
		_pos{ pos[0] + matrix._pos[0], pos[1] + matrix._pos[1] },
		_size{ size } {}

	const double& operator[](const std::array<size_t, 2> idx) const noexcept {
		return _storage[_orig_size * (idx[0] + _pos[0]) + idx[1] + _pos[1]];
	}

	size_t rows() const noexcept { return _size[0]; }

	size_t cols() const noexcept { return _size[1]; }

	MatrixViewIterator begin() const noexcept { return MatrixViewIterator(_storage, 0, _orig_size, _size[1]); }

	MatrixViewIterator end() const noexcept {
		auto result =
			MatrixViewIterator(_storage + (_size[0] - 1) * _orig_size + _size[1] - 1, _orig_size, _orig_size, _size[1]);
		result++;
		return result;
	}

	Matrix operator*(const MatrixView& other) const noexcept { auto result = Matrix(); }

	operator bool() const noexcept { return _size[0] != 0 and _size[1] != 0; }

	void print() const noexcept {
		auto i = begin();
		const auto it_end = end();

		while (i != it_end) {
			std::cout << (*i) << '\t';
			if (i++) {
				std::cout << '\n';
			}
		}
	}

private:

	const double* _storage;
	const size_t _orig_size;
	const std::array<size_t, 2> _pos;
	const std::array<size_t, 2> _size;
};

MatrixView Matrix::as_view() const noexcept {
	return MatrixView(*this, { 0, 0 }, { size(), size() });
}

Matrix& Matrix::set_at(const std::array<size_t, 2> idx, const MatrixView& other) noexcept {
	for (size_t row = idx[0]; row != idx[0] + other.rows() and row != _size; ++row) {
		for (size_t column = idx[1]; column != idx[1] + other.cols() and column != _size; ++column) {
			(*this)[{ row, column }] = other[{ row - idx[0], column - idx[1] }];
		}
	}
	return *this;
}

Matrix& Matrix::add_at(const std::array<size_t, 2> idx, const MatrixView& other) noexcept {
	for (size_t row = idx[0]; row != idx[0] + other.rows() and row != _size; ++row) {
		for (size_t column = idx[1]; column != idx[1] + other.cols() and column != _size; ++column) {
			auto& from_other = other[{ row - idx[0], column - idx[1] }];
			op((*this)[{ row, column }] += from_other);
		}
	}
	return *this;
}

Matrix& Matrix::sub_at(const std::array<size_t, 2> idx, const MatrixView& other) noexcept {
	for (size_t row = idx[0]; row != idx[0] + other.rows() and row != _size; ++row) {
		for (size_t column = idx[1]; column != idx[1] + other.cols() and column != _size; ++column) {
			auto& from_other = other[{ row - idx[0], column - idx[1] }];
			op((*this)[{ row, column }] -= from_other);
		}
	}
	return *this;
}

/// @brief Rozdziela zadaną macierz na 4 mniejsze
/// @brief ich rozmiar jest potęgą 2
// std::array<Matrix, 4> split_matrix(const Matrix& A) noexcept {
//	// rozszerza macierz do rozmiaru potęgi 2
//
//	// rozmiar podmacierzy po rozszerzeniu
//	const auto subsize = std::bit_ceil(A.size()) >> 1;
//	std::array<Matrix, 4> result{
//		Matrix(subsize),
//		Matrix(subsize),
//		Matrix(subsize),
//		Matrix(subsize),
//	};
//
//	// kopiowanie wartości z A
//	result[0].set_at({ 0, 0 }, A, { 0, 0 });
//	result[1].set_at({ 0, 0 }, A, { 0, subsize });
//	result[2].set_at({ 0, 0 }, A, { subsize, 0 });
//	result[3].set_at({ 0, 0 }, A, { subsize, subsize });
//
//	return result;
// }

std::array<MatrixView, 4> split_matrix(const MatrixView& A) noexcept {
	// rozszerza macierz do rozmiaru potęgi 2

	const auto subsize = std::bit_ceil(std::max(A.rows(), A.cols())) >> 1;
	std::array<MatrixView, 4> result{
		MatrixView(A, { 0, 0 }, { subsize, subsize }),
		MatrixView(A, { 0, subsize }, { subsize, A.cols() - subsize }),
		MatrixView(A, { subsize, 0 }, { A.rows() - subsize, subsize }),
		MatrixView(A, { subsize, subsize }, { A.rows() - subsize, A.cols() - subsize }),
	};

	return result;
}

Matrix binet_recursive2(const MatrixView& A, const MatrixView& B) noexcept {
	if (A.rows() == A.cols()) {
		const auto subsize = std::bit_ceil(A.rows()) >> 1;
		if (A.rows() == 1) {
			return op(Matrix({ { (A[{ 0, 0 }]) * (B[{ 0, 0 }]) } }));
		}

		auto&& [A_11, A_12, A_21, A_22] = split_matrix(A);
		auto&& [B_11, B_12, B_21, B_22] = split_matrix(B);

		auto C = Matrix(A.rows());
		// const auto subsize = A_11.size();

		{
			auto C_11 = binet_recursive(A_11, B_11);
			C_11 += binet_recursive(A_12, B_21);
			C.set_at({ 0, 0 }, C_11);
		}
		{
			auto C_12 = binet_recursive(A_11, B_12);
			C_12 += binet_recursive(A_12, B_22);
			C.set_at({ 0, subsize }, C_12);
		}
		{
			auto C_21 = binet_recursive(A_21, B_11);
			C_21 += binet_recursive(A_22, B_21);
			C.set_at({ subsize, 0 }, C_21);
		}
		{
			auto C_22 = binet_recursive(A_21, B_12);
			C_22 += binet_recursive(A_22, B_22);
			C.set_at({ subsize, subsize }, C_22);
		}

		return C;
	} else {
		if (A.cols() == 1) {
			auto&& [A11, A12, A21, A22] = split_matrix(A);
			auto&& [B11, B12, B21, B22] = split_matrix(B);

			auto C = Matrix()
		}
	}
}

Matrix binet_recursive(const Matrix& A, const Matrix& B) noexcept {}

// Matrix strassen_recursive(const Matrix& A, const Matrix& B) noexcept {
//	ASSERT(A.size() == B.size());
//	if (A.size() == 1) {
//		return op(Matrix({ { (A[{ 0, 0 }]) * (B[{ 0, 0 }]) } }));
//	}
//
//	auto&& [A_11, A_12, A_21, A_22] = split_matrix(A);
//	auto&& [B_11, B_12, B_21, B_22] = split_matrix(B);
//
//	auto M_1 = strassen_recursive(A_11 + A_22, B_11 + B_22);
//	auto M_2 = strassen_recursive(A_21 + A_22, B_11);
//	auto M_3 = strassen_recursive(A_11, B_12 - B_22);
//	auto M_4 = strassen_recursive(A_22, B_21 - B_11);
//	auto M_5 = strassen_recursive(A_11 + A_12, B_22);
//	auto M_6 = strassen_recursive(A_21 - A_11, B_11 + B_12);
//	auto M_7 = strassen_recursive(A_12 - A_22, B_21 + B_22);
//
//	auto C = Matrix(A.size());
//	const auto subsize = A_11.size();
//
//	{
//		auto C_11 = M_1 + M_4 - M_5 + M_7;
//		C.set_at({ 0, 0 }, C_11);
//	}
//	{
//		auto C_12 = M_3 + M_5;
//		C.set_at({ 0, subsize }, C_12);
//	}
//	{
//		auto C_21 = M_2 + M_4;
//		C.set_at({ subsize, 0 }, C_21);
//	}
//	{
//		auto C_22 = M_1 - M_2 + M_3 + M_6;
//		C.set_at({ subsize, subsize }, C_22);
//	}
//
//	return C;
// }

int main() {
	auto mat1 = Matrix{
		{ 1, 2, 3 },
		{ 1, 2, 3 },
		{ 1, 2, 3 }
	};
	auto mat2 = Matrix{
		{ 4, 5, 6 },
		{ 4, 5, 6 },
		{ 4, 5, 6 }
	};

	auto view = MatrixView(mat1, { 0, 0 }, { 2, 3 });

	auto subviews = split_matrix(view);

	for (auto&& sv : subviews) {
		std::cout << (bool)sv << '\n';
	}

	/*mat1 = Matrix{
		{ 1, 2 },
		{ 1, 2 }
	};
	mat2 = Matrix{
		{ 4, 5 },
		{ 4, 5 }
	};*/

	// const size_t N = 256;

	// auto file = std::ofstream("times.txt");
	// file << "N\ttime(ms)\n";

	// for (size_t i = 2; i != 150 + 1; ++i) {
	//	const size_t N = i;
	//	//(size_t)1 << i;

	//	auto mat1 = Matrix::random(N, 0);
	//	auto mat2 = Matrix::random(N, 1);

	//	operation_counting::reset();
	//	{
	//		auto start = clk::now();
	//		auto mat3 = binet_recursive(mat1, mat2, mat1.size());
	//		auto end = clk::now();

	//		/*std::cout << std::format(
	//			"BINET:\n+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
	//			operation_counting::counter[0],
	//			operation_counting::counter[1],
	//			operation_counting::counter[2],
	//			operation_counting::counter[3],
	//			chr::duration_cast<chr::milliseconds>(end - start)
	//		);*/
	//		file << std::format("{}\t{}\n", N, chr::duration_cast<chr::milliseconds>(end - start).count());
	//		operation_counting::reset();
	//	}
	//	/*{
	//		auto start = clk::now();
	//		auto mat3 = strassen_recursive(mat1, mat2);
	//		auto end = clk::now();

	//		std::cout << std::format(
	//			"STRASSEN:\n+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
	//			operation_counting::counter[0],
	//			operation_counting::counter[1],
	//			operation_counting::counter[2],
	//			operation_counting::counter[3],
	//			chr::duration_cast<chr::milliseconds>(end - start)
	//		);
	//		operation_counting::reset();
	//	}*/
	//	// std::cout << '\n';
	//}
}
