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
#include <iostream>
#include <random>

namespace chr = std::chrono;
using clk = chr::high_resolution_clock;

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

	/// @brief Kopiuje dane z innej macierzy
	/// @param idx - miejsce gdzie skopiować
	/// @param other - macierz z której kopiować
	/// @param subidx - miejsce skąd kopiować
	Matrix& set_at(
		const std::array<size_t, 2> idx,
		const Matrix& other,
		const std::array<size_t, 2> subidx = { 0, 0 }
	) noexcept {
		for (size_t row = idx[0]; row != idx[0] + other.size() - subidx[0] and row != _size; ++row) {
			for (size_t column = idx[1]; column != idx[1] + other.size() - subidx[1] and column != _size; ++column) {
				(*this)[{ row, column }] = other[{ row - idx[0] + subidx[0], column - idx[1] + subidx[1] }];
			}
		}
		return *this;
	}

	size_t size() const noexcept { return _size; }

	size_t rows() const noexcept { return _size; }

	size_t columns() const noexcept { return _size; }

	static Matrix random(
		size_t size,
		std::random_device::result_type seed = 0,
		double a = 0.00000001,
		double b = 1.0
	) noexcept {
		auto mt = std::mt19937(seed);
		auto random_T = std::uniform_real_distribution<double>(a, b);

		auto result = Matrix(size);
		for (auto&& v : result) {
			v = random_T(mt);
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

	const size_t _size;
	std::vector<double> _storage{};
};

/// @brief Rozdziela zadaną macierz na 4 mniejsze
/// @brief ich rozmiar jest potęgą 2
std::array<Matrix, 4> split_matrix(const Matrix& A) noexcept {
	// rozszerza macierz do rozmiaru potęgi 2

	// rozmiar podmacierzy po rozszerzeniu
	const auto subsize = std::bit_ceil(A.size()) >> 1;
	std::array<Matrix, 4> result{
		Matrix(subsize),
		Matrix(subsize),
		Matrix(subsize),
		Matrix(subsize),
	};

	// kopiowanie wartości z A
	result[0].set_at({ 0, 0 }, A, { 0, 0 });
	result[1].set_at({ 0, 0 }, A, { 0, subsize });
	result[2].set_at({ 0, 0 }, A, { subsize, 0 });
	result[3].set_at({ 0, 0 }, A, { subsize, subsize });

	return result;
}

Matrix binet_recursive(const Matrix& A, const Matrix& B) noexcept {
	if (A.size() == 1) {
		return op(Matrix({ { (A[{ 0, 0 }]) * (B[{ 0, 0 }]) } }));
	}

	auto&& [A_11, A_12, A_21, A_22] = split_matrix(A);
	auto&& [B_11, B_12, B_21, B_22] = split_matrix(B);

	auto C = Matrix(A.size());
	const auto subsize = A_11.size();

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
}

Matrix strassen_recursive(const Matrix& A, const Matrix& B) noexcept {
	if (A.size() == 1) {
		return op(Matrix({ { (A[{ 0, 0 }]) * (B[{ 0, 0 }]) } }));
	}

	auto&& [A_11, A_12, A_21, A_22] = split_matrix(A);
	auto&& [B_11, B_12, B_21, B_22] = split_matrix(B);

	auto M_1 = strassen_recursive(A_11 + A_22, B_11 + B_22);
	auto M_2 = strassen_recursive(A_21 + A_22, B_11);
	auto M_3 = strassen_recursive(A_11, B_12 - B_22);
	auto M_4 = strassen_recursive(A_22, B_21 - B_11);
	auto M_5 = strassen_recursive(A_11 + A_12, B_22);
	auto M_6 = strassen_recursive(A_21 - A_11, B_11 + B_12);
	auto M_7 = strassen_recursive(A_12 - A_22, B_21 + B_22);

	auto C = Matrix(A.size());
	const auto subsize = A_11.size();

	{
		auto C_11 = M_1 + M_4 - M_5 + M_7;
		C.set_at({ 0, 0 }, C_11);
	}
	{
		auto C_12 = M_3 + M_5;
		C.set_at({ 0, subsize }, C_12);
	}
	{
		auto C_21 = M_2 + M_4;
		C.set_at({ subsize, 0 }, C_21);
	}
	{
		auto C_22 = M_1 - M_2 + M_3 + M_6;
		C.set_at({ subsize, subsize }, C_22);
	}

	return C;
}

int main() {
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

	/*mat1 = Matrix{
		{ 1, 2 },
		{ 1, 2 }
	};
	mat2 = Matrix{
		{ 4, 5 },
		{ 4, 5 }
	};*/

	// należy założyć że N zwiększy się do najbliższej >= potęgi 2
	// const size_t N = 256;

	for (size_t i = 1; i != 9 + 1; ++i) {
		const size_t N = (size_t)1 << i;

		auto mat1 = Matrix::random(N, 0);
		auto mat2 = Matrix::random(N, 1);

		{
			auto start = clk::now();
			auto mat3 = binet_recursive(mat1, mat2);
			auto end = clk::now();

			std::cout << std::format("{}\n", (end - start).count());

			/*std::cout << std::format(
				"+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
				operation_counting::counter[0],
				operation_counting::counter[1],
				operation_counting::counter[2],
				operation_counting::counter[3],
				chr::duration_cast<chr::milliseconds>(end - start)
			);*/
		}
		operation_counting::reset();
		{
			auto start = clk::now();
			auto mat3 = strassen_recursive(mat1, mat2);
			auto end = clk::now();

			std::cout << std::format("{}\n", (end - start).count());

			/*std::cout << std::format(
				"+: {}\n-: {}\n*: {}\n/: {}\ntime: {}\n",
				operation_counting::counter[0],
				operation_counting::counter[1],
				operation_counting::counter[2],
				operation_counting::counter[3],
				chr::duration_cast<chr::milliseconds>(end - start)
			);*/
		}
		std::cout << '\n';
	}
}
