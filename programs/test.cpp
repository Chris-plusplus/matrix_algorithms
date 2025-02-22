#include <iostream>
#include <format>
#include <array>
#include <string_view>
#include <chrono>
#include <cstdlib>

namespace operation_counting {

enum class operation {
	add,
	sub,
	mul,
	div,
	unknown,
};

size_t counter[(size_t)operation::unknown]{};

consteval std::array<size_t, 4> parse_operation(std::string_view str) noexcept {
	std::array<size_t, 4> result{};
	for (auto&& c : str) {
		if (c == '+') {
			++result[(size_t)operation::add];
		}
		else if (c == '-') {
			++result[(size_t)operation::sub];
		}
		else if (c == '*') {
			++result[(size_t)operation::mul];
		}
		else if (c == '/') {
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

#define op_impl(...) (operation_counting::increment_counter<operation_counting::parse_operation(#__VA_ARGS__)>(), __VA_ARGS__)

#define op(...) op_impl(__VA_ARGS__)
//#define op(...) __VA_ARGS__

}

namespace chr = std::chrono;
using clk = chr::high_resolution_clock;

int main() {
	double x = rand();
	double y = rand();

	for (size_t j = 0; j != 20; ++j) {
		auto start = clk::now();
		for (size_t i = 0; i != 10'000'000; ++i) {
			op(x = x + y);
			op(x = x * y + x);
			op(x = x / y + x);
			op(x = x - y);
		}
		auto end = clk::now();
		std::cout << std::format("10'000'000 * 4 ops took: {}\n", chr::duration_cast<chr::milliseconds>(end - start));
	}

	std::cout << std::format("x = {}\ny = {}\n", x, y);
}