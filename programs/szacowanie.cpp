#include <array>
#include <format>
#include <iostream>
#include <map>
#include <ranges>

int main() {
	auto map = std::map<unsigned int, std::array<double, 4>>();

	size_t size;
	std::cout << "sample size: ";
	std::cin >> size;

	for (size_t i = 0; i != size; ++i) {
		unsigned int N;
		double t1;
		double t2;
		double t3;
		double t4;

		std::cin >> N >> t1 >> t2 >> t3 >> t4;
		map[N] = { t1, t2, t3, t4 };
	}

	auto get = [&map](unsigned int N) -> std::array<double, 4> {
		if (map.contains(N)) {
			return map[N];
		} else {
			auto ub = map.upper_bound(N);
			auto ubm1 = std::prev(ub);

			return { (ub->second[0] + ubm1->second[0]) / 2,
					 (ub->second[1] + ubm1->second[1]) / 2,
					 (ub->second[2] + ubm1->second[2]) / 2,
					 (ub->second[3] + ubm1->second[3]) / 2 };
		}
	};

	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	double sum4 = 0;

	unsigned int counter = 0;
	for (auto&& [key, val] : std::views::reverse(map)) {
		if (key > 5) {
			++counter;
			auto got = get(key / 2);
			sum1 += double(val[0] / got[0]);
			sum2 += double(val[1] / got[1]);
			sum3 += double(val[2] / got[2]);
			sum4 += double(val[3] / got[3]);
		}
	}
	sum1 /= counter;
	sum2 /= counter;
	sum3 /= counter;
	sum4 /= counter;

	std::cout << sum1 << '\n';
	std::cout << sum2 << '\n';
	std::cout << sum3 << '\n';
	std::cout << sum4 << '\n';
}
