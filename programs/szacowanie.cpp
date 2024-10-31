#include <format>
#include <iostream>
#include <map>
#include <ranges>

int main() {
	auto map = std::map<unsigned int, std::pair<unsigned int, unsigned int>>();

	size_t size;
	std::cout << "sample size: ";
	std::cin >> size;

	for (size_t i = 0; i != size; ++i) {
		unsigned int N;
		unsigned int t1;
		unsigned int t2;

		std::cin >> N >> t1 >> t2;
		map[N] = { t1, t2 };
	}

	auto get = [&map](unsigned int N) -> std::pair<unsigned int, unsigned int> {
		if (map.contains(N)) {
			return map[N];
		} else {
			auto ub = map.upper_bound(N);
			auto ubm1 = std::prev(ub);

			return { (ub->second.first + ubm1->second.first) / 2, (ub->second.second + ubm1->second.second) / 2 };
		}
	};

	double sum1 = 0;
	double sum2 = 0;

	unsigned int counter = 0;
	for (auto&& [key, val] : std::views::reverse(map)) {
		if (key > 50) {
			++counter;
			auto got = get(key / 2);
			sum1 += double(val.first / got.first);
			sum2 += double(val.second / got.second);
		}
	}
	sum1 /= counter;
	sum2 /= counter;

	std::cout << sum1 << '\n';
	std::cout << sum2 << '\n';
}
