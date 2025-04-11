#include <iostream>
#include <vector>


int main() {
	std::vector<int> vec1;
	std::vector<int> vec2;
	std::vector<std::vector<int>> vec3;

	vec2.push_back(1);
	vec2.push_back(2);

	for (int v : vec2) {
		std::cout << v << "\n";
	}

	vec3.push_back(vec2);
	
	for (int v : vec3) {
		std::cout << v << "\n";
	}
}
