#include <json.hpp>
#include <iostream>
#include <fstream>

// テスト用コード

using json = nlohmann::json;

int main(){
	json molinfo;
	std::ifstream ifs("input.json");
	std::cout << "OK" << std::endl;
	if (ifs.fail()) {
		std::cerr << "File do not exist.\n";
		return 1;
	}
	ifs >> molinfo;

	const int N_atoms = molinfo["resID"].size();

	for (int i=0; i<N_atoms; ++i) {
		std::cout << molinfo["resID"][i] << " ";
	}
	std::cout << std::endl;
	return 0;
}

