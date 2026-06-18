#include "include/rebuild/nuclear_data.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "include/rebuild/constants.h"

namespace brill {

namespace {

std::string &AssetsPathStorage() {
	static std::string assets_path = "assets";
	return assets_path;
}

std::filesystem::path AssetsFilePath(const std::string &filename) {
	if (filename.empty()) {
		throw std::runtime_error("Assets filename cannot be empty.");
	}
	std::filesystem::path path = std::filesystem::path(AssetsPathStorage()) / filename;
	if (!std::filesystem::is_regular_file(path)) {
		throw std::runtime_error("Assets file not found: " + path.string());
	}
	return path;
}

} // namespace

void SetAssetsPath(const std::string &path) {
	AssetsPathStorage() = path.empty() ? "assets" : path;
}

NuclearData GetNuclear(int z, int a, int q) {
	const std::string amdc_filename = "amdc_ion_2020.txt";
	std::ifstream fin(AssetsFilePath(amdc_filename));
	if (!fin.good()) {
		throw std::runtime_error("amdc_ion_2020.txt not found");
	}
	std::string line;
	std::getline(fin, line);
	NuclearData data{-1, -1, -1.0, ""};
	while (fin.good() && (data.z != z || data.a != a)) {
		fin >> data;
	}
	fin.close();
	if (data.z != z || data.a != a) {
		throw std::runtime_error("Search NuclearData failed");
	}
	if (q > 0) {
		data.mass += q * electron_mass / atomic_mass_unit;
	}
	return data;
}

double GetMass(int z, int a, int q) {
	return GetNuclear(z, a, q).mass * atomic_mass_unit;
}

double GetMassInUnit(int z, int a, int q) {
	return GetNuclear(z, a, q).mass;
}

} // namespace brill
