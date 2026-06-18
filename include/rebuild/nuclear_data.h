#ifndef __BRILL_NUCLEAR_DATA_H__
#define __BRILL_NUCLEAR_DATA_H__

#include <iostream>
#include <string>

namespace brill {

struct NuclearData {
	int z;
	int a;
	double mass;
	std::string name;
};

inline std::istream &operator>>(std::istream &is, NuclearData &nuclear) {
	is >> nuclear.z >> nuclear.a >> nuclear.name >> nuclear.mass;
	return is;
}

inline std::ostream &operator<<(std::ostream &os, const NuclearData &nuclear) {
	os << nuclear.z << " " << nuclear.a << " " << nuclear.name << " " << nuclear.mass;
	return os;
}

void SetAssetsPath(const std::string &path);
NuclearData GetNuclear(int z, int a, int q = 0);
double GetMass(int z, int a, int q = 0);
double GetMassInUnit(int z, int a, int q = 0);

} // namespace brill

#endif
