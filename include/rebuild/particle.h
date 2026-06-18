#ifndef __BRILL_REBUILD_PARTICLE_H__
#define __BRILL_REBUILD_PARTICLE_H__

#include <array>
#include <cmath>
#include <initializer_list>
#include <string>
#include <vector>

#include <Math/Vector3D.h>
#include <catima/catima.h>

namespace brill {

class Particle {
public:
	Particle(
		int z,
		int a,
		double kinetic = 0.0,
		const ROOT::Math::XYZVector &direction = ROOT::Math::XYZVector(0, 0, 1),
		int q = 0,
		double excitation = 0.0
	);

	Particle(const Particle &other) = default;
	~Particle() = default;

	Particle& SetMass0(double mass);
	inline double Mass0() const { return mass_; }
	inline double Mass() const { return mass_ + excitation_; }
	Particle& LostKineticEnergy(catima::Material material);
	Particle& AddKineticEnergy(double energy);
	Particle& SetKineticEnergy(double energy);
	inline double KineticEnergy() const { return kinetic_; }
	inline double Energy() const { return energy_; }
	Particle& SetExcitationEnergy(double energy);
	inline double ExcitationEnergy() const { return excitation_; }
	Particle& SetDirection(const ROOT::Math::XYZVector &direction);
	inline ROOT::Math::XYZVector Direction() const { return direction_; }
	inline double Polar() const { return direction_.Theta(); }
	inline double Azimuthal() const { return direction_.Phi(); }
	Particle& SetMomentum(const ROOT::Math::XYZVector &momentum);
	Particle& SetMomentum(double momentum);
	inline double Momentum() const { return momentum_; }
	inline ROOT::Math::XYZVector MomentumVector() const { return momentum_ * direction_; }
	Particle operator+(const Particle &other) const;
	Particle operator-(const Particle &other) const;

private:
	int z_;
	int a_;
	int q_;
	double mass_;
	double excitation_;
	double energy_;
	double kinetic_;
	double momentum_;
	ROOT::Math::XYZVector direction_;
};

inline double MomentumFromKinetic(double mass, double kinetic) {
	return sqrt((2.0 * mass + kinetic) * kinetic);
}

inline double KineticFromMomentum(double mass, double momentum) {
	return sqrt(pow(momentum, 2.0) + pow(mass, 2.0)) - mass;
}

catima::Material GasMaterial(
	std::initializer_list<std::array<double, 3>> compound,
	double pressure,
	double length
);

catima::Material SolidMaterial(
	std::initializer_list<std::array<double, 3>> compound,
	double thickness
);

using ParticleList = std::vector<Particle>;

void Scatter(
	const Particle &beam,
	const Particle &target,
	const ROOT::Math::XYZVector &scatter_direction,
	Particle &fragment0,
	Particle &fragment1
);

void Breakup(
	const Particle &parent,
	const ROOT::Math::XYZVector &breakup_direction,
	Particle &fragment0,
	Particle &fragment1
);

ROOT::Math::XYZVector Rotate(
	const ROOT::Math::XYZVector &parent,
	const ROOT::Math::XYZVector &child
);

ROOT::Math::XYZVector DirectionVector(
	double theta,
	double phi
);

} // namespace brill

#endif
