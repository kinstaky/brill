#include "include/rebuild/particle.h"

#include <cmath>

#include "include/rebuild/constants.h"
#include "include/rebuild/nuclear_data.h"

namespace brill {

Particle::Particle(
	int z,
	int a,
	double kinetic,
	const ROOT::Math::XYZVector &direction,
	int q,
	double excitation
) :
	z_(z),
	a_(a),
	q_(q),
	excitation_(excitation),
	kinetic_(kinetic),
	direction_(direction.Unit())
{
	mass_ = GetMass(z_, a_);
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
}

Particle& Particle::SetMass0(double mass) {
	mass_ = mass;
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}

Particle& Particle::LostKineticEnergy(catima::Material material) {
	material.thickness(material.thickness() / cos(direction_.Theta()));
	double mass_in_u = (mass_ + excitation_) / atomic_mass_unit;
	catima::Projectile projectile(mass_in_u, z_, 0, kinetic_ / mass_in_u);
	kinetic_ = catima::energy_out(projectile, material) * mass_in_u;
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}

Particle& Particle::AddKineticEnergy(double energy) {
	kinetic_ += energy;
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}

Particle& Particle::SetKineticEnergy(double energy) {
	kinetic_ = energy;
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}

Particle& Particle::SetExcitationEnergy(double excitation) {
	excitation_ = excitation;
	energy_ = mass_ + excitation_ + kinetic_;
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}

Particle& Particle::SetDirection(const ROOT::Math::XYZVector &direction) {
	direction_ = direction.Unit();
	return *this;
}

Particle& Particle::SetMomentum(const ROOT::Math::XYZVector &momentum) {
	direction_ = momentum.Unit();
	momentum_ = momentum.R();
	kinetic_ = KineticFromMomentum(mass_ + excitation_, momentum_);
	energy_ = mass_ + excitation_ + kinetic_;
	return *this;
}

Particle& Particle::SetMomentum(double momentum) {
	momentum_ = momentum;
	kinetic_ = KineticFromMomentum(mass_ + excitation_, momentum_);
	energy_ = mass_ + excitation_ + kinetic_;
	return *this;
}

Particle Particle::operator+(const Particle &other) const {
	Particle result(z_ + other.z_, a_ + other.a_);
	ROOT::Math::XYZVector momentum_vector = MomentumVector() + other.MomentumVector();
	result.direction_ = momentum_vector.Unit();
	result.momentum_ = momentum_vector.R();
	result.energy_ = Energy() + other.Energy();
	double mass = sqrt(pow(result.energy_, 2.0) - pow(result.momentum_, 2.0));
	result.excitation_ = mass - result.mass_;
	result.kinetic_ = result.energy_ - mass;
	return result;
}

Particle Particle::operator-(const Particle &other) const {
	Particle result(z_ - other.z_, a_ - other.a_);
	ROOT::Math::XYZVector momentum_vector = MomentumVector() - other.MomentumVector();
	result.direction_ = momentum_vector.Unit();
	result.momentum_ = momentum_vector.R();
	result.energy_ = Energy() - other.Energy();
	double mass = sqrt(pow(result.energy_, 2.0) - pow(result.momentum_, 2.0));
	result.excitation_ = mass - result.mass_;
	result.kinetic_ = result.energy_ - mass;
	return result;
}

catima::Material SolidMaterial(
	std::initializer_list<std::array<double, 3>> compound,
	double thickness
) {
	catima::Material material;
	for (auto &element : compound) {
		material.add_element(
			GetMassInUnit(element[1], element[0]),
			element[1],
			element[2]
		);
	}
	material.thickness(thickness * 1e-3);
	return material;
}

catima::Material GasMaterial(
	std::initializer_list<std::array<double, 3>> compound,
	double pressure,
	double length
) {
	catima::Material material;
	int molar_mass = 0;
	for (auto &element : compound) {
		material.add_element(
			GetMassInUnit(element[1], element[0]),
			element[1],
			element[2]
		);
		molar_mass += int(element[0] * element[2]);
	}
	pressure /= 0.0075006;
	material.density(molar_mass * pressure / (gas_constant * 293.15) * 1e-6);
	material.thickness_cm(length * 100.0);
	return material;
}

void Scatter(
	const Particle &beam,
	const Particle &target,
	const ROOT::Math::XYZVector &scatter_direction,
	Particle &fragment0,
	Particle &fragment1
) {
	double polar = scatter_direction.Theta();
	double azimuth = scatter_direction.Phi();
	double beta_mass_center = beam.Momentum() / (beam.Energy() + target.Mass());
	double gamma_mass_center = 1.0 / sqrt(1.0 - beta_mass_center * beta_mass_center);
	double reaction_energy = sqrt(
		(beam.Energy() + target.Mass() + beam.Momentum())
		* (beam.Energy() + target.Mass() - beam.Momentum())
	);
	double exit_momentum_center =
		sqrt(
			(reaction_energy - fragment0.Mass() - fragment1.Mass())
			* (reaction_energy - fragment0.Mass() + fragment1.Mass())
			* (reaction_energy + fragment0.Mass() - fragment1.Mass())
			* (reaction_energy + fragment0.Mass() + fragment1.Mass())
		) / (2.0 * reaction_energy);
	double exit_momentum_center_parallel = exit_momentum_center * cos(polar);
	double exit_momentum_center_vertical = exit_momentum_center * sin(polar);
	double exit_energy_center = sqrt(
		exit_momentum_center * exit_momentum_center
		+ fragment0.Mass() * fragment0.Mass()
	);
	double exit_energy = gamma_mass_center * exit_energy_center
		+ gamma_mass_center * beta_mass_center * exit_momentum_center_parallel;
	double exit_momentum_parallel = gamma_mass_center * exit_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * exit_energy_center;
	double exit_momentum_vertical = exit_momentum_center_vertical;
	double exit_angle = fabs(atan(exit_momentum_vertical / exit_momentum_parallel));
	exit_angle = exit_momentum_parallel > 0 ? exit_angle : pi - exit_angle;
	fragment0.SetKineticEnergy(exit_energy - fragment0.Mass());
	fragment0.SetDirection(Rotate(
		beam.Direction(),
		ROOT::Math::XYZVector(
			sin(exit_angle) * cos(azimuth),
			sin(exit_angle) * sin(azimuth),
			cos(exit_angle)
		)
	));

	double recoil_momentum_center = exit_momentum_center;
	double recoil_momentum_center_parallel = -recoil_momentum_center * cos(polar);
	double recoil_momentum_center_vertical = -recoil_momentum_center * sin(polar);
	double recoil_energy_center = sqrt(
		recoil_momentum_center * recoil_momentum_center
		+ fragment1.Mass() * fragment1.Mass()
	);
	double recoil_energy = gamma_mass_center * recoil_energy_center
		+ gamma_mass_center * beta_mass_center * recoil_momentum_center_parallel;
	double recoil_momentum_parallel =
		gamma_mass_center * recoil_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * recoil_energy_center;
	double recoil_momentum_vertical = recoil_momentum_center_vertical;
	double recoil_angle = fabs(atan(recoil_momentum_vertical / recoil_momentum_parallel));
	recoil_angle = recoil_momentum_parallel > 0 ? recoil_angle : pi - recoil_angle;
	fragment1.SetKineticEnergy(recoil_energy - fragment1.Mass());
	fragment1.SetDirection(Rotate(
		beam.Direction(),
		ROOT::Math::XYZVector(
			sin(recoil_angle) * cos(azimuth - pi),
			sin(recoil_angle) * sin(azimuth - pi),
			cos(recoil_angle)
		)
	));
}

void Breakup(
	const Particle &parent,
	const ROOT::Math::XYZVector &breakup_direction,
	Particle &fragment0,
	Particle &fragment1
) {
	double polar = breakup_direction.Theta();
	double azimuth = breakup_direction.Phi();
	double beta_center = parent.Momentum() / parent.Energy();
	double gamma_center = 1.0 / sqrt(1.0 - pow(beta_center, 2.0));
	double parent_energy_center = parent.Mass();
	double fragment_momentum =
		sqrt(
			(parent_energy_center - fragment0.Mass() - fragment1.Mass())
			* (parent_energy_center - fragment0.Mass() + fragment1.Mass())
			* (parent_energy_center + fragment0.Mass() - fragment1.Mass())
			* (parent_energy_center + fragment0.Mass() + fragment1.Mass())
		) / (2.0 * parent.Mass());

	double fragment0_momentum_center_parallel = fragment_momentum * cos(polar);
	double fragment0_momentum_center_vertical = fragment_momentum * sin(polar);
	double fragment0_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment0.Mass(), 2.0)
	);
	double fragment0_energy = gamma_center * fragment0_energy_center
		+ gamma_center * beta_center * fragment0_momentum_center_parallel;
	double fragment0_momentum_parallel =
		gamma_center * fragment0_momentum_center_parallel
		+ gamma_center * beta_center * fragment0_energy_center;
	double fragment0_momentum_vertical = fragment0_momentum_center_vertical;
	double fragment0_angle = fabs(atan(fragment0_momentum_vertical / fragment0_momentum_parallel));
	fragment0_angle = fragment0_momentum_parallel > 0 ? fragment0_angle : pi - fragment0_angle;
	fragment0.SetKineticEnergy(fragment0_energy - fragment0.Mass());
	fragment0.SetDirection(Rotate(
		parent.Direction(),
		ROOT::Math::XYZVector(
			sin(fragment0_angle) * cos(azimuth),
			sin(fragment0_angle) * sin(azimuth),
			cos(fragment0_angle)
		)
	));

	double fragment1_momentum_center_parallel = -fragment_momentum * cos(polar);
	double fragment1_momentum_center_vertical = -fragment_momentum * sin(polar);
	double fragment1_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment1.Mass(), 2.0)
	);
	double fragment1_energy = gamma_center * fragment1_energy_center
		+ gamma_center * beta_center * fragment1_momentum_center_parallel;
	double fragment1_momentum_parallel =
		gamma_center * fragment1_momentum_center_parallel
		+ gamma_center * beta_center * fragment1_energy_center;
	double fragment1_momentum_vertical = fragment1_momentum_center_vertical;
	double fragment1_angle = fabs(atan(fragment1_momentum_vertical / fragment1_momentum_parallel));
	fragment1_angle = fragment1_momentum_parallel > 0 ? fragment1_angle : pi - fragment1_angle;
	fragment1.SetKineticEnergy(fragment1_energy - fragment1.Mass());
	fragment1.SetDirection(Rotate(
		parent.Direction(),
		ROOT::Math::XYZVector(
			sin(fragment1_angle) * cos(azimuth - pi),
			sin(fragment1_angle) * sin(azimuth - pi),
			cos(fragment1_angle)
		)
	));
}

ROOT::Math::XYZVector Rotate(
	const ROOT::Math::XYZVector &parent,
	const ROOT::Math::XYZVector &child
) {
	double x =
		+sin(child.Theta()) * cos(child.Phi()) * cos(parent.Theta()) * cos(parent.Phi())
		-sin(child.Theta()) * sin(child.Phi()) * sin(parent.Phi())
		+cos(child.Theta()) * sin(parent.Theta()) * cos(parent.Phi());
	double y =
		+sin(child.Theta()) * cos(child.Phi()) * cos(parent.Theta()) * sin(parent.Phi())
		+sin(child.Theta()) * sin(child.Phi()) * cos(parent.Phi())
		+cos(child.Theta()) * sin(parent.Theta()) * sin(parent.Phi());
	double z =
		-sin(child.Theta()) * cos(child.Phi()) * sin(parent.Theta())
		+cos(child.Theta()) * cos(parent.Theta());

	double theta = fabs(atan(sqrt(pow(x, 2.0) + pow(y, 2.0)) / z));
	theta = z > 0 ? theta : pi - theta;
	double phi = atan(y / x);
	if (x < 0) {
		phi = y > 0 ? phi + pi : phi - pi;
	}

	return ROOT::Math::XYZVector(
		sin(theta) * cos(phi),
		sin(theta) * sin(phi),
		cos(theta)
	);
}

ROOT::Math::XYZVector DirectionVector(double theta, double phi) {
	return ROOT::Math::XYZVector(
		sin(theta) * cos(phi),
		sin(theta) * sin(phi),
		cos(theta)
	);
}

} // namespace brill
