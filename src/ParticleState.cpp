#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

#include "HepPID/ParticleIDMethods.hh"

#include <cstdlib>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir, double pmass, double mcharge): id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.), mcharge(0.)
{
	setId(id, pmass, mcharge);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy); // prevent negative energies
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId, double new_pmass, double new_mcharge) {
	id = newId;
	if (isNucleus(id)) {
		pmass = nuclearMass(id);
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} 
	else if (HepPID::isDyon(id)) {
		setMass(new_pmass);
		setMcharge(new_mcharge);
		charge = HepPID::charge(id) * eplus; //returns positive for 411xyz0, negative for 412xyz0. Need to flip sign if id is negative
		if (id < 0)
			charge *= -1;
	}
	else {
		if (abs(id) == 11)
			pmass = mass_electron;
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

void ParticleState::setMass(double new_pmass) {
	pmass = new_pmass * kilogram;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getMcharge() const {
	return mcharge;
}

void ParticleState::setMcharge(double new_mcharge) {
	mcharge = fabs(new_mcharge * ampere * meter);
	if (id < 0)
		mcharge *= -1; //anti-monopole
}

double ParticleState::getLorentzFactor() const {
	return 1 + energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	energy = (lf - 1) * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light * sqrt(1 - 1 / (1 + energy / pmass / c_squared) / (1 + energy / pmass / c_squared));
}

Vector3d ParticleState::getMomentum() const {
	return direction * sqrt( (energy + pmass * c_squared) * (energy + pmass * c_squared) - (pmass * c_squared) * (pmass * c_squared) ) / c_light;
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction << " (direction), ";
	ss << "q = " << charge << " C, ";
	ss << "m = " << pmass / (gigaelectronvolt / c_squared) << " GeV/c^2, ";
	ss << "gD = " << mcharge / gD << " gD";
	return ss.str();
}

} // namespace crpropa
