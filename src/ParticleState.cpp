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
		charge = HepPID::charge(id) * eplus;
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
	mcharge = new_mcharge * ampere * meter;
	if (id < 0)
		mcharge *= -1; //anti-monopole
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	energy = lf * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction;
	return ss.str();
}

} // namespace crpropa
