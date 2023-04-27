#include "Monopole.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleMass.h"

#include <cstdlib>
#include <sstream>

namespace crpropa {
/*
MParticlState
*/

MParticleState::MParticleState(int id, double E, Vector3d pos, Vector3d dir, double pmass, double mcharge)
{
	setId(id, pmass, mcharge);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

void MParticleState::setId(int newId, double new_pmass, double new_mcharge) {
	ParticleState::setId(newId);	//Sets the charge
	if (isDyon(newId)) {
		setMass(new_pmass);
		setMcharge(new_mcharge);
		//setCharge(HepPID::charge(newId) * eplus); //returns positive for 411xyz0, negative for 412xyz0. Need to flip sign if id is negative
	}
}

double MParticleState::getMass() const {
	return mass;
}

void MParticleState::setMass(double new_pmass) {
	mass = new_pmass * kilogram;
}

double MParticleState::getMcharge() const {
	return mcharge;
}

void MParticleState::setMcharge(double new_mcharge) {
	mcharge = fabs(new_mcharge * ampere * meter);
	if (getId() < 0)
		mcharge *= -1; //anti-monopole
}

double MParticleState::getLorentzFactor() const {
	return 1 + getEnergy() / (getMass() * c_squared);
}

void MParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	setEnergy((lf - 1) * getMass() * c_squared);
}

Vector3d MParticleState::getVelocity() const {
	return getDirection() * c_light * sqrt(1 - 1 / (1 + getEnergy() / getMass() / c_squared) / (1 + getEnergy() / getMass() / c_squared));
}

Vector3d MParticleState::getMomentum() const {
	return getDirection() * sqrt( (getEnergy() + getMass() * c_squared) * (getEnergy() + getMass() * c_squared) - (getMass() * c_squared) * (getMass() * c_squared) ) / c_light;
}

std::string MParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << getId() << ", ";
	ss << "E = " << getEnergy() / EeV << " EeV, ";
	ss << "x = " << getPosition() / Mpc << " Mpc, ";
	ss << "p = " << getDirection() << " (direction), ";
	ss << "q = " << getCharge() << " C, ";
	ss << "m = " << mass / (gigaelectronvolt / c_squared) << " GeV/c^2, ";
	ss << "gD = " << mcharge / gD << " gD";
	return ss.str();
}

// ----------------------------------------------------------------------------
/*
MCandidate
*/

MCandidate::MCandidate(int id, double E, Vector3d pos, Vector3d dir, double pmass, double mcharge, double z, double weight, std::string tagOrigin) {
	Candidate::setRedshift(z);
	Candidate::setTrajectoryLength(0);
	Candidate::setWeight(weight);
	Candidate::setCurrentStep(0);
	Candidate::setNextStep(0);
	Candidate:setActive(true);
	parent = 0;
	Candidate::setTagOrigin(tagOrigin);
	MParticleState Mstate(id, E, pos, dir, pmass, mcharge);
	source = Mstate;
	created = Mstate;
	previous = Mstate;
	current = Mstate;
	Msource = Mstate;
	Mcreated = Mstate;
	Mprevious = Mstate;
	Mcurrent = Mstate;

#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber()++);}
#elif defined(__GNUC__)
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber() + 1);}
#else
		#pragma omp critical
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber()++);}
#endif
}


MCandidate::MCandidate(const MParticleState &Mstate) :
		Msource(Mstate), Mcreated(Mstate), Mcurrent(Mstate), Mprevious(Mstate){
	
	source = Mstate;
	created = Mstate;
	previous = Mstate;
	current = Mstate;
		
	Candidate::setRedshift(0);
	Candidate::setTrajectoryLength(0);
	Candidate::setCurrentStep(0);
	Candidate::setNextStep(0);
	Candidate:setActive(true);
	parent = 0;
	Candidate::setTagOrigin("PRIM");

#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber()++);}
#elif defined(__GNUC__)
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber() + 1);}
#else
		#pragma omp critical
		{Candidate::setSerialNumber(Candidate::getNextSerialNumber()++);}
#endif
}

void MCandidate::setStepRadiation(double radiation)  {
	stepRadiation = radiation;
}

double MCandidate::getStepRadiation() const {
	return stepRadiation;
}

// Source ---------------------------------------------------------------------
void MSource::add(SourceFeature* property) {
	features.push_back(property);
}

ref_ptr<Candidate> MSource::getCandidate() const {
	ref_ptr<Candidate> candidate = static_cast<Candidate*>(new MCandidate());
	for (int i = 0; i < features.size(); i++)
		(*features[i]).prepareCandidate(*candidate);
	return candidate;
}

std::string MSource::getDescription() const {
	std::stringstream ss;
	ss << "Cosmic ray source\n";
	for (int i = 0; i < features.size(); i++)
		ss << "    " << features[i]->getDescription();
	return ss.str();
}


// SourceParticleMonopole ----------------------------------------------------------------------------
SourceParticleMonopole::SourceParticleMonopole(int id, double pmass, double mcharge) :
		id(id), pmass(pmass), mcharge(mcharge) {
	setDescription();
}

void SourceParticleMonopole::prepareCandidate(Candidate& candidate) const {
	MCandidate& Mcandidate = *MCandidate::convertToMCandidate(&candidate);
	ParticleState &source = Mcandidate.source;
	MParticleState &Msource = Mcandidate.Msource;
	
	Msource.setPosition(source.getPosition());
	Msource.setDirection(source.getDirection());
	Msource.setEnergy(source.getEnergy());
	
	Msource.setId(id);
	Msource.setMass(pmass);
	Msource.setMcharge(mcharge);
	
	Mcandidate.created = Msource;
	Mcandidate.current = Msource;
	Mcandidate.previous = Msource;
	Mcandidate.Mcreated = Msource;
	Mcandidate.Mcurrent = Msource;
	Mcandidate.Mprevious = Msource;
}

void SourceParticleMonopole::setDescription() {
	std::stringstream ss;
	ss << "SourceParticleMonopole: " << id << ", ";
	ss << "Mass = " << pmass/(gigaelectronvolt/c_squared)<< " GeV/c^2, ";
	if (id>0) {
		ss << "Magnetic Charge = " << mcharge*ampere*meter/gD << " gD\n";
	}
	else {
		ss << "Magnetic Charge = -" << mcharge*ampere*meter/gD << " gD\n";
	}
	description = ss.str();
}

// ----------------------------------------------------------------------------
/*
MonopoleSimulationModule
*/
void MonopoleSimulationModule::process(Candidate* c) const {
	//Safetycheck before Downconvert
	//if (!isDyon(candidate_base->current.getId())) {
	//	return;
	//}
	
	//Downconvert
	MCandidate *candidate = MCandidate::convertToMCandidate(c);
	MParticleState &current = candidate->Mcurrent;
	
	//Copy current to Mcurrent
	current.setEnergy(candidate -> current.getEnergy());
	current.setDirection(candidate ->current.getDirection());
	
	// save the new previous particle state
	candidate->previous = current;
	candidate->Mprevious = current;
	
	//Call the virtual function for normal processing
	Mprocess(candidate, current);
	
	//Copy Mcurrent to current
	candidate -> current = current;
}

} // namespace crpropa
