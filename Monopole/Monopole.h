#ifndef MONOPOLE_H
#define MONOPOLE_H

#include "crpropa/Vector3.h"
#include "crpropa/ParticleState.h"
#include "crpropa/Referenced.h"
#include "crpropa/AssocVector.h"
#include "crpropa/Variant.h"
#include "crpropa/Candidate.h"
#include "crpropa/Source.h"
#include "crpropa/Module.h"

#include <vector>
#include <map>
#include <sstream>
#include <stdint.h>

namespace crpropa {
static const double gD =  3.291059798e-9*ampere * meter;
/**
 @class MParticleState
 @brief State of the particle: ID, energy, position, direction; modified for monopoles

 The MParticleState defines the state of a cosmic ray, which
 is travelling at any velocity
 The cosmic ray state is defined by particle ID, energy and position and
 direction vector.
 For faster lookup mass and charge of the particle are stored as members.
 */
class MParticleState: public ParticleState {
private:
	double mcharge; ///<particle magnetic charge
	double mass; ///<particle mass to override ParticleState pmass

public:
	/** Constructor for a particle state.
	 @param id			id of the particle following the PDG numbering scheme
	 @param energy		kinetic energy of the particle [in Joules]
	 @param position	vector containing the coordinates of the particle [in meters]
	 @param direction	vector containing the direction of motion of the particle
	 @param pmass		mass of particle [in kg]
	 @param mcharge		magnetic charge of particle [in A*m]
	 @param charge		electric charge of particle [in C]
	 */
	MParticleState(int id = 0, double energy = 0,
			Vector3d position = Vector3d(0, 0, 0),
			Vector3d direction = Vector3d(-1, 0, 0),
			double pmass = 0, double mcharge = 0);

	/** Set particle ID.
	 This follows the PDG numbering scheme:
	  https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
	 @param newId		id to be assigned to the particle 
	 @param pmass		mass to be assigned to the particle if it is a dyon [in kilograms]
	 @param mcharge		magnetic charge to be assigned to the particle if it is a dyon [in A*m]
	 */
	void setId(int newId, double new_pmass = 0, double new_mcharge = 0);
	/** Get particle ID
	 @returns Particle ID (in PDG format).
	 */
	std::string getDescription() const;

	// ======== Helper methods ========

	/** Get magnetic charge of the particle.
	 @returns Magnetic charge of the particle [in A*m]
	 */
	double getMcharge() const;
	/** Set magnetic charge of the particle.
	 @param new_mcharge		new magnetic charge to be assigned to the particle [in A*m]
	 */
	void setMcharge(double new_mcharge); 
	/** Get mass of the particle.
	 @returns Mass of the particle [in kg]
	 */
	double getMass() const;
	/** Set mass of the particle.
	 @param new_pmass		new mass to be assigned to the particle [in kilograms]
	 */
	void setMass(double new_pmass);

	/** Set Lorentz factor and modify the particle's energy accordingly.
	 @param gamma		Lorentz factor
	 */
	void setLorentzFactor(double gamma);
	/** Get Lorentz factor
	 @returns Lorentz factor of particle
	 */
	double getLorentzFactor() const;

	/** Get velocity: direction times the speed of light.
	 @returns Velocity of particle [m/s]
	 */
	Vector3d getVelocity() const;
	/** Get momentum: direction times energy divided by the speed of light 
	 @returns The momentum [kg m/s]
	*/
	Vector3d getMomentum() const;
}; // MParticleState

/**
 @class Candidate Candidate.h include/crpropa/Candidate.h
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray and the simulation itself.
 */
class MCandidate: public Candidate {
public:
	MParticleState Msource; /**< Particle state at the source */
	MParticleState Mcreated; /**< Particle state of parent particle at the time of creation */
	MParticleState Mcurrent; /**< Current particle state */
	MParticleState Mprevious; /**< Particle state at the end of the previous step */
	
private:
	double stepRadiation; /**<Electromagnetic radiation lost at current step */
	
public:
	MCandidate(
		int id = 0,
		double energy = 0,
		Vector3d position = Vector3d(0, 0, 0),
		Vector3d direction = Vector3d(-1, 0, 0),
		double pmass = 0,
		double mcharge = 0,
		double z = 0,
		double weight = 1., 
		std::string tagOrigin = "PRIM"
	);

	/**
	 Creates a candidate, initializing the Candidate::source, Candidate::created,
	 Candidate::previous and Candidate::current state with the argument.
	 */
	MCandidate(const MParticleState &Mstate);
	
	//Helper functions to store and retrieve the radiative losses at each step for debugging and verification
	void setStepRadiation(double radiation);
	double getStepRadiation() const;
	
	//Downconvert pointer of type *Candidate to pointer of type *MCandidate
	static MCandidate* convertToMCandidate(Candidate *candidate) {return static_cast<MCandidate*>(candidate);}  
}; //MCandidate

/**
 @class MSource
 @brief General source of particles compatible with monopoles

 This class is a container for source features.
 The source prepares a new candidate by passing it to all its source features
 to be modified accordingly.
 */
class MSource: public SourceInterface {
	std::vector<ref_ptr<SourceFeature> > features;
public:
	void add(SourceFeature* feature);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
}; //MSource

/**
 @class SourceParticleMonopole
 @brief Magnetic Monopole at the source

 This feature assigns a monopole to the source.
 Monopoles are identified following the PDG numbering scheme:
   https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
 Monopoles need their mass and magnetic charge defined explicitly as they aren't encoded in the id
 */
class SourceParticleMonopole: public SourceFeature {
	int id;
	double pmass;
	double mcharge;
public:
	/** Constructor for a source with a sign
	 @param id		id of the particle following the PDG numbering scheme
	*/
	SourceParticleMonopole(int id, double pmass = 100*gigaelectronvolt/c_squared, double mcharge = 1*gD);
	void prepareCandidate(Candidate& candidate) const override;
	void setDescription();
}; // SourceParticleMonopole

/**
 @class MonopoleSimulationBase
 @brief Abstract class that implements pre/postprocessing of Candidate for monopole simulation modules
 */
class MonopoleSimulationModule: public Module {
public:
 	void process (Candidate* c) const override;
 	virtual void Mprocess (MCandidate *candidate, MParticleState& current) const = 0; 
};//MonopoleSimulationModule

} // namespace crpropa

#endif

