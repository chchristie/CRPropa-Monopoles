#include "crpropa/magneticField/MagneticField.h"
#include "kiss/logger.h"
#include "Monopole.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class MonopoleRadiation
 @brief Radiation of magnetically charged particles in magnetic fields.

 This module simulates the continuous energy loss of magnetically charged particles in magnetic fields, c.f. Jackson.
 The magnetic field is specified either by a MagneticField or by a RMS field strength value.
 The module limits the next step size to ensure a fractional energy loss dE/E < limit (default = 0.1).
 Optionally, photons above a threshold (default E > 10^7 eV) are created as secondary particles.
 Note that the large number of secondary photons per propagation can cause memory problems.
 To mitigate this, use thinning. However, this still does not solve the problem completely.
 For this reason, a break-condition stops tracking secondary photons and reweights the current ones. 
 */
class MonopoleRadiation: public MonopoleSimulationModule {
private:
	ref_ptr<MagneticField> field; ///< MagneticField instance
	double Brms; ///< Brms value in case no MagneticField is specified
	double limit; ///< fraction of energy loss length to limit the next step
	double thinning; ///< thinning parameter for weighted-sampling (maximum 1, minimum 0)
	bool havePhotons; ///< flag for production of secondary photons
	int maximumSamples; ///< maximum number of samples of photons (break condition; defaults to 100; 0 or <0 means no sampling)
	double secondaryThreshold; ///< threshold energy for secondary photons
	std::vector<double> tabx; ///< tabulated fraction E_photon/E_critical from 10^-6 to 10^2 in 801 log-spaced steps
	std::vector<double> tabCDF; ///< tabulated CDF of radiation spectrum
	std::string interactionTag = "SYN";

public:
	/** Constructor
	 @param field			magnetic field object
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param nSamples		number of photons to be sampled and added as candidates
	 @param limit			step size limit as fraction of mean free path
	 */
	MonopoleRadiation(ref_ptr<MagneticField> field, bool havePhotons = false, double thinning = 0, int nSamples = 0, double limit = 0.1);
	/** Constructor
	 @param field			RMS of the magnetic field (if magnetic-field object not provided)
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param nSamples		number of photons to be sampled and added as candidates
	 @param limit			step size limit as fraction of mean free path
	 */
	MonopoleRadiation(double Brms = 0, bool havePhotons = false, double thinning = 0, int nSamples = 0, double limit = 0.1);
	
	void setField(ref_ptr<MagneticField> field);
	void setBrms(double Brms);	
	void setHavePhotons(bool havePhotons);
	void setThinning(double thinning);
	void setLimit(double limit);
	/** Set the maximum number of photons that will be allowed to be added as candidates. 
	 This choice depends on the problem at hand. It must be such that all relevant physics is captured with the sample. Weights are added accordingly and the column 'weight' must be added to output.
	 @param nmax	maximum number of photons to be sampled
	 */
	void setMaximumSamples(int nmax);
	/** Photons above the secondary energy threshold are added as candidates.
	 This may lead to a quick increase in memory.
	 @param threshold	energy threshold above which photons will be added [in Joules]
	 */
	void setSecondaryThreshold(double threshold);	
	void setInteractionTag(std::string tag);
	ref_ptr<MagneticField> getField();

	double getBrms();
	bool getHavePhotons();
	double getThinning();
	double getLimit();
	int getMaximumSamples();
	double getSecondaryThreshold() const;
	std::string getInteractionTag() const;

	void initSpectrum();
	/** Propagates the particle. Is called once per iteration.
	 * @param candidate	 The Candidate is a passive object, that holds the information about the state of the cosmic ray and the simulation itself. 
	   @param current	Current is a reference to the Mcurrent member of candidate*/
	void Mprocess(MCandidate *candidate, MParticleState& current) const override;
	std::string getDescription() const;
	
	/** Get magnetic field vector at current candidate position
	 * @param pos   current position of the candidate
	 * @param z	 current redshift is needed to calculate the magnetic field
	 * @return	  magnetic field vector at the position pos 
	 */
	Vector3d getFieldAtPosition(Vector3d pos, double z) const;
};
/** @}*/

} // namespace crpropa

