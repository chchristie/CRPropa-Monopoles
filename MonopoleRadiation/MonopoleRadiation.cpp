#include "MonopoleRadiation.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

MonopoleRadiation::MonopoleRadiation(ref_ptr<MagneticField> field, bool havePhotons, double thinning, int nSamples, double limit) {
	setField(field);
	setBrms(0);
	//initSpectrum(); //Only necessary for photon secondaries
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
}

MonopoleRadiation::MonopoleRadiation(double Brms, bool havePhotons, double thinning, int nSamples, double limit) {
	setBrms(Brms);
	//initSpectrum(); //Only necessary for photon secondaries
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
}

void MonopoleRadiation::setField(ref_ptr<MagneticField> f) {
	this->field = f;
}

ref_ptr<MagneticField> MonopoleRadiation::getField() {
	return field;
}

void MonopoleRadiation::setBrms(double Brms) {
	this->Brms = Brms;
}

double MonopoleRadiation::getBrms() {
	return Brms;
}

void MonopoleRadiation::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

bool MonopoleRadiation::getHavePhotons() {
	return havePhotons;
}

void MonopoleRadiation::setThinning(double thinning) {
	this->thinning = thinning;
}

double MonopoleRadiation::getThinning() {
	return thinning;
}

void MonopoleRadiation::setLimit(double limit) {
	this->limit = limit;
}

double MonopoleRadiation::getLimit() {
	return limit;
}

void MonopoleRadiation::setMaximumSamples(int nmax) {
	maximumSamples = nmax;
}

int MonopoleRadiation::getMaximumSamples() {
	return maximumSamples;
}

void MonopoleRadiation::setSecondaryThreshold(double threshold) {
	secondaryThreshold = threshold;
}

double MonopoleRadiation::getSecondaryThreshold() const {
	return secondaryThreshold;
}

void MonopoleRadiation::initSpectrum() {
	std::string filename = getDataPath("MonopoleRadiation/spectrum.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("MonopoleRadiation: could not open file " + filename);

	// clear previously loaded interaction rates
	tabx.clear();
	tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabx.push_back(pow(10, a));
				tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void MonopoleRadiation::process(Candidate *candidate) const{
	double mcharge = fabs(candidate->current.getMcharge());
	if (mcharge == 0)
		return; // only charged particles

	// calculate gyroradius, evaluated at the current position
	double z = candidate->getRedshift();
	double B;
	if (field.valid()) {
		Vector3d Bvec = field->getField(candidate->current.getPosition(), z);
		B = Bvec.cross(candidate->current.getDirection()).getR();
	} else {
		B = sqrt(2. / 3) * Brms; // average perpendicular field component
	}
	B *= pow(1 + z, 2); // cosmological scaling

	//Get helper values
	double lf = candidate->current.getLorentzFactor();
	double step = candidate->getCurrentStep() / (1 + z); // step size in local frame
	Vector3d beta = candidate->current.getVelocity() / c_light;
	Vector3d dbeta = (candidate->current.getVelocity() - candidate->previous.getVelocity()) / step;	//d(beta)/dt = vf-vi/c*t = vf-vi/step
	
	// calculate energy loss
	double dEdx = mu0 / 6 / M_PI * pow(lf, 6) * pow(mcharge / c_light, 2) * (dbeta.dot(dbeta) - beta.cross(dbeta).dot(beta.cross(dbeta))); // Jackson p. 666 (14.26)
	double dE = step * dEdx;
	candidate->setStepRadiation(dE);

	// apply energy loss and limit next step
	double E = candidate->current.getEnergy();
	candidate->current.setEnergy(E - dE);
	candidate->limitNextStep(limit * E / dEdx);

	// optionally add secondary photons
	if (not(havePhotons))
		return;

	// Everything below this has not been altered for magnetically charged particles
	
	/*double Ecrit = 3. / 4 * h_planck / M_PI * c_light * pow(lf, 3) / Rg;
	if (14 * Ecrit < secondaryThreshold)
		return;

	// draw photons up to the total energy loss
	// if maximumSamples is reached before that, compensate the total energy afterwards
	Random &random = Random::instance();
	double dE0 = dE;
	std::vector<double> energies;
	int counter = 0;
	while (dE > 0) {
		// draw random value between 0 and maximum of corresponding cdf
		// choose bin of s where cdf(x) = cdf_rand -> x_rand
		size_t i = random.randBin(tabCDF); // draw random bin (upper bin boundary returned)
		double binWidth = (tabx[i] - tabx[i-1]);
		double x = tabx[i-1] + random.rand() * binWidth; // draw random x uniformly distributed in bin
		double Ephoton = x * Ecrit;

		// if the remaining energy is not sufficient check for random accepting
		if (Ephoton > dE) {
			if (random.rand() > (dE / Ephoton))
				break; // not accepted
		}

		// only activate the "per-step" sampling if maximumSamples is explicitly set.
		if (maximumSamples > 0) {
			if (counter >= maximumSamples) 
				break;			
		}

		// store energies in array
		energies.push_back(Ephoton);

		// energy loss
		dE -= Ephoton;

		// counter for sampling break condition;
		counter++;
	}

	// while loop before gave total energy which is just a fraction of the required
	double w1 = 1;
	if (maximumSamples > 0 && dE > 0)
		w1 = 1. / (1. - dE / dE0); 

	// loop over sampled photons and attribute weights accordingly
	for (int i = 0; i < energies.size(); i++) {
		double Ephoton = energies[i];
		double f = Ephoton / (E - dE0);
		double w = w1 / pow(f, thinning);

		// thinning procedure: accepts only a few random secondaries
		if (random.rand() < pow(f, thinning)) {
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			if (Ephoton > secondaryThreshold) // create only photons with energies above threshold
				candidate->addSecondary(22, Ephoton, pos, w, interactionTag);
		}
	} */
}

std::string MonopoleRadiation::getDescription() const {
	std::stringstream s;
	s << "Monopole radiation";
	if (field.valid())
		s << " for specified magnetic field";
	else
		s << " for Brms = " << Brms / nG << " nG";
	if (havePhotons)
		s << ", photons E > " << secondaryThreshold / eV << " eV";
	else
		s << ", no photons";
	if (maximumSamples > 0)
		s << "maximum number of photon samples: " << maximumSamples;
	if (thinning > 0)
		s << "thinning parameter: " << thinning; 
	return s.str();
}

void MonopoleRadiation::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string MonopoleRadiation::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
