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

void MonopoleRadiation::Mprocess(MCandidate *candidate, MParticleState &current) const {
	double mcharge = fabs(current.getMcharge());
	if (mcharge == 0)
		return; // only charged particles

	// get magnetic field
	Vector3d pos = current.getPosition();
	double z = candidate->getRedshift();
	Vector3d B = getFieldAtPosition(pos, z);

	//Get helper values
	double lf = current.getLorentzFactor();
	double step = candidate->getCurrentStep(); // step size in local frame
	Vector3d v = current.getVelocity();
	Vector3d F = mcharge*B + current.getCharge()*v.cross(B); //Force on particle
	double m = current.getMass();
	
	//Vector3d a = F_mag / m * (dir / cos * (1/ pow(lf, 3)  -1 / lf) + Fdir / lf);
	//Vector3d a = 1 / m / lf * (F - F.dot(v)*v / c_squared);	
	//if (F.isParallelTo(v, 0.0001)) a = F / m / pow(lf, 3); //Limit floating point errors in linear accelerator
	//else a = 1 / m / lf * (F - F.dot(v)*v / c_squared);
	
	// calculate energy loss
	//double P = mu0 / 6 / M_PI * pow(mcharge/ m, 2) * pow(1/c_light, 3) * (F.dot(F) * pow(lf, 2) - pow(p.dot(F) / m, 2)/c_squared); // Jackson p. 666 (14.26)
	//double P = mu0 * pow(mcharge, 2) * pow (lf, 6) / 6 / M_PI / c_light * (a.dot(a) / c_squared - v.cross(a).dot(v.cross(a)) / c_squared / c_squared); 
	double A = mu0 / 6 / M_PI * pow(mcharge / m, 2) * pow(1/c_light, 3);
	double P = A * (F.getR2() + pow(lf, 2) * v.cross(F).getR2() / c_squared);
	double dE = P * step / c_light;
	candidate->setStepRadiation(dE);

	// apply energy loss and limit next step
	double E = current.getEnergy();
	current.setEnergy(E - dE);
	candidate->limitNextStep(limit * E * c_light / P);

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

Vector3d MonopoleRadiation::getFieldAtPosition(Vector3d pos, double z) const {
	Vector3d B(0, 0, 0);
	try {
		// check if field is valid and use the field vector at the
		// position pos with the redshift z
		if (field.valid())
			B = field->getField(pos, z);
		else B = Brms;
	} catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "MonopoleRadiation: Exception in MonopoleRadiation::getFieldAtPosition.\n"
				<< e.what();
	}	
	return B;
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
