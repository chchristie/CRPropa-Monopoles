#include "crpropa/module/ElasticScattering.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {

const double ElasticScattering::lgmin = 6.;  // minimum log10(Lorentz-factor)
const double ElasticScattering::lgmax = 14.; // maximum log10(Lorentz-factor)
const size_t ElasticScattering::nlg = 201;   // number of Lorentz-factor steps
const double ElasticScattering::epsmin = log10(2 * eV) + 3;    // log10 minimum photon background energy in nucleus rest frame for elastic scattering
const double ElasticScattering::epsmax = log10(2 * eV) + 8.12; // log10 maximum photon background energy in nucleus rest frame for elastic scattering
const size_t ElasticScattering::neps = 513; // number of photon background energies in nucleus rest frame

ElasticScattering::ElasticScattering(ref_ptr<PhotonField> f) {
	setPhotonField(f);
}

void ElasticScattering::setPhotonField(ref_ptr<PhotonField> photonField) {
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("ElasticScattering: " + fname);
	initRate(getDataPath("ElasticScattering/rate_" + fname.substr(0,3) + ".txt"));
	initCDF(getDataPath("ElasticScattering/cdf_" + fname.substr(0,3) + ".txt"));
}

void ElasticScattering::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("ElasticScattering: could not open file " + filename);

	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() == '#') {
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		double r;
		infile >> r;
		if (!infile)
			break;
		tabRate.push_back(r / Mpc);
	}

	infile.close();
}

void ElasticScattering::initCDF(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("ElasticScattering: could not open file " + filename);

	tabCDF.clear();
	std::string line;
	double a;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;

		std::stringstream lineStream(line);
		lineStream >> a;

		std::vector<double> cdf(neps);
		for (size_t i = 0; i < neps; i++) {
			lineStream >> a;
			cdf[i] = a;
		}
		tabCDF.push_back(cdf);
	}

	infile.close();
}

void ElasticScattering::process(Candidate *candidate) const {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();

	if (not isNucleus(id))
		return;

	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));
	if ((lg < lgmin) or (lg > lgmax))
		return;

	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	double step = candidate->getCurrentStep();
	while (step > 0) {

		double rate = interpolateEquidistant(lg, lgmin, lgmax, tabRate);
		rate *= Z * N / double(A);  // TRK scaling
		rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);  // cosmological scaling

		// check for interaction
		Random &random = Random::instance();
		double randDist = -log(random.rand()) / rate;
		if (step < randDist)
			return;

		// draw random background photon energy from CDF
		size_t i = floor((lg - lgmin) / (lgmax - lgmin) * (nlg - 1)); // index of closest gamma tabulation point
		size_t j = random.randBin(tabCDF[i]) - 1; // index of next lower tabulated eps value
		double binWidth = (epsmax - epsmin) / (neps - 1); // logarithmic bin width
		double eps = pow(10, epsmin + (j + random.rand()) * binWidth);

		// boost to lab frame
		double cosTheta = 2 * random.rand() - 1;
		double E = eps * candidate->current.getLorentzFactor() * (1. - cosTheta);

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		candidate->addSecondary(22, E, pos, 1., interactionTag);

		// repeat with remaining step
		step -= randDist;
	}
}

void ElasticScattering::setInteractionTag(std::string tag) {
	this -> interactionTag = tag;
}

std::string ElasticScattering::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
