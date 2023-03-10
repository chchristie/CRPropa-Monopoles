#include "crpropa/module/EMCascade.h"
#include "crpropa/Cosmology.h"
#include "crpropa/Units.h"

#include "dint/DintEMCascade.h"
#include "kiss/logger.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

namespace crpropa {

EMCascade::EMCascade() : nE(170), logEmin(7), logEmax(24), dlogE(0.1) {
	KISS_LOG_WARNING << "EMCascade is deprecated and is no longer supported. Please use the EM* (EMPairProduction, EMInverseComptonScattering, ...) modules instead.\n";
	setDistanceBinning(1000 * Mpc, 1000);
}

void EMCascade::setDistanceBinning(double Dmax, int nD) {
	this->Dmax = Dmax;
	this->nD = nD;
	this->dD = Dmax / nD;
	init();
}

void EMCascade::init() {
	photonHist.reserve(nD * nE);
	photonHist.assign(nD * nE, 0);
	electronHist.reserve(nD * nE);
	electronHist.assign(nD * nE, 0);
	positronHist.reserve(nD * nE);
	positronHist.assign(nD * nE, 0);
}

std::string EMCascade::getDescription() const {
	std::stringstream s;
	s << "EMCascade";
	return s.str();
}

void EMCascade::process(Candidate *candidate) const {
	int id = candidate->current.getId();
	if ((id != 22) and (id != 11) and (id != -11))
		return;

	candidate->setActive(false);

	double logE = log10(candidate->current.getEnergy() / eV);
	double D = candidate->current.getPosition().getR();  // distance to (0,0,0)

	if ((logE < logEmin) or (logE > logEmax))
		return;
	if (D > Dmax)
		return;

	int iE = (logE - logEmin) / dlogE;
	int iD = D / dD;
	int i = (iD * nE) + iE;

#pragma omp critical
	{
	if (id == 22)
		photonHist[i] += 1;
	else if (id == 11)
		electronHist[i] += 1;
	else
		positronHist[i] += 1;
	}
}

void EMCascade::save(const std::string &filename) {
	std::ofstream outfile(filename.c_str());
	if (!outfile) {
		std::stringstream s;
		s << "EMCascade: could not open " << filename;
		throw std::runtime_error(s.str());
	}
	outfile << "# D/Mpc log10(E/eV) nPhotons nElectrons nPositrons\n";
	for (int i = 0; i < (nD * nE); i++) {
		div_t divresult = div(i, nE);
		double D = (divresult.quot + 0.5) * dD / Mpc;
		double logE = logEmin + (divresult.rem + 0.5) * dlogE;
		outfile << D << "\t";
		outfile << logE << "\t";
		outfile << photonHist[i] << "\t";
		outfile << electronHist[i] << "\t";
		outfile << positronHist[i] << "\n";
	}
	outfile.close();
}

void EMCascade::load(const std::string &filename) {
	std::ifstream infile(filename.c_str());
	if (!infile) {
		std::stringstream s;
		s << "EMCascade: could not open " << filename;
		throw std::runtime_error(s.str());
	}

	infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header
	double D, lE, h1, h2, h3;
	for (int i = 0; i < (nD * nE); i++) {
		infile >> D >> lE >> h1 >> h2 >> h3;
		if (!infile.good())
			throw std::runtime_error("EMCascde: error reading file");
		photonHist[i] += h1;
		electronHist[i] += h2;
		positronHist[i] += h3;
	}
	infile.close();
}

void EMCascade::runCascade(const std::string &filename, int IRBFlag,
		int RadioFlag, double Bfield, double cutCascade) {

	// set up DINT
	std::string dataPath = getDataPath("dint");
	double B = Bfield / gauss;
	double h = H0() * Mpc / 1000;
	DintEMCascade dint(IRBFlag, RadioFlag, dataPath, B, h, omegaM(), omegaL());

	Spectrum inputSpectrum, outputSpectrum;
	NewSpectrum(&inputSpectrum, nE);
	NewSpectrum(&outputSpectrum, nE);
	InitializeSpectrum(&outputSpectrum);

	// step-wise cascade calculation
	for (int iD = nD - 1; iD >= 0; iD--) {
		// make output of previous step the new input and reset output
		SetSpectrum(&inputSpectrum, &outputSpectrum);
		InitializeSpectrum(&outputSpectrum);

		// add new particles to input
		double count = 0;
		for (int iE = 0; iE < nE; iE++) {
			int i = (iD * nE) + iE;
			inputSpectrum.spectrum[PHOTON][iE] += photonHist[i];
			inputSpectrum.spectrum[ELECTRON][iE] += electronHist[i];
			inputSpectrum.spectrum[POSITRON][iE] += positronHist[i];
			count += inputSpectrum.spectrum[PHOTON][iE];
			count += inputSpectrum.spectrum[ELECTRON][iE];
			count += inputSpectrum.spectrum[POSITRON][iE];
		}
		// skip step if input spectrum empty
		if (count == 0)
			continue;

		// start and stop distance [Mpc,light travel] from bin center to bin center
		double D1 = comoving2LightTravelDistance( (iD + 0.5) * dD );
		double D0 = comoving2LightTravelDistance( std::max((iD - 0.5) * dD, 0.) );

		// propagate distance step
		dint.propagate(D1/Mpc, D0/Mpc, &inputSpectrum, &outputSpectrum, cutCascade);
	}

	// write output
	std::ofstream outfile(filename.c_str());
	if (!outfile) {
		std::stringstream s;
		s << "EMCascade: could not open " << filename;
		throw std::runtime_error(s.str());
	}
	outfile << "# log10(E/eV) photons electrons positrons\n";
	for (int iE = 0; iE < nE; iE++) {
		outfile << std::setw(5) << logEmin + (iE + 0.5) * dlogE;
		for (int s = 0; s < 3; s++)
			outfile << std::setw(13) << outputSpectrum.spectrum[s][iE];
		outfile << "\n";
	}
	outfile.close();

	// clear the histogram
	photonHist.assign(nD * nE, 0);
	electronHist.assign(nD * nE, 0);
	positronHist.assign(nD * nE, 0);

	DeleteSpectrum(&outputSpectrum);
	DeleteSpectrum(&inputSpectrum);
}

} // namespace crpropa
