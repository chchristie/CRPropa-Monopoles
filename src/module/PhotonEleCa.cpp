#include "crpropa/module/PhotonEleCa.h"
#include "crpropa/Units.h"

#include "EleCa/Propagation.h"
#include "EleCa/Particle.h"
#include "EleCa/Common.h"
#include "kiss/logger.h"

#include <vector>

namespace crpropa {

PhotonEleCa::PhotonEleCa(const std::string background,
		const std::string &outputFilename) :
		propagation(new eleca::Propagation), saveOnlyPhotonEnergies(false) {
	KISS_LOG_WARNING << "EleCa propagation is deprecated and is no longer supported. Please use the EM* (EMPairProduction, EMInverseComptonScattering, ...) modules instead.\n";
	propagation->ReadTables(getDataPath("EleCa/eleca.dat"));
	propagation->InitBkgArray(background);
	output.open(outputFilename.c_str());
}

PhotonEleCa::~PhotonEleCa() {
}

void PhotonEleCa::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return; // do nothing if not a photon

	double z = candidate->getRedshift();
	if (z == 0)
		z = eleca::Mpc2z(
				(candidate->current.getPosition() - observer).getR() / Mpc);
	eleca::Particle p0(candidate->current.getId(),
			candidate->current.getEnergy() / eV, z);
	std::vector<eleca::Particle> ParticleAtMatrix;
	std::vector<eleca::Particle> ParticleAtGround;
	ParticleAtMatrix.push_back(p0);

	while (ParticleAtMatrix.size() > 0) {

		eleca::Particle p1 = ParticleAtMatrix.back();
		ParticleAtMatrix.pop_back();

		if (p1.IsGood()) {
			propagation->Propagate(p1, ParticleAtMatrix, ParticleAtGround);
		}
	}

#pragma omp critical
	{
		if (saveOnlyPhotonEnergies) {
			for (int i = 0; i < ParticleAtGround.size(); ++i) {
				eleca::Particle &p = ParticleAtGround[i];
				if (p.GetType() != 22)
					continue;
				output << p.GetEnergy() << "\n";
			}
		} else {
			propagation->WriteOutput(output, p0, ParticleAtGround);
		}
	}

	candidate->setActive(false);
	return;
}

void PhotonEleCa::setObserver(const Vector3d &position) {
	observer = position;
}

void PhotonEleCa::setSaveOnlyPhotonEnergies(bool photonsOnly) {
	saveOnlyPhotonEnergies = photonsOnly;
}

std::string PhotonEleCa::getDescription() const {
	std::stringstream s;
	s << "PhotonEleCa";
	return s.str();
}

} // namespace crpropa
