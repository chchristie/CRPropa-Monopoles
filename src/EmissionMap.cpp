#include "crpropa/EmissionMap.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"

#include "kiss/logger.h"

#include <fstream>

namespace crpropa {

CylindricalProjectionMap::CylindricalProjectionMap() : nPhi(360), nTheta(180), dirty(false), pdf(nPhi* nTheta, 0), cdf(nPhi* nTheta, 0) {
	sPhi = 2. * M_PI / nPhi;
	sTheta = 2. / nTheta;
}

CylindricalProjectionMap::CylindricalProjectionMap(size_t nPhi, size_t nTheta) : nPhi(nPhi), nTheta(nTheta), dirty(false), pdf(nPhi* nTheta, 0), cdf(nPhi* nTheta, 0) {
	sPhi = 2 * M_PI / nPhi;
	sTheta = 2. / nTheta;
}

void CylindricalProjectionMap::fillBin(const Vector3d& direction, double weight) {
	size_t bin = binFromDirection(direction);
	fillBin(bin, weight);
}

void CylindricalProjectionMap::fillBin(size_t bin, double weight) {
	pdf[bin] += weight;
	dirty = true;
}

Vector3d CylindricalProjectionMap::drawDirection() const {
	if (dirty)
		updateCdf();

	size_t bin = Random::instance().randBin(cdf);

	return directionFromBin(bin);
}

bool CylindricalProjectionMap::checkDirection(const Vector3d &direction) const {
	size_t bin = binFromDirection(direction);
	return pdf[bin];
}


const std::vector<double>& CylindricalProjectionMap::getPdf() const {
	return pdf;
}

std::vector<double>& CylindricalProjectionMap::getPdf() {
	return pdf;
}

const std::vector<double>& CylindricalProjectionMap::getCdf() const {
	return cdf;
}

size_t CylindricalProjectionMap::getNPhi() {
	return nPhi;
}

size_t CylindricalProjectionMap::getNTheta() {
	return nTheta;
}

/*
 * Cylindrical Coordinates
 * iPhi -> [0, 2*pi]
 * iTheta -> [0, 2]
 *
 * Spherical Coordinates
 * phi -> [-pi, pi]
 * theta -> [0, pi]
 */
size_t CylindricalProjectionMap::binFromDirection(const Vector3d& direction) const {
	// convert to cylindrical
	double phi = direction.getPhi() + M_PI;
	double theta = sin(M_PI_2 - direction.getTheta()) + 1;

	// to indices
	size_t iPhi = phi / sPhi;
	size_t iTheta = theta / sTheta;

	// interleave
	size_t bin =  iTheta * nPhi + iPhi;
	return bin;
}

Vector3d CylindricalProjectionMap::directionFromBin(size_t bin) const {
	// deinterleave
	double iPhi = bin % nPhi;
	double iTheta = bin / nPhi;

	// any where in the bin
	iPhi += Random::instance().rand();
	iTheta += Random::instance().rand();

	// cylindrical Coordinates
	double phi = iPhi * sPhi;
	double theta = iTheta * sTheta;

	// sphericala Coordinates
	phi = phi - M_PI;
	theta = M_PI_2 - asin(theta - 1);

	Vector3d v;
	v.setRThetaPhi(1.0, theta, phi);
	return v;
}

void CylindricalProjectionMap::updateCdf() const {
	if (dirty) {
		cdf[0] = pdf[0];
		for (size_t i = 1; i < pdf.size(); i++) {
			cdf[i] = cdf[i-1] + pdf[i];
		}
		dirty = false;
	}
}

EmissionMap::EmissionMap() : minEnergy(0.0001 * EeV), maxEnergy(10000 * EeV),
	nEnergy(8*2), nPhi(360), nTheta(180) {
	logStep = log10(maxEnergy / minEnergy) / nEnergy;
}

EmissionMap::EmissionMap(size_t nPhi, size_t nTheta, size_t nEnergy) : minEnergy(0.0001 * EeV), maxEnergy(10000 * EeV),
	nEnergy(nEnergy), nPhi(nPhi), nTheta(nTheta) {
	logStep = log10(maxEnergy / minEnergy) / nEnergy;
}

EmissionMap::EmissionMap(size_t nPhi, size_t nTheta, size_t nEnergy, double minEnergy, double maxEnergy) : minEnergy(minEnergy), maxEnergy(maxEnergy), nEnergy(nEnergy), nPhi(nPhi), nTheta(nTheta) {
	logStep = log10(maxEnergy / minEnergy) / nEnergy;
}

double EmissionMap::energyFromBin(size_t bin) const {
	return pow(10, log10(minEnergy) + logStep * bin);
}

size_t EmissionMap::binFromEnergy(double energy) const {
	return log10(energy / minEnergy) / logStep;
}

void EmissionMap::fillMap(int pid, double energy, const Vector3d& direction, double weight) {
	getMap(pid, energy)->fillBin(direction, weight);
}

void EmissionMap::fillMap(const ParticleState& state, double weight) {
	fillMap(state.getId(), state.getEnergy(), state.getDirection(), weight);
}

EmissionMap::map_t &EmissionMap::getMaps() {
	return maps;
}

const EmissionMap::map_t &EmissionMap::getMaps() const {
	return maps;
}

bool EmissionMap::drawDirection(int pid, double energy, Vector3d& direction) const {
	key_t key(pid, binFromEnergy(energy));
	map_t::const_iterator i = maps.find(key);

	if (i == maps.end() || !i->second.valid()) {
		return false;
	} else {
		direction = i->second->drawDirection();
		return true;
	}
}

bool EmissionMap::drawDirection(const ParticleState& state, Vector3d& direction) const {
	return drawDirection(state.getId(), state.getEnergy(), direction);
}

bool EmissionMap::checkDirection(int pid, double energy, const Vector3d& direction) const {
	key_t key(pid, binFromEnergy(energy));
	map_t::const_iterator i = maps.find(key);

	if (i == maps.end() || !i->second.valid()) {
		return false;
	} else {
		return i->second->checkDirection(direction);
	}
}

bool EmissionMap::checkDirection(const ParticleState& state) const {
	return checkDirection(state.getId(), state.getEnergy(), state.getDirection());
}

bool EmissionMap::hasMap(int pid, double energy) {
    key_t key(pid, binFromEnergy(energy));
    map_t::iterator i = maps.find(key);
    if (i == maps.end() || !i->second.valid())
		return false;
	else
		return true;
}

ref_ptr<CylindricalProjectionMap> EmissionMap::getMap(int pid, double energy) {
	key_t key(pid, binFromEnergy(energy));
	map_t::iterator i = maps.find(key);
	if (i == maps.end() || !i->second.valid()) {
		ref_ptr<CylindricalProjectionMap> cpm = new CylindricalProjectionMap(nPhi, nTheta);
		maps[key] = cpm;
		return cpm;
	} else {
		return i->second;
	}
}

void EmissionMap::save(const std::string &filename) {
	std::ofstream out(filename.c_str());
	out.imbue(std::locale("C"));

	for (map_t::iterator i = maps.begin(); i != maps.end(); i++) {
		if (!i->second.valid())
			continue;
		out << i->first.first << " " << i->first.second << " " << energyFromBin(i->first.second) << " ";
		out << i->second->getNPhi() << " " << i->second->getNTheta();
		const std::vector<double> &pdf = i->second->getPdf();
		for (size_t i = 0; i < pdf.size(); i++)
			out << " " << pdf[i];
		out << std::endl;
	}
}

void EmissionMap::merge(const EmissionMap *other) {
	if (other == 0)
		return;
	map_t::const_iterator i = other->getMaps().begin();
	map_t::const_iterator end = other->getMaps().end();
	for(;i != end; i++) {
		if (!i->second.valid())
			continue;

		std::vector<double> &otherpdf = i->second->getPdf();
		ref_ptr<CylindricalProjectionMap> cpm = getMap(i->first.first, i->first.second);

		if (otherpdf.size() != cpm->getPdf().size()) {
			throw std::runtime_error("PDF size mismatch!");
			break;
		}

		for (size_t k = 0; k < otherpdf.size(); k++) {
			cpm->fillBin(k, otherpdf[k]);
		}
	}
}

void EmissionMap::merge(const std::string &filename) {
	EmissionMap em;
	em.load(filename);
	merge(&em);
}

void EmissionMap::load(const std::string &filename) {
	std::ifstream in(filename.c_str());
	in.imbue(std::locale("C"));

	while(in.good()) {
		key_t key;
		double tmp;
		size_t nPhi_, nTheta_;
		in >> key.first >> key.second >> tmp;
		in >> nPhi_ >> nTheta_;

		if (!in.good()) {
			KISS_LOG_ERROR << "Invalid line: " << key.first << " " << key.second << " " << nPhi_ << " " << nTheta_;
			break;
		}

		if (nPhi != nPhi_) {
			KISS_LOG_WARNING << "nPhi mismatch: " << nPhi << " " << nPhi_;
		}
		if (nTheta != nTheta_) {
			KISS_LOG_WARNING << "nTheta mismatch: " << nTheta << " " << nTheta_;
		}

		ref_ptr<CylindricalProjectionMap> cpm = new CylindricalProjectionMap(nPhi_, nTheta_);
		std::vector<double> &pdf = cpm->getPdf();
		for (size_t i = 0; i < pdf.size(); i++)
			in >> pdf[i];

		if (in.good()) {
			maps[key] = cpm;
			// std::cout << "added " << key.first << " " << key.second << std::endl;
		} else {
			KISS_LOG_WARNING << "Invalid data in line: " << key.first << " " << key.second << " " << nPhi_ << " " << nTheta_ << "\n";
		}
	}

}

} // namespace crpropa
