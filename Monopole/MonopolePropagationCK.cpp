#include "MonopolePropagationCK.h"

#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {

// Cash-Karp coefficients
const double cash_karp_a[] = {
	0., 0., 0., 0., 0., 0.,
	1. / 5., 0., 0., 0., 0., 0.,
	3. / 40., 9. / 40., 0., 0., 0., 0.,
	3. / 10., -9. / 10., 6. / 5., 0., 0., 0.,
	-11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0.,
	1631. / 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253. / 4096., 0.
};

const double cash_karp_b[] = {
	37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.
};

const double cash_karp_bs[] = {
	2825. / 27648., 0., 18575. / 48384., 13525. / 55296., 277. / 14336., 1. / 4.
};

void MonopolePropagationCK::tryStep(const Y &y, Y &out, Y &error, double h,
		MParticleState &particle, double z) const {
	std::vector<Y> k;
	k.reserve(6);

	out = y;
	error = Y(0);

	// calculate the sum of b_i * k_i
	for (size_t i = 0; i < 6; i++) {

		Y y_n = y;
		for (size_t j = 0; j < i; j++)
			y_n += k[j] * a[i * 6 + j] * h;

		// update k_i
		k[i] = dYdt(y_n, particle, z);

		out += k[i] * b[i] * h;
		error += k[i] * (b[i] - bs[i]) * h;
	}
}

MonopolePropagationCK::Y MonopolePropagationCK::dYdt(const Y &y, MParticleState &p, double z) const {
	// Derivative of position is velocity
	Vector3d velocity = p.getVelocity();
	
	// get B field at particle position
	Vector3d B = getFieldAtPosition(y.x, z);

	// Lorentz force: du/dt = dp/dt = F = g*B + q*vxB
	Vector3d dudt = p.getMcharge() * B + p.getCharge() * velocity.cross(B);
	return Y(velocity, dudt);
}

MonopolePropagationCK::MonopolePropagationCK(ref_ptr<MagneticField> field, double tolerance,
		double minStep, double maxStep) :
		minStep(0) {
	setField(field);
	setTolerance(tolerance);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);

	// load Cash-Karp coefficients
	a.assign(cash_karp_a, cash_karp_a + 36);
	b.assign(cash_karp_b, cash_karp_b + 6);
	bs.assign(cash_karp_bs, cash_karp_bs + 6);
}

void MonopolePropagationCK::Mprocess(MCandidate *candidate, MParticleState &current) const {
	Y yIn(current.getPosition(), current.getMomentum());
	double step = maxStep;

	// rectilinear propagation for neutral particles
	if (current.getCharge() == 0 && current.getMcharge() == 0) {
		step = clip(candidate->getNextStep(), minStep, maxStep);
		current.setPosition(yIn.x + current.getVelocity() * step / c_light);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	Y yOut, yErr;
	double newStep = step;
	double z = candidate->getRedshift();


	// if minStep is the same as maxStep the adaptive algorithm with its error
	// estimation is not needed and the computation time can be saved:
	if (minStep == maxStep){
		tryStep(yIn, yOut, yErr, step / c_light, current, z);
	} else {
		step = clip(candidate->getNextStep(), minStep, maxStep);
		newStep = step;
		double r = 42;  // arbitrary value

		// try performing step until the target error (tolerance) or the minimum/maximum step size has been reached
		while (true) {
			tryStep(yIn, yOut, yErr, step / c_light, current, z);
			r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
			if (r > 1) {  // large direction error relative to tolerance, try to decrease step size
				if (step == minStep)  // already minimum step size
					break;
				else {
					newStep = step * 0.95 * pow(r, -0.2);
					newStep = std::max(newStep, 0.1 * step); // limit step size decrease
					newStep = std::max(newStep, minStep); // limit step size to minStep
					step = newStep;
				}
			} else {  // small direction error relative to tolerance, try to increase step size
				if (step != maxStep) {  // only update once if maximum step size yet not reached
					newStep = step * 0.95 * pow(r, -0.2);
					newStep = std::min(newStep, 5 * step); // limit step size increase
					newStep = std::min(newStep, maxStep); // limit step size to maxStep
				}
				break;
			}
		}
	}

	current.setPosition(yOut.x);
	double m = current.getMass();
	double E = sqrt(pow(m*c_squared, 2) + pow(yOut.u.getR()*c_light, 2)) - m*c_squared; 
	current.setEnergy(E);
	current.setDirection(yOut.u.getUnitVector());
	
	candidate->setCurrentStep(step);
	candidate->setNextStep(newStep);
}

void MonopolePropagationCK::setField(ref_ptr<MagneticField> f) {
	field = f;
}

ref_ptr<MagneticField> MonopolePropagationCK::getField() const {
	return field;
}

Vector3d MonopolePropagationCK::getFieldAtPosition(Vector3d pos, double z) const {
	Vector3d B(0, 0, 0);
	try {
		// check if field is valid and use the field vector at the
		// position pos with the redshift z
		if (field.valid())
			B = field->getField(pos, z);
	} catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "MonopolePropagationCK: Exception in MonopolePropagationCK::getFieldAtPosition.\n"
				<< e.what();
	}	
	return B;
}

void MonopolePropagationCK::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"PropagationCK: target error not in range 0-1");
	tolerance = tol;
}

void MonopolePropagationCK::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("PropagationCK: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("PropagationCK: minStep > maxStep");
	minStep = min;
}

void MonopolePropagationCK::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("PropagationCK: maxStep < minStep");
	maxStep = max;
}

double MonopolePropagationCK::getTolerance() const {
	return tolerance;
}

double MonopolePropagationCK::getMinimumStep() const {
	return minStep;
}

double MonopolePropagationCK::getMaximumStep() const {
	return maxStep;
}

std::string MonopolePropagationCK::getDescription() const {
	std::stringstream s;
	s << "Propagation of dyons in magnetic fields using the Cash-Karp method.";
	s << " Target error: " << tolerance;
	s << ", Minimum Step: " << minStep / kpc << " kpc";
	s << ", Maximum Step: " << maxStep / kpc << " kpc";
	return s.str();
}

} // namespace crpropa
