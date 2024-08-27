#include "include/spectrum/asymmetric_voigt.h"

#include <cmath>
#include <complex>

#include <RooMath.h>

namespace ribll {

AsymmetricVoigtian::AsymmetricVoigtian(
	const char *name,
	const char *title,
	RooAbsReal &x,
	RooAbsReal &mean,
	RooAbsReal &g,
	RooAbsReal &sigma,
	double threshold
)
: RooAbsPdf(name, title)
, x_("x", "Dependent", this, x)
, mean_("mean", "Mean", this, mean)
, g_("width", "Breit-Wigner width factor", this, g)
, sigma_("sigma", "Resolution", this, sigma)
, threshold_(threshold) {}


AsymmetricVoigtian::AsymmetricVoigtian(
	const AsymmetricVoigtian &other,
	const char *name
)
: RooAbsPdf(other, name)
, x_("x", this, other.x_)
, mean_("mean", this, other.mean_)
, g_("width", this, other.g_)
, sigma_("sigma", this, other.sigma_)
, threshold_(other.threshold_) {}


double AsymmetricVoigtian::evaluate() const {
	// return 0 value for under threshold
	if (x_ < threshold_) return 0.0;

	double s = sigma_;
	double w = g_ * sqrt(x_ - threshold_);
	double coef = -0.5 / (s * s);
	double arg = x_ - mean_;

	// return constant for zero width and sigma
	if (s == 0.0 && w == 0.0) return 1.0;
	// Breit-Wigner for zero sigma
	if (s == 0.0) return 1.0 / (arg*arg + 0.25*w*w);
	// Gauss for zero width
	if (w == 0.0) return exp(coef*arg*arg);

	// actual Voigtian for non-trivial width and sigma
	double c = 1.0 / (sqrt(2.0) * s);
	double a = 0.5 * c * w;
	double u = c * arg;
	std::complex<double> z(u, a);
	std::complex<double> v(0.0);

	v = RooMath::faddeeva(z);

	return c * v.real();
}


};