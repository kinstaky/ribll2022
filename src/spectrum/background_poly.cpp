#include "include/spectrum/background_poly.h"

namespace ribll {

BackgroundPoly::BackgroundPoly(
	const char *name,
	const char *title,
	RooAbsReal &x,
	RooAbsReal &p1,
	RooAbsReal &p0
)
: RooAbsPdf(name, title)
, x_("x", "Dependent", this, x)
, p1_("p1", "p1", this, p1)
, p0_("p0", "p0", this, p0) {}


BackgroundPoly::BackgroundPoly(
	const BackgroundPoly &other,
	const char *name
)
: RooAbsPdf(other, name)
, x_("x", this, other.x_)
, p1_("p1", this, other.p1_)
, p0_("p0", this, other.p0_) {}


double BackgroundPoly::evaluate() const {
	double x = x_;
	double result = -x*x + p1_*x + p0_;
	if (result <= 0.0) return 0.0;
	return result;
}


}		// namespace ribll