#include "include/spectrum/background_poly.h"

namespace ribll {

BackgroundPoly::BackgroundPoly(
	const char *name,
	const char *title,
	RooAbsReal &x,
	RooAbsReal &root0,
	RooAbsReal &root1
)
: RooAbsPdf(name, title)
, x_("x", "Dependent", this, x)
, root0_("root0", "root0", this, root0)
, root1_("root1", "root1", this, root1) {}


BackgroundPoly::BackgroundPoly(
	const BackgroundPoly &other,
	const char *name
)
: RooAbsPdf(other, name)
, x_("x", this, other.x_)
, root0_("root0", this, other.root0_)
, root1_("root1", this, other.root1_) {}


double BackgroundPoly::evaluate() const {
	double x = x_;
	double result = -(x-root0_) * (x-root1_);
	if (result <= 0.0) return 0.0;
	return result;
}


}		// namespace ribll