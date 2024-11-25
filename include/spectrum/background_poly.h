#ifndef __BACKGROUND_POLY_H__
#define __BACKGROUND_POLY_H__

#include <RooAbsPdf.h>
#include <RooRealProxy.h>

namespace ribll {

class BackgroundPoly : public RooAbsPdf {
public:

	/// @brief empty constructor
	///
	BackgroundPoly() {}


	/// @brief normal constructor
	/// @param[in] name P.D.F. name
	/// @param[in] title P.D.F. title
	/// @param[in] x observable
	/// @param[in] root0 the first root
	/// @param[in] root1 the second root
	///
	BackgroundPoly(
		const char *name,
		const char *title,
		RooAbsReal &x,
		RooAbsReal &root0,
		RooAbsReal &root1
	);


	/// @brief copy constructor
	/// @param[in] other another object to copy
	/// @param[in] name P.D.F. name
	///
	BackgroundPoly(
		const BackgroundPoly &other,
		const char *name = nullptr
	);


	/// @brief clone
	/// @param[in] new_name new name
	/// @returns pointer to new instance
	///
	TObject *clone(const char *new_name) const override {
		return new BackgroundPoly(*this, new_name);
	}

protected:
	RooRealProxy x_;
	RooRealProxy root0_;
	RooRealProxy root1_;

	double evaluate() const override;
};

};

#endif	// __BACKGROUND_POLY_H__