#ifndef __ASYMMETRIC_VOIGT_H__
#define __ASYMMETRIC_VOIGT_H__

#include <RooAbsPdf.h>
#include <RooRealProxy.h>

namespace ribll {

class AsymmetricVoigtian : public RooAbsPdf {
public:
	/// @brief empty constructor
	///
	AsymmetricVoigtian() {}


	/// @brief normal constructor
	/// @param[in] name name of this P.D.F.
	/// @param[in] title P.D.F. title
	/// @param[in] x observable
	/// @param[in] mean mean of the distribution
	/// @param[in] g width parameter of the distribution, Gamma = g*sqrt(E-Er)
	/// @param[in] sigma resolution width
	/// @param[in] threshold threshold to breakup
	/// @param[in] efficiency efficiency graph
	///
	AsymmetricVoigtian(
		const char *name,
		const char *title,
		RooAbsReal &x,
		RooAbsReal &mean,
		RooAbsReal &g,
		RooAbsReal &sigma,
		double threshold,
		TGraph *efficiency = nullptr
	);


	/// @brief copy constructor
	/// @param[in] other the other voigt to copy
	/// @param[in] name P.D.F. name
	///
	AsymmetricVoigtian(
		const AsymmetricVoigtian &other,
		const char *name = nullptr
	);


	/// @brief clone
	/// @param[in] new_name new name
	/// @returns pointer to new instance
	///
	TObject *clone(const char *new_name) const override {
		return new AsymmetricVoigtian(*this, new_name);
	}

protected:
	RooRealProxy x_;
	RooRealProxy mean_;
	RooRealProxy g_;
	RooRealProxy sigma_;
	double threshold_;
	TGraph *efficiency_;

	double evaluate() const override;
};

}	// namespace ribll

#endif	// __ASYMMETRIC_VOIGT_H__