#ifndef __T0D1_H__
#define __T0D1_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d1 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d1(unsigned int run, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~T0d1() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	inline size_t FrontStrip() const override {
		return 64;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const override {
		return 64;
	}


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	int NormalizeFilter(int iteration);

protected:

	//-------------------------------------------------------------------------
	//								geometry
	//-------------------------------------------------------------------------

	/// @brief calculate the position from strip index
	/// @param[in] front_strip front strip
	/// @param[in] back_strip front strip
	/// @returns vector point to the position
	///
	virtual ROOT::Math::XYZVector CalculatePosition(
		double front_strip,
		double back_strip
	) const override;


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(
		TChain *chain,
		int iteration
	) override;
};



}		// namespace ribll

#endif		// __T0D1_H__