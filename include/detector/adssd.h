#ifndef __ADSSD_H__
#define __ADSSD_H__

#include <string>

#include <Math/Vector3D.h>
#include <TMath.h>

#include "include/detector/dssd.h"
#include "include/event/dssd_event.h"

namespace ribll {

class Adssd : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	/// @param[in] tag trigger tag
	///
	Adssd(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~Adssd() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	virtual inline size_t FrontStrip() const {
		return 16;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const {
		return 8;
	}

	//-------------------------------------------------------------------------
	//								merge
	//-------------------------------------------------------------------------

	/// @brief merge adjacent event in the same side and merge events of two sides
	/// @param[in] energy_diff tolerant energy relateive difference
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Merge(double energy_diff) override;

protected:
	//-------------------------------------------------------------------------
	//								geometry
	//-------------------------------------------------------------------------

	/// @brief calculate the position from strip index
	/// @param[in] front_strip front strip
	/// @param[in] back_strip back strip
	/// @returns vector class point to position
	///
	ROOT::Math::Polar3DVector CalculatePosition(
		unsigned short front_strip,
		unsigned short back_strip
	) const;


	ROOT::Math::Polar3DVector center_;
	std::pair<double, double> radius_range_;
	std::pair<double, double> phi_range_;
};

}	// namespace ribll

#endif 		// __ADSSD_H__