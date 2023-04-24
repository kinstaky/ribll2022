#ifndef __T0D2_H__
#define __T0D2_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d2 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d2(unsigned int run, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~T0d2() = default;

protected:
	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(TChain *chain, int iteration) override;


	/// @brief check whether energy is suitable for fitting
	/// @param[in] side side to normalize
	/// @param[in] event fundamental event
	/// @returns true if pass check, false not pass
	///
	virtual bool NormEnergyCheck(
		size_t side,
		const DssdNormalizeEvent &event
	) const override;
};

}		// namespace ribll

#endif		// __T0D2_H__