#ifndef __T0D3_H__
#define __T0D3_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d3 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d3(unsigned int run, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~T0d3() = default;


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	int NormalizeFilter(int iteration);

protected:
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

#endif		// __T0D3_H__