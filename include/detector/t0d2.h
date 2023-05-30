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


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief filter normalize events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	int NormalizeFilter(int iteration);


	//-------------------------------------------------------------------------
	//								time
	//-------------------------------------------------------------------------

	/// @brief analyze time
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AnalyzeTime() override;


	// /// @brief check whether time curve is appropriate
	// /// @param[in] condition 0-single hit, 1-double hit small, 2-double hit big
	// /// @param[in] side front (0) or back (1)
	// /// @param[in] strip strip number
	// /// @param[in] energy normalized energy
	// /// @param[in] time normalized time
	// /// @returns true if pass time check, false otherwise
	// ///
	// virtual bool CheckTime(
	// 	int condition,
	// 	size_t side,
	// 	unsigned short strip,
	// 	double energy,
	// 	double time
	// ) override;

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

#endif		// __T0D2_H__