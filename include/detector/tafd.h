#ifndef __TAF_H__
#define __TAF_H__

#include "include/detector/adssd.h"

namespace ribll {

class Tafd : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	/// @param[in] tag trigger tag
	///
	Tafd(unsigned int run, unsigned int index, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~Tafd() = default;


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief extract trigger with detector events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ExtractTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate() override;

protected:

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

private:
	// detector index, 0 to 5
	int index_;
};

}		// namespace ribll

#endif		// __TAF_H__