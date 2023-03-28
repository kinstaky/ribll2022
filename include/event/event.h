#ifndef __EVENT_H__
#define __EVENT_H__

#include "include/defs.h"

#include <TTree.h>

namespace ribll {

class Event {
public:
	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) = 0;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) = 0;
};

}

#endif		// __EVENT_H__