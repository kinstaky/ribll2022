#include "include/detector/adssd.h"

#include <iostream>


namespace ribll {

Adssd::Adssd(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Dssd(run, name, tag) {
}


int Adssd::MatchTrigger(double, double) {
	return Detector::VmeMatchTrigger<DssdFundamentalEvent>();
}


// int Adssd::Merge() {
// 		// input file name
// 	TString fundamental_file_name;
// 	fundamental_file_name.Form(
// 		"%s%s%s-fundamental"
// 	)
	
// 	DssdFundamentalEvent fundamental_event;

// }

}		// namespace ribll