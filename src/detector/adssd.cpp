#include "include/detector/adssd.h"

namespace ribll {

Adssd::Adssd(unsigned int run, const std::string &name)
: Dssd(run, name) {
}


Taf::Taf(unsigned int run, unsigned int index)
: Adssd(run, "taf"+std::to_string(index)) {

}


Tab::Tab(unsigned int run, unsigned int index)
: Adssd(run, "tab"+std::to_string(index)) {
}

}