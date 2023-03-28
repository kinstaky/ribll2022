#include "include/detector/tab.h"

namespace ribll {

Tab::Tab(unsigned int run, unsigned int index, const std::string &tag)
: Adssd(run, "tab"+std::to_string(index), tag) {
}

}		// namespace ribll