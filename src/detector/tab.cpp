#include "include/detector/tab.h"

namespace ribll {

Tab::Tab(unsigned int run, unsigned int index)
: Adssd(run, "tab"+std::to_string(index)) {
}

}		// namespace ribll