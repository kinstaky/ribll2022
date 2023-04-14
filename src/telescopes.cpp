#include "include/telescopes.h"

namespace ribll {

std::shared_ptr<Telescope> CreateTelescope(
	const std::string &name,
	unsigned int run,
	const std::string &tag
) {
	if (name == "t0") {
		return std::make_shared<T0>(run, tag);
	} else if (name.size() == 4 && name.substr(0, 3) == "taf") {
		unsigned int index = name[3] - '0';
		if (index <= 5) {
			return std::make_shared<Taf>(run, index, tag);
		}
	}

	std::cerr << "Error: Create telescope " << name
		<< " with tag " << tag << " failed.\n";
	return nullptr;
}

}		// namespace ribll