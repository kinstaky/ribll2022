#include <iostream>

#include "include/energy_loss.h"

using namespace ribll::el;

int main() {
	std::vector<ProjectileMaterial> list = {
		{"1H", "Si"},
		{"1H", "Al"},
		{"1H", "Mylar"},
		{"2H", "Si"},
		{"2H", "Al"},
		{"2H", "Mylar"},
		{"3H", "Si"},
		{"3H", "Al"},
		{"3H", "Mylar"},
		{"3He", "Si"},
		{"3He", "Al"},
		{"3He", "Mylar"},
		{"4He", "Si"},
		{"4He", "Al"},
		{"4He", "Mylar"}
	};
	if (Initialize(list)) {
		std::cerr << "Error: Initialize energy loss calculator failed.\n";
		return -1;
	}
	return 0;
}