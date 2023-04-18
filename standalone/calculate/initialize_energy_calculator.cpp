#include <iostream>

#include "include/defs.h"
#include "include/calculator/range_energy_calculator.h"
#include "include/calculator/delta_energy_calculator.h"

using namespace ribll::elc;

int main() {
	const std::vector<ProjectileMaterial> list{
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
		{"4He", "Mylar"},
		{"7Li", "Si"},
		{"10Be", "Si"},
		{"14C", "Si"}
	};
	if (RangeEnergyCalculator::Initialize(list)) {
		std::cerr << "Error: Initialize range-energy calculator failed.\n";
		return -1;
	}

	const std::vector<std::string> t0_projectiles{
		"1H",
		"2H",
		"4He",
		"7Li",
		"10Be",
		"14C"
	};
	if (DeltaEnergyCalculator::Initialize(
		"t0", ribll::t0_thickness, t0_projectiles
	)) {
		std::cerr << "Error: Initialize delta-energy calculator failed.\n";
		return -1;
	}
	return 0;
}