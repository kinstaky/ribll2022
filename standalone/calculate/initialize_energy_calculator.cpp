#include <iostream>

#include "include/defs.h"
#include "include/calculator/lost_energy_calculator.h"
#include "include/calculator/range_energy_calculator.h"
#include "include/calculator/delta_energy_calculator.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll::elc;

int main() {
	const std::vector<std::string> loss_projectiles {
		"4He", "10Be", "14C"
	};

	if (LostEnergyCalculator::Initialize(loss_projectiles)) {
		std::cerr << "Error: Initialize loss-energy calculator failed.\n";
		return -1;
	}


	const std::vector<ProjectileMaterial> list {
		{"1H", "Si"},
		{"1H", "Al"},
		{"1H", "Mylar"},
		{"1H", "CD2"},
		{"2H", "Si"},
		{"2H", "Al"},
		{"2H", "Mylar"},
		{"2H", "CD2"},
		{"3H", "Si"},
		{"3H", "Al"},
		{"3H", "Mylar"},
		{"3He", "Si"},
		{"3He", "Al"},
		{"3He", "Mylar"},
		{"4He", "Si"},
		{"4He", "Al"},
		{"4He", "Mylar"},
		{"4He", "CD2"},
		{"6He", "Si"},
		{"6Li", "Si"},
		{"7Li", "Si"},
		{"7Be", "Si"},
		{"9Be", "Si"},
		{"10Be", "Si"},
		{"10Be", "CD2"},
		{"10B", "Si"},
		{"11B", "Si"},
		{"12B", "Si"},
		{"13B", "Si"},
		{"12C", "Si"},
		{"13C", "Si"},
		{"14C", "Si"},
		{"14C", "CD2"}
	};
	if (RangeEnergyCalculator::Initialize(list)) {
		std::cerr << "Error: Initialize range-energy calculator failed.\n";
		return -1;
	}

	const std::vector<std::string> t0_projectiles{
		"1H",
		"2H",
		"3H",
		"3He",
		"4He",
		"6He",
		"6Li",
		"7Li",
		"9Be",
		"10Be",
		"14C"
	};
	if (DeltaEnergyCalculator::Initialize(
		"t0", ribll::t0_thickness, t0_projectiles
	)) {
		std::cerr << "Error: Initialize delta-energy calculator failed.\n";
		return -1;
	}

	const std::vector<std::string> t0_projectiles2 = {"4He", "10Be"};
	if (D2EnergyCalculator::Initialize(
		ribll::t0_thickness, t0_projectiles2
	)) {
		std::cerr << "Error: Initialize d2-energy calculator failed.\n";
		return -1;
	}
	return 0;
}