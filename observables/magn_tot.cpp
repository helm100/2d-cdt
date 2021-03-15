// Copyright 2018 Joren Brunekreef and Andrzej Görlich
#include "magn_tot.hpp"

void TotMagn::process() {
	std::string tmp; // = std::to_string(Universe::magnTot);
	
	for (auto m : Universe::magnTot) {
		tmp += std::to_string(m);
		tmp += " ";
	}
	tmp.pop_back();

	output = tmp;
}
