// Copyright 2018 Joren Brunekreef and Andrzej Görlich
#ifndef MAGN_TOT_HPP_
#define MAGN_TOT_HPP_

#include "../observable.hpp"
#include "../universe.hpp"

class TotMagn : public Observable {
public:
	TotMagn(std::string id) : Observable(id) { name = "magn_tot"; }

	void process();
};
#endif  // MAGN_TOT_HPP_
