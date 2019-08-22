// Copyright 2018 Joren Brunekreef and Andrzej Görlich
#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <random>
#include "universe.hpp"
#include "observable.hpp"

class Simulation {
public:
	static int sweepSize;
	static double lambda;
	static double gsq;

	static void start(int sweeps, int sweepSize_, double lambda_, int targetVolume_ = 0);

	static void addObservable(Observable& o) {
		observables.push_back(&o);
	}

private:
	static std::default_random_engine rng;

	static int targetVolume;
	static double constexpr epsilon = 0.004;

	static std::vector<Observable*> observables;

	static void sweep();

	static bool moveAdd();
	static bool moveDelete();
	static bool moveFlip();

	static void prepare();
};
#endif  // SIMULATION_HPP_