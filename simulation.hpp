// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#pragma once

#include <random>
#include <vector>
#include <algorithm>
#include "universe.hpp"
#include "observable.hpp"

class Simulation {
public:
	static double lambda;
	static double isingJ;

	static void start(int sweeps, double lambda_, int targetVolume_, int seed = 0, double J = 1);

	static void addObservable(Observable& o) {
		observables.push_back(&o);
	}

	static bool pinch;

	static std::array<int, 3> moveFreqs;
	static int attemptMove();

private:
	static std::default_random_engine rng;

	static int targetVolume;
	static double epsilon;
	static bool measuring;

	static std::vector<Observable*> observables;
	static std::vector<Triangle::Label> spinCluster;
	static double avgClusterSize;

	static int dot(std::array<int,Triangle::nSpins> s1, std::array<int,Triangle::nSpins> s2) {
		int res = 0;
		for (int i = 0; i < Triangle::nSpins; i++) res += s1[i]*s2[i];
		return res;
	}

	static void sweep();

	static bool moveAdd();
	static bool moveDelete();
	static bool moveFlip();
	static bool moveSpinFlip();
	static double moveSpinClusterFlip();

	static void prepare();

	static void tune();
	static void grow();
	static void thermalize();
	static void sweepSpin();
};
