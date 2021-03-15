// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#include "simulation.hpp"
#include <iostream>
#include <chrono>
using namespace std::chrono;

std::default_random_engine Simulation::rng(0);  // TODO(JorenB): seed properly
int Simulation::targetVolume = 0;
double Simulation::lambda = 0;
double Simulation::epsilon = 0.02;//0.02
double Simulation::isingJ = 1.0;
std::vector<Observable*> Simulation::observables;
std::vector<Triangle::Label> Simulation::spinCluster (0,0);

std::array<int, 3> Simulation::moveFreqs = {1, 1, 1};

void Simulation::start(int measurements, double lambda_, int targetVolume_, int seed, double J) {
	targetVolume = targetVolume_;
	lambda = lambda_; // (Triangle::nSpins+1)*0.6931; //
	isingJ = J;

	for (auto o : observables) {
		o->clear();
	}

	rng.seed(seed);

	std::cout << "Number of spins: " << Triangle::nSpins << std::endl;

	if (!Universe::imported) {
		tune();
		// grow();
		// thermalize();
		// std::string exportFn = "geom/geometry-V"+std::to_string(targetVolume)+"-sl"+std::to_string(Universe::nSlices)+"-J"+std::to_string(int(isingJ*1000))+"-s"+std::to_string(seed)+".txt";
		// Universe::exportGeometry(exportFn, lambda);
	} else {
		thermalize();
	}

	// for (int i = 0; i < 20; i++) {
	// 	moveSpinFlip();
	// 	prepare();
	// 	for (int i = 0; i < 8; i++) std::cout << Universe::magnTot[i] << " ";
	// 	std::cout << std::endl;
	// }

	/*
	std::array<int,8> aSpin = Universe::trianglesAll.pick()->getSpins();
	std::array<int,8> oSpin = Universe::trianglesAll.pick()->getSpins();
	int dotP = dot(aSpin,oSpin);
	std::cout << "Spins 1:" << std::endl;
	for (auto s : aSpin) std::cout << s << " ";
	std::cout << std::endl << "Spins 2:" << std::endl;
	for (auto s : oSpin) std::cout << s << " ";
	std::cout << std::endl << "dot product:" << std::endl << dotP << std::endl;
	*/

	for (int i = 0; i < measurements; i++) {

		for (int i = 0; i < Triangle::nSpins*50; i++) moveSpinClusterFlip();
		sweep();

		Universe::check();
		printf("m%d, ", i);
		fflush(stdout);
	}
	printf("\b\b \n");

}

int Simulation::attemptMove() {
	std::array<int, 3> cumFreqs = {0, 0, 0};
	int freqTotal = 0;
	int prevCumFreq = 0;
	for (int i = 0; i < moveFreqs.size(); i++) {
		freqTotal += moveFreqs[i];
		cumFreqs[i] = prevCumFreq + moveFreqs[i];
		prevCumFreq = cumFreqs[i];
	}

	std::uniform_int_distribution<> moveGen(0, freqTotal-1);
	std::uniform_int_distribution<> binGen(0, 1);

	int move = moveGen(rng);

	if (move < cumFreqs[0]) {

		if (binGen(rng) == 0) {
			if (moveAdd()) return 1;
		} else {
			if (moveDelete()) return 2;
		}
		
	} else if (move < cumFreqs[1]) {
		if (moveFlip()) return 3;
	} else if (cumFreqs[1] <= move) {
		for (int i = 0; i < Triangle::nSpins; i++) moveSpinFlip();
		return 4;
	}

	return 0;
}

void Simulation::sweep() {
	std::uniform_int_distribution<> uniform_int(0, 3);

	std::array<int, 5> moves = {0, 0, 0, 0, 0};
	for (int i = 0; i < 500 * targetVolume; i++) {
		moves[attemptMove()]++;
	}

	// std::cout << "mvs " ;
	// for (auto m : moves) std::cout << m << "|" ;
	// std::cout << std::endl;

	do {
		attemptMove();
	} while (Triangle::size() != targetVolume);

	prepare();
	for (auto o : observables) {
		o->measure();
	}
}

bool Simulation::moveAdd() {
	double n0 = Vertex::size();
	// double n2 = Triangle::size();
	double n0_four = Universe::verticesFour.size();

	// double edS = exp(- 2 * lambda);
	// double rg = n2 / (n2 + 2.0);
	// double ar = edS * rg;

	double ar = n0 / (n0_four + 1.0) * exp(- 2 * lambda);
	if (targetVolume > 0) {
		double expesp = exp(2 * epsilon);
		ar *= Triangle::size() < targetVolume ? expesp : 1 / expesp;
	}

	Triangle::Label t = Universe::trianglesAll.pick();

	if (Universe::sphere) {
		if (t->time == 0) return false;
	}
	
	auto s = t->getSpins();
	auto sr = t->getTriangleRight()->getSpins();
	auto sc = t->getTriangleCenter()->getSpins();
	auto scr = t->getTriangleCenter()->getTriangleRight()->getSpins();

	std::uniform_int_distribution<> binGen(0,1);
	std::array<int,Triangle::nSpins> s2;
	for (int i = 0; i < Triangle::nSpins; i++) s2[i] = 2*binGen(rng)-1;
	std::array<int,Triangle::nSpins> sc2;
	for (int i = 0; i < Triangle::nSpins; i++) sc2[i] = 2*binGen(rng)-1;

	double dE = -isingJ*(dot(sr,s2)-dot(sr,s)+dot(scr,sc2)-dot(scr,sc)+dot(s,s2)+dot(sc,sc2)+dot(s2,sc2));

	double arsf = pow(2,2*Triangle::nSpins)*exp(-dE); 
	ar *= arsf;

	if (ar < 1.0) {
		std::uniform_real_distribution<> uniform(0.0, 1.0);
		double r = uniform(rng);
		if (r > ar) return false;
	}

	// int newS2 [8] = {s2,1,1,1,1,1,1,1};
	// int newSc2 [8] = {sc2,1,1,1,1,1,1,1}; 

	Universe::insertVertex(t,s2,sc2);
	return true;
}

				// \sl /s \sr /
				// /scl\sc/scr\

bool Simulation::moveDelete() {
	if (Universe::verticesFour.size() == 0) return false;

	double n0 = Vertex::size();
	// double n2 = Triangle::size();
	double n0_four = Universe::verticesFour.size();

	// double edS = exp(2 * lambda);
	// double rg = n2 / (n2 - 2.0);
	// double ar = edS * rg;

	double ar = n0_four / (n0 - 1.0) * exp(2 * lambda);
	if (targetVolume > 0) {
		double expesp = exp(2 * epsilon);
		ar *= Triangle::size() < targetVolume ? 1 / expesp : expesp;
	}



	Vertex::Label v = Universe::verticesFour.pick();
	// auto t = Universe::trianglesAll.pick();
	// auto v = t->getVertexLeft();
	// if (v->nUp + v->nDown != 4) return false;

	if (Universe::sliceSizes[v->time] < 4) return false;  // reject moves that shrink slices below size 3

	auto sr = v->getTriangleRight()->getSpins();
	auto sl = v->getTriangleLeft()->getSpins();
	auto scr = v->getTriangleRight()->getTriangleCenter()->getSpins();
	auto scl = v->getTriangleLeft()->getTriangleCenter()->getSpins();
	auto srr = v->getTriangleRight()->getTriangleRight()->getSpins();
	auto scrr = v->getTriangleRight()->getTriangleCenter()->getTriangleRight()->getSpins();

	double dE = -isingJ*(dot(srr,sl)-dot(srr,sr)+dot(scrr,scl)-dot(scrr,scr)-dot(sr,scr)-dot(sl,sr)-dot(scl,scr));

	double arsf = pow(2,-2*Triangle::nSpins)*exp(-dE);
	ar *= arsf;

	if (ar < 1.0) {
	std::uniform_real_distribution<> uniform(0.0, 1.0);
	double r = uniform(rng);
	if (r > ar) return false;
	}

	Universe::removeVertex(v);

	return true;
}

bool Simulation::moveFlip() {
	if (Universe::trianglesFlip.size() == 0) return false;

	auto t = Universe::trianglesFlip.pick();

	int wa = Universe::trianglesFlip.size();
	int wb = wa;
	if (t->type == t->getTriangleLeft()->type) {
		wb++;
	} else {
		wb--;
	}

	if (t->getTriangleRight()->type == t->getTriangleRight()->getTriangleRight()->type) {
	   	wb++;
	} else {
		wb--;
	}

	double ar = 1.0*wa/wb;

	if (ar < 1.0) {
		std::uniform_real_distribution<> uniform(0.0, 1.0);
		double r = uniform(rng);
		if (r > ar) return false;
	}

	Triangle::Label tr = t->getTriangleRight();
	Triangle::Label tc = t->getTriangleCenter();
	Triangle::Label trc = tr->getTriangleCenter();
	auto s = t->getSpins();
	auto sr = tr->getSpins();
	auto sc = tc->getSpins();
	auto src = trc->getSpins();
	
	if ( 1 == 1 ) { // s != sr If the spins are unequal, determine the acceptance ratio. Else, flipping will not have effect on the energy.
		double dE = -isingJ*(dot(s,src)+dot(sr,sc)-dot(s,sc)-dot(sr,src));
		if ( dE > 0 ) {
			double ar = exp(-dE);
			std::uniform_real_distribution<> uniform(0.0, 1.0);
			double r = uniform(rng);
			if (r > ar) return false;
		}
	}

	Universe::flipLink(t);

	return true;
}

bool Simulation::moveSpinFlip() {
	Triangle::Label t = Universe::trianglesAll.pick();
	std::uniform_int_distribution<> posGen(0, Triangle::nSpins-1);
	int pos = posGen(rng);

	double curSpin =  t->getSpin(pos);

	int sl = t->getTriangleLeft()->getSpin(pos);
	int sr = t->getTriangleRight()->getSpin(pos);
	int sc = t->getTriangleCenter()->getSpin(pos);
	int neighborSum = sl + sr + sc;

	double dE = isingJ*-2.0*curSpin*double(neighborSum) ;//= E_old - E_new

	if (dE < 0) {
		std::uniform_real_distribution<> uniform(0.0, 1.0);
		double r = uniform(rng);
		if (r > exp(dE)) return false;
	}

	// std::cout << pos;

	Universe::flipSpin(t,pos);
	return true;
}

double Simulation::moveSpinClusterFlip() {
	assert(spinCluster.empty());
	std::uniform_int_distribution<> posGen(0, Triangle::nSpins-1);
	int pos = posGen(rng);

	double r;
	double addProb = 1-exp(-2*isingJ);
	std::uniform_real_distribution<> uniform(0.0, 1.0);

	std::vector<Triangle::Label> newNeighbors;
	newNeighbors.push_back(Universe::trianglesAll.pick());

	auto selectNeighbors = [&] (Triangle::Label t) {
		Triangle::Label neighbors [3] = {t->getTriangleLeft(),t->getTriangleRight(),t->getTriangleCenter()};
		for (auto tn : neighbors) {
			if (tn->getSpin(pos) == t->getSpin(pos) && !tn->marked) {
				r = uniform(rng);
				if (r < addProb) {
					newNeighbors.push_back(tn);
					tn->marked = 1;
				}
			}
		}
	};

	do {
		Triangle::Label addedTriangle = newNeighbors.back();
		addedTriangle->marked = 1;
		spinCluster.push_back(addedTriangle);
		newNeighbors.pop_back();
		selectNeighbors(addedTriangle);
	} while (!newNeighbors.empty());
	

	for (auto t : spinCluster) {
		t->marked = 0;
		Universe::flipSpin(t,pos);
	}

	// std::cout << pos;

	double clusterSize = spinCluster.size();
	assert(clusterSize <= Universe::trianglesAll.size());
	spinCluster.clear();

	return clusterSize;
}

void Simulation::prepare() {
	Universe::updateVertexData();
	Universe::updateTriangleData();
	Universe::updateLinkData();
}

void Simulation::tune() {
	fflush(stdout);
	std::vector<int> volumes;
	double epsilonAim = 0.05; //0.02
	epsilon = epsilonAim; //0.02
	int ConsecutiveSmallEps = 0;

	// //initialize to expected value
	// double lambstr[10] = {0.692,0.696,0.708,0.727,0.754,0.788,0.832,0.886,0.950,1.024};
	// if ( isingJ < 1 ) {
	// 	lambda = (Triangle::nSpins+1)*lambstr[ (int) (10*isingJ) ];
	// } else {
	// 	lambda = 0.69 + Triangle::nSpins*(3*isingJ)/2;
	// }
	printf("Chose initial lambda %f,\n",lambda);

	bool done = false;
	int tuneSteps = 800; //50
	// std::array<int, 5> moves = {0, 0, 0, 0, 0};
	for (int k = 0; k < tuneSteps && !done; k++) {
		for (int i = 0; i < targetVolume; i++) {
			for (int j = 0; j < 100; j++) attemptMove(); //moves[attemptMove()]++;

			volumes.push_back(Triangle::size());
		}

		double avg = 0.0;
		for (auto v : volumes) avg += static_cast<double>(v);
		avg /= volumes.size();

		double sd = 0.0;
		for (auto v : volumes) sd += (static_cast<double>(v) - avg)*(static_cast<double>(v) - avg);
		sd /= volumes.size();

		if ((targetVolume - avg)*(targetVolume - avg) < 2*sd) {
			epsilon *= 0.7;
			if (epsilon < epsilonAim) {
				ConsecutiveSmallEps += 1;
				epsilon = epsilonAim;
				lambda -= 0.003 * (targetVolume - avg) / sqrt(sd);
			}
		} else if ((targetVolume - avg)*(targetVolume - avg) > 8*sd) {
			ConsecutiveSmallEps = 0;
			epsilon *= 1.2;
			if (epsilon > 5.0) epsilon = 5.0;
		} else if ((targetVolume - avg)*(targetVolume - avg) < 0.04*targetVolume*targetVolume) {
			ConsecutiveSmallEps = 0;
			lambda += 0.6*(avg - targetVolume)/abs((avg-targetVolume)) * epsilon;
		}
		volumes.clear();
		// if (k >= tuneSteps && abs(avg-targetVolume) < 0.1*targetVolume && epsilon < 0.021) done = true;
		if ( (k >= tuneSteps-1 || ConsecutiveSmallEps > 30 ) && abs(avg-targetVolume) < 0.1*targetVolume && epsilon <= epsilonAim) {
			done = true;
			printf("tuned to %f after %d steps.\n", lambda, k);
		}
	}

}

void Simulation::grow() {
	int growSteps = 0;
	printf("growing...\n");
	do {
		for (int i = 0; i < targetVolume; i++) attemptMove();
		fflush(stdout);
		growSteps++;
		if (growSteps > 150) {printf("Doesn't grow properly...\n"); exit(1);}
	} while (Triangle::size() < targetVolume);
	printf("grown in %d sweeps\n", growSteps);
}

void Simulation::thermalize() {
	int thermSteps = 0;
	printf("thermalizing...\n");
	fflush(stdout);
	double coordBound = log(2 * targetVolume) / static_cast<double>(log(2));
	int maxUp, maxDown;
	do { //thermalize spacetime
		for (int i = 0; i < 100 * targetVolume; i++) attemptMove();
		fflush(stdout);

		prepare();
		maxUp = 0;
		maxDown = 0;
		for (auto v : Universe::vertices) {
			int nup = 0, ndown = 0;
			for (auto vn : Universe::vertexNeighbors.at(v)) {
				if (vn->time > v->time || (v->time == Universe::nSlices-1 && vn->time == 0)) nup++;
				if (vn->time < v->time || (v->time == 0 && vn->time == Universe::nSlices-1)) ndown++;
			}
			if (nup > maxUp) maxUp = nup;
			if (ndown > maxDown) maxDown = ndown;
		}

		thermSteps++;
		if (Triangle::size() > 5*targetVolume) {printf("Isn't tuned properly...\n"); exit(1);}
	} while (maxUp > coordBound || maxDown > coordBound);

	int spinThermSteps = 0;
	if (isingJ > 0.69 && !Universe::imported) {
		spinThermSteps = 25;
		for (int i = 0; i < 50000; i++)	moveSpinClusterFlip();
		for (int i = 0; i < 20*500*targetVolume; i++) attemptMove();
	} else {
		spinThermSteps = 5;
		for (int i = 0; i < 20000; i++)	moveSpinClusterFlip();
		for (int i = 0; i < 3*500*targetVolume; i++) attemptMove();
	} //spin thermalization

	printf("thermalized in %d sweeps\n", thermSteps+spinThermSteps);
}