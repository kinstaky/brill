#ifndef __PARTCILE_EVENT_H__
#define __PARTICLE_EVENT_H__

#include <string>
#include <TTree.h>

namespace brill {

struct ParticleEvent {
	int num;
	int charge[8];
	int mass[8];
	double energy[8];
	double time[8];
	double x[8];
	double y[8];
	double z[8];
	double px[8];
	double py[8];
	double pz[8];
	bool stop[8];
	int last[8];
};

void SetupInput(TTree *tree, ParticleEvent &event, const std::string &prefix="");
void SetupOutput(TTree *tree, const ParticleEvent &event);
void Reset(ParticleEvent &event);
}

#endif