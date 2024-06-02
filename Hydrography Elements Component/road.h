#pragma once
#include "Config.h"
#include"BasicFeature.h"
class Road : public BasicFeature
{
public:
	// flow functions
	Water_API Road(Config* param);
	Water_API void init()override;
	Water_API void update()override;
	Water_API void finalize()override;
	// inner functions
	Water_API void LoadRoadDepth(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
private:
	Config* parameter;
};

