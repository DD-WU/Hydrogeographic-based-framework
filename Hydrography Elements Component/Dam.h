#pragma once
#include "Config.h"
#include"BasicFeature.h"
class Dam : public BasicFeature
{
public :
	// flow functions
	Water_API Dam(Config* param);
	Water_API void init()override;
	Water_API void update()override;
	Water_API void finalize()override;
	// inner functions
	Water_API void LoadWeir(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void DamQ(States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API double CalcWeirQx(int i, int j, Pars* Parptr, Arrays* Arrptr, Solver* Solverptr, States* Statesptr);
	Water_API double CalcWeirQy(int i, int j, Pars* Parptr, Arrays* Arrptr, Solver* Solverptr, States* Statesptr);
	// data exchange functions
private:
	Config* parameter;
};

