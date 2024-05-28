#pragma once
#include"Config.h"
#include"BasicFeature.h"
class Buildings : public BasicFeature
{
public:
	//flow functions
	Water_API Buildings(Config* param);
	Water_API void init()override;
	Water_API void update()override;
	Water_API void finalize()override;
	//inner functions
	Water_API void LoadDSM(Fnames*, States*, Pars*, Arrays*, int*);
	Water_API void LoadBuildingHeight(Fnames*, States*, Pars*, Arrays*, int*);
	Water_API void LoadARF(Fnames*, States*, Pars*, Arrays*, int*);
	Water_API void LoadRoofStartWaterDepth(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void LoadRoofManningsn(Fnames* Fnameptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void BuildingsQ(States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API double CalcARFQx(int i, int j, States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr, double* TSptr);
	Water_API double CalcARFQy(int i, int j, States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr, double* TSptr);

	//data exchange functions
	Water_API void set_inflow_lisflood(int j, double value);
	Water_API int get_gird_x(int j);
	Water_API int get_gird_y(int j);
	Water_API void set_waterdepth(int row, int col, double t);
	Water_API double get_waterdepth(int row, int col);
	Water_API double* get_BCVar(int j);
	Water_API double get_dx();
	Water_API map<string, int> get_building_pipe_index();
private:
	Config* paramater;
	map<char*, char*> building_pipe_index;
};


