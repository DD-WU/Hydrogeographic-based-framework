#pragma once
#include"Config.h"
#include"BasicFeature.h"
#ifdef Water_API//宏重定义解决办法
#undef Water_API
#endif

#if defined _WIN32
#define Water_API __declspec(dllexport)
/* Calling convention, stdcall in windows, cdecl in the rest of the world */
#define CALLCONV __stdcall
#else
#define Water_API
#define CALLCONV
#endif


#define MAXSTRINGLEN 1024
#define MAXDIMS 6
struct riverBoudaryIndex
{
	int i = -9999;
	int j = -9999;
};
class River : public BasicFeature
{
public:
	// flow functions
	Water_API River(Config*);
	Water_API void init();
	Water_API void update();
	Water_API void finalize();
	
	//inner functions
	Water_API void LoadRiverNetwork(Fnames* Fnameptr,States* Statesptr,Pars* Parptr,vector<ChannelSegmentType>* ChannelSegmentsVecPtr,
		Arrays* Arrptr,vector<QID7_Store>* QID7_Vec_Ptr,vector<int>* RiversIndexVecPtr,int* verbose);
	Water_API void LoadRiver(Fnames* Fnameptr, States* Statesptr,Pars* Parptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr,
		Arrays* Arrptr, vector<QID7_Store>* QID7_Vec_Ptr, vector<int>* RiversIndexVecPtr, int* verbose);
	Water_API void Load_river_pipe_index(Fnames* Fnameptr, States* Statesptr);
	Water_API void UpdateChannelsVector(States*, ChannelSegmentType*, vector<QID7_Store>*, QID7_Store*, int*); // CCS
	Water_API void ChannelQ_Diff(double deltaT, States*, Pars*, Solver*, BoundCs*, ChannelSegmentType*, Arrays*, vector<int>*, int*);
	Water_API void ChannelQ_Diff1(double deltaT, States*, Pars*, Solver*, BoundCs*, ChannelSegmentType*, Arrays*, vector<int>*, int*);
	Water_API void calcJ(double* x, double* xn, double** J, double dt, ChannelSegmentType* csp, int chseg, int HoutFREE);
	Water_API void bandec(double** a, int n, int m1, int m2, double** al, int indx[], double& d);
	Water_API void banbks(double** a, int n, int m1, int m2, double** al, int indx[], double b[]);
	Water_API void SWAP(double& a, double& b);
	Water_API double BankQ(int chani, ChannelSegmentType*, Pars*, Arrays*);
	Water_API void calcF(double* x, double* xn, double* f, double dt, ChannelSegmentType* csp, 
		Pars* Parptr, Arrays* Arrptr, double Qin, int chseg, double WSout, int HoutFREE, Solver* Solverptr, int low);
	Water_API int signR(double a);
	Water_API void SetChannelStartH(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, vector<int>*, int*);
	Water_API void SetChannelStartHfromQ(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, Solver*, vector<int>*, int*);
	Water_API void CalcChannelStartQ(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, vector<int>*, int*);
	Water_API double CalcEnergySlope(double n, double w, double h, double Q);

	//Data exchange functions
	Water_API char* get_river_name(int i);
	Water_API int get_river_chsz(int i);
	Water_API int get_river_count();
	Water_API int get_river_index_size();
	Water_API riverBoudaryIndex* get_river_index(string rivername); 
	Water_API map< string, riverBoudaryIndex*> get_river_pipe_index();
	Water_API double get_river_width(riverBoudaryIndex river);
	Water_API double get_river_elevation(riverBoudaryIndex river);
	Water_API double get_river_H(riverBoudaryIndex river);
	Water_API riverBoudaryIndex* get_boundary(char* PSNAME);
	Water_API void set_river_inflow(double depth, riverBoudaryIndex name);

private:
	Config* parameter;
	map< string,riverBoudaryIndex*> river_pipe_index;
};

