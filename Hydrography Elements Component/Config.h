#pragma once
#include <time.h>
#include <iostream>
#include "BasicExcel.hpp"
using namespace YExcel;

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

// define basic constants
#define ON 1
#define OFF 0
#define CHKINTERVAL 1.0 // default checkpoint interval in hours runtime
#define NULLVAL -9999.0 // MT: define ascii file NULL value as constant

// limits
#define MAXPI 20000

// MT: disable visual C++ depreciation warnings so we can see other warnings
#pragma warning( disable : 4996)
#include"tool.h"

struct Arrays {
	/*! DEM, Water height, Flow in x-direction and Flow in y-direction */
	double* DEM; // Digital elevation model
	double* DSM; // 当ARF模式时，我们用到dsm和dem两个东西
	int* RiverFlag;
	double* ARF;// Digital elevation model
	double* H,*temp_H, * Qx, * Qy, * Qxold, * Qyold, * U, * V, * Hflowx, * Hflowy;
	/* TRENT additions  */
	double* HU, * HV, * RSHU, * LSHU, * RSHV, * LSHV, * BSHU, * TSHU, * BSHV, * TSHV, * FHx, * FHUx, * FHVx, * FHy, * FHUy, * FHVy;
	/* Flow direction map for Rainfall*/
	int* FlowDir; // CCS: added to hold DEM flow direction map for routing shallow rainfall flow 13/03/2012
	double* Route_dH, * RouteInt; // CCS: added to record routing scheme dH and interval

	/* ---------------- */
	double* maxH, * maxHtm, * initHtm, * totalHtm;
	double* Manningsn, * SGCManningsn;
	double* paerial, * pbound;
	double* Weir_hc, * Weir_Cd, * Weir_m, * Weir_w;
	int* Weir_Identx, * Weir_Identy, * Weir_Fixdir, * Weir_Typ;
	double* evap, * rain;
	int* ChanMask, * SegMask;
	double* TRecx, * TRecy; // MT: add to record TStep
	double* LimQx, * LimQy; // MT: add to record Qlimits
	double* Vx, * Vy, * maxVx, * maxVy, * Vc; // JCN: added to record velocity
	double* maxVc, * maxVcH, * maxHaz; // JCN added to calculate hazard
	double* SGCwidth, * SGCz, * QxSGold, * QySGold, * SGCbfH, * SGCVol, * SGCdVol, * SGCbfV, * SGCc, * SGCFlowWidth, * SGCdx, * SGCcat_area; // JCN added to store widths and depths
	double* SGCQin; //JMH
	double* dx, * dy, * dA; // CCS added for lat long data
	int* SGCgroup;
};
//-------------------------------------------
// Files
struct Files {
	FILE* mass_fp;
	FILE* stage_fp;
	FILE* vel_fp;
	FILE* gau_fp;
};

struct Fnames {
//building file name
	char ARFname[80];
	char DSMname[80];
	char BuildingStartWaterName[80];
	char BuildingManningn[80];
	char BuildingheightName[80];
//building file name
//Lake file name
    char Lakename[80];
	char Lakestartfilename[80];
//Lake file name
//Road file name
    char Roadname[80];
//Road file name
//river
	char River2dFlagName[80];
	char demfilename[80];
	char startfilename[80];
	char chanfilename[80];
	//char Qfilename[80]; Not used! JCN
	char resrootname[255];
	char dirrootname[255];
	char qfilename[80];
	char nfilename[80];
	char SGCnfilename[80];
	char porfilename[80];
	char rivername[80];
	char river_pipe_index_name[80];
	char land_pipe_index_name[80];
	char land_river_index_name[80];
	char bcifilename[80];
	char bdyfilename[80];
	char weirfilename[80];
	char opfilename[80];
	char stagefilename[80];
	char ascheaderfilename[80];
	char multiriverfilename[80]; // CCS
	char checkpointfilename[255]; // used to write a checkpoint file (note it is placed in the input file dir not the results dir)
	char loadCheckpointFilename[255]; // explicit checkpoint file to start run (specify with -loadcheck option, defaults to checkpointfilename)
	char evapfilename[255];
	char rainfilename[255];
	char logfilename[255];
	char SGCwidthfilename[255]; // JN sub grid channel widths
	char SGCbankfilename[255]; // JN sub grid channel bank elevations
	char SGCbedfilename[250];  // JN sub grid channel bed elevation
	char SGCcat_areafilename[250]; // JN sub grid channel accumulation area
	char SGCchangroupfilename[250];
	char SGCchanpramsfilename[250];
	char gaugefilename[255];
};
//-------------------------------------------
// Boundary Conditions
struct BoundCs {
	int* xpi;
	int* ypi;
	int* PS_Ident;
	int   numPS;
	double* PS_Val;
	double* PS_qold;
	double* PS_qSGold;
	char* PS_Name;
	int* BC_Ident;
	double* BC_Val;
	double** BCVarlist;
	char* BC_Name;
	double Qpoint;
	double Qin;
	double Qout;
	double QChanOut;
	double VolInMT; // added by JCN stores volume in over mass inteval
	double VolOutMT; // added by JCN stores volume out over mass inteval
};
//-------------------------------------------
// Stage
struct Stage {
	int Nstages, Ngauges;
	double* stage_loc_x, * stage_loc_y;
	double* gauge_loc_x, * gauge_loc_y, * gauge_dist;
	int* stage_grid_x, * stage_grid_y, * vel_grid_xy, * stage_check;
	int* gauge_grid_xy, * gauge_grid_x, * gauge_grid_y;
	int* gauge_dir, * gauge_cells;
};

// SGC parameters
struct SGCprams {
	int NSGCprams;
	int* SGCchantype;
	double SGCbetahmin;
	double* SGCp, * SGCr, * SGCs, * SGCn, * SGCm, * SGCa;
	double* SGCgamma, * SGCbeta1, * SGCbeta2, * SGCbeta3, * SGCbeta4, * SGCbeta5;
};

//-------------------------------------------
// Simulation States
struct States {
	int ChannelPresent;
	int river2d_couple;
	int river_couple;
	int land_couple;
	int buildingHeight;
	int TribsPresent;
	int NCFS;
	int save_depth;
	int save_elev;
	int out_dir;
	int single_op;
	int multi_op;
	int calc_area;
	int calc_meandepth;
	int calc_volume;
	int save_stages;
	int adaptive_ts;
	int acceleration; // PB: Flag to switch to acceleration formulation
	int qlim; // TJF: Flag for qlim version
	int debugmode;
	int save_Qs;
	int calc_infiltration;
	int call_gzip;
	int alt_ascheader;
	int checkpoint;
	int checkfile;
	int calc_evap;
	int rainfall; // TJF: added for time varying, spatially uniform rainfall
	int routing; // CCS: added for routing routine.
	int reset_timeinit;
	int profileoutput;
	int porosity;
	int weirs;
	int save_Ts;   // MT: added flag to output adaptive timestep
	int save_QLs;  // MT: added flag to output Qlimits
	int diffusive; // MT: added flag to indicate wish to use diffusive channel solver instead of default kinematic
	int startq;    // MT: added flag to indicate wish to use start flow to calculate initial water depths throughout channel
	int logfile;   // MT: added flag to record logfile
	int startfile; // MT: added flag to note use of startfile
	int start_ch_h; // MT: added flag to note use of starting H in channel
	int comp_out; // TJF: added to make computational output information optional
	int chainagecalc; // MT: added so user can switch off grid independent river chainage calculation
	int mint_hk; // JN: added to request maxH, maxHtm totalHtm and initHtm be calulated at the mass interval
	int Roe; // JN/IV: added to use Roe solver
	int killsim; // MDW: added to flag kill of simulation after specified run time
	int dhoverw; // TJF: added as a switch for dhlin (ON - dhlin set by command line/parfile; OFF - dhlin prescribed by gradient 0.0002 Cunge et al. 1980)
	int drychecking; //JN Option to turn DryCheck off
	int voutput; // exports velocity esimates based on Q's of Roe velocity (JCN)
	int steadycheck; // MDW: added flag to check for model steady state
	int hazard; // JN additional module for calculating hazards
	int startq2d; // JN: initalises inertial model with uniform flow for Qold
	int Roe_slow; // JN: ghost cell version of Roe solver
	int multiplerivers; // CCS multiple river switch
	int SGC; // JN sub gird channels r
	int SGCbed; // JN sub grid bed elevation file to override hydraulic geometry
	int SGCcat_area; // JN sub grid channel accumulated area to override hydraulic geometry based on width
	int SGCchangroup; // turns on distributed channe groups
	int SGCchanprams; // parameters for distributed channel groups
	int binary_out; // JN binary raster output
	int gsection; // JN virtual gauge sections
	int binarystartfile; // JN load a binary start file
	int startelev; // used to use an elevation file for the startfile
	int latlong; // CCS: added for lat-long coordinate systems
	int SGCbfh_mode; // JCN switches model to use parameter p as bank full depth
	int SGCA_mode; // JCN switches model to use parameter p as bank full Area
	int dist_routing; // JCN turnes on spatially distributed routing velocity
	int SGCvoutput; // JCN Turns on sub-grid channel velocity output
	//building related
	int ARFFlag;
	int BuildingFlag;
	int BuildingStartWaterFlag;
	int BuildingManningFlag;
	//building related
};

//-------------------------------------------
// Model Parameters
struct Pars {
	int xsz, ysz;
	double dx, dx_sqrt;
	double dy, dA;
	double FPn;
	double tlx, tly, blx, bly;
	double SaveInt, MassInt;
	double SaveTotal, MassTotal;
	int SaveNo;
	int op_multinum;
	double* op_multisteps;
	int* op_multiswitch;
	double op;
	double InfilRate;
	double InfilLoss, EvapLoss, RainLoss; // previous mass interval loss
	double InfilTotalLoss, EvapTotalLoss, RainTotalLoss; // cumulative loss
	double checkfreq, nextcheck;
	double reset_timeinit_time;
	double maxelev, zlev; // Water depth dependent porosity
	int zsz; // Water depth dependent porosity
	int Por_Ident;
	double dAPor;
	char** ascheader;
	double ch_start_h; // starting depth of channel flow. default to 2m or read from par file.
	double killsim_time; // time to kill simulation
	double steadyQdiff, steadyQtol, steadyInt, steadyTotal; // used for checking steady-state
	double SGC_p; // sub grid channel width depth exponent
	double SGC_r; // sub grid channel width depth mutiplier
	int SGCchan_type; // JCN bank slop for trapazoidal channel
	double SGC_s, SGC_2, SGC_n; // JCN trapazodal channel slope dependent constant
	double* SGCprams; // pointer to table of SGC parameters
	double Routing_Speed, RouteInt; // CCS variables controlling routing speed in rainfall routing routine
	double RouteSfThresh; // CCS water surface slope at which routing scheme takes over from shallow water eqn is SGC mode (when routing==ON).
	double SGC_m, SGC_a; // allows a meander coefficient to be set for the sub-grid model, default 1, allows channel upstream area to be set, defaul -1;
	double min_dx, min_dy, min_dx_dy; // CCS added to hold minimum values of dx and dy when using lat-long projected grids.
};

// Solver settings
struct Solver {
	double ARF;
	double t;
	double g;
	double divg;
	double cfl;
	int    ts_multiple; // channel timestep multiple for running 1D decoupled from 2D
	long   Nit, itCount;
	double Sim_Time;
	double InitTstep; // Maximum timestep
	double Tstep;  // Adapting timestep
	double MinTstep;  // Stores minimum timestep during simulation
	double SolverAccuracy;
	int dynsw; // Switch for full dynamic steady state or diffusive steady state
	double Impfactor;
	double Hds;
	double vol1, vol2;
	double Qerror;
	double Verror;
	double FArea; // Store flooded area
	double DepthThresh, MomentumThresh, MaxHflow;
	double dhlin;
	double htol;
	double Qlimfact; // MT added to allow user to relax Qlimit
	double itrn_time;
	double itrn_time_now;
	double SGCtmpTstep; // JCN added to enable time step calculatin in UpdateH for SGC method
	time_t time_start;
	time_t time_finish;
	time_t time_check;
	double theta; //GAMA added for q-centred scheme
	int fricSolver2D; //GAMA: Solves the friction term using the vectorial (2D) scheme
};

//-------------------------------------------
// ChannelSegmentType
struct ChannelSegmentType {
	double* Chandx;
	double* Shalf;
	double* Chainage;
	double* ChanQ; // only for recording Q for output in profile
	double* A;
	double* NewA;
	double* ChanWidth;
	double* ChanN;
	int* ChanX;
	int* ChanY;
	double* Q_Val;
	double** QVarlist;
	int* Q_Ident;
	char* Q_Name;
	double* BankZ;
	int chsz;
	int Next_Segment;
	int Next_Segment_Loc;
	int N_Channel_Segments;
	double JunctionH; // allows recording of H data for dummy junction point of tributary, without overwriting main channel info
	double JunctionDEM; // allows recording of DEM data for dummy junction point of tributary, without overwriting main channel info
};

//-------------------------------------------
/* QID7_Store // CCS for temp storage of trib boundary condition info when (Q_Ident_tmp[i]==7) in LoadRiver. A vector containing these
   structures is built in LoadRiver function and used in UpdateChannelsVector function. */
struct QID7_Store {
	int trib;
	int Next_Segment_Loc;
	int chseg;
	int RiverID;
};
//-------------------------------------------


class Excel
{
public:
	BasicExcelWorksheet* read_excel(const char* excel_name, const char* sheet_name);
	void excel_close();
	BasicExcelWorksheet* sheet = NULL;
private:
	BasicExcel* xls;
};
class Config
{
public:
	// Work flow function
	Water_API Config(const char* config, const char* sheet);
	Water_API void init();
	Water_API void update_time();
	Water_API void finalize();

	// inner function
	Water_API void readParamFile_lisflood(const char* file, const char* sheet);
	Water_API void LoadDEM(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void LoadRain(Fnames* Fnameptr, Arrays* Arrptr, int* verbose);
	Water_API void LoadEvap(Fnames*, Arrays*, int*);
	Water_API void BCs(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr);
	Water_API void DryCheck(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API void FPInfiltration(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API void Evaporation(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API void Rainfall(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr);
	Water_API int fexist(const char* filename);
	Water_API int init_lisflood();
	Water_API double DomainVol(States* Statesptr, Pars* Parptr, ChannelSegmentType* ChannelSegments, Arrays* Arrptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr);
	Water_API void write_regular_output(Fnames* Fnameptr, Solver* Solverptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr);
	Water_API void fileoutput(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr);
	Water_API void write_ascfile(char* root, int SaveNumber, char* extension, double* data, double* dem, int outflag, States* Statesptr, Pars* Parptr);
	Water_API void write_binrasterfile(char* root, int SaveNumber, char* extension, double* data, double* dem, int outflag, States* Statesptr, Pars* Parptr);
	Water_API void UpdateH(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr);
	Water_API void UpdateV(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr);	
	Water_API void BoundaryFlux(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr);
	Water_API double PorArea(int i, int j, Pars* Parptr, Arrays* Arrptr);

	// Exchange function
	Water_API double get_current_time();
	Water_API double get_sim_time();
	Water_API int get_PSNum();
	Water_API int* get_BCPSNum();
	Water_API char* get_PSName();
	Water_API void set_verbose();
	Water_API double get_save_time();

	// LISFLOOD-FP Structure
	States* Statesptr;
	Pars* Parptr;
	Fnames* Fnameptr;
	Solver* Solverptr;
	BoundCs* BCptr;
	Stage* Stageptr;
	ChannelSegmentType* CSTypePtr;
	vector<int>* RiversIndexVecPtr;
	int* RiversIndexPtr;
	vector<ChannelSegmentType>* ChannelSegmentsVecPtr;
	ChannelSegmentType* ChannelSegments;
	Arrays* Arrptr;
	char* para[3];//for river2d
	int* verbose;
	char t1[80];
	char tmp_sys_com[255];
	Files Fps;
	double Previous_t;      // previous time channel was calculated
	int steadyCount;
	int tstep_counter;   // start at -1 so that in first run through we calculate river
	double tstep_channel; // channel timestep
	const char* sheet_name;
	const char* config_name;
	bool init_flag;
	Excel* excel;
};
#include"pipe_dependency\swmm5.h"
#include"pipe_dependency\bmi_swmm.h"
class Config_pipe
{
public:
	Water_API Config_pipe(char* config_file, char* sheet);
	Water_API void init();
	Water_API void update_time();
	Water_API void finalize();

	//inner functions
	Water_API void readParamFile_SWMM( char* file,  char* sheet);
	
	//data exchange functions
	Water_API double get_current_time();
	Water_API double get_Elapsed_Time();
	Water_API double get_Total_Duration();
	Water_API double get_New_Routing_Time();
	Water_API void  set_inflow(double* a);
	Water_API void get_overflow(double* b);
	Water_API char* get_pipe_name();
	Water_API int get_pipe_number();
	Water_API double get_pipe_elevation(int j);
	Water_API void setOutletDepth(int j, double z);
	Water_API bool isOutlet(int j);
	Water_API double get_Link_Flow(int j);
	Water_API int get_upstream_Link(char* ps_name);
	Water_API double get_position_node_x(int j);
	Water_API double get_position_node_y(int j);
	Water_API bool isRiverCouple();
	Water_API bool isLandCouple();
	Water_API char* get_river_pipe_index_name();
	Water_API char* get_land_pipe_index_name();
private:
	Excel* excel;
	const char* sheet_name;
	const char* config_name;
	char* inpfile, * outfile, * rptfile;
	int*  OpenFlag;           // TRUE if a project has been opened
	int* StartedFlag;        // TRUE if a simulation has been started
	int* ResultsFlag;      // TRUE if output to be saved to binary file
	int* Exception_Count;       // number of exceptions handled
	int* Runoff;             // TRUE if runoff is computed
	int* Routing;            // TRUE if flow routing is computed
	double* current_time;
	double* Elapsed_Time;
	int rivercouple=0;
	int landcouple = 0;
	char river_pipe_index_name[80];
	char land_pipe_index_name[80];
};


