
#include <cstdio>
#include <string>
#include <sstream>
//#include <string>
//#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "bmi_swmm.h"

extern "C" {
    #include "headers.h"
    #include "swmm5.h"
}


using namespace std;


enum VARIABLE_LABEL {
    H,
    DEM,
    Qx,
    Qy,
    rain,
    dx,
    dx_sqrt,
    dy,
    dA,
    SGCwidth,
    SGCQin,
    SGCz,
    blx,
    bly,
    BC_numPS,
    BC_xpi,
    BC_ypi,
    QxSGold,
    QySGold,

    //lmy start
    Sim_Time,
    killsim,
    itrn_time_now,
    killsim_time,
    steadycheck,
    t,
    steadyTotal,
    //steadyCount,
    // lmy end
    // if variable is not found
    elapsed_time,
    NOVAR
};

VARIABLE_LABEL variable2label(std::string const& variable) {
    if (variable == "H") return H;
    if (variable == "DEM") return DEM;
    if (variable == "Qx") return Qx;
    if (variable == "Qy") return Qy;
    if (variable == "rain") return rain;
    if (variable == "dx") return dx;
    if (variable == "dx_sqrt") return dx_sqrt;
    if (variable == "dy") return dy;
    if (variable == "dA") return dA;
    if (variable == "SGCwidth") return SGCwidth;
    if (variable == "SGCQin") return SGCQin;
    if (variable == "SGCz") return SGCz;
    if (variable == "blx") return blx;
    if (variable == "bly") return bly;
    if (variable == "BC.numPS") return BC_numPS;
    if (variable == "BC.xpi") return BC_xpi;
    if (variable == "BC.ypi") return BC_ypi;
    if (variable == "QxSGold") return QxSGold;
    if (variable == "QySGold") return QySGold;
    //lmy s
    if (variable == "Sim_Time") return Sim_Time;
    if (variable == "killsim") return killsim;
    if (variable == "itrn_time_now") return itrn_time_now;
    if (variable == "killsim_time") return killsim_time;
    if (variable == "steadycheck") return steadycheck;
    if (variable == "t") return t;
    if (variable == "steadyTotal") return steadyTotal;

    if (variable == "elapsed_time") return elapsed_time;
    //if (variable == "steadyCount") return steadyCount;
    //lmy e
    //if (variable == "FArea") return FArea;
    else return NOVAR;
};

/* Store callback */
//Logger logger = NULL;

/* Logger function */
//void _log(Level_swmm level, std::string msg);

extern "C" int init(int, const char* []);
extern "C" int init_iterateq();

//功能函数，string转double
double stringToNum(const string& str)
{
    istringstream iss(str);
    double num;
    iss >> num;
    return num;
}

void trim(string & s)
{
    int index = 0;
    if (!s.empty())
    {
        while ((index = s.find(' ', index)) != string::npos)
        {
            s.erase(index, 1);
        }
    }    
}


int BMI_API  initialize_swmm(vector<string> config_file)
{
    int result = swmm_open((char*)config_file[0].c_str(), (char*)config_file[1].c_str(), (char*)config_file[2].c_str());
    swmm_start(1);
    return result;
}

int BMI_API update_swmm(double dt)
{
    // set timestep and do:
    double oldStep = dt;
    //if (dt != -1) {
    //    Solverptr->Tstep = dt;
    //}
    //iterateq_step();
    //// restore dt to default
    //Solverptr->Tstep = oldStep;
    swmm_step(&oldStep);
    dt = oldStep;
    ElapsedTime = oldStep;
    return 1;
}

int BMI_API  finalize_swmm()
{
    swmm_end();
    swmm_report();
    swmm_close();
    return 0;
}
double BMI_API get_current_time_swmm(double *dt)
{
    if (dt != nullptr)
    {
        *dt = ElapsedTime;
    }
    return ElapsedTime;
}

double BMI_API get_TotalDuration()
{
    return TotalDuration;
}

double BMI_API get_NewRoutingTime()
{
    return NewRoutingTime;
}

__declspec(dllexport) char* __stdcall get_pipe_name_swmm(char* inpfile)
{
    FILE* fp_swmm;
    char line[100] = "";
    char* PS_name = new char[5000 * 80];
    float  xc, yc;
    int p_number = 0;
    int i, j = 0;
    int* numps = NULL;
    fp_swmm = fopen(inpfile, "r");
    if (fp_swmm == NULL) { printf("no inp files for swmm5\n"); exit(0); }
    while (fgets(line, 100, fp_swmm))
    {
        line[strlen(line) - 1] = 0;

        if (strcmp("[COORDINATES]", line) == 0)                    //[COORDINATES]  nodes  outfalls 
        {
            fgets(line, 100, fp_swmm);
            fgets(line, 100, fp_swmm);
            for (i = 0; i < 10000; i++)     //  NOT p_number
            {
                fscanf(fp_swmm, " %s %f  %f ", PS_name + i * 80, &xc, &yc);
                if (strcmp("[VERTICES]", PS_name + i * 80) == 0)  break;          // may changed according to inp file         
            }
        }
    }
    return PS_name;
}
int BMI_API  get_pipe_number_swmm()
{
    return Nobjects[NODE];
}
double BMI_API get_pipe_elevation_swmm(int j)
{
    return Node[j].invertElev * UCF(LENGTH);
}
void BMI_API setOutletDepth_swmm(int j, double z) {
    if (Node[j].type != OUTFALL) return;
    double   x, y;                     // x,y values in table
    double   stage;                    // water elevation at outfall (ft)
    int      k;                        // table index
    int      i = Node[j].subIndex;     // outfall index
    if (Outfall[i].type!= TIDAL_OUTFALL)
    {
        return;
    }
    k = Outfall[i].tideCurve;
    table_getFirstEntry(&Curve[k], &x, &y);
    Curve[k].thisEntry->y = z * UCF(LENGTH);
    while (table_getNextEntry(&Curve[k], &x, &y))
    {
        Curve[k].thisEntry->y = z * UCF(LENGTH);

    }
    table_tseriesInit(&Curve[k]);
}
 double __stdcall get_Link_Flow_swmm(int j) {
     return Link[j].newFlow;
}
int __stdcall get_upstream_Link_swmm(char*ps_name)
{
    for (int i = 0; i < Nobjects[LINK]; i++)
    {
        int n1 = Link[i].node1;
        int n2 = Link[i].node2;
        char* id1 = Node[n1].ID;
        char* id2 = Node[n2].ID;
        if (strcmp(id1, ps_name) == 0 || strcmp(id2, ps_name) == 0)
        {
            return i;
        }
    }
    return -1;
}
int BMI_API get_Node_Id_swmm(char* ps_name)
{
    for (int i = 0; i < Nobjects[NODE]; i++)
    {
        char* id1 = Node[i].ID;
        if (strcmp(id1, ps_name) == 0) {
            return i;
        }
    }
    return -1;
}
bool __stdcall isOutlet_swmm(int j)
{
    if (Node[j].type != OUTFALL) return false;
    return true;
}
double __stdcall get_position_node_x_swmm(int j)
{
    return Node[j].X;
}
double __stdcall get_position_node_y_swmm(int j)
{
    return Node[j].Y;
}
void BMI_API set_inflow_swmm(double* a) {
    int     j, p;
    double  q, w;
    TExtInflow* inflow;
    DateTime currentDate = getDateTime(NewRoutingTime);    
    TTableEntry* P;
    int type;
    int tseries = -1;
    // --- for each node with a defined external inflow
    for (j = 0; j < Nobjects[NODE]; j++)
    {
        inflow = Node[j].extInflow;
        if (!inflow) continue;

        // --- get flow inflow
        q = 0.0;
        while (inflow)
        {
            type = inflow->type;
            tseries = inflow->tSeries;
            if (inflow->type == FLOW_INFLOW)
            {
                double x = 0;
                double y = 0;
                table_getFirstEntry(&Tseries[tseries], &x, &y);
                Tseries[tseries].thisEntry->y = a[j];
                while (table_getNextEntry(&Tseries[tseries], &x, &y))
                {
                    Tseries[tseries].thisEntry->y = a[j];
                    
                }
                table_tseriesInit(&Tseries[tseries]);
                q = inflow_getExtInflow(inflow, currentDate);
                break;
            }
            else inflow = inflow->next;
        }
        if (fabs(q) < FLOW_TOL) q = 0.0;

        // --- add flow inflow to node's lateral inflow
        Node[j].newLatFlow += q;
        /*massbal_addInflowFlow(EXTERNAL_INFLOW, q);*/

        // --- add on any inflow (i.e., reverse flow) through an outfall
        if (Node[j].type == OUTFALL && Node[j].oldNetInflow < 0.0)
        {
            q = q - Node[j].oldNetInflow;
        }    
    }
}

void BMI_API get_overflow_swmm( double* b)
{
    int j;
    for (j = 0; j < Nobjects[NODE]; j++)
    {
        if (Node[j].type == OUTFALL) continue;
        b[j] = Node[j].overflow * UCF(FLOW);
    }
}
