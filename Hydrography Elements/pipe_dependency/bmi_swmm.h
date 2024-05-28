/* -*- c-file-style: "stroustrup" -*- */
/* Please use the stroustrup coding standard: */


#ifndef BMI_SWMM_API_H
#define BMI_SWMM_API_H

#define BMI_API_VERSION_MAJOR 1
#define BMI_API_VERSION_MINOR 0

#if defined _WIN32
#define BMI_API __declspec(dllexport) __stdcall
/* Calling convention, stdcall in windows, cdecl in the rest of the world */
#define CALLCONV __stdcall
#else
#define BMI_API
#define CALLCONV
#endif


#define MAXSTRINGLEN 1024
#define MAXDIMS 6
#include <stddef.h>
#include <string>
#include <vector>
#include <map>
using std::vector;
using std::string;
using std::map;

#ifdef __cplusplus
extern "C" {
#endif


    /* process interface.*/
    int BMI_API  initialize_swmm(vector<string> config_file);
    int BMI_API  update_swmm(double dt);
    int BMI_API  finalize_swmm();


    /* data interface */
    double BMI_API  get_current_time_swmm(double *t);
    double BMI_API  get_TotalDuration();
    double BMI_API  get_NewRoutingTime();
    void BMI_API set_inflow_swmm(double* a);
    void BMI_API get_overflow_swmm(double* b);
    __declspec(dllexport) char* __stdcall get_pipe_name_swmm(char* inpfile);
   int  BMI_API get_pipe_number_swmm();
   double BMI_API get_pipe_elevation_swmm(int j);
    void BMI_API setOutletDepth_swmm(int j, double z);
    bool BMI_API isOutlet_swmm(int j);
    double BMI_API get_Link_Flow_swmm(int j);
    int BMI_API get_upstream_Link_swmm(char* ps_name);
    int BMI_API get_Node_Id_swmm(char* ps_name);
    double BMI_API get_position_node_x_swmm(int j);
    double BMI_API get_position_node_y_swmm(int j);
#ifdef __cplusplus
}
#endif


#endif
