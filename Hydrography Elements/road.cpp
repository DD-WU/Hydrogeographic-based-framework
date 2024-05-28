#include "road.h"

Road::Road(Config* param)
{
    this->parameter = param;
}

void Road::init()
{
    LoadRoadDepth(parameter->Fnameptr, parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->verbose);
}

void Road::update()
{
}

void Road::finalize()
{
}

void Road::LoadRoadDepth(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j, gr;
    char dum[800];
    double no_data_value = -9999;

    fp = fopen(Fnameptr->Roadname, "r");
    if (fp == NULL)
    {
        if (*verbose == ON) printf("\nWARNING: Unable to load startfile:\t%s\t\n", Fnameptr->Roadname);
        return;
    }

    printf("Loading initial depths:\t%s\t", Fnameptr->Roadname);


    for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
    fscanf(fp, "%s %lf", dum, &no_data_value);

    for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
    {
        double temp = 0;
        fscanf(fp, "%lf", &temp);
        if (temp == no_data_value) {

        }
        else
        {
            if (Arrptr->DEM[i + j * Parptr->xsz] > 0) {//这个主要是因为fenhu这个案例道路数据不太准，其他区域把这个if删了就可以了
                Arrptr->DEM[i + j * Parptr->xsz] += temp;
            }
        }
    }
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");

    return;
}
