#include "Lake.h"

Lake::Lake(Config* param)
{
	this->parameter = param;
}

void Lake::init()
{
    LoadLakeStartWaterDepth(parameter->Fnameptr, parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->verbose);
}

void Lake::update()
{
}

void Lake::finalize()
{
}

void Lake::LoadLakeStartWaterDepth(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j, gr;
    char dum[800];
    double no_data_value = -9999;

    fp = fopen(Fnameptr->Lakestartfilename, "r");
    if (fp == NULL)
    {
        if (*verbose == ON) printf("\nWARNING: Unable to load startfile:\t%s\t\n", Fnameptr->Lakestartfilename);
        return;
    }

    if (*verbose == ON) printf("Loading initial depths:\t%s\t", Fnameptr->Lakestartfilename);


    for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
    fscanf(fp, "%s %lf", dum, &no_data_value);

    for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
    {
        double temp=0;
        fscanf(fp, "%lf", &temp);
        if (temp == no_data_value) {

        }
        else
        {
            Arrptr->H[i + j * Parptr->xsz] = temp;
        }
    }
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");

    return;
}
