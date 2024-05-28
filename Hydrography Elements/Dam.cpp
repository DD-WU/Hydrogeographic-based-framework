#include "Dam.h"

Dam::Dam(Config* param)
{
    this->parameter = param;
    parameter-> Statesptr->weirs = ON;
}

void Dam::init()
{
    LoadWeir(parameter->Fnameptr, parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->verbose);
}

void Dam::update()
{
    DamQ(parameter->Statesptr, parameter->Parptr, parameter->Solverptr, parameter->Arrptr);
}

void Dam::finalize()
{
}

void Dam::LoadWeir(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    //One-directional flow functionality added by Matt Wilson 13 Feb 2504

    FILE* fp;
    int nw, i, j, xi, yi, p0, p1;
    double x, y, z0, z1;
    char char_tmp[10], buff[80];
    char tag_w1[] = "W";
    char tag_w2[] = "w";
    char tag_e1[] = "E";
    char tag_e2[] = "e";
    char tag_s1[] = "S";
    char tag_s2[] = "s";
    char tag_n1[] = "N";
    char tag_n2[] = "n";
    //tags for one-directional flow (culverts)
    char tag_wf1[] = "WF";
    char tag_wf2[] = "wf";
    char tag_ef1[] = "EF";
    char tag_ef2[] = "ef";
    char tag_sf1[] = "SF";
    char tag_sf2[] = "sf";
    char tag_nf1[] = "NF";
    char tag_nf2[] = "nf";
    //tags for bridge/culvert
    char tag_wb1[] = "WB";
    char tag_wb2[] = "wb";
    char tag_eb1[] = "EB";
    char tag_eb2[] = "eb";
    char tag_sb1[] = "SB";
    char tag_sb2[] = "sb";
    char tag_nb1[] = "NB";
    char tag_nb2[] = "nb";

    fp = fopen(Fnameptr->weirfilename, "r");
    if (fp == NULL)
    {
        if (*verbose == ON) printf("Weirs off\n");
        return;
    }
    if (*verbose == ON) printf("Loading weir information:\t%s\n", Fnameptr->weirfilename);

    Statesptr->weirs = ON;

    j = 0;
    do { buff[j] = fgetc(fp); } while (buff[j++] != '\n');
    buff[j] = '\0';
    sscanf(buff, "%i", &nw);

    Arrptr->Weir_hc = new double[nw];
    Arrptr->Weir_Cd = new double[nw];
    Arrptr->Weir_m = new double[nw];
    Arrptr->Weir_w = new double[nw];
    Arrptr->Weir_Typ = new int[nw]; // type of structure... weir = 0, bridge = 1;

    Arrptr->Weir_Fixdir = new int[nw];   // Fixed flow directions
    Arrptr->Weir_Identx = new int[(Parptr->xsz + 1) * (Parptr->ysz + 1)];
    Arrptr->Weir_Identy = new int[(Parptr->xsz + 1) * (Parptr->ysz + 1)];

    // Defalut to -1 for no weir link, and 0 for fixed flow direction
    for (i = 0; i <= Parptr->xsz; i++) for (j = 0; j <= Parptr->ysz; j++)
    {
        Arrptr->Weir_Identx[i + j * (Parptr->xsz + 1)] = -1;
        Arrptr->Weir_Identy[i + j * (Parptr->xsz + 1)] = -1;
    }

    for (i = 0; i < nw; i++)
    {
        j = 0;	   // load buffer until EOL
        do { buff[j] = fgetc(fp); } while (buff[j++] != '\n');
        buff[j] = '\0';  // Finish off string

        if (sscanf(buff, "%lf %lf %s %lf %lf %lf %lf",
            &x, &y, char_tmp, Arrptr->Weir_Cd + i, Arrptr->Weir_hc + i, Arrptr->Weir_m + i, Arrptr->Weir_w + i) != 7)
            Arrptr->Weir_w[i] = Parptr->dx;

        xi = (int)((x - Parptr->blx) / Parptr->dx);
        yi = (int)((Parptr->tly - y) / Parptr->dy);

        // Unfixed flow direction = 0.
        if (strcmp(char_tmp, tag_w1) == 0 || strcmp(char_tmp, tag_w2) == 0)
        {
            Arrptr->Weir_Identx[xi + 1 + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_e1) == 0 || strcmp(char_tmp, tag_e2) == 0)
        {
            Arrptr->Weir_Identx[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_s1) == 0 || strcmp(char_tmp, tag_s2) == 0)
        {
            Arrptr->Weir_Identy[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_n1) == 0 || strcmp(char_tmp, tag_n2) == 0)
        {
            Arrptr->Weir_Identy[xi + (yi + 1) * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 0;
        }
        // Control tags for one-directional flow (culverts)
        // Fixed flow directions: N = 1, E = 2, S = 3, W = 4.
        if (strcmp(char_tmp, tag_wf1) == 0 || strcmp(char_tmp, tag_wf2) == 0)
        {
            Arrptr->Weir_Identx[xi + 1 + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 4;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_ef1) == 0 || strcmp(char_tmp, tag_ef2) == 0)
        {
            Arrptr->Weir_Identx[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 2;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_sf1) == 0 || strcmp(char_tmp, tag_sf2) == 0)
        {
            Arrptr->Weir_Identy[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 3;
            Arrptr->Weir_Typ[i] = 0;
        }
        if (strcmp(char_tmp, tag_nf1) == 0 || strcmp(char_tmp, tag_nf2) == 0)
        {
            Arrptr->Weir_Identy[xi + (yi + 1) * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 1;
            Arrptr->Weir_Typ[i] = 0;
        }
        // control tags for bridge
        if (strcmp(char_tmp, tag_wb1) == 0 || strcmp(char_tmp, tag_wb2) == 0)
        {
            Arrptr->Weir_Identx[xi + 1 + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 1;
        }
        if (strcmp(char_tmp, tag_eb1) == 0 || strcmp(char_tmp, tag_eb2) == 0)
        {
            Arrptr->Weir_Identx[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 1;
        }
        if (strcmp(char_tmp, tag_sb1) == 0 || strcmp(char_tmp, tag_sb2) == 0)
        {
            Arrptr->Weir_Identy[xi + yi * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 1;
        }
        if (strcmp(char_tmp, tag_nb1) == 0 || strcmp(char_tmp, tag_nb2) == 0)
        {
            Arrptr->Weir_Identy[xi + (yi + 1) * (Parptr->xsz + 1)] = i;
            Arrptr->Weir_Fixdir[i] = 0;
            Arrptr->Weir_Typ[i] = 1;
        }

        // now a check to make sure that Arrptr->Weir_hc is greater than the ground elevation,
        // this is especally important for the SGC model where bed elevations may change.

        // first index the cell of the weir and get z0
        p0 = xi + yi * Parptr->xsz;
        if (Statesptr->SGC == ON && Arrptr->SGCwidth[p0] > 0.0) z0 = Arrptr->SGCz[p0];
        else z0 = Arrptr->DEM[p0];
        // Then index the cell it flows from and get the elevation
        if (strcmp(char_tmp, tag_w1) == 0 || strcmp(char_tmp, tag_w2) == 0 || strcmp(char_tmp, tag_wf1) == 0 || strcmp(char_tmp, tag_wf2) == 0 || strcmp(char_tmp, tag_wb1) == 0 || strcmp(char_tmp, tag_wb2) == 0)
        {
            p1 = xi + 1 + yi * Parptr->xsz;
        }
        if (strcmp(char_tmp, tag_e1) == 0 || strcmp(char_tmp, tag_e2) == 0 || strcmp(char_tmp, tag_ef1) == 0 || strcmp(char_tmp, tag_ef2) == 0 || strcmp(char_tmp, tag_eb1) == 0 || strcmp(char_tmp, tag_eb2) == 0)
        {
            p1 = xi - 1 + yi * Parptr->xsz;
        }
        if (strcmp(char_tmp, tag_s1) == 0 || strcmp(char_tmp, tag_s2) == 0 || strcmp(char_tmp, tag_sf1) == 0 || strcmp(char_tmp, tag_sf2) == 0 || strcmp(char_tmp, tag_sb1) == 0 || strcmp(char_tmp, tag_sb2) == 0)
        {
            p1 = xi + (yi - 1) * Parptr->xsz;
        }
        if (strcmp(char_tmp, tag_n1) == 0 || strcmp(char_tmp, tag_n2) == 0 || strcmp(char_tmp, tag_nf1) == 0 || strcmp(char_tmp, tag_nf2) == 0 || strcmp(char_tmp, tag_nb1) == 0 || strcmp(char_tmp, tag_nb2) == 0)
        {
            p1 = xi + (yi + 1) * Parptr->xsz;
        }
        if (Statesptr->SGC == ON && Arrptr->SGCwidth[p1] > 0.0) z1 = Arrptr->SGCz[p1];
        else z1 = Arrptr->DEM[p1];

        // now work out if either of the elevations (z0,z1) are above the weir crest hight.
        if (Arrptr->Weir_hc[i] < z0 || Arrptr->Weir_hc[i] < z1)
        {
            if (*verbose == ON)
            {
                if (Arrptr->Weir_Typ[i] == 0)
                {
                    printf("WARNING: Weir crest height is below DEM\n");
                    // for sub-grid model increase the crest height
                    //if(Statesptr->SGC==ON)
                    //{
                    Arrptr->Weir_hc[i] = Tool::getmax(z0, z1);
                    printf("Weir number %i crest height increased to %.3f m\n", i, Arrptr->Weir_hc[i]);
                    //}
                }
                //else
                //{
                //printf("WARNING: Bridge soffit height is below DEM converted to weir!!\n");
                //Arrptr->Weir_Typ[i] = 0;
                //}
            }
        }
        // need check for bridge less than or equal to subgrid width ?.......
    }

    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");
    return;
}

void Dam::DamQ(States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr) {
    int i, j;
    double h0, h1, ThreadTS, TmpTstep;
    double* hptr0, * hptr1, * qptr, * TSptr;
    int* wiptr;
    TmpTstep = Solverptr->Tstep;

    // Calculate Qx
#pragma omp parallel for private( i, h0, h1, hptr0,qptr,wiptr,TSptr,ThreadTS) 
    for (j = 0; j < Parptr->ysz; j++)
    {
        hptr0 = Arrptr->H + j * Parptr->xsz;
        qptr = Arrptr->Qx + j * (Parptr->xsz + 1) + 1;
        wiptr = Arrptr->Weir_Identx + j * (Parptr->xsz + 1) + 1;
        // initialise thread time step for openMP
        ThreadTS = TmpTstep;
        TSptr = &ThreadTS;
        for (i = 0; i < Parptr->xsz - 1; i++)
        {
            h0 = *hptr0;
            h1 = *(hptr0 + 1);
            if (h0 > Solverptr->DepthThresh || h1 > Solverptr->DepthThresh)
            {
                if (*wiptr != -1) {
                    *qptr = 0.0;
                    *qptr = CalcWeirQx(i, j, Parptr, Arrptr, Solverptr, Statesptr); // timestep update needed here TSptr}
                }
            }
            qptr++;
            hptr0++;
            wiptr++;
        }
    }
    // Calculate Qy
//#pragma omp section
#pragma omp parallel for private( i, h0, h1, hptr0,hptr1,qptr,wiptr,TSptr,ThreadTS)
    for (j = 0; j < Parptr->ysz - 1; j++)
    {
        hptr0 = Arrptr->H + j * Parptr->xsz;
        hptr1 = Arrptr->H + (j + 1) * Parptr->xsz;
        qptr = Arrptr->Qy + (j + 1) * (Parptr->xsz + 1);
        wiptr = Arrptr->Weir_Identy + (j + 1) * (Parptr->xsz + 1);
        // initialise thread time step for openMP
        ThreadTS = TmpTstep;
        TSptr = &ThreadTS;
        for (i = 0; i < Parptr->xsz; i++)
        {
            h0 = *hptr0;
            h1 = *hptr1;

            if (h0 > Solverptr->DepthThresh || h1 > Solverptr->DepthThresh)
            {
                if (Statesptr->weirs == ON && *wiptr != -1) {
                    *qptr = 0.0;
                    *qptr = CalcWeirQy(i, j, Parptr, Arrptr, Solverptr, Statesptr); // timestep update needed here TSptr
                }

            }
            hptr0++;
            hptr1++;
            qptr++;
            wiptr++;
        }
    }
};

double Dam::CalcWeirQx(int i, int j, Pars* Parptr, Arrays* Arrptr, Solver* Solverptr, States* Statesptr)
{
    double z0, z1, h0, h1, Q, hu, hd;
    int p0, p1, pq0, weir_id;
    double g;  // gravity accel

    g = Solverptr->g;

    p0 = i + j * Parptr->xsz;
    p1 = i + 1 + j * Parptr->xsz;
    pq0 = i + j * (Parptr->xsz + 1) + 1;

    z0 = Arrptr->DEM[p0];
    z1 = Arrptr->DEM[p1];
    h0 = Arrptr->H[p0];
    h1 = Arrptr->H[p1];
    weir_id = Arrptr->Weir_Identx[i + 1 + j * (Parptr->xsz + 1)];

    //Weir equation implementation altered by Rich Dawson 12 Jan 2004/3 Feb 2004.
    //One-directional flow functionality added by Matt Wilson 13 Feb 2004.
    //Bridge equation added by Jeff Neal 1 Sep 2012.

    Q = 0.0;

    if (Arrptr->Weir_Typ[weir_id] == 0) // simulate a weir
    {
        if ((z0 + h0) > (z1 + h1))		// Flow in +x direction
        {
            if ((h0 + z0) > Arrptr->Weir_hc[weir_id] && h0 > 0) // check depth is above weir and that the cell is wet
            {
                if (Arrptr->Weir_Fixdir[weir_id] == 0 || Arrptr->Weir_Fixdir[weir_id] == 2) // check for one-directional flow (culvert)
                {
                    hu = h0 + z0 - Arrptr->Weir_hc[weir_id]; // upstream head
                    hd = h1 + z1 - Arrptr->Weir_hc[weir_id]; // downstream head
                    if ((hd / hu) < Arrptr->Weir_m[weir_id]) Q = Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * pow(hu, (1.5)); // Free flow
                    else								Q = Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * hu * (sqrt(hu - hd)) / sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
                }
            }
        }

        if ((z0 + h0) < (z1 + h1))		// Flow in -x direction
        {
            if ((h1 + z1) > Arrptr->Weir_hc[weir_id] && h1 > 0) // check depth is above weir and that the cell is wet
            {
                if (Arrptr->Weir_Fixdir[weir_id] == 0 || Arrptr->Weir_Fixdir[weir_id] == 4) // check for one-directional flow (culvert)
                {
                    hu = h1 + z1 - Arrptr->Weir_hc[weir_id]; // upstream head
                    hd = h0 + z0 - Arrptr->Weir_hc[weir_id]; // downstream head
                    if ((hd / hu) < Arrptr->Weir_m[weir_id]) Q = -Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * pow(hu, (1.5)); // Free flow
                    else								Q = -Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * hu * (sqrt(hu - hd)) / sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
                }
            }
        }
    }

    return(Q);
}

double Dam::CalcWeirQy(int i, int j, Pars* Parptr, Arrays* Arrptr, Solver* Solverptr, States* Statesptr)
{
    double z0, z1, h0, h1, Q, hu, hd;
    int p0, p1, pq0, weir_id;

    double g;  // gravity accel


    g = Solverptr->g;

    pq0 = i + (j + 1) * (Parptr->xsz + 1);

    p0 = i + j * Parptr->xsz;
    p1 = i + (j + 1) * Parptr->xsz;
    z0 = Arrptr->DEM[p0];
    z1 = Arrptr->DEM[p1];
    h0 = Arrptr->H[p0];
    h1 = Arrptr->H[p1];
    weir_id = Arrptr->Weir_Identy[i + (j + 1) * (Parptr->xsz + 1)];
    //Weir equation implementation altered by Rich Dawson 12 Jan 2004/3 Feb 2004.
    //One-directional flow functionality added by Matt Wilson 13 Feb 2004.

    Q = 0.0;
    if (Arrptr->Weir_Typ[weir_id] == 0)
    {
        if ((z0 + h0) > (z1 + h1))		// Flow in +y direction
        {
            if ((h0 + z0) > Arrptr->Weir_hc[weir_id] && h0 > 0) // check depth is above weir and that the cell is wet
            {
                if (Arrptr->Weir_Fixdir[weir_id] == 0 || Arrptr->Weir_Fixdir[weir_id] == 3)  // check for one-directional flow (culvert)
                {
                    hu = h0 + z0 - Arrptr->Weir_hc[weir_id]; // upstream head
                    hd = h1 + z1 - Arrptr->Weir_hc[weir_id]; // downstream head
                    if ((hd / hu) < Arrptr->Weir_m[weir_id]) Q = Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * pow(hu, (1.5)); // Free flow
                    else								Q = Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * (hu) * (sqrt(hu - hd)) / sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
                }
            }
        }
        if ((z0 + h0) < (z1 + h1))			// Flow in -y direction
        {
            if ((h1 + z1) > Arrptr->Weir_hc[weir_id] && h1 > 0) // check depth is above weir and that the cell is wet
            {
                if (Arrptr->Weir_Fixdir[weir_id] == 0 || Arrptr->Weir_Fixdir[weir_id] == 1)  // check for one-directional flow (culvert)
                {
                    hu = h1 + z1 - Arrptr->Weir_hc[weir_id]; // upstram head
                    hd = h0 + z0 - Arrptr->Weir_hc[weir_id]; // downstream head

                    if ((hd / hu) < Arrptr->Weir_m[weir_id]) Q = -Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * pow(hu, (1.5)); // Free flow
                    else								Q = -Arrptr->Weir_Cd[weir_id] * Arrptr->Weir_w[weir_id] * (hu) * (sqrt(hu - hd)) / sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
                }
            }
        }
    }

    return(Q);
}
