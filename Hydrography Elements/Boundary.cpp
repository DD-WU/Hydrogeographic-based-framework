#include "Boundary.h"
Boundary::Boundary(const char* config, const char* sheet)
{
    parameter = new Config(config, sheet);

    if (parameter)
    {
        cout << "Creating Boundary module..." << endl;

        cout << "Using Qlim formulation for floodplain flow" << endl;
    }
}
Boundary::Boundary(Config* lisf)
{
    parameter = lisf;
    if (parameter)
    {
        cout << "   Creating Boundary module..." << endl;
        cout << "   " << endl;
    }
}
void Boundary::init() {
    LoadBCs(
        parameter->Fnameptr, 
        parameter->Statesptr, 
        parameter->Parptr, 
        parameter->BCptr, 
        parameter->Arrptr, 
        parameter->verbose);
    LoadBCVar(
        parameter->Fnameptr, 
        parameter->Statesptr, 
        parameter->Parptr, 
        parameter->BCptr, 
        parameter->CSTypePtr, 
        parameter->Arrptr,
        parameter->ChannelSegmentsVecPtr, 
        parameter->verbose);
};
//-----------------------------------------------------------------------------------
// LOADS FILE GIVING IDENTIFIERS FOR EACH BOUNDARY CELL FROM .bci FILE
// (e.g. HFIX, QVAR etc)
void Boundary::LoadBCs(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, BoundCs* BCptr, Arrays* Arrptr, int* verbose)
{
    int numBCs, i, j, BCi1, BCi2, tmpi;
    double start, finish;
    FILE* fp;
    char buff[800], buff2[800], buff3[800], side;
    double BC_tmp;
    int pi = -1, maxpi = MAXPI;  // increase max from 10 to 20K (MT)
    double px, py;

    int* new_xpi, * new_ypi, * new_PS_Ident;
    double* new_PS_Val, * new_PS_qold, * new_PS_qSGold;
    char* new_PS_Name;

    // POINT SOURCE STUFF
    BCptr->xpi = new int[maxpi];
    BCptr->ypi = new int[maxpi];
    for (i = 0; i < maxpi; i++) BCptr->xpi[i] = BCptr->ypi[i] = -1;
    BCptr->PS_Ident = new int[maxpi];
    BCptr->PS_Val = new double[maxpi];
    BCptr->PS_Name = new char[maxpi * 80];
    for (i = 0; i < maxpi; i++) { BCptr->PS_Ident[i] = 0; BCptr->PS_Val[i] = -1.0; BCptr->PS_Name[i] = '\0'; }
    BCptr->numPS = -1;

    // BOUNDARY CONDITION STUFF
    numBCs = 2 * Parptr->xsz + 2 * Parptr->ysz;
    BCptr->BC_Ident = new int[numBCs];
    BCptr->BC_Val = new double[numBCs];
    BCptr->BC_Name = new char[numBCs * 80];
    for (i = 0; i < numBCs; i++) { BCptr->BC_Ident[i] = 0; BCptr->BC_Val[i] = -1.0; BCptr->BC_Name[i] = '\0'; }

    fp = fopen(Fnameptr->bcifilename, "rb");
    if (fp == NULL) return;
    if (*verbose == ON) printf("Loading boundary condition IDs:\t%s\n", Fnameptr->bcifilename);

    while (!feof(fp))
    {
        BCi1 = BCi2 = -1; side = '\0';
        // Read NSEW and location, and determine start/finish of BC-->(BCi1,BCi2)
        fscanf(fp, "%s", buff);
        if (feof(fp)) break;
        if (buff[0] == 'N')
        {
            fscanf(fp, "%lf%lf", &start, &finish);

            if (start < Parptr->blx) start = Parptr->blx;
            if (start > Parptr->blx + Parptr->xsz * Parptr->dx) start = Parptr->blx + Parptr->xsz * Parptr->dx;
            if (finish < Parptr->blx) finish = Parptr->blx;
            if (finish > Parptr->blx + Parptr->xsz * Parptr->dx) finish = Parptr->blx + Parptr->xsz * Parptr->dx;

            BCi1 = (int)((start - Parptr->blx) / Parptr->dy);
            BCi2 = (int)((finish - Parptr->blx) / Parptr->dy);
            if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
            BCi2--; side = 'N';
        }

        if (buff[0] == 'W')
        {
            fscanf(fp, "%lf%lf", &start, &finish);

            if (start < Parptr->bly) start = Parptr->bly;
            if (start > Parptr->tly) start = Parptr->tly;
            if (finish < Parptr->bly) finish = Parptr->bly;
            if (finish > Parptr->tly) finish = Parptr->tly;

            BCi1 = (int)(2 * Parptr->xsz + 2 * Parptr->ysz - (Parptr->tly - start) / Parptr->dy);
            BCi2 = (int)(2 * Parptr->xsz + 2 * Parptr->ysz - (Parptr->tly - finish) / Parptr->dy);
            if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
            BCi2--; side = 'W';
        }

        if (buff[0] == 'S')
        {
            fscanf(fp, "%lf%lf", &start, &finish);

            if (start < Parptr->blx) start = Parptr->blx;
            if (start > Parptr->blx + Parptr->xsz * Parptr->dx) start = Parptr->blx + Parptr->xsz * Parptr->dx;
            if (finish < Parptr->blx) finish = Parptr->blx;
            if (finish > Parptr->blx + Parptr->xsz * Parptr->dx) finish = Parptr->blx + Parptr->xsz * Parptr->dx;

            BCi1 = (int)(2 * Parptr->xsz + Parptr->ysz - (start - Parptr->blx) / Parptr->dy);
            BCi2 = (int)(2 * Parptr->xsz + Parptr->ysz - (finish - Parptr->blx) / Parptr->dy);

            if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
            BCi2--; side = 'S';
        }

        if (buff[0] == 'E')
        {
            fscanf(fp, "%lf%lf", &start, &finish);

            if (start < Parptr->bly) start = Parptr->bly;
            if (start > Parptr->tly) start = Parptr->tly;
            if (finish < Parptr->bly) finish = Parptr->bly;
            if (finish > Parptr->tly) finish = Parptr->tly;

            BCi1 = (int)(Parptr->xsz + (Parptr->tly - start) / Parptr->dy);
            BCi2 = (int)(Parptr->xsz + (Parptr->tly - finish) / Parptr->dy);
            if (BCi1 > BCi2) { tmpi = BCi1; BCi1 = BCi2; BCi2 = tmpi; }
            BCi2--; side = 'E';
        }

        // Read locations of point sources or point free
        if (buff[0] == 'P' || buff[0] == 'F')
        {
            pi++;
            fscanf(fp, "%lf%lf", &px, &py);
            BCptr->xpi[pi] = (int)((px - Parptr->blx) / Parptr->dx);
            BCptr->ypi[pi] = (int)((Parptr->tly - py) / Parptr->dy);
        }

        // Read free boundary condition locations
        // load buffer until EOL
        j = 0;
        do { buff2[j] = fgetc(fp); } while (buff2[j++] != '\n' && !feof(fp));
        buff2[j - 1] = '\0';               // Finish off string
        // get buff so you know boundary type
        BC_tmp = -1;
        sscanf(buff2, "%s%lf", buff, &BC_tmp);

        // If a FREE surface boundary condition
        if ((!strcmp(buff, "FREE") || !strcmp(buff, "free")) && BCi1 > -1)
        {
            // if BC_tmp is -1 there is no slope specifed... use local slope from elevation model (origional mehod)
            if (BC_tmp < -0.999) // ie -1 (done like this as double)
            {
                if (*verbose == ON) printf("FREE on %c side start %lf end %lf\n", side, start, finish);
            }
            else
            {
                if (*verbose == ON) printf("FREE on %c side start %lf end %lf using slope %.5f\n", side, start, finish, BC_tmp);
            }
            for (i = BCi1; i <= BCi2; i++)
            {
                BCptr->BC_Ident[i] = 1;
                // store floodplain slope in BCptr->BC_Val[i]
                if (Statesptr->adaptive_ts == ON || Statesptr->qlim == ON)
                {
                    if (BC_tmp < -0.999) BCptr->BC_Val[i] = BC_tmp; // make BC_Val equal to -1 to use local water surface slope (origional lisflood)
                    else BCptr->BC_Val[i] = sqrt(BC_tmp); // sqrt of user specified slope for diffusive version (jcn)
                }
                else
                {
                    BCptr->BC_Val[i] = BC_tmp; // user specified slope or -1(for local slope) for accelleration or any other version (jcn)
                }
            }
        }
        // Read in fixed values of H for boundary
        else if ((!strcmp(buff, "HFIX") || !strcmp(buff, "hfix")) && BCi1 > -1)
        {
            sscanf(buff2, "%s%lf", buff, &BC_tmp);
            if (*verbose == ON) printf("HFIX at %lf on %c side start %lf end %lf\n", BC_tmp, side, start, finish);

            for (i = BCi1; i <= BCi2; i++)
            {
                BCptr->BC_Val[i] = BC_tmp;
                BCptr->BC_Ident[i] = 2;
            }
        }
        // Read in fixed values if Q for boundary
        else if ((!strcmp(buff, "QFIX") || !strcmp(buff, "qfix")) && BCi1 > -1)
        {
            sscanf(buff2, "%s%lf", buff, &BC_tmp);
            if (*verbose == ON) printf("QFIX at %lf on %c side start %lf end %lf\n", BC_tmp, side, start, finish);

            for (i = BCi1; i <= BCi2; i++)
            {
                BCptr->BC_Val[i] = BC_tmp;
                BCptr->BC_Ident[i] = 4;
            }
        }

        //	Read boundary names for varying values
        else if ((!strcmp(buff, "QVAR") || !strcmp(buff, "qvar")) && BCi1 > -1)
        {
            //fscanf(fp,"%s",buff);
            sscanf(buff2, "%s%s", buff3, buff);
            if (*verbose == ON)
                printf("QVAR from bdy file %s on %c side start %lf end %lf\n",
                    buff, side, start, finish);
            for (i = BCi1; i <= BCi2; i++)
            {
                strcpy(BCptr->BC_Name + i * 80, buff);
                BCptr->BC_Ident[i] = 5;
            }
        }
        else if ((!strcmp(buff, "HVAR") || !strcmp(buff, "hvar")) && BCi1 > -1)
        {
            //fscanf(fp,"%s",buff);
            sscanf(buff2, "%s%s", buff3, buff);
            if (*verbose == ON)
                printf("HVAR from bdy file %s on %c side start %lf end %lf\n", buff, side, start, finish);
            for (i = BCi1; i <= BCi2; i++)
            {
                strcpy(BCptr->BC_Name + i * 80, buff);
                BCptr->BC_Ident[i] = 3;
            }
        }
        // Fixed/Varying values/names for point sources
        // Note these need to come after the boundary conditions in the code else both will get implemented!
        // I'm not convinced &&BCptr->xpi[pi]>-1&&pi>-1 is strict enough (JCN)
        else if ((!strcmp(buff, "HFIX") || !strcmp(buff, "hfix")) && BCptr->xpi[pi] > -1 && pi > -1)
        {
            //fscanf(fp,"%lf",&BC_tmp);
            sscanf(buff2, "%s%s", buff3, buff);
            if (*verbose == ON) printf("HFIX at point [%lf,%lf] %lf\n", px, py, BC_tmp);
            BCptr->PS_Val[pi] = BC_tmp;
            BCptr->PS_Ident[pi] = 2;
        }
        else if ((!strcmp(buff, "QFIX") || !strcmp(buff, "qfix")) && BCptr->xpi[pi] > -1 && pi > -1)
        {
            //fscanf(fp,"%lf",&BC_tmp);
            if (*verbose == ON) printf("QFIX at point [%lf,%lf] %lf\n", px, py, BC_tmp);
            BCptr->PS_Val[pi] = BC_tmp;
            BCptr->PS_Ident[pi] = 4;
        }
        else if ((!strcmp(buff, "QVAR") || !strcmp(buff, "qvar")) && BCptr->xpi[pi] > -1 && pi > -1)
        {
            //fscanf(fp,"%s",buff);
            sscanf(buff2, "%s%s", buff3, buff);
            if (*verbose == ON)
                printf("QVAR at point [%lf,%lf] %s\n", px, py, buff);
            strcpy(BCptr->PS_Name + pi * 80, buff);
            BCptr->PS_Ident[pi] = 5;
        }
        else if ((!strcmp(buff, "HVAR") || !strcmp(buff, "hvar")) && BCptr->xpi[pi] > -1 && pi > -1)
        {
            //fscanf(fp,"%s",buff);
            sscanf(buff2, "%s%s", buff3, buff);
            if (*verbose == ON) printf("HVAR at point [%lf,%lf] %s\n", px, py, buff);
            strcpy(BCptr->PS_Name + pi * 80, buff);
            BCptr->PS_Ident[pi] = 3;
        }
        else if ((!strcmp(buff, "FREE") || !strcmp(buff, "free")) && BCptr->xpi[pi] > -1 && pi > -1 && Statesptr->SGC == ON)
        {
            // point FREE boundary for internal SGC boundaries
            if (*verbose == ON) printf("FREE at point [%lf,%lf] with slope %lf\n", px, py, BC_tmp);
            BCptr->PS_Val[pi] = BC_tmp;
            BCptr->PS_Ident[pi] = 6;
        }
        else
        {
            if (*verbose == ON) printf("WARNING: Incorrect boundary condition in .bci file\n");
        }

    }

    if (pi > -1)
    {
        pi++;
        new_xpi = new int[pi];
        new_ypi = new int[pi];
        new_PS_Ident = new int[pi];
        new_PS_Val = new double[pi];
        new_PS_qold = new double[pi]();
        new_PS_qSGold = new double[pi]();
        new_PS_Name = new char[pi * 80];

        for (i = 0; i < pi; i++)
        {
            new_xpi[i] = BCptr->xpi[i];
            new_ypi[i] = BCptr->ypi[i];
            new_PS_Ident[i] = BCptr->PS_Ident[i];
            new_PS_Val[i] = BCptr->PS_Val[i];
            for (j = 0; j < 80; j++) new_PS_Name[i * 80 + j] = BCptr->PS_Name[i * 80 + j];
        }
        BCptr->xpi = new_xpi;
        BCptr->ypi = new_ypi;
        BCptr->PS_Ident = new_PS_Ident;
        BCptr->PS_Val = new_PS_Val;
        BCptr->PS_qold = new_PS_qold;
        BCptr->PS_qSGold = new_PS_qSGold;
        BCptr->PS_Name = new_PS_Name;

        BCptr->numPS = pi;
    }

    if (*verbose == ON) printf("Done.\n\n");
    //  LoadBCVar();

    fclose(fp);

    return;
}
//-----------------------------------------------------------------------------
// LOAD TIME VARYING BOUNDARY CONDITIONS FROM .bdy FILE
void Boundary::LoadBCVar(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, BoundCs* BCptr, ChannelSegmentType* ChannelSegments, Arrays* Arrptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr, int* verbose)
{
    FILE* fp;
    int i, j, nbdy = 0, ndata, numBCs, chseg;
    char buff[255], units[80];
    double** peterpointer;
    bool flag = false; int temp;
    numBCs = 2 * Parptr->xsz + 2 * Parptr->ysz;
    BCptr->BCVarlist = new double* [numBCs];
    for (i = 0; i < numBCs; i++) BCptr->BCVarlist[i] = NULL;

    if (Statesptr->ChannelPresent == ON) {
        for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
        {
            ChannelSegments[chseg].QVarlist = new double* [ChannelSegments[chseg].chsz];//ChannelSegments[chseg].chsz
            for (i = 0; i < ChannelSegments[chseg].chsz; i++) ChannelSegments[chseg].QVarlist[i] = NULL;
        }
    }

    fp = fopen(Fnameptr->bdyfilename, "r");
    if (fp == NULL) return;
    if (*verbose == ON) printf("Loading time varying boundary conditions:\t%s\n", Fnameptr->bdyfilename);

    while (!feof(fp))
    {
        j = 0;  // skip 1st comment line

        if (nbdy == 0) {
            do {
                buff[j] = fgetc(fp);
                if (feof(fp)) break;
            } while (buff[j++] != '\n');

        }

        fscanf(fp, "%s", buff);
        if (feof(fp) || buff[0] == '\n') break;

        // Check through list (2d domain) of boundary names, if a match found, assign
        // pixel to this boundary and set peterpointer to point to
        // relevant structure. If none found set to NULL, boundary condition data is
        // unassigned but file is still read through to get to next set of data.
        peterpointer = NULL;
        for (i = 0; i < numBCs; i++)
        {
            if (!strcmp(buff, (BCptr->BC_Name + i * 80)))
            {
                BCptr->BC_Val[i] = nbdy;
                peterpointer = BCptr->BCVarlist;
            }
        }

        // Check through list (river channel) of boundary names, if a match found, assign
        // channel node to this boundary and set peterpointer to point to
        // relevant structure.
        if (Statesptr->ChannelPresent == ON)
        {
            for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
            {
                for (i = 0; i < ChannelSegments[chseg].chsz; i++)
                {
                    if (!strcmp(buff, (ChannelSegments[chseg].Q_Name + i * 80)))
                    {
                        ChannelSegments[chseg].Q_Val[i] = nbdy;
                        peterpointer = ChannelSegments[chseg].QVarlist;
                    }
                }
            }
        }

        // same for point sources
        for (i = 0; i < BCptr->numPS; i++)
        {
            if (!strcmp(buff, (BCptr->PS_Name + i * 80)))
            {
                BCptr->PS_Val[i] = nbdy;
                peterpointer = BCptr->BCVarlist;
            }
        }

        // Set up arrays, 1 element per time step for BCVarlist, 1 per data
        // point for Vargiven and Tgiven
        fscanf(fp, "%i%s", &ndata, units);
        if (peterpointer != NULL)
        {
            int* x, * y;
            peterpointer[nbdy] = new double[ndata * 2 + 2];
            peterpointer[nbdy][ndata * 2 + 1] = -1;

            // Go through the motions even if peterpointer==NULL
            for (i = 0; i < ndata; i++) fscanf(fp, "%lf%lf", peterpointer[nbdy] + i * 2, peterpointer[nbdy] + i * 2 + 1);
            if (!strcmp(units, "hours")) for (i = 0; i < ndata; i++) *(peterpointer[nbdy] + i * 2 + 1) *= 3600;
            else if (!strcmp(units, "days")) for (i = 0; i < ndata; i++) *(peterpointer[nbdy] + i * 2 + 1) *= (3600 * 24);

            fgetc(fp);		// Get end of line character
        }

        //    for(i=0;i<ndata*2+2;i++){
        //      printf("\nLoadBCVar: i=%d, peterpointer[0]+i*2=%lf, peterpointer[nbdy]+i*2+1=%lf",i,*(peterpointer[0]+i*2),*(peterpointer[0]+i*2+1));
        //    }

        if (peterpointer == NULL)
        {
            printf("WARNING: bdy %s is unreferenced - data ignored.\n", buff);
            continue;
        }

        nbdy++;

        if (*verbose == ON) printf("bdy %s read.\n", buff);
    }

    // Check for bdy names not found in bdy file
    if (Statesptr->ChannelPresent == ON) {
        for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) for (i = 0; i < ChannelSegments[chseg].chsz; i++) if (ChannelSegments[chseg].Q_Ident[i] == 5 && ChannelSegments[chseg].Q_Val[i] < 0) // CCS
        {
            printf("WARNING: bdy %s in river file not found in bdy file - ignored.\n", ChannelSegments[chseg].Q_Name + i * 80);
            ChannelSegments[chseg].Q_Ident[i] = 0;
        }
    }
    for (i = 0; i < numBCs; i++) if ((BCptr->BC_Ident[i] == 5 || BCptr->BC_Ident[i] == 3) && BCptr->BC_Val[i] < 0)
    {
        printf("WARNING: bdy %s in bci file not found in bdy file - ignored.\n", BCptr->BC_Name + i * 80);
        for (j = 0; j < numBCs; j++) if (!strcmp(BCptr->BC_Name + i * 80, BCptr->BC_Name + j * 80)) BCptr->BC_Ident[j] = 0;
    }

    for (i = 0; i < BCptr->numPS; i++) if ((BCptr->PS_Ident[i] == 5 || BCptr->PS_Ident[i] == 3) && BCptr->PS_Val[i] < 0)
    {
        printf("WARNING: bdy %s in bci file not found in bdy file - ignored.\n", BCptr->BC_Name + i * 80);
        for (j = 0; j < BCptr->numPS; j++) if (!strcmp(BCptr->PS_Name + i * 80, BCptr->PS_Name + j * 80)) BCptr->PS_Ident[j] = 0;
    }

    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");
    return;
}

