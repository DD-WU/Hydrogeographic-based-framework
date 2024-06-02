#include "River.h"
#include"tool.h"
#include <sstream>
Water_API River::River(Config* param)
{
	parameter = param;
	parameter->Statesptr->ChannelPresent = ON;
    if (parameter)
    {
        cout << "   Creating river module..." << endl;
        cout << "   " << endl;
        if (param->Statesptr->diffusive == ON) {
            cout << "   Use diffusive module...\n" << endl;
            cout << "   " << endl;
        }
        cout << "please notes that river's update function must be in front of surface according to lisflood-fp" << endl;
    }
    vector<QID7_Store> QID7; //#CCS A temporary store for some ChannelSegments variables that cannot be written to the correct location during LoadRiver().
    vector<QID7_Store>* QID7_Vec_Ptr; // CCS
    QID7_Vec_Ptr = &QID7; // CCS
    if (parameter->Arrptr->DEM == NULL) {
        cout << "The discretion of river depend on the surface, however, we have not detected your DEM " << endl;
        abort();
    }
    else
    {
        LoadRiverNetwork(parameter->Fnameptr, parameter->Statesptr, parameter->Parptr, parameter->ChannelSegmentsVecPtr, parameter->Arrptr, QID7_Vec_Ptr, parameter->RiversIndexVecPtr, parameter->verbose); // CCS
    }
    parameter->CSTypePtr = &(parameter->ChannelSegmentsVecPtr->front()); // CCS has to be defined after LoadRiverNetwork has completed.
    parameter->RiversIndexPtr = &(parameter->RiversIndexVecPtr->front());  // CCS has to be defined after LoadRiverNetwork has completed.
    if (QID7.size() != 0) // CCS If there are any tribs then we need to copy the terms from the temp store to the correct place.
    {
        QID7_Store* QID7Ptr = &QID7[0]; // CCS
        UpdateChannelsVector(parameter->Statesptr, parameter->CSTypePtr, QID7_Vec_Ptr, QID7Ptr, parameter->RiversIndexPtr); // CCS
    }
    // apply different starting methods for channel

}

void River::init()
{

    if (parameter->Statesptr->ChannelPresent == ON)
    {
        // calc initial steady state flows down channel
        CalcChannelStartQ(parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->CSTypePtr, parameter->RiversIndexVecPtr, parameter->RiversIndexPtr);

        if (parameter->Statesptr->startfile == ON)
        {
            // start file is specified. Do nothing, as starting H values for channel already read in from the startfile.
            SetChannelStartHfromQ(parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->CSTypePtr, parameter->Solverptr, parameter->RiversIndexVecPtr, parameter->RiversIndexPtr);
        }
        else if (parameter->Statesptr->startq == ON)
        {
            // Kinematic: Uses the kinematic initial solution to calculate H from Q
            // Diffusive: Uses diffusive steady state initial solution (default) or can use full dynamic steady state
            // initial if turned on using -dynsw on command line or "ch_dynamic" in the parameter file

            // use the flows to calculate a starting H
            SetChannelStartHfromQ(parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->CSTypePtr, parameter->Solverptr, parameter->RiversIndexVecPtr, parameter->RiversIndexPtr);
        }
        else
        {
            // set channel start H to default or user defined H
            SetChannelStartH(parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->CSTypePtr, parameter->RiversIndexVecPtr, parameter->RiversIndexPtr);
        }
        if (parameter->Statesptr->river_couple == 1)
        {
            Load_river_pipe_index(parameter->Fnameptr, parameter->Statesptr);
        }
    }
    for (int j = 0; j < parameter->Parptr->ysz; j++) {
        double* temp_hptr, * temp_hptr1;
        temp_hptr = parameter->Arrptr->H + j * parameter->Parptr->xsz;
        temp_hptr1 = parameter->Arrptr->temp_H + j * parameter->Parptr->xsz;
        for (int i = 0; i < parameter->Parptr->xsz; i++) {
            if (*temp_hptr != 0) {
                *temp_hptr1 = *temp_hptr;
            }
            temp_hptr++;
            temp_hptr1++;
        }
    }
}

void River::update()
{
    parameter->ChannelSegments = &(parameter->ChannelSegmentsVecPtr->front());
    if (parameter->Statesptr->ChannelPresent == ON)
    {

       parameter-> tstep_channel = parameter->Solverptr->t - parameter->Previous_t; // calculate river timestep (deltaT since last river calc)
        //This place is not really well written with "if" judgment, if you want to add functions, you can try to use "else" to update it, what's more, you can also try to reconstruct this part. By Wang.
        if (parameter->Statesptr->diffusive == ON)
            ChannelQ_Diff1(parameter->tstep_channel,
                parameter->Statesptr,
                parameter->Parptr,
                parameter->Solverptr,
                parameter->BCptr,
                parameter->ChannelSegments,
                parameter->Arrptr,
                parameter->RiversIndexVecPtr,
                parameter->RiversIndexPtr);
        else cout << "to be written,"<<endl;
        // end part
        parameter->tstep_counter = 0;  // set timestep counter to zero
        parameter->Previous_t = parameter->Solverptr->t; // record previous timestep

    }
}

void River::finalize()
{
    //to do
}

//-----------------------------------------------------------------------------
// Set start water depths for channels
void River::SetChannelStartH(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, vector<int>* RiversIndexVecPtr, int* RiversIndexPtr)
{

    int i, chseg, pi, pj;
    int nriv, high, low; // CCS For mulitple river loop
    ChannelSegmentType* csp;

    for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS Multiple river loop
    {
        high = RiversIndexPtr[nriv] - 1;
        if (nriv == 0)
        {
            low = 0;
        }
        else
        {
            low = RiversIndexPtr[nriv - 1];
        }

        for (chseg = low; chseg <= high; chseg++) // CCS
        {
            csp = ChannelSegments + chseg;

            // Setup starting water depth
            for (i = 0; i < csp->chsz; i++)
            {
                pi = csp->ChanX[i];
                pj = csp->ChanY[i];

                // set all water depths to start value (either default or defined in par file)
                if (chseg != low && i == csp->chsz - 1) // trib and last point - ie dummy junction node // CCS
                {
                    csp->JunctionH = Parptr->ch_start_h;
                }
                else // all points on main channel and all but last one on tribs
                {
                    Arrptr->H[pi + pj * Parptr->xsz] = Parptr->ch_start_h;
                    Arrptr->H[pi + pj * Parptr->xsz] = Parptr->ch_start_h;
                }
            }
        }
    }

}

//-----------------------------------------------------------------------------
// calculate initial channel flows assuming steady state and no fp flow
void River::CalcChannelStartQ(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, vector<int>* RiversIndexVecPtr, int* RiversIndexPtr)
{

    int i, chseg, pi, pj;
    int nriv, high, low; // CCS For multiple river loop
    double Qstart;
    ChannelSegmentType* csp;

    // reverse loop to get trib flows to pass on to d/s channel
    for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS
    {
        high = RiversIndexPtr[nriv] - 1;
        if (nriv == 0)
        {
            low = 0;
        }
        else
        {
            low = RiversIndexPtr[nriv - 1];
        }

        for (chseg = high; chseg >= low; chseg--) // CCS
        {
            csp = ChannelSegments + chseg;

            //set Q to zero at start
            Qstart = 0;

            // loop through and pickup any inflows adding them up as we go downstream and calculating the flow depth as we go
            for (i = 0; i < csp->chsz; i++)
            {
                pi = csp->ChanX[i];
                pj = csp->ChanY[i];

                if (csp->Q_Ident[i] == 4)
                {
                    // fixed flow 
                    Qstart += csp->Q_Val[i];
                }
                if (csp->Q_Ident[i] == 5)
                {
                    // interpolate from hydrograph
                    Qstart += Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[i]], 0);
                }
                if (csp->Q_Ident[i] == 7)
                {
                    // tributary connection
                    Qstart += csp->Q_Val[i];
                }

                csp->ChanQ[i] = Qstart;

            }

            if (chseg != low) // only for tribs record flow in d/s channel // CCS
            {
                ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc] = Qstart;
            }
        }
    }

    return;
}
//-----------------------------------------------------------------------------
// Set start water depths for channels based on initial channel flows
void River::SetChannelStartHfromQ(States* Statesptr, Pars* Parptr, Arrays* Arrptr, ChannelSegmentType* ChannelSegments, Solver* Solverptr, vector<int>* RiversIndexVecPtr, int* RiversIndexPtr)
{

    int i, chseg, pi, pj, iter, pid, pjd;
    double divg, Elevus, Elevds, dh, err_d, Jus = 0.0, Jds = 0.0, Eus = 0.0, Eds = 0.0, dx, eps;
    double WSbc;  // value used to hold last channel element flow h (d/s boundary condition)
    int HoutFREE = OFF; // flag to indicate calc h BC from slope on  the fly later.
    int nodes;  // temporary variable to hold number of cells/nodes in segment
    int nriv, low, high; // CCS For multiple channel loop
    ChannelSegmentType* csp;
    divg = Solverptr->divg;

    // Diffusive channel solver start H from Q
    if (Statesptr->diffusive == ON)
    {
        for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS
        {
            high = RiversIndexPtr[nriv] - 1;
            if (nriv == 0)
            {
                low = 0;
            }
            else
            {
                low = RiversIndexPtr[nriv - 1];
            }
            for (chseg = low; chseg <= high; chseg++) // CCS
            {
                csp = ChannelSegments + chseg;
                // loop through and calc H from previously worked out Q (reverse loop)
                for (i = csp->chsz - 1; i >= 0; i--)
                {
                    pi = csp->ChanX[i];
                    pj = csp->ChanY[i];
                    // setup temporary variable - just because it is easier to read code
                    nodes = csp->chsz;

                    if (chseg != low && i == csp->chsz - 1) // trib and last point - ie dummy junction node // CCS
                    {
                        //csp->JunctionH=CalcA(csp->ChanN[i],csp->Shalf[i],csp->ChanWidth[i],csp->ChanQ[i])/csp->ChanWidth[i]; // calc new area from this flow and calc h
                        csp->JunctionH = Arrptr->H[pi + pj * Parptr->xsz];
                    }
                    else if (chseg == low && i == csp->chsz - 1) // last point in main channel // CCS
                    {
                        if (csp->Q_Ident[nodes - 1] == 1) //free boundary, calc h from slope ## This BC needs some stability work
                        {
                            // set flag so we can calc h from the flow and slope "on the fly" in CalcF() function
                            WSbc = 0; // set dummy value as will be calculated later
                            HoutFREE = ON;
                        }
                        else if (csp->Q_Ident[nodes - 1] == 2) // fixed H out 
                        {
                            // note - we will need to subtract DEM Elev to get h from water elevation entered
                            WSbc = csp->Q_Val[nodes - 1];
                        }
                        else if (csp->Q_Ident[nodes - 1] == 3) //interpolate H from stage hydrograph
                        {
                            // note - we will need to subtract DEM Elev to get h from water elevation entered
                            WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], Solverptr->t);
                        }
                        else if (csp->Q_Ident[nodes - 1] == 8) //interpolate H from rating curve ## This BC needs some stability work
                        {
                            // set flag to 2 so we can calc h from the stage discharge curve "on the fly" in CalcF() function
                            HoutFREE = 2;
                            WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], csp->Q_Val[nodes - 1]);
                        }
                        else
                        {
                            printf("WARNING: Main Channel has no d/s BC\n");
                            // already checked in input function, and set to FREE by default, but leave this here just in case of any odd bugs
                        }
                        if (HoutFREE == ON) // only if HoutFREE flagged ON
                        {
                            // calc h "on the fly" from slope for FREE BC
                            if (csp->Q_Val[csp->chsz - 1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
                            {
                                Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = Tool::CalcA(csp->ChanN[nodes - 1], csp->Shalf[nodes - 1], csp->ChanWidth[nodes - 1], csp->ChanQ[nodes - 1])
                                    / csp->ChanWidth[nodes - 1];
                            }
                            else // else use user supplied slope
                            {
                                Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = Tool::CalcA(csp->ChanN[nodes - 1], sqrt(csp->Q_Val[nodes - 1]), csp->ChanWidth[nodes - 1], csp->ChanQ[nodes - 1])
                                    / csp->ChanWidth[nodes - 1];
                            }
                        }
                        else if (HoutFREE == 2) // HoutFREE flagged for stage discharge curve
                        {
                            // subtract dem to get water depth
                            Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = WSbc - Arrptr->DEM[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz];
                            // This BC needs some stability work
                        }
                        else // HoutFREE not flagged
                        {
                            // subtract dem to get water depth
                            if (chseg == low) Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = WSbc - Arrptr->DEM[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz]; // CCS
                            else         Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = WSbc - csp->JunctionDEM; //trib uses dummy junction node
                        }
                    }

                    /*
                    Full dynamic steady state solution to incorporate dys/dx and u/g du/dx terms to improve stability of
                    initial water surface profile
                    */
                    else // all points on main channel and all but last one on trib
                    {
                        iter = 0;
                        pid = csp->ChanX[i + 1];
                        pjd = csp->ChanY[i + 1];
                        Arrptr->H[pi + pj * Parptr->xsz] = Arrptr->H[pid + pjd * Parptr->xsz];
                        Elevds = Arrptr->DEM[pid + pjd * Parptr->xsz] + Arrptr->H[pid + pjd * Parptr->xsz];
                        Jds = CalcEnergySlope(csp->ChanN[i + 1], csp->ChanWidth[i + 1], Arrptr->H[pid + pjd * Parptr->xsz], csp->ChanQ[i + 1]);
                        Eds = Elevds + pow(csp->ChanQ[i + 1] / (csp->ChanWidth[i + 1] * Arrptr->H[pid + pjd * Parptr->xsz]), 2) * divg * Solverptr->dynsw;
                        dx = (csp->Chainage[i + 1] - csp->Chainage[i]);
                        dh = 0.00001;
                        eps = 1.5 * dh;
                        err_d = -2.0 * dh;
                        do {
                            iter = iter + 1;
                            Arrptr->H[pi + pj * Parptr->xsz] = Arrptr->H[pi + pj * Parptr->xsz] + dh * signR(err_d);
                            Elevus = Arrptr->DEM[pi + pj * Parptr->xsz] + Arrptr->H[pi + pj * Parptr->xsz];
                            Jus = CalcEnergySlope(csp->ChanN[i], csp->ChanWidth[i], Arrptr->H[pi + pj * Parptr->xsz], csp->ChanQ[i]);
                            Eus = Elevus + pow(csp->ChanQ[i] / (csp->ChanWidth[i] * Arrptr->H[pi + pj * Parptr->xsz]), 2) * divg * Solverptr->dynsw;
                            err_d = 0.5 * (Jus + Jds) * dx + Eds - Eus;
                        } while (fabs(err_d) > eps);
                    }
                }
            }
        }
    }

    // Kinematic channel solver start H from Q
    else
    {
        for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS Multiple channel loop
        {
            high = RiversIndexPtr[nriv] - 1;
            if (nriv == 0)
            {
                low = 0;
            }
            else
            {
                low = RiversIndexPtr[nriv - 1];
            }

            // reverse loop 
            for (chseg = high; chseg >= low; chseg--) // CCS
            {
                csp = ChannelSegments + chseg;

                // loop through and calc H from previously worked out Q
                for (i = 0; i < csp->chsz; i++)
                {
                    pi = csp->ChanX[i];
                    pj = csp->ChanY[i];

                    if (chseg != low && i == csp->chsz - 1) // trib and last point - ie dummy junction node // CCS
                    {
                        csp->JunctionH = Tool::CalcA(csp->ChanN[i], csp->Shalf[i], csp->ChanWidth[i], csp->ChanQ[i]) / csp->ChanWidth[i]; // calc new area from this flow and calc h
                    }
                    else // all points on main channel and all but last one on tribs
                    {
                        Arrptr->H[pi + pj * Parptr->xsz] = Tool::CalcA(csp->ChanN[i], csp->Shalf[i], csp->ChanWidth[i], csp->ChanQ[i]) / csp->ChanWidth[i]; // calc new area from this flow and calc h
                    }
                }
            }
        }

    }



    //write_profile("FREE",-1,".profile",Statesptr,ChannelSegments,Arrptr,Parptr);

    return;
}

void River::LoadRiverNetwork(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr, Arrays* Arrptr, vector<QID7_Store>* QID7_Vec_Ptr, vector<int>* RiversIndexVecPtr, int* verbose)
{
    FILE* rfp;  // local file pointer
    int i, n;
    int tmp_size;

    if (Statesptr->multiplerivers == 0)
    {
        LoadRiver(Fnameptr, Statesptr, Parptr, ChannelSegmentsVecPtr, Arrptr, QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver once.
        tmp_size = ChannelSegmentsVecPtr->size();
        RiversIndexVecPtr->push_back(tmp_size); // CCS
    }
    else if (Statesptr->multiplerivers == 1)
    {
        rfp = fopen(Fnameptr->multiriverfilename, "r");
        fscanf(rfp, "%i", &n);
        if (*verbose == ON) printf("Loading %i Rivers\n\n", n);

        for (i = 0; i < n; i++)
        {
            fscanf(rfp, "%s", Fnameptr->rivername); //Scan the next .river filename from the .rivers file and asign it to Fnameptr->rivername.
            LoadRiver(Fnameptr, Statesptr, Parptr, ChannelSegmentsVecPtr, Arrptr, QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver in loop.
            tmp_size = ChannelSegmentsVecPtr->size();
            RiversIndexVecPtr->push_back(tmp_size); /* Builds the RiversIndex vector so we know where one river stops and the next starts within
                                                       the ChannelSegments vector. */
        }

        if (*verbose == ON) printf("%i Rivers Loaded Successfully.\n\n", n);
    }

    return;
}
//-----------------------------------------------------------------------------
// LOAD RIVER DATA FROM FILE
void River::LoadRiver(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr, Arrays* Arrptr, vector<QID7_Store>* QID7_Vec_Ptr, vector<int>* RiversIndexVecPtr, int* verbose)
{
    FILE* fp;  // local file pointer
    int npoints, * xpi, * ypi;
    int* trib; // temp xs array to record any trib connections
    double* xp, * yp, * wp, * np, * hp, * cp, * rp, * ap, total_length = 0.0, * qp;
    char buff[800], buff2[800], buff3[800];
    int i, j, pi, pj, ni, nj, oldpi, oldpj, i1, i2;
    double tmp1, tmp2, tmp3;
    double grad; // temporary gradient calculation variable eventually stored in csp->Shalf
    char* Q_Name_tmp;
    char buffer[80];
    int* Q_Ident_tmp;
    int count, chseg, SegOutNo, tmp_int;

    // MSH: csp is a utility pointer,

    fp = fopen(Fnameptr->rivername, "r");

    if (fp == NULL) return;

    Statesptr->ChannelPresent = ON;

    if (*verbose == ON) printf("Loading channel information:\t%s\n", Fnameptr->rivername);

    fscanf(fp, "%s", buffer);
    if (!strcmp(buffer, "Tribs") || !strcmp(buffer, "tribs") || !strcmp(buffer, "TRIBS"))
    {
        Statesptr->TribsPresent = ON;

        // MSH: Since we haven't allocated the memory yet, we can't assign the number of channel segments in
        // the first element of the ChannelSegments array - so read into a temp variable
        fscanf(fp, "%i", &tmp_int);
        if (*verbose == ON) printf("%i Channel Segments\n", tmp_int);
    }
    else
    {
        rewind(fp);
        tmp_int = 1;
        if (*verbose == ON) printf("%i Channel Segment\n", tmp_int);
    }

    for (chseg = 0; chseg < tmp_int; chseg++) // CCS Slight reorganisation of old LoadRiver function but fundamentally unchanged.
    {
        ChannelSegmentType tmp_chan; // CCS
        tmp_chan.Next_Segment = tmp_chan.Next_Segment_Loc = -1; // CCS from old code
        tmp_chan.N_Channel_Segments = tmp_int; // CCS
        ChannelSegmentType* csp = &tmp_chan; // CCS

        fscanf(fp, "%i", &npoints);
        if (*verbose == ON) printf("%i points in channel segment %i\n", npoints, chseg);

        //setup local temporary arrays, Note the () at the end ensures all elements are initialised to zero
        xp = new double[npoints]();
        yp = new double[npoints]();
        wp = new double[npoints]();
        np = new double[npoints]();
        hp = new double[npoints]();
        cp = new double[npoints]();
        qp = new double[npoints]();
        rp = new double[npoints](); // chainage ratio between entered cross sections and cell chainage
        ap = new double[npoints](); // actual entered xs chainage
        trib = new int[npoints]();
        Q_Name_tmp = new char[npoints * 80]();
        Q_Ident_tmp = new int[npoints]();

        for (i = 0; i < npoints; i++)
        {
            // set default qp and trib value (all the other arrays are pre zeroed with new command and end brackets () )
            qp[i] = -1;
            trib[i] = -1;

            fscanf(fp, "%lf %lf", xp + i, yp + i); // Load x,y values.

            // load buffer until EOL
            j = 0;
            do { buff[j] = fgetc(fp); } while (buff[j++] != '\n');
            buff[j] = '\0';									// Finish off string
            if (sscanf(buff, "%lf%lf%lf", &tmp1, &tmp2, &tmp3) == 3)
            {        										// Only store values if 3 reads successful
                wp[i] = tmp1;
                np[i] = tmp2;
                hp[i] = tmp3;
                if (*verbose == ON)
                    printf("Xsec %4i\tw=%8.3f n=%5.3f z=%6.3f\n", i, wp[i], np[i], hp[i]);
            }

            if (sscanf(buff, "%lf%lf%lf%s%s", &tmp1, &tmp1, &tmp1, buff2, buff3) == 5  // 4+5th item found must be Q_Ident
                || sscanf(buff, "%s%s", buff2, buff3) == 2) // OR No channel info - just Q in
            {
                if (!strcmp(buff2, "QFIX") || !strcmp(buff2, "qfix"))
                {
                    Q_Ident_tmp[i] = 4;
                    sscanf(buff3, "%lf", qp + i);
                    if (*verbose == ON) printf("Xsec %4i\tQFIX at %7.2f\n", i, qp[i]);
                }
                if (!strcmp(buff2, "QVAR") || !strcmp(buff2, "qvar"))
                {
                    Q_Ident_tmp[i] = 5;
                    strcpy(Q_Name_tmp + i * 80, buff3);
                    if (*verbose == ON) printf("Xsec %4i\tQVAR from bdy file %s\n", i, Q_Name_tmp + i * 80);
                }
                if (!strcmp(buff2, "QOUT") || !strcmp(buff2, "qout"))
                {
                    Q_Ident_tmp[i] = 6;
                    sscanf(buff3, "%i", &SegOutNo);
                    if (*verbose == ON) printf("Xsec %4i\tQOUT from segment %i discharges into segment %i\n", i, chseg, SegOutNo);
                }
                if (!strcmp(buff2, "TRIB") || !strcmp(buff2, "trib"))
                {
                    Q_Ident_tmp[i] = 7;
                    sscanf(buff3, "%i", &SegOutNo);
                    if (*verbose == ON) printf("Xsec %4i\tQin from segment %i\n", i, SegOutNo);
                    trib[i] = SegOutNo;
                }
                if (!strcmp(buff2, "FREE") || !strcmp(buff2, "free"))
                    // normal depth based on slope, if -1 then use last channel segment slope else use slope supplied
                    // NOT fully working for Diffusive - stability issues
                {
                    Q_Ident_tmp[i] = 1;
                    sscanf(buff3, "%lf", qp + i);
                    if (qp[i] < -0.999) // ie -1 (done like this as double)
                    {
                        if (*verbose == ON) printf("Xsec %4i\tFREE using end slope\n", i);
                    }
                    else
                    {
                        if (*verbose == ON) printf("Xsec %4i\tFREE using slope %7.4f\n", i, qp[i]);
                    }
                }
                if (!strcmp(buff2, "HFIX") || !strcmp(buff2, "hfix"))
                {
                    Q_Ident_tmp[i] = 2;
                    sscanf(buff3, "%lf", qp + i);
                    if (*verbose == ON) printf("Xsec %4i\tHFIX at %7.2f\n", i, qp[i]);
                }
                if (!strcmp(buff2, "HVAR") || !strcmp(buff2, "hvar"))
                {
                    Q_Ident_tmp[i] = 3;
                    strcpy(Q_Name_tmp + i * 80, buff3);
                    if (*verbose == ON) printf("Xsec %4i\tHVAR from bdy file %s\n", i, Q_Name_tmp + i * 80);
                }
                if (!strcmp(buff2, "RATE") || !strcmp(buff2, "rate"))
                    // NOT fully working for Diffusive - stability issues
                {
                    Q_Ident_tmp[i] = 8;
                    strcpy(Q_Name_tmp + i * 80, buff3);
                    if (*verbose == ON) printf("Xsec %4i\tRATE from bdy file %s\n", i, Q_Name_tmp + i * 80);
                }
            }
        }

        if (*verbose == ON) printf("Channel data read for segment %i - interpolating values.\n", chseg);



        // Estimate number of channel pixels - to ensure we allocate enough memory for temporary arrays xpi,ypi
        // total length divided by cell size and then double this value.
        ap[0] = 0;
        for (i = 0; i < npoints - 1; i++)
        {
            // calc straight line chainage between entered cross sections - used for cell independent chainage calcs. Add up chainage
            ap[i + 1] = ap[i] + sqrt((xp[i + 1] - xp[i]) * (xp[i + 1] - xp[i]) + (yp[i + 1] - yp[i]) * (yp[i + 1] - yp[i]));
        }
        xpi = new int[int(2.0 * ap[npoints - 1] / Parptr->dx)]();
        ypi = new int[int(2.0 * ap[npoints - 1] / Parptr->dx)]();



        // Insert channel into DEM grid
        oldpi = (int)((xp[0] - Parptr->tlx) / Parptr->dx);
        oldpj = (int)((Parptr->tly - yp[0]) / Parptr->dy);
        total_length = 0.0;

        count = 0;
        for (i = 1; i < npoints; i++)
        {
            pi = (int)((xp[i] - Parptr->tlx) / Parptr->dx);
            pj = (int)((Parptr->tly - yp[i]) / Parptr->dy);
            for (j = 0; j <= 1000; j++)					// Take very small steps and insert channel
            {															// whenever x,y position changes
                ni = oldpi + ((pi - oldpi) * j / 1000);
                nj = oldpj + ((pj - oldpj) * j / 1000);
                if (ni >= 0 && ni < Parptr->xsz && nj >= 0 && nj < Parptr->ysz) // check it stays within DEM
                {
                    if (count == 0)								// Always insert first point
                    {
                        xpi[count] = ni;
                        ypi[count] = nj;
                        Arrptr->ChanMask[ni + nj * Parptr->xsz] = 1; // mark mask with value of 1 - will renumber later in order
                        count++;
                    }
                    else if (ni != xpi[count - 1] || nj != ypi[count - 1]) // if grid location changes
                    {
                        if (Arrptr->ChanMask[ni + nj * Parptr->xsz] == -1)   // channel mask not set
                        {
                            xpi[count] = ni;
                            ypi[count] = nj;
                            Arrptr->ChanMask[ni + nj * Parptr->xsz] = 1; // mark mask with value of 1 - will renumber later in order
                            total_length += Parptr->dx * sqrt(pow(ni - xpi[count - 1], (2.0)) + pow(nj - ypi[count - 1], (2.0)));
                            count++;
                        }
                        else
                        {
                            // channel mask set, so likely that it is trib junction point.
                            // NOTE, cannot have crossing channels !!!
                            // DO NOT mark mask with value of 1 - as this is the junction
                            // of the trib with main channel so is already marked for main channel
                            xpi[count] = ni;
                            ypi[count] = nj;
                            total_length += Parptr->dx * sqrt(pow((double)(ni - xpi[count - 1]), 2) + pow((double)(nj - ypi[count - 1]), 2));
                            count++;
                        }
                    }
                }
            }
            oldpi = pi;
            oldpj = pj;
            cp[i] = total_length;
            rp[i] = (cp[i] - cp[i - 1]) / (ap[i] - ap[i - 1]);
        }
        csp->chsz = count;

        if (count == 0) printf("\nWARNING: no overlap with DEM cells for channel %i.\n", chseg);

        // Set up other channel rasters and fill in values
        csp->ChanX = new int[csp->chsz]();
        csp->ChanY = new int[csp->chsz]();
        csp->Chandx = new double[csp->chsz]();
        csp->Shalf = new double[csp->chsz]();
        csp->A = new double[csp->chsz]();
        csp->NewA = new double[csp->chsz]();
        csp->Chainage = new double[csp->chsz]();
        csp->ChanQ = new double[csp->chsz](); // only used to record Q values for output in profile - not used in calc
        csp->ChanWidth = new double[csp->chsz]();
        csp->ChanN = new double[csp->chsz]();
        csp->Q_Val = new double[csp->chsz]();
        csp->BankZ = new double[csp->chsz]();
        csp->Q_Name = new char[csp->chsz * 80]();
        csp->Q_Ident = new int[csp->chsz]();

        for (i = 0; i < csp->chsz; i++)
        {
            csp->ChanX[i] = xpi[i];
            csp->ChanY[i] = ypi[i];
        }


        // Find chainage and dx along channel
        for (i = 0; i < csp->chsz - 1; i++)
            csp->Chandx[i] = sqrt(pow(Parptr->dx * (csp->ChanX[i] - csp->ChanX[i + 1]), 2) +
                pow(Parptr->dy * (csp->ChanY[i] - csp->ChanY[i + 1]), 2));

        // assume dx for last cell is same as last segment
        csp->Chandx[csp->chsz - 1] = sqrt(pow(Parptr->dx * (csp->ChanX[csp->chsz - 1] - csp->ChanX[csp->chsz - 2]), 2) +
            pow(Parptr->dy * (csp->ChanY[csp->chsz - 1] - csp->ChanY[csp->chsz - 2]), 2));
        csp->Chainage[0] = 0;

        // add up dx to get chainage
        for (i = 1; i < csp->chsz; i++)csp->Chainage[i] = csp->Chainage[i - 1] + csp->Chandx[i - 1];

        // Fill in channel mask
        for (i = 0; i < csp->chsz; i++)
        {
            pi = csp->ChanX[i];
            pj = csp->ChanY[i];

            // renumber channel mask
            Arrptr->ChanMask[pi + pj * Parptr->xsz] = i;
            // set bank level
            csp->BankZ[i] = Arrptr->DEM[pi + pj * Parptr->xsz];
            // mark segment mask
            Arrptr->SegMask[pi + pj * Parptr->xsz] = chseg;
        }

        // Adjust chainage calcs so that it is independent of cell size. ####
        if (Statesptr->chainagecalc == ON)
        {
            if (*verbose == ON) printf("Cell size independent channel chainage calculations are ON.\n");
            // adjust each dx by ratio calculated previously
            for (i = 0; i < npoints; i++)
            {
                for (j = 0; j < csp->chsz; j++)
                {
                    if (csp->Chainage[j] >= cp[i] && csp->Chainage[j] < (cp[i + 1]))
                    {
                        csp->Chandx[j] = csp->Chandx[j] / rp[i + 1]; // adjust by ratio calculated previously
                    }
                }
            }
            // calc last seg dx as same as penultimate
            csp->Chandx[csp->chsz - 1] = csp->Chandx[csp->chsz - 2];

            // add up chainage again
            for (i = 1; i < csp->chsz; i++)csp->Chainage[i] = csp->Chainage[i - 1] + csp->Chandx[i - 1];

            // also adjust original xs chainage as this was based on centre cell distance
            for (i = 1; i < npoints; i++) cp[i] = ap[i];
        }
        else if (*verbose == ON) printf("Cell size independent channel chainage calculations are OFF.\n");




        // Interpolate width, Mannings n
        i1 = 0; i2 = 0;
        while (i2 < npoints - 1)
        {
            for (i2 = i1 + 1; i2 < npoints; i2++)
            {
                if (wp[i2] > 0.0) break; // loop until next non zero value found
            }
            for (i = 0; i < csp->chsz; i++)
            {
                if (csp->Chainage[i] >= cp[i1] && csp->Chainage[i] <= (cp[i2] + 0.1)) // add 0.1m to cp[] to get last point
                {
                    pi = csp->ChanX[i];
                    pj = csp->ChanY[i];
                    csp->ChanWidth[i] = wp[i1] + (wp[i2] - wp[i1]) * (csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
                    csp->ChanN[i] = np[i1] + (np[i2] - np[i1]) * (csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
                    if (i == csp->chsz - 1 && chseg != 0)
                    {
                        // special case where trib junction with main channel. Do not overwrite main channel bed elevation but
                        // otherwise record channel info for trib end point in dummy node (allows channel BC link for diffusive
                        // and correct slope calc for kinematic).
                        csp->JunctionDEM = hp[i1] + (hp[i2] - hp[i1]) * (csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
                    }
                    else
                    {
                        Arrptr->DEM[pi + pj * Parptr->xsz] = hp[i1] + (hp[i2] - hp[i1]) * (csp->Chainage[i] - cp[i1]) / (cp[i2] - cp[i1]);
                    }
                }
            }
            i1 = i2;
        }


        // Interpolate slope
        for (i = 0; i < csp->chsz; i++)
        {
            // find gradient of segment
            if (chseg == 0) // main channel
            {
                if (i == csp->chsz - 1) // last point
                {
                    // special case for last point as we can only know the slope of the segment behind it
                    grad = (Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i - 1] + csp->ChanY[i - 1] * Parptr->xsz])
                        / csp->Chandx[i - 1];
                }
                else // all other points
                {
                    grad = (Arrptr->DEM[csp->ChanX[i + 1] + csp->ChanY[i + 1] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
                        / csp->Chandx[i];
                }
            }
            else // trib special case at end due to junction
            {
                if (i == csp->chsz - 2) // last but one point
                {
                    // tribs use dummy node for last point - ie junction ##
                    grad = (csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
                        / csp->Chandx[i];
                }
                else if (i == csp->chsz - 1) // last point
                {
                    // tribs use dummy node for last point - ie junction ##
                    grad = (csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i - 1] + csp->ChanY[i - 1] * Parptr->xsz])
                        / csp->Chandx[i - 1];
                }
                else // all other points as main channel
                {
                    grad = (Arrptr->DEM[csp->ChanX[i + 1] + csp->ChanY[i + 1] * Parptr->xsz] - Arrptr->DEM[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz])
                        / csp->Chandx[i];
                }
            }


            // check if slope is positive or negative
            if (grad >= 0)
            {
                if (Statesptr->diffusive == ON)
                {
                    // only keep uphill slopes for diffusive
                    csp->Shalf[i] = -1 * sqrt(fabs(grad));
                }
                else
                {
                    // for kinematic we just pretend it is a downhill
                    csp->Shalf[i] = sqrt(fabs(grad));
                    // also warn user, in case they don't know
                    if (*verbose == ON) printf("\nWARNING: Kinematic solver BUT uphill slope at point %i for channel %i.\n", i - 1, chseg);
                }
            }
            else
            {
                csp->Shalf[i] = sqrt(fabs(grad));
            }
        }

        // Fill in Q boundary conditions from _tmp arrays
        for (i = 0; i < npoints; i++) // loop through the input cross-sections
        {
            if (Q_Ident_tmp[i] != 0)  // if there is a boundary condition
            {
                for (j = 0; j < csp->chsz; j++) // loop through the cross-sections mapped onto the dem space
                {
                    if ((cp[i] + 0.1) >= csp->Chainage[j] && (cp[i] + 0.1) < (csp->Chainage[j] + csp->Chandx[j]))
                        // make sure we only apply the BC to one point
                    {
                        csp->Q_Ident[j] = Q_Ident_tmp[i]; // copy type across
                        csp->Q_Val[j] = qp[i];            // copy value across
                        if (Q_Ident_tmp[i] == 3 || Q_Ident_tmp[i] == 5 || Q_Ident_tmp[i] == 8)	// check if BC type has name
                        {
                            strcpy(csp->Q_Name + j * 80, Q_Name_tmp + i * 80); // copy name across
                        }
                        if (Q_Ident_tmp[i] == 7) // only for tribs
                        {
                            /*
                              record in the trib data the location and segment it links to
                              ChannelSegments[trib[i]].Next_Segment_Loc=j; // CCS
                              ChannelSegments[trib[i]].Next_Segment=chseg; // CCS

                              ^^ The above code is left commented out to show why the following QID7_Store vector is needed.
                              As these terms need to be written to an instance of ChannelSegmentType not pointed to by csp
                              (ChannelSegments[trib[i]]), we store them in the QID7 vector and use the UpdateChannelVector function
                              to move the contents to the correct place after LoadRiver has finished. // CCS
                            */

                            QID7_Store QID7_tmp; // CCS create temp instance of QID7_Store and populate struct:
                            QID7_tmp.chseg = chseg;
                            QID7_tmp.Next_Segment_Loc = j;
                            QID7_tmp.trib = trib[i];
                            if (Statesptr->multiplerivers == ON)
                            {
                                QID7_tmp.RiverID = RiversIndexVecPtr->size(); /*#CCS# this allows us to keep track of which river we are in when later using QID7.
                                                                              1st river ID will be 0, 2nd will be 1, etc.*/
                            }
                            QID7_Vec_Ptr->push_back(QID7_tmp); // CCS push_back temp instance into external vector
                        }
                    }
                }
            }
        }

        // release memory used for temporary variables
        delete[] xp;
        delete[] yp;
        delete[] wp;
        delete[] np;
        delete[] hp;
        delete[] xpi;
        delete[] ypi;
        delete[] cp;
        delete[] Q_Name_tmp;

        if (Statesptr->diffusive == 1 && csp->Q_Ident[csp->chsz - 1] == 0)
        {
            if (*verbose == ON) printf("\nWARNING: Channel %i has no d/s BC using FREE with channel slope\n", chseg); // warn user that for diffusive no BC is set so using free
            csp->Q_Ident[csp->chsz - 1] = 1;  // copy type across
            csp->Q_Val[csp->chsz - 1] = -1;   // copy value across
        }

        ChannelSegmentsVecPtr->push_back(tmp_chan);

    } // end of channel segment loop


    fclose(fp);
    if (*verbose == ON) printf("Done.\n\n");

    return;
}
//加载索引
void River::Load_river_pipe_index(Fnames* Fnameptr, States* Statesptr)
{
    ifstream file(Fnameptr->river_pipe_index_name);
    if (!file.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return ;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        // 使用空格分隔键和值
        if (iss >> key >> value) {
            // 在这里你可以对key和value进行处理，比如存储到map中
            riverBoudaryIndex* temp = get_boundary((char*)key.c_str());
            river_pipe_index[key] = temp;
        }
    }
    file.close();
}

//----------------------------------------------------------------------------
// Calculate energy slope using mannings, slope^1/2, channel width, water depth
// and initial Q
double River::CalcEnergySlope(double n, double w, double h, double Q)
{
    double s;
    s = pow((Q * n) / (w * pow(h, (5. / 3.))), 2);
    return(s);
}

void River::UpdateChannelsVector(States* Statesptr, ChannelSegmentType* CSTypePtr, vector<QID7_Store>* QID7_Vec_Ptr, QID7_Store* QID7Ptr, int* RiversIndexPtr)
{
    int vecsize, i, n;
    if (Statesptr->ChannelPresent == OFF) return;
    vecsize = QID7_Vec_Ptr->size();
    for (i = 0; i < vecsize; i++)
    {
        if (Statesptr->multiplerivers == OFF)
        {
            n = QID7Ptr[i].trib;
            CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg;
            CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
        }
        else if (Statesptr->multiplerivers == ON)
        {
            if (QID7Ptr[i].RiverID == 0)
            {
                n = QID7Ptr[i].trib;
                CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg;
                CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
            }
            else
            {
                n = QID7Ptr[i].trib + RiversIndexPtr[QID7Ptr[i].RiverID - 1]; //Make sure we write to the correct CSTypePtr index when not in first river.
                CSTypePtr[n].Next_Segment = QID7Ptr[i].chseg + RiversIndexPtr[QID7Ptr[i].RiverID - 1]; //Again, when not in first river we need to ensure correct ID of next segment.
                CSTypePtr[n].Next_Segment_Loc = QID7Ptr[i].Next_Segment_Loc;
            }
        }
    }

    return;
}

void River::ChannelQ_Diff(double deltaT, States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, ChannelSegmentType* ChannelSegments, Arrays* Arrptr, vector<int>* RiversIndexVecPtr, int* RiversIndexPtr)
{
    // main arrays for diffusive
    double* x, * xn, * f, ** J, ** JinvL;
    int* indx;
    double d, res;

    int i;		// loop counter
    int chseg;	// channel segment loop counter
    int nodes;	// temporary variable to hold number of cells/nodes in segment
    double Qbc;	// flow value used to hold first channel element flow area at start (u/s boundary condition)
    double WSbc;  // value used to hold last channel element flow h (d/s boundary condition)
    int HoutFREE = OFF; // flag to indicate calc h BC from slope on  the fly later.
    double qc;    //temporary flow value used to hold channel segment outflow
    ChannelSegmentType* csp;  //local pointer to channel segment
    int exitcrit = 1;	// variable for determing the exit criterion for the NR solver; see manual for details (TJF)

    // new diffusive solver declarations
    int maxiterations = 200;		// maximum number of iterations for solver
    int itcount = 0;				  // iteration counter

    int nriv, low, high; // CCS For multiple river loop
    double qc_total; // CCS

    double* qc_store = new double[RiversIndexVecPtr->size()]; //CCS preallocate array for storing parallel channel solver output 

#pragma omp parallel for private(high, low, chseg, itcount, HoutFREE, csp, nodes, x, xn, f, indx, J, i, JinvL, Qbc, WSbc, res, qc) //parrallelised by CCS
    for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS 
    {
        high = RiversIndexPtr[nriv] - 1;
        if (nriv == 0)
        {
            low = 0;
        }
        else
        {
            low = RiversIndexPtr[nriv - 1];
        }
        // main loop for channel segments
        for (chseg = low; chseg <= high; chseg++) // CCS
        {
            // set iteration counter to zero for each channel
            itcount = 0;
            // set default flags
            HoutFREE = OFF;
            // set up local pointer to this segment
            csp = ChannelSegments + chseg;
            // setup temporary variable - just because it is easier to read code
            nodes = csp->chsz;

            // allocate memory for local solver arrays
            x = new double[nodes * 2]();
            xn = new double[nodes * 2]();
            f = new double[nodes * 2]();
            indx = new int[nodes * 2]();
            J = new double* [nodes * 2]();
            for (i = 0; i < 2 * nodes; i++) J[i] = new double[5]();
            JinvL = new double* [2 * nodes]();
            for (i = 0; i < 2 * nodes; i++) JinvL[i] = new double[5]();

            // calculate Areas from water height and channel width
            for (i = 0; i < nodes; i++)
            {
                if (chseg != low && i == nodes - 1) { 
                    x[2 * i] = xn[2 * i] = csp->JunctionH * csp->ChanWidth[i]; 
                } // trib and last point - ie dummy junction node // CCS
                else {
                    x[2 * i] = xn[2 * i] = Arrptr->H[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] * csp->ChanWidth[i];
                    if ( x[2 * i] == 0)
                    {
                        x[2 * i] = xn[2 * i] =0.01* csp->ChanWidth[i];//如果这里是0，之后calcJ会有inf，lisflood的一维渠道模块写的是真不行，。
                    }
                    if (isnan(x[2 * i])) {
                        printf("");
                    }
                } // all other points
            }

            // fill in Q values
            for (i = 0; i < nodes; i++)
            {
                // use flows from last iteration as start point for this one.
                // for first timestep, these are already filled in by CalcChannelStartQ()
                x[2 * i + 1] = xn[2 * i + 1] = csp->ChanQ[i];
            }

            // Fix first channel element with inflow area, ie u/s BC
            if (csp->Q_Ident[0] == 4) Qbc = csp->Q_Val[0];  // fixed Q inflow
            if (csp->Q_Ident[0] == 5) Qbc =Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[0]], Solverptr->t); //interpolate Q from hydrograph

            // work out d/s boundary condition ie d/s H 
            if (chseg > low) //trib // CCS
            {
                // get H from d/s channel - add DEM value to convert to elevation - to transfer as d/s BC to trib (other bc types are elevations). 
                // Also covers the case when d/s trib bed elevation is different from main channel.
                WSbc = Arrptr->H[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc] + ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc] * Parptr->xsz]
                    + Arrptr->DEM[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc] + ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc] * Parptr->xsz];
            }
            else  // main channel d/s BC
            {
                if (csp->Q_Ident[nodes - 1] == 1) //free boundary, calc h from slope ## This BC needs some stability work
                {
                    // set flag so we can calc h from the flow and slope "on the fly" in CalcF() function
                    WSbc = 0; // set dummy value as will be calculated later
                    HoutFREE = ON;
                }
                else if (csp->Q_Ident[nodes - 1] == 2) // fixed H out 
                {
                    // note - we will need to subtract DEM Elev to get h from water elevation entered
                    WSbc = csp->Q_Val[nodes - 1];
                }
                else if (csp->Q_Ident[nodes - 1] == 3) //interpolate H from stage hydrograph
                {
                    // note - we will need to subtract DEM Elev to get h from water elevation entered
                    WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], Solverptr->t);
                }
                else if (csp->Q_Ident[nodes - 1] == 8) //interpolate H from rating curve ## This BC needs some stability work
                {
                    // set flag to 2 so we can calc h from the stage discharge curve "on the fly" in CalcF() function
                    HoutFREE = 2;
                    WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], x[2 * nodes - 3]);
                }
                else
                {
                    printf("WARNING: Main Channel has no d/s BC\n");
                    // already checked in input function, and set to FREE by default, but leave this here just in case of any odd bugs
                }
            }

            // Use implicit newton raphson scheme to solve 
            // Utilises LU decomposition - Crout's method to solve the banded diagonal matrix of linear equations
            do
            {
                res = 0.0; // reset the res incase we use the max rather than the norm
                // setup function matrix
                calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low);
                for (int i = 0; i < nodes * 2; i++)
                {
                    if (isnan(f[i])) {
                        printf("");
                    }
                }
                // setup Jacobian matrix - Q determinant of function
                calcJ(x, xn, J, deltaT, csp, chseg, HoutFREE);
                for (int i = 0; i < nodes * 2; i++)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        if (isnan(J[i][j])|| isinf(J[i][j])) {
                            printf("");
                        }
                    }
                }
                bandec(J, 2 * nodes, 2, 2, JinvL, indx, d);
                for (int i = 0; i < nodes * 2; i++)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        if (isnan(J[i][j])) {
                            printf("");
                        }
                    }
                }
                banbks(J, 2 * nodes, 2, 2, JinvL, indx, f); // On input, f is function vector from CalcF. On output, f is solution.
                //mprove(); // Routine for iterative improvement - not currently implemented (TJF).

                for (i = 0; i < 2 * nodes; i++) xn[i] -= f[i];
                for (i = 0; i < nodes; i++) if (xn[2 * i] < 1.0) xn[2 * i] = 1.0;

                switch (exitcrit) // Switch for different exit criteria to improve channel model solution. See manual for details.
                {
                case 1:
                    res = Tool::norm(f, 2 * nodes); // norm of x (the solution)
                    break;
                case 2:
                    for (i = 0; i < 2 * nodes; i++) res = Tool::getmax(fabs(f[i]), res); // max of x (the solution)
                    break;
                case 3:
                    calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low); // recalculate f(x) for exit criteria
                    res = Tool::norm(f, 2 * nodes);  // norm of f(x) (the function vector)
                    break;
                case 4:
                    calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low); // recalculate f(x) for exit criteria
                    for (i = 0; i < 2 * nodes; i++) res = Tool::getmax(fabs(f[i]), res); // max of f(x) (the function vector)
                    break;
                }
            } while (res > Solverptr->SolverAccuracy && itcount++ < maxiterations);
            if (itcount >= maxiterations) printf("WARNING: Max iterations exceeded diffusive channel at t=%.3f in channel %i.\n", Solverptr->t, chseg); //Warning for iterations

            // update H based on new areas calculated, divided by channel width
            for (i = 0; i < nodes; i++)
            {
                if (chseg != low && i == nodes - 1)	csp->JunctionH = xn[2 * i] / csp->ChanWidth[i];  // trib and last point - ie dummy junction node // CCS
                else						Arrptr->H[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] = xn[2 * i] / csp->ChanWidth[i]; // all other points

                csp->ChanQ[i] = xn[2 * i + 1]; // record Q for profile output
            }

            // get outflow from last channel segment
            qc = xn[2 * (nodes - 1) + 1];

            if (chseg > low) // trib #CCS#
            {
                // set the QVal of the next segment to the outflow from the end of this segment
                ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc] = qc;
            }
            else // main channel
            {
                // store the main channel outflow CCS
                qc_store[nriv] = qc;
                // update downstream water depth
                Solverptr->Hds = Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz];
            }

            // clean up memory
            for (i = 0; i < 2 * nodes; i++) delete[] J[i];
            for (i = 0; i < 2 * nodes; i++) delete[] JinvL[i];
            delete[] J;
            delete[] JinvL;
            delete[] x;
            delete[] xn;
            delete[] f;
            delete[] indx;

        } // end of main river channel segment loop
    }
    //Calculate QChanOut from all rivers CCS
    qc_total = 0;
    for (int n = 0; n < (int)RiversIndexVecPtr->size(); n++)
    {
        qc_total = qc_total + qc_store[n];
    }
    delete[] qc_store; // delete qc store array

    BCptr->QChanOut = qc_total;

    return;
}
;
void River::ChannelQ_Diff1(double deltaT, States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, ChannelSegmentType* ChannelSegments, Arrays* Arrptr, vector<int>* RiversIndexVecPtr, int* RiversIndexPtr)
{
    // main arrays for diffusive
    double* x, * xn, * f, ** J, ** JinvL;
    int* indx;
    double d, res;

    int i;		// loop counter
    int chseg;	// channel segment loop counter
    int nodes;	// temporary variable to hold number of cells/nodes in segment
    double Qbc;	// flow value used to hold first channel element flow area at start (u/s boundary condition)
    double WSbc;  // value used to hold last channel element flow h (d/s boundary condition)
    int HoutFREE = OFF; // flag to indicate calc h BC from slope on  the fly later.
    double qc;    //temporary flow value used to hold channel segment outflow
    ChannelSegmentType* csp;  //local pointer to channel segment
    int exitcrit = 1;	// variable for determing the exit criterion for the NR solver; see manual for details (TJF)

    // new diffusive solver declarations
    int maxiterations = 200;		// maximum number of iterations for solver
    int itcount = 0;				  // iteration counter

    int nriv, low, high; // CCS For multiple river loop
    double qc_total; // CCS

    double* qc_store = new double[RiversIndexVecPtr->size()]; //CCS preallocate array for storing parallel channel solver output 

#pragma omp parallel for private(high, low, chseg, itcount, HoutFREE, csp, nodes, x, xn, f, indx, J, i, JinvL, Qbc, WSbc, res, qc) //parrallelised by CCS
    for (nriv = 0; nriv < (int)RiversIndexVecPtr->size(); nriv++) // CCS 
    {
        high = RiversIndexPtr[nriv] - 1;
        if (nriv == 0)
        {
            low = 0;
        }
        else
        {
            low = RiversIndexPtr[nriv - 1];
        }
        // main loop for channel segments
        for (chseg = low; chseg <= high; chseg++) // CCS
        {
            // set iteration counter to zero for each channel
            itcount = 0;
            // set default flags
            HoutFREE = OFF;
            // set up local pointer to this segment
            csp = ChannelSegments + chseg;
            // setup temporary variable - just because it is easier to read code
            nodes = csp->chsz;

            // allocate memory for local solver arrays
            x = new double[nodes * 2]();
            xn = new double[nodes * 2]();
            f = new double[nodes * 2]();
            indx = new int[nodes * 2]();
            J = new double* [nodes * 2]();
            for (i = 0; i < 2 * nodes; i++) J[i] = new double[5]();
            JinvL = new double* [2 * nodes]();
            for (i = 0; i < 2 * nodes; i++) JinvL[i] = new double[5]();

            // calculate Areas from water height and channel width
            for (i = 0; i < nodes; i++)
            {
                if (chseg != low && i == nodes - 1) {
                    x[2 * i] = xn[2 * i] = csp->JunctionH * csp->ChanWidth[i];
                } // trib and last point - ie dummy junction node // CCS
                else {
                    x[2 * i] = xn[2 * i] = Arrptr->H[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] * csp->ChanWidth[i];
                    if (x[2 * i] == 0)
                    {
                        x[2 * i] = xn[2 * i] = 0.01 * csp->ChanWidth[i];//如果这里是0，之后calcJ会有inf，lisflood的一维渠道模块写的是真不行，。
                    }
                    if (isnan(x[2 * i])) {
                        printf("");
                    }
                } // all other points
            }

            // fill in Q values
            for (i = 0; i < nodes; i++)
            {
                // use flows from last iteration as start point for this one.
                // for first timestep, these are already filled in by CalcChannelStartQ()
                x[2 * i + 1] = xn[2 * i + 1] = csp->ChanQ[i];
            }

            // Fix first channel element with inflow area, ie u/s BC
            if (csp->Q_Ident[0] == 4) Qbc = csp->Q_Val[0];  // fixed Q inflow
            if (csp->Q_Ident[0] == 5) Qbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[0]], Solverptr->t); //interpolate Q from hydrograph

            // work out d/s boundary condition ie d/s H 
            if (chseg > low) //trib // CCS
            {
                // get H from d/s channel - add DEM value to convert to elevation - to transfer as d/s BC to trib (other bc types are elevations). 
                // Also covers the case when d/s trib bed elevation is different from main channel.
                WSbc = Arrptr->H[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc] + ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc] * Parptr->xsz]
                    + Arrptr->DEM[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc] + ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc] * Parptr->xsz];
            }
            else  // main channel d/s BC
            {
                if (csp->Q_Ident[nodes - 1] == 1) //free boundary, calc h from slope ## This BC needs some stability work
                {
                    // set flag so we can calc h from the flow and slope "on the fly" in CalcF() function
                    WSbc = 0; // set dummy value as will be calculated later
                    HoutFREE = ON;
                }
                else if (csp->Q_Ident[nodes - 1] == 2) // fixed H out 
                {
                    // note - we will need to subtract DEM Elev to get h from water elevation entered
                    WSbc = csp->Q_Val[nodes - 1];
                }
                else if (csp->Q_Ident[nodes - 1] == 3) //interpolate H from stage hydrograph
                {
                    // note - we will need to subtract DEM Elev to get h from water elevation entered
                    WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], Solverptr->t);
                }
                else if (csp->Q_Ident[nodes - 1] == 8) //interpolate H from rating curve ## This BC needs some stability work
                {
                    // set flag to 2 so we can calc h from the stage discharge curve "on the fly" in CalcF() function
                    HoutFREE = 2;
                    WSbc = Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes - 1]], x[2 * nodes - 3]);
                }
                else
                {
                    printf("WARNING: Main Channel has no d/s BC\n");
                    // already checked in input function, and set to FREE by default, but leave this here just in case of any odd bugs
                }
            }

            // Use implicit newton raphson scheme to solve 
            // Utilises LU decomposition - Crout's method to solve the banded diagonal matrix of linear equations
            do
            {
                res = 0.0; // reset the res incase we use the max rather than the norm
                // setup function matrix
                calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low);
                for (int i = 0; i < nodes * 2; i++)
                {
                    if (isnan(f[i])) {
                        printf("");
                    }
                }
                // setup Jacobian matrix - Q determinant of function
                calcJ(x, xn, J, deltaT, csp, chseg, HoutFREE);
                for (int i = 0; i < nodes * 2; i++)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        if (isnan(J[i][j]) || isinf(J[i][j])) {
                            printf("");
                        }
                    }
                }
                bandec(J, 2 * nodes, 2, 2, JinvL, indx, d);
                for (int i = 0; i < nodes * 2; i++)
                {
                    for (int j = 0; j < 5; j++)
                    {
                        if (isnan(J[i][j])) {
                            printf("");
                        }
                    }
                }
                banbks(J, 2 * nodes, 2, 2, JinvL, indx, f); // On input, f is function vector from CalcF. On output, f is solution.
                //mprove(); // Routine for iterative improvement - not currently implemented (TJF).

                for (i = 0; i < 2 * nodes; i++) xn[i] -= f[i];
                for (i = 0; i < nodes; i++) if (xn[2 * i] < 1.0) xn[2 * i] = 1.0;

                switch (exitcrit) // Switch for different exit criteria to improve channel model solution. See manual for details.
                {
                case 1:
                    res = Tool::norm(f, 2 * nodes); // norm of x (the solution)
                    break;
                case 2:
                    for (i = 0; i < 2 * nodes; i++) res = Tool::getmax(fabs(f[i]), res); // max of x (the solution)
                    break;
                case 3:
                    calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low); // recalculate f(x) for exit criteria
                    res = Tool::norm(f, 2 * nodes);  // norm of f(x) (the function vector)
                    break;
                case 4:
                    calcF(x, xn, f, deltaT, csp, Parptr, Arrptr, Qbc, chseg, WSbc, HoutFREE, Solverptr, low); // recalculate f(x) for exit criteria
                    for (i = 0; i < 2 * nodes; i++) res = Tool::getmax(fabs(f[i]), res); // max of f(x) (the function vector)
                    break;
                }
            } while (res > Solverptr->SolverAccuracy && itcount++ < maxiterations);
            if (itcount >= maxiterations) printf("WARNING: Max iterations exceeded diffusive channel at t=%.3f in channel %i.\n", Solverptr->t, chseg); //Warning for iterations

            // update H based on new areas calculated, divided by channel width
            for (i = 0; i < nodes; i++)
            {
                if (chseg != low && i == nodes - 1)	csp->JunctionH = xn[2 * i] / csp->ChanWidth[i];  // trib and last point - ie dummy junction node // CCS
                else						Arrptr->temp_H[csp->ChanX[i] + csp->ChanY[i] * Parptr->xsz] = xn[2 * i] / csp->ChanWidth[i]; // all other points

                csp->ChanQ[i] = xn[2 * i + 1]; // record Q for profile output
            }

            // get outflow from last channel segment
            qc = xn[2 * (nodes - 1) + 1];

            if (chseg > low) // trib #CCS#
            {
                // set the QVal of the next segment to the outflow from the end of this segment
                ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc] = qc;
            }
            else // main channel
            {
                // store the main channel outflow CCS
                qc_store[nriv] = qc;
                // update downstream water depth
                Solverptr->Hds = Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz];
            }

            // clean up memory
            for (i = 0; i < 2 * nodes; i++) delete[] J[i];
            for (i = 0; i < 2 * nodes; i++) delete[] JinvL[i];
            delete[] J;
            delete[] JinvL;
            delete[] x;
            delete[] xn;
            delete[] f;
            delete[] indx;

        } // end of main river channel segment loop
    }
    //Calculate QChanOut from all rivers CCS
    qc_total = 0;
    for (int n = 0; n < (int)RiversIndexVecPtr->size(); n++)
    {
        qc_total = qc_total + qc_store[n];
    }
    delete[] qc_store; // delete qc store array

    BCptr->QChanOut = qc_total;

    return;
}
;
void River::calcJ(double* x, double* xn, double** J, double dt, ChannelSegmentType* csp, int chseg, int HoutFREE)
{
    double a, q, dx;
    int i;
    int n = csp->chsz;
    double th = 1;
    double w, mann;
    double n2w43th, fu3;

    fu3 = 5. / 3.; // five divided by three for later use

    for (i = 0; i < 5; i++) J[0][i] = 0;
    J[0][3] = 1;

    for (i = 0; i < 2 * n - 3; i += 2)
    {
        dx = csp->Chandx[i / 2];

        w = csp->ChanWidth[i / 2];
        n2w43th = csp->ChanN[i / 2] * csp->ChanN[i / 2] * pow(w, (4. / 3.)) * th;

        // Setting the odd rows of the compact form of the Jacobian [Continuity Equation]
        J[i + 1][0] = 0;
        J[i + 1][1] = 0.5 / dt;
        J[i + 1][2] = -th / dx;
        J[i + 1][3] = J[i + 1][1];
        J[i + 1][4] = -J[i + 1][2];

        // Setting the even rows of the compact form of the Jacobian [Momentum Equation]
        a = th * (xn[i] + xn[i + 2]) / 2 + (1 - th) * (x[i] + x[i + 2]) / 2;
        q = th * (xn[i + 1] + xn[i + 3]) / 2 + (1 - th) * (x[i + 1] + x[i + 3]) / 2;

        J[i + 2][0] = th / (w * dx) + 10 * n2w43th * fabs(q) * q / (6 * pow(a, (13. / 3.)));//这个如果xn可以为0，a,可以是0，J[i + 2][0] 就可以是inf
        if (isinf(J[i + 2][0]))
        {
            printf("    ");
        }
        J[i + 2][1] = -n2w43th * fabs(q) / pow(a, (10. / 3.));
        J[i + 2][2] = -J[i + 2][0];
        J[i + 2][3] = J[i + 2][1];
        J[i + 2][4] = 0;
    }

    // if FREE boundary, reset the final lines of the Jacobian as the derivatives of Manning's equation rather than
    // the diffusive form of the de St. Venant.
    if (HoutFREE == ON) // only if HoutFREE flagged ON
    {
        for (i = 0; i < 5; i++) J[2 * n - 2][i] = 0; // reset to 0
        if (csp->Q_Val[csp->chsz - 1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
        {
            // Combining Manning's and diffusive continuity and momentum for final section
            // derivative of combination wrt Q
            mann = (pow(xn[2 * n - 2], (5. / 3.)) * csp->Shalf[n - 1]) / (csp->ChanN[n - 1] * pow(csp->ChanWidth[n - 1], (2. / 3.)));
            J[2 * n - 2][2] = -pow(csp->ChanN[n - 1], 2.0) * pow(csp->ChanWidth[n - 2], (4. / 3.)) * pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-10. / 3.)) * (mann + xn[2 * n - 3]) / 2;

            // Combining Manning's and diffusive continuity and momentum for final section
            // derivative of combination wrt A
            J[2 * n - 2][3] = -(1 / (csp->Chandx[n - 1] * csp->ChanWidth[n - 1]));
            J[2 * n - 2][3] -= pow(csp->ChanN[n - 1], 2.0) * pow(csp->ChanWidth[n - 2], (4. / 3.)) * pow((mann + xn[2 * n - 3]) / 2, 2.0) * (10. / 6.) * pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-13. / 3.));
            J[2 * n - 2][3] += pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-10. / 3.)) * ((mann + xn[2 * n - 3]) / 2) * -fu3 * pow(xn[2 * n - 2], (2. / 3.)) * (csp->Shalf[n - 1] / (csp->ChanN[n - 1] * pow(csp->ChanWidth[n - 1], (2. / 3.))));

            // Old version
            //J[2*n-2][2]+=(csp->ChanN[n-1]*pow((1-alpha)*csp->ChanWidth[n-2]+alpha*csp->ChanWidth[n-1],(2./3.)))/csp->Shalf[n-1]; // derivative of Mannings wrt Q
            //J[2*n-2][3]=-fu3*(pow((1-alpha)*xn[2*n-4]+alpha*xn[2*n-2],(2./3.))); // derivative of Mannings wrt A
        }
        else // user supplied slope
        {
            // Combining Manning's and diffusive continuity and momentum for final section
            // derivative of combination wrt Q
            mann = (pow(xn[2 * n - 2], (5. / 3.)) * sqrt(csp->Q_Val[n - 1])) / (csp->ChanN[n - 1] * pow(csp->ChanWidth[n - 1], (2. / 3.)));
            J[2 * n - 2][2] = -pow(csp->ChanN[n - 1], 2.0) * pow(csp->ChanWidth[n - 2], (4. / 3.)) * pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-10. / 3.)) * (mann + xn[2 * n - 3]) / 2;

            // Combining Manning's and diffusive continuity and momentum for final section
            // derivative of combination wrt A
            J[2 * n - 2][3] = -(1 / (csp->Chandx[n - 1] * csp->ChanWidth[n - 1]));
            J[2 * n - 2][3] -= pow(csp->ChanN[n - 1], 2.0) * pow(csp->ChanWidth[n - 2], (4. / 3.)) * pow((mann + xn[2 * n - 3]) / 2, 2.0) * (10. / 6.) * pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-13. / 3.));
            J[2 * n - 2][3] += pow((xn[2 * n - 2] + xn[2 * n - 4]) / 2, (-10. / 3.)) * ((mann + xn[2 * n - 3]) / 2) * -fu3 * pow(xn[2 * n - 2], (2. / 3.)) * (sqrt(csp->Q_Val[n - 1]) / (csp->ChanN[n - 1] * pow(csp->ChanWidth[n - 1], (2. / 3.))));

            //J[2*n-2][2]+=(csp->ChanN[n-1]*pow((1-alpha)*csp->ChanWidth[n-2]+alpha*csp->ChanWidth[n-1],(2./3.)))/sqrt(csp->Q_Val[n-1]); // derivative of Mannings wrt Q
            //J[2*n-2][3]+=-fu3*(pow((1-alpha)*xn[2*n-4]+alpha*xn[2*n-2],(2./3.))); // derivative of Mannings wrt A
        }
    }

    for (i = 0; i < 5; i++) J[2 * n - 1][i] = 0;
    J[2 * n - 1][1] = 1;

    return;
}
//====================================================================================================
// bandec(), banbks(), SWAP() and mprove() below are all from Chapter 2.4 in the book 
// "Numerical Recipes in C" p51-54. These allow solution of band diagonal 
// linear systems by LU decomposition (Crout's method). Rewritten for matrix 
// indexes starting at zero rather than 1.
//---------------------------------------------------------------------------

#define TINY 1.e-20

void River::bandec(double** a, int n, int m1, int m2, double** al, int indx[], double& d)
// Given an n x n band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal rows,
// compactly stored in the array a[1..n][1..m1+m2+1]. The diagonal elements are in a[1..n][m1+1]. 
// Subdiagonal elements are in a[j..n][1..m1] (with j > 1 appropriate to the number of elements 
// on each subdiagonal). Superdiagonal elements are in a[1..j][m1+2..m1+m2+1] with j < n appropriate 
// to the number of elements on each superdiagonal. This routine constructs an LU decomposition of 
// a rowwise permutation of A. The upper triangular matrix replaces a, while the lower triangular 
// matrix is returned in al[1..n][1..m1]. indx[1..n] is an output vector which records the row 
// permutation effected by the partial pivoting; d is output as +-1 depending on whether the 
// number of row interchanges was even or odd, respectively. This routine is used in combination 
// with banbks to solve band-diagonal sets of equations.
{
    int i, j, k, l;
    int mm;
    double dum;

    mm = m1 + m2+ 1 ;/**/
    l = m1;
    for (i = 1; i <= m1; i++) // Rearrange the storage a bit.
    {
        for (j = m1 + 2 - i; j <= mm; j++) {
            a[i - 1][j - l - 1] = a[i - 1][j - 1];
            if (isnan(a[i - 1][j - l - 1]) || isinf(a[i - 1][j - l - 1]))
            {
                printf("    ");
            }
        }
        l--;
        for (j = mm - l; j <= mm; j++) a[i - 1][j - 1] = 0.0;
    }

    d = 1.0;
    l = m1;
    for (k = 1; k <= n; k++) // For each row...
    {
        dum = a[k - 1][0];
        i = k;
        if (l < n) l++;
        for (j = k + 1; j <= l; j++) // Find the pivot element.
        {
            if (fabs(a[j - 1][1 - 1]) > fabs(dum))
            {
                dum = a[j - 1][0];
                i = j;
            }
        }
        indx[k - 1] = i;
        if (dum == 0.0) a[k - 1][0] = TINY;

        //  Matrix is algorithmically singular, but proceed anyway with
        //  TINY pivot (desirable in some applications).

        if (i != k) // Interchange rows.
        {
            d = -(d);
            for (j = 1; j <= mm; j++) SWAP(a[k - 1][j - 1], a[i - 1][j - 1]);
        }
        for (i = k + 1; i <= l; i++) // Do the elimination.
        {
            dum = a[i - 1][0] / a[k - 1][0];
            if (isnan(dum) || isinf(dum))
            {
                printf("    ");
            }
            al[k - 1][i - k - 1] = dum;
            for (j = 2; j <= mm; j++) {
                a[i - 1][j - 2] = a[i - 1][j - 1] - dum * a[k - 1][j - 1];
                if (isnan(a[i - 1][j - 2])||isinf(a[i - 1][j - 2]))
                {
                    printf("    ");
                }
            }
            a[i - 1][mm - 1] = 0.0;
        }
    }
    return;
}

//---------------------------------------------------------------------------
// swap a and b
void River::SWAP(double& a, double& b)
{
    double dum = a;
    a = b;
    b = dum;
    return;
}

//---------------------------------------------------------------------------
void River::calcF(double* x, double* xn, double* f, double dt, ChannelSegmentType* csp, Pars* Parptr,
    Arrays* Arrptr, double Qin, int chseg, double WSout, int HoutFREE, Solver* Solverptr, int low)
{
    double q, a, dx;
    int i;

    int nodes = csp->chsz;
    double th = 1.0;
    double w, w0, w1;
    double w43, mann;
    double nSq;
    double Hout; //h calculated at downstream boundary condition

    for (i = 0; i < nodes - 1; i++)
    {

        dx = csp->Chandx[i];

        // need to use both channel widths
        w0 = csp->ChanWidth[i];
        w1 = csp->ChanWidth[i + 1];
        w = (w0 + w1) / 2.0;
        w43 = pow(w0, (4. / 3.));

        nSq = 1. / ((1. / csp->ChanN[i] + 1. / csp->ChanN[i + 1]) / 2.); // calc inverse average for 2 cross-sections
        nSq = nSq * nSq; // square for later use

        f[2 * i + 1] = (xn[2 * i] + xn[2 * i + 2] - x[2 * i] - x[2 * i + 2]) / (2. * dt);
        f[2 * i + 1] += th * (xn[2 * i + 3] - xn[2 * i + 1]) / dx;
        f[2 * i + 1] += (1 - th) * (x[2 * i + 3] - x[2 * i + 1]) / dx;  // if th=1, use fully implicit method

        // Add overbank flows - only if not in startup mode
        f[2 * i + 1] += BankQ(i + 1, csp, Parptr, Arrptr) / dx;

        // check to see if there are any inflows at the node and add them
        //
        if (i != 0)
            // make sure it is not start of channel inflow BC
        {
            if (csp->Q_Ident[i + 1] == 4)
            {
                // fixed flow 
                f[2 * i + 1] -= csp->Q_Val[i + 1] / dx;
            }
            if (csp->Q_Ident[i + 1] == 5)
            {
                // interpolate from hydrograph
                f[2 * i + 1] -= Tool::InterpBC(csp->QVarlist[(int)csp->Q_Val[i + 1]], Solverptr->t) / dx;
            }
            if (csp->Q_Ident[i + 1] == 7)
            {
                // tributary connection
                f[2 * i + 1] -= csp->Q_Val[i + 1] / dx;
            }
        }

        f[2 * i + 2] = pow(csp->Shalf[i], 2.0);
        if (csp->Shalf[i] < 0) f[2 * i + 2] = -f[2 * i + 2];

        // new variable width code
        f[2 * i + 2] -= (((th / w1) * xn[2 * i + 2]) - ((th / w0) * xn[2 * i])) / dx;
        f[2 * i + 2] -= ((((1 - th) / w1) * x[2 * i + 2]) - (((1 - th) / w0) * x[2 * i])) / dx; // if th=1, use fully implicit method

        a = th * (xn[2 * i] + xn[2 * i + 2]) / 2;
        a += (1 - th) * (x[2 * i] + x[2 * i + 2]) / 2;
        q = th * (xn[2 * i + 1] + xn[2 * i + 3]) / 2;
        q += (1 - th) * (x[2 * i + 1] + x[2 * i + 3]) / 2;
        f[2 * i + 2] -= nSq * w43 * q * fabs(q) * pow(a, (-10. / 3.0));
    }

    // boundary conditions
    // Qdiff of first element - Qin
    f[0] = xn[1] - Qin;


    if (HoutFREE == ON) // only if HoutFREE flagged ON
    {
        // calc h "on the fly" from slope for FREE BC
        if (csp->Q_Val[csp->chsz - 1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
        {
            // use a combination of previous iteration (xn) and previous timestep (x) for estimate of Q - set by th
            // if th=1, use fully implicit method  
            Hout = Tool::CalcA(csp->ChanN[nodes - 1], csp->Shalf[nodes - 1], csp->ChanWidth[nodes - 1], ((1 - th) * x[2 * nodes - 1]) + th * xn[2 * nodes - 1])
                / csp->ChanWidth[nodes - 1];

            // Recalculate f(x) at last section by substituting Manning's into diffusive form of momentum equation 
            f[2 * nodes - 2] = pow(csp->Shalf[nodes - 1], 2.0);
            f[2 * nodes - 2] -= (((th / csp->ChanWidth[nodes - 1]) * xn[2 * nodes - 2]) - ((th / csp->ChanWidth[nodes - 2]) * xn[2 * nodes - 4])) / csp->Chandx[nodes - 2];
            a = th * (xn[2 * nodes - 4] + xn[2 * nodes - 2]) / 2;
            mann = (pow(xn[2 * nodes - 2], (5. / 3.)) * csp->Shalf[nodes - 1]) / (csp->ChanN[nodes - 1] * pow(csp->ChanWidth[nodes - 1], (2. / 3.)));
            q = th * (xn[2 * nodes - 3] + mann) / 2;
            nSq = 1. / ((1. / csp->ChanN[nodes - 2] + 1. / csp->ChanN[nodes - 1]) / 2.); // calc inverse average for 2 cross-sections
            nSq = nSq * nSq;
            w43 = pow(csp->ChanWidth[nodes - 2], (4. / 3.));
            f[2 * nodes - 2] -= nSq * w43 * q * fabs(q) * pow(a, (-10. / 3.0));
            //f[2*nodes-2]=((1-alpha)*xn[2*nodes-3]+alpha*xn[2*nodes-1]*csp->ChanN[nodes-1]*pow((1-alpha)*csp->ChanWidth[nodes-2]+alpha*csp->ChanWidth[nodes-1],(2./3.)))/csp->Shalf[nodes-1]-pow((1-alpha)*xn[2*nodes-4]+alpha*xn[2*nodes-2],(5./3.));
        }
        else // else use user supplied slope
        {
            // use a combination of previous iteration (xn) and previous timestep (x) for estimate of Q - set by th
            // if th=1, use fully implicit method 
            Hout = Tool::CalcA(csp->ChanN[nodes - 1], sqrt(csp->Q_Val[nodes - 1]), csp->ChanWidth[nodes - 1], ((1 - th) * x[2 * nodes - 1]) + th * xn[2 * nodes - 1])
                / csp->ChanWidth[nodes - 1];

            // Recalculate f(x) at last section by substituting Manning's into diffusive form of momentum equation 
            f[2 * nodes - 2] = csp->Q_Val[nodes - 1];
            f[2 * nodes - 2] -= (((th / csp->ChanWidth[nodes - 1]) * xn[2 * nodes - 2]) - ((th / csp->ChanWidth[nodes - 2]) * xn[2 * nodes - 4])) / csp->Chandx[nodes - 2];
            a = th * (xn[2 * nodes - 4] + xn[2 * nodes - 2]) / 2;
            mann = (pow(xn[2 * nodes - 2], (5. / 3.)) * sqrt(csp->Q_Val[nodes - 1])) / (csp->ChanN[nodes - 1] * pow(csp->ChanWidth[nodes - 1], (2. / 3.)));
            q = th * (xn[2 * nodes - 3] + mann) / 2;
            nSq = 1. / ((1. / csp->ChanN[nodes - 2] + 1. / csp->ChanN[nodes - 1]) / 2.); // calc inverse average for 2 cross-sections
            nSq = nSq * nSq;
            w43 = pow(csp->ChanWidth[nodes - 2], (4. / 3.));
            f[2 * nodes - 2] -= nSq * w43 * q * fabs(q) * pow(a, (-10. / 3.0));
            //f[2*nodes-2]=((1-alpha)*xn[2*nodes-3]+alpha*xn[2*nodes-1]*csp->ChanN[nodes-1]*pow((1-alpha)*csp->ChanWidth[nodes-2]+alpha*csp->ChanWidth[nodes-1],(2./3.)))/sqrt(csp->Q_Val[nodes-1])-pow((1-alpha)*xn[2*nodes-4]+alpha*xn[2*nodes-2],(5./3.));
        }
    }
    else if (HoutFREE == 2) // HoutFREE flagged for stage discharge curve
    {
        // subtract dem to get water depth
        Hout = WSout - Arrptr->DEM[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz];
        // This BC needs some stability work
    }
    else // HoutFREE not flagged
    {
        // subtract dem to get water depth
        if (chseg == low) Hout = WSout - Arrptr->DEM[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz]; // CCS

        else         Hout = WSout - csp->JunctionDEM; //trib uses dummy junction node
    }

    // Set global variables to the value of Hout calculated here
    if (chseg != 0) csp->JunctionH = Hout;
    else Arrptr->H[csp->ChanX[nodes - 1] + csp->ChanY[nodes - 1] * Parptr->xsz] = Hout;

    // Area diff of last element - bc area
    f[2 * nodes - 1] = xn[2 * nodes - 2] - (Hout)*csp->ChanWidth[nodes - 1];

    return;
}
//---------------------------------------------------------------------------
// CALCULATE BANK FLOWS FOR ChannelQ()

double River::BankQ(int chani, ChannelSegmentType* ChannelSegments, Pars* Parptr, Arrays* Arrptr)
{
    double qbank0, q0, q1, q2, q3;
    int pi0, pj0;

    pi0 = ChannelSegments->ChanX[chani];
    pj0 = ChannelSegments->ChanY[chani];

    // Find Qbank0 from 4-neighbours, masking out channel flows
    q0 = Arrptr->Qx[pi0 + 1 + pj0 * (Parptr->xsz + 1)];

    q1 = -Arrptr->Qx[pi0 + pj0 * (Parptr->xsz + 1)];

    q2 = Arrptr->Qy[pi0 + (pj0 + 1) * (Parptr->xsz + 1)];

    q3 = -Arrptr->Qy[pi0 + pj0 * (Parptr->xsz + 1)];

    qbank0 = q0 + q1 + q2 + q3;

    return(qbank0);
}
//---------------------------------------------------------------------------

void River::banbks(double** a, int n, int m1, int m2, double** al, int indx[], double b[])
// Given the arrays a, al, and indx as returned from bandec, and given a right-hand side vector
// b[1..n], solves the band diagonal linear equations A . x = b. The solution vector x overwrites
// b[1..n]. The other input arrays are not modified, and can be left in place for successive calls
// with different right-hand sides.
{
    int i, k, l;
    int mm;
    double dum;

    mm = m1 + m2 + 1;
    l = m1;
    for (k = 1; k <= n; k++) // Forward substitution, unscrambling the permuted rows
    {                  // as we go.
        i = indx[k - 1];
        if (i != k) SWAP(b[k - 1], b[i - 1]);
        if (l < n) l++;
        for (i = k + 1; i <= l; i++)
        {
            b[i - 1] -= al[k - 1][i - k - 1] * b[k - 1];
            if (isnan(b[i - 1]))
            {
                printf("");
            }
        }
    }

    l = 1;
    for (i = n; i >= 1; i--) // Backsubstitution.
    {
        dum = b[i - 1];
        for (k = 2; k <= l; k++) dum -= a[i - 1][k - 1] * b[k + i - 2];
        b[i - 1] = dum / a[i - 1][0];
        if (l < mm) l++;
    }
}
//----------------------------------------------------------------------------
// Return the sign of a number

int River::signR(double a)
{
    int b;

    if (a >= 0) b = 1;
    else b = -1;

    return(b);
}

char* River::get_river_name(int i) {
    ChannelSegmentType* tmp_chan = &( parameter->ChannelSegmentsVecPtr->front()) + i;
    return tmp_chan->Q_Name;
}
int River::get_river_chsz(int i) {
    ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + i;
    return tmp_chan->chsz;
}
int River::get_river_count() {
    return parameter->ChannelSegmentsVecPtr->size();
}
int River::get_river_index_size()
{
    return river_pipe_index.size();
}
riverBoudaryIndex* River::get_river_index(string rivername)
{
    return river_pipe_index[rivername];
}
map<string, riverBoudaryIndex*> River::get_river_pipe_index()
{
    return river_pipe_index;
}
double River::get_river_width(riverBoudaryIndex river) {
    ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + river.i;
    return tmp_chan->ChanWidth[river.j];
}
double River::get_river_elevation(riverBoudaryIndex river) {
    ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + river.i;
    int pxi = tmp_chan->ChanX[river.j];
    int pyi = tmp_chan->ChanY[river.j];
    return parameter->Arrptr->DEM[pxi + pyi * (parameter->Parptr->xsz)];
}
double River::get_river_H(riverBoudaryIndex river) {
    ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + river.i;
    int pxi = tmp_chan->ChanX[river.j];
    int pyi = tmp_chan->ChanY[river.j];
    return parameter->Arrptr->H[pxi + pyi * (parameter->Parptr->xsz)];
}
riverBoudaryIndex* River::get_boundary(char* PSNAME) {
    riverBoudaryIndex* temp = NULL;
    for (int i = 0; i < parameter->ChannelSegmentsVecPtr->size(); i++)
    {
        ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + i;
        for (int j = 0; j < tmp_chan->chsz; j++) {
            if (strcmp((tmp_chan->Q_Name + j * 80), PSNAME) == 0) {
                temp = new riverBoudaryIndex();
                temp->i = i;
                temp->j = j;
            }
        }
    }
    return temp;
};
void  River::set_river_inflow(double depth, riverBoudaryIndex name) {
    if (fabs(depth - (-1.0)) < 0.01) {
        depth = -0.99;
    }
    ChannelSegmentType* tmp_chan = &(parameter->ChannelSegmentsVecPtr->front()) + name.i;
    int nbdy = tmp_chan->Q_Val[name.j];
    double* data = tmp_chan->QVarlist[nbdy];
    int k = 0;
    while (data[k] != -1)//change Q
    {
        if (k % 2 == 0) {
            data[k] = depth;// overflow  into LISFLOOD unit m^2/s
        }
        k++;
    }
}