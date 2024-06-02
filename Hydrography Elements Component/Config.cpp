#include "Config.h"
#include"pipe_dependency\bmi_swmm.h"
#include "pipe_dependency\headers.h"
Config::Config(const char* config, const char* sheet)
{
    //xls name
    sheet_name=sheet;
    config_name=config;
    init_flag = false;
    RiversIndexVecPtr = new vector<int>(); // CCS
    // lisflood structure
	Arrptr = new Arrays();
	Fnameptr = new Fnames();
	Statesptr = new States();
	Parptr = new Pars();
	Solverptr = new Solver();
	BCptr = new BoundCs();
	Stageptr = new Stage();
    ChannelSegmentsVecPtr = new vector<ChannelSegmentType>();

	verbose = new int;
	*verbose = OFF;
    Solverptr->DepthThresh = 1e-3;
    Parptr->FPn = 0.06;
    Parptr->dx = 10.0;
    Parptr->dy = 10.0;
    Parptr->dA = 100.0;
    Parptr->tlx = 0.0;
    Parptr->tly = 0.0;
    Parptr->blx = 0.0;
    Parptr->bly = 0.0;
    Parptr->FPn = 0.06;
    Parptr->SaveInt = 1000.0;
    Parptr->SaveTotal = 0.0;
    Parptr->MassInt = 100.0;
    Parptr->MassTotal = 0.0;
    Parptr->SaveNo = 0;
    Parptr->op = 100.0;
    Parptr->InfilLoss = 0.0;
    Parptr->EvapLoss = 0.0;
    Parptr->RainLoss = 0.0;
    Parptr->InfilTotalLoss = 0.0;
    Parptr->EvapTotalLoss = 0.0;
    Parptr->RainTotalLoss = 0.0;
    Parptr->ch_start_h = 2.0; // default start water depth for channel
    Parptr->steadyQtol = 0.0005; // tolerance for steady-state definition
    Parptr->steadyInt = 1800.0; // interval at which to assess steady-state
    Parptr->steadyTotal = 0.0;
    Parptr->min_dx = 10.0; // CCS Holds min_dx value (needed for variable dimension lat-long grids)
    Parptr->min_dy = 10.0; // CCS Holds min_dy value (needed for variable dimension lat-long grids)
    Parptr->min_dx_dy = 10.0; // CCS Holds min of min_dx and min_dy values (needed for variable dimension lat-long grids)
    
    // Define initial values for boundary conditions
    BCptr->Qin = 0.0;
    BCptr->Qout = 0.0;
    BCptr->VolInMT = 0.0;
    BCptr->VolOutMT = 0.0;
    // Define initial values for solver settings
    Solverptr->Sim_Time = 3600.0;
    Solverptr->InitTstep = 10.0;		// Maximum timestep
    Solverptr->Nit = 360;
    Solverptr->itCount = 0;
    Solverptr->t = 0.0;
    Solverptr->g = 9.8065500000000;
    Solverptr->divg = (1 / (2 * Solverptr->g));
    Solverptr->SolverAccuracy = 1e-4;
    Solverptr->dynsw = 0; // Switch for full dynamic steady state (1) or diffusive steady state (0)
    Solverptr->DepthThresh = 1e-3;
    Solverptr->MomentumThresh = 1e-2;
    Solverptr->MaxHflow = 10.0;
    Solverptr->Hds = 0.0;
    Solverptr->Qerror = 0.0;
    Solverptr->Verror = 0.0;
    Solverptr->dhlin = 0.01;
    Solverptr->htol = 1.0;
    Solverptr->Qlimfact = 1.0;
    Solverptr->itrn_time = 0.0;
    Solverptr->itrn_time_now = 0.0;
    Solverptr->ts_multiple = 1;  // default to x1 timestep decouple multiple
    Solverptr->theta = 1.0; // GAMA (for q-centred numerical scheme), 1.0= semi-implicit version (Bates et al 2010);
    Solverptr->fricSolver2D = ON; //GAMA: uses the 2D friction scheme as default
    
    // Define initial values for arrays
    Arrptr->Manningsn = NULL;

    // Define default values for SimStates instance of States
    Statesptr->diffusive = OFF;	// CCS added default state
    Statesptr->ChannelPresent = OFF;
    Statesptr->TribsPresent = ON;
    Statesptr->NCFS = ON;
    Statesptr->save_depth = ON;
    Statesptr->save_elev = ON;
    Statesptr->out_dir = OFF;
    Statesptr->single_op = OFF;
    Statesptr->multi_op = OFF;
    Statesptr->calc_area = OFF;
    Statesptr->calc_meandepth = OFF;
    Statesptr->calc_volume = OFF;
    Statesptr->save_stages = OFF;
    Statesptr->adaptive_ts = OFF;
    Statesptr->qlim = OFF; //TJF: Switch for qlim version, default is OFF
    Statesptr->acceleration = OFF; //PB: Switch for acceleration version, default is OFF
    Statesptr->debugmode = OFF;
    Statesptr->save_Qs = OFF;
    Statesptr->calc_infiltration = OFF;
    Statesptr->call_gzip = OFF;
    Statesptr->alt_ascheader = OFF;
    Statesptr->checkpoint = OFF;
    Statesptr->checkfile = OFF;
    Statesptr->calc_evap = OFF;
    Statesptr->routing = OFF; //CCS: Switch for rainfall routing routine
    Statesptr->rainfall = OFF;
    Statesptr->reset_timeinit = OFF;
    Statesptr->profileoutput = OFF;
    Statesptr->porosity = OFF;
    Statesptr->weirs = OFF;
    Statesptr->save_Ts = OFF;
    Statesptr->save_QLs = OFF;
    Statesptr->startq = OFF;
    Statesptr->logfile = OFF;
    Statesptr->startfile = OFF;
    Statesptr->start_ch_h = OFF;
    Statesptr->comp_out = OFF;
    Statesptr->chainagecalc = ON;
    Statesptr->mint_hk = OFF;
    Statesptr->Roe = OFF;
    Statesptr->killsim = OFF;
    Statesptr->dhoverw = OFF;
    Statesptr->drychecking = ON;
    Statesptr->voutput = OFF;
    Statesptr->steadycheck = OFF;
    Statesptr->hazard = OFF;
    Statesptr->startq2d = OFF;
    Statesptr->Roe_slow = OFF;
    Statesptr->multiplerivers = OFF;
    Statesptr->SGC = OFF;
    Statesptr->SGCbed = OFF;
    Statesptr->SGCcat_area = OFF;
    Statesptr->SGCchangroup = OFF;
    Statesptr->SGCchanprams = OFF;
    Statesptr->SGCbfh_mode = OFF;
    Statesptr->SGCA_mode = OFF;
    Statesptr->binary_out = OFF;
    Statesptr->gsection = OFF;
    Statesptr->binarystartfile = OFF;
    Statesptr->startelev = OFF;
    Statesptr->latlong = OFF;
    Statesptr->dist_routing = OFF;
    Statesptr->SGCvoutput = OFF; // switch for sub-grid channel velocity output
    readParamFile_lisflood(config_name, sheet_name);
}

void Config::init()
{
    if (!init_flag) {
        init_lisflood();
        init_flag = true;
    }
}

void Config::readParamFile_lisflood(const char* config_file, const char* sheet)
{
    if (*verbose == ON) cout<<"Loading parameters... "<<endl;
    excel=new Excel();
    excel->read_excel(config_file, sheet);
    if (excel->sheet) {//如果表存在
        size_t maxRows = excel->sheet->GetTotalRows();//获取工作表总行数
        size_t maxCols = excel->sheet->GetTotalCols();//总列数
        for (size_t r = 0; r < maxRows; ++r)
        {
            for (size_t c = 0; c < maxCols; ++c)
            {
                BasicExcelCell* cell = excel->sheet->Cell(r, c);
                if (cell->Type() == BasicExcelCell::STRING)//找到标识符（字符串类型）
                {
                    //读数据
                    string temp = cell->GetString();
                    if (temp == "initial_tstep") Solverptr->InitTstep = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "resroot") strcpy(Fnameptr->resrootname,excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "dirroot") {
                        strcpy(Fnameptr->dirrootname, excel->sheet->Cell(r, c + 1)->GetString());
                        Statesptr->out_dir = ON;
                    };
                    if (temp == "sim_time") Solverptr->Sim_Time = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "saveint") Parptr->SaveInt = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "massint") Parptr->MassInt = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "fpfric") Parptr->FPn = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "startq") Statesptr->startq = ON;
                    if (temp == "ch_start_h")
                    {
                        Parptr->ch_start_h = excel->sheet->Cell(r, c + 1)->GetDouble();
                        Statesptr->start_ch_h = ON;
                    }
                    if (temp == "DEMfile") strcpy(Fnameptr->demfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "DSMfile") {
                        Statesptr->buildingHeight = 0;
                        strcpy(Fnameptr->DSMname, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "ARFfile") {
                        Statesptr->ARFFlag = ON;
                        strcpy(Fnameptr->ARFname, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp=="Damfile") strcpy(Fnameptr->weirfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "BuildingStartWaterName") {
                        Statesptr->BuildingStartWaterFlag = ON;
                        strcpy(Fnameptr->BuildingStartWaterName, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "BuildingManningn") {
                        Statesptr->BuildingManningFlag = ON;
                        strcpy(Fnameptr->BuildingManningn, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "Buildingheight") {
                        Statesptr->buildingHeight = 1;
                        strcpy(Fnameptr->BuildingheightName, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "Building_Flag") {
                        Statesptr->buildingcouple = 1;
                        strcpy(Fnameptr->BuildingFlag, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "manningfile") strcpy(Fnameptr->nfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "ARF") Solverptr->ARF = excel->sheet->Cell(r, c + 1)->GetDouble();
                    if (temp == "qoutput") Statesptr->save_Qs = ON;
                    if (temp == "voutput") Statesptr->voutput = ON;
                    if (temp == "diffusive") Statesptr->diffusive = ON;
                    if (temp == "evaporation")
                    {
                        Statesptr->calc_evap = ON;
                        strcpy(Fnameptr->evapfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "rainfall") // Enable rainfall CCS
                    {
                        Statesptr->rainfall = ON;
                        strcpy(Fnameptr->rainfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "multiriverfile") {
                        strcpy(Fnameptr->multiriverfilename, excel->sheet->Cell(r, c + 1)->GetString());
                        Statesptr->multiplerivers = ON;
                        if (*verbose == ON) printf("\nMultiple river mode selected\n");
                    }
                    if (temp == "riverfile") strcpy(Fnameptr->rivername, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "river_pipe_index") { 
                        strcpy(Fnameptr->river_pipe_index_name, excel->sheet->Cell(r, c + 1)->GetString());
                        Statesptr->river_couple = 1;
                    }
                    if (temp == "land_pipe_index") {
                        strcpy(Fnameptr->land_pipe_index_name, excel->sheet->Cell(r, c + 1)->GetString());
                        Statesptr->land_couple = 1;
                    }
                    if (temp == "Land_river_index") {
                        strcpy(Fnameptr->land_river_index_name, excel->sheet->Cell(r, c + 1)->GetString());
                        Statesptr->river2d_couple = 1;
                    }
                    if (temp == "Lakefile") { 
                        Statesptr->startfile = 1;
                        strcpy(Fnameptr->Lakestartfilename, excel->sheet->Cell(r, c + 1)->GetString()); 
                    }
                    if (temp == "Roadfile") strcpy(Fnameptr->Roadname, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "bcifile") strcpy(Fnameptr->bcifilename, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "bdyfile") strcpy(Fnameptr->bdyfilename, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "riverDirRoot") {
                        const char* source =excel->sheet->Cell(r, c + 1)->GetString();
                        para[0] = new char[strlen(source) + 1];
                        strcpy(para[0], source);
                    }
                    if (temp == "riverResRoot") {
                        const char* source = excel->sheet->Cell(r, c + 1)->GetString();
                        para[1] = new char[strlen(source) + 1];
                        strcpy(para[1], source);
                    }
                    if (temp == "RiverName") {
                        const char* source = excel->sheet->Cell(r, c + 1)->GetString();
                        para[2] = new char[strlen(source) + 1];
                        strcpy(para[2], source);
                    }
                    if (temp == "RiverFlag") {
                        strcpy(Fnameptr->River2dFlagName, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                }
            }
        }
        excel->excel_close();
    }
}

int Config::fexist(const char* filename)
{
    struct stat buffer;
    if (stat(filename, &buffer)) return 0;
    return 1;
}

int Config::init_lisflood()
{
    int result = 0;
    int i;
    char strtmp[256];

    if (Solverptr->t == 0)
    {
        Solverptr->Tstep = Solverptr->InitTstep;
        Solverptr->MinTstep = Solverptr->InitTstep;
    }
    if (Statesptr->out_dir == ON)
    {
        char tmp_sys_com[255]; // temporary string to hold system command
        if (fexist(Fnameptr->dirrootname) == 0) // check if it doesn't exist
        {
            //create output folder
            sprintf(tmp_sys_com, "%s%s", "mkdir ", Fnameptr->dirrootname);
            system(tmp_sys_com);
        }
    }
    //update the resroot to include the folder information
    strcpy(strtmp, Fnameptr->resrootname);
#ifdef unix
    sprintf(Fnameptr->resrootname, "%s/%s", Fnameptr->dirrootname, strtmp);
#elif __APPLE__
    sprintf(Fnameptr->resrootname, "%s/%s", Fnameptr->dirrootname, strtmp);
#else
    sprintf(Fnameptr->resrootname, "%s\\%s", Fnameptr->dirrootname, strtmp);
#endif
    time_t ts = time(0);
    tm timeS = *localtime(&ts);
    if (*verbose) {
        printf("\nStart Date: %d/%d/%d \n", timeS.tm_mday, timeS.tm_mon + 1, timeS.tm_year + 1900);
        printf("Start Time: %d:%d:%d \n\n", timeS.tm_hour, timeS.tm_min, timeS.tm_sec);
    }
    sprintf(t1, "%s%s", Fnameptr->resrootname, ".mass");
    Fps.mass_fp = fopen(t1, "w");
    if (Fps.mass_fp != NULL)
    {
        if (Solverptr->t == 0) fprintf(Fps.mass_fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-Inf+Evap\n");
        else
        {
            // make a note in the mass file that this is a restart point - user can then edit the overlap out if they want a continuous mass file record.
            fprintf(Fps.mass_fp, "####################################################### Checkpoint restart ########################################################\n");
            fprintf(Fps.mass_fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-Inf+Evap\n");
            fflush(Fps.mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
        }
    }
    else
    {
        if (*verbose == ON)
        {
            printf("Unable to open mass balance file: %s", t1);
            result = 1;
            return result;
        }
    }

    //stage output file
    if (Statesptr->save_stages == ON) {
        sprintf(t1, "%s%s", Fnameptr->resrootname, ".stage");
        if (Statesptr->checkpoint == ON && Solverptr->t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
            Fps.stage_fp = fopen(t1, "a");
        }
        else {
            Fps.stage_fp = fopen(t1, "w");
        }
        if (Fps.stage_fp != NULL)
        {
            if (Solverptr->t == 0.0 || Statesptr->checkpoint == OFF)
            {
                fprintf(Fps.stage_fp, "Stage output, depth (m). Stage locations from: %s\n\n", Fnameptr->stagefilename);
                fprintf(Fps.stage_fp, "Stage information (stage,x,y,elev):\n");
                for (i = 0; i < Stageptr->Nstages; i++)
                {
                    if (Statesptr->SGC == ON && Arrptr->SGCwidth[Stageptr->stage_grid_x[i] + Stageptr->stage_grid_y[i] * Parptr->xsz] > 0) // if a SUB GRID channel is present export the channel bed elevation)
                    {
                        if (Stageptr->stage_check[i] == 1) fprintf(Fps.stage_fp, "%d\t%.4f\t%.4f\t%.4f\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i], Arrptr->SGCz[Stageptr->stage_grid_x[i] + Stageptr->stage_grid_y[i] * Parptr->xsz]);
                        else fprintf(Fps.stage_fp, "%d\t%.4f\t%.4f\tn/a\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i]);
                    }
                    else
                    {
                        if (Stageptr->stage_check[i] == 1) fprintf(Fps.stage_fp, "%d\t%.4f\t%.4f\t%.4f\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i], Arrptr->DEM[Stageptr->stage_grid_x[i] + Stageptr->stage_grid_y[i] * Parptr->xsz]);
                        else fprintf(Fps.stage_fp, "%d\t%.4f\t%.4f\tn/a\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i]);
                    }
                }
                fprintf(Fps.stage_fp, "\nOutput, depths:\n");
                fprintf(Fps.stage_fp, "Time; stages 1 to %d\n", Stageptr->Nstages);
            }
            else
            {
                fprintf(Fps.stage_fp, "####################################################### Checkpoint restart ########################################################\n");
                fflush(Fps.stage_fp);
            }
        }
        else
        {
            if (*verbose == ON) printf("Unable to open stage output file: %s", t1);
            Statesptr->save_stages = OFF;
        }

    }
    //velocity output file
    if (Statesptr->save_stages == ON && Statesptr->voutput == ON)
    {
        sprintf(t1, "%s%s", Fnameptr->resrootname, ".velocity");
        if (Statesptr->checkpoint == ON && Solverptr->t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
            Fps.vel_fp = fopen(t1, "a");
        }
        else {
            Fps.vel_fp = fopen(t1, "w");
        }
        if (Fps.vel_fp != NULL) {
            if (Solverptr->t == 0) {
                fprintf(Fps.vel_fp, "Velocity output, velocity (ms-1). Velocity locations from: %s\n\n", Fnameptr->stagefilename);
                fprintf(Fps.vel_fp, "Stage information (stage,x,y,elev):\n");
                for (i = 0; i < Stageptr->Nstages; i++) {
                    if (Stageptr->stage_check[i] == 1) fprintf(Fps.vel_fp, "%d\t%.4f\t%.4f\t%.4f\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i], Arrptr->DEM[Stageptr->stage_grid_x[i] + Stageptr->stage_grid_y[i] * Parptr->xsz]);
                    else fprintf(Fps.vel_fp, "%d\t%.4f\t%.4f\tn/a\n", i + 1, Stageptr->stage_loc_x[i], Stageptr->stage_loc_y[i]);
                }
                fprintf(Fps.vel_fp, "\nOutput, depths:\n");
                fprintf(Fps.vel_fp, "Time; velocities 1 to %d\n", Stageptr->Nstages);
            }
            else {
                fprintf(Fps.vel_fp, "####################################################### Checkpoint restart ########################################################\n");
                fflush(Fps.vel_fp);
            }
        }
        else {
            if (*verbose == ON) printf("Unable to open velocity output file: %s", t1);
            Statesptr->save_stages = OFF;
        }
    }

    //velocity output file
    if (Statesptr->gsection == ON)
    {
        sprintf(t1, "%s%s", Fnameptr->resrootname, ".discharge");
        if (Statesptr->checkpoint == ON && Solverptr->t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
            Fps.gau_fp = fopen(t1, "a");
        }
        else {
            Fps.gau_fp = fopen(t1, "w");
        }
        if (Fps.gau_fp != NULL) {
            if (Solverptr->t == 0) {
                fprintf(Fps.gau_fp, "Discharge output, discharge (m3s-1). Discharge locations from: %s\n\n", Fnameptr->gaugefilename);
                fprintf(Fps.gau_fp, "Time; discharge 1 to %d\n", Stageptr->Ngauges);
            }
            else {
                fprintf(Fps.gau_fp, "####################################################### Checkpoint restart ########################################################\n");
                fflush(Fps.gau_fp);
            }
        }
        else {
            if (*verbose == ON) printf("Unable to open discharge output file: %s", t1);
            Statesptr->gsection = OFF;
        }
    }
    //start simulation
    time(&Solverptr->time_start);
    LoadDEM(
        Fnameptr,
        Statesptr,
        Parptr,
        Arrptr,
        verbose);
    if (Statesptr->calc_evap == ON)
        LoadEvap(
            Fnameptr,
            Arrptr,
            verbose);
    if (Statesptr->rainfall == ON)
        LoadRain(
            Fnameptr,
            Arrptr,
            verbose);
    Previous_t = Solverptr->t - Solverptr->Tstep;

    steadyCount = 0;
    tstep_counter = -1;   // start at -1 so that in first run through we calculate river

    tstep_channel = 0; // channel timestep
}

double Config::DomainVol(States* Statesptr, Pars* Parptr, ChannelSegmentType* ChannelSegments, Arrays* Arrptr, vector<ChannelSegmentType>* ChannelSegmentsVecPtr)
{
    double vol = 0.0;
    int i, j, pi, pj, chseg, pH0, p0;
    double HeightOverBank, por0, dAPor;
    ChannelSegmentType* csp;

#pragma omp parallel for private( j,pH0,por0,dAPor,p0) reduction( + : vol )
    for (i = 0; i < Parptr->xsz; i++) for (j = 0; j < Parptr->ysz; j++)
    {
        p0 = i + j * Parptr->xsz;
        if (Statesptr->porosity == ON)
        {
            if (Parptr->Por_Ident == 1 || Parptr->Por_Ident == 3)
            {
                por0 = Arrptr->paerial[p0];
                dAPor = Parptr->dA * por0;
                if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * dAPor;
            }
            else if (Parptr->Por_Ident == 2 || Parptr->Por_Ident == 4)
            {
                pH0 = (int)(Arrptr->H[p0] / Parptr->zlev);
                if (pH0 > (Parptr->maxelev / Parptr->zlev)) pH0 = (int)(Parptr->maxelev / Parptr->zlev);
                por0 = Arrptr->paerial[i + j * Parptr->xsz + pH0 * Parptr->xsz * Parptr->ysz];
                dAPor = Parptr->dA * por0;
                if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * dAPor;
            }
        }
        else if (Statesptr->SGC == ON)
        {
            // calculate add on volume
            vol += Arrptr->SGCVol[p0];
        }
        else
        {
            if (Arrptr->ChanMask[p0] == -1) vol += Arrptr->H[p0] * Parptr->dA;
        }
    }

    // Calculate Channel Volume if present
    if (Statesptr->ChannelPresent == ON)
    {
        for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++)
        {
            csp = ChannelSegments + chseg;
            for (i = 0; i < csp->chsz; i++)
            {
                pi = csp->ChanX[i];
                pj = csp->ChanY[i];
                vol += Arrptr->H[pi + pj * Parptr->xsz] * csp->ChanWidth[i] * csp->Chandx[i];

                HeightOverBank = Arrptr->DEM[pi + pj * Parptr->xsz] + Arrptr->H[pi + pj * Parptr->xsz] - csp->BankZ[i];

                if (Statesptr->NCFS && HeightOverBank > 0 && csp->ChanWidth[i] < Parptr->dx)
                    vol += csp->Chandx[i] * (Parptr->dx - csp->ChanWidth[i]) * HeightOverBank;
            }
        }
    }

    return(vol);
}

void Config::finalize()
{
    time(&Solverptr->time_finish);

    // get system time and echo for user
    if (*verbose == ON) {
        time_t tf = time(0);
        tm timeF = *localtime(&tf);
        printf("\nFinish Date: %d/%d/%d \n", timeF.tm_mday, timeF.tm_mon + 1, timeF.tm_year + 1900);
        printf("Finish Time: %d:%d:%d \n\n", timeF.tm_hour, timeF.tm_min, timeF.tm_sec);
    }
    //iteration time
    Solverptr->itrn_time = Solverptr->itrn_time + difftime(Solverptr->time_finish, Solverptr->time_start);
    if (*verbose == ON) printf("\n  Total computation time: %.2lf mins\n\n", (Solverptr->itrn_time / 60.0));

    if (Statesptr->logfile == ON)
    {
        freopen("CON", "w", stdout);
        printf("\nLisflood run finished see log file for run details");
    }

    if (Statesptr->save_stages == ON) fclose(Fps.stage_fp);
    sprintf(t1, "%s%s", Fnameptr->resrootname, ".stage");
    if (Statesptr->call_gzip == ON) {
        sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", t1);
        system(tmp_sys_com);
    }

    //fclose(Fps.mass_fp);
    sprintf(t1, "%s%s", Fnameptr->resrootname, ".mass");
    if (Statesptr->call_gzip == ON) {
        sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", t1);
        system(tmp_sys_com);
    }
    delete Arrptr;
    delete Fnameptr;
    delete Statesptr;
    delete Parptr;
    delete Solverptr;
    delete BCptr;
    delete Stageptr;
    delete verbose;
}

void Config::update_time()
{
    int i, j;
    size_t ptr;
    Files* Fptr = &Fps;
    BCs(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
    if (Statesptr->drychecking == ON) DryCheck(Parptr, Solverptr, Arrptr);
    // Infiltration, evaporation and rainfall routines after time step update (TJF)11
    if (Statesptr->calc_infiltration == ON)
        FPInfiltration(Parptr,Solverptr, Arrptr);
    if (Statesptr->calc_evap == ON)
        Evaporation(Parptr, Solverptr, Arrptr);
    if (Statesptr->rainfall == ON &&Statesptr->routing == OFF)
        Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
    UpdateH(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
    BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
    // Update t with final Tstep calculated
    if (Solverptr->t > 0.0) Solverptr->MinTstep = Tool::getmin(Solverptr->MinTstep, Solverptr->Tstep);
    Solverptr->t += Solverptr->Tstep;
    Solverptr->itCount += 1;
    // Update maxH, maxHtm, totalHtm and initHtm at the mass interval if mint_hk is specifed in the .par file OR at every time step if not
    if (Solverptr->t >= Parptr->MassTotal || Statesptr->mint_hk == OFF)
    {
#pragma omp parallel for private (j, ptr)
        for (i = 0; i < Parptr->xsz; i++) for (j = 0; j < Parptr->ysz; j++)
        {
            ptr = i + j * Parptr->xsz;
            if ((Arrptr->initHtm[ptr] == (NULLVAL)) && (Arrptr->H[ptr] > 0.01)) Arrptr->initHtm[ptr] = Solverptr->t / 3600.0;
            if (Arrptr->H[ptr] > 0.01) Arrptr->totalHtm[ptr] += Solverptr->Tstep / 3600.0;
            // Update maximum water depths, and time of maximum (in hours)
            if (Arrptr->H[ptr] > Arrptr->maxH[ptr])
            {
                Arrptr->maxH[ptr] = Arrptr->H[ptr];
                Arrptr->maxHtm[ptr] = Solverptr->t / 3600.0;
            }
        }
    }
    // Calculate mass balance error
    if (Solverptr->t >= Parptr->MassTotal)
    {
        Solverptr->vol2 = DomainVol(Statesptr, Parptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr); // CCS

        // calc losses for this mass interval
        double loss = (Parptr->InfilTotalLoss - Parptr->InfilLoss) + (Parptr->EvapTotalLoss - Parptr->EvapLoss) - (Parptr->RainTotalLoss - Parptr->RainLoss);

        //Solverptr->Qerror=BCptr->Qin-BCptr->Qout-(Solverptr->vol2+loss-Solverptr->vol1)/Parptr->MassInt;
        // New version using VolInMT and VolOutMT
        // volume error
        Solverptr->Verror = BCptr->VolInMT - BCptr->VolOutMT - (Solverptr->vol2 + loss - Solverptr->vol1);
        // Q error
        Solverptr->Qerror = Solverptr->Verror / Parptr->MassInt;
        // reset to 0.0
        BCptr->VolInMT = 0.0;
        BCptr->VolOutMT = 0.0;

        // record cumulative loss for next time.
        Parptr->InfilLoss = Parptr->InfilTotalLoss;
        Parptr->EvapLoss = Parptr->EvapTotalLoss;
        Parptr->RainLoss = Parptr->RainTotalLoss;

        // Calculate flood area
        double FloodArea = 0.0;
        double dA = Parptr->dA;
#pragma omp parallel for private(j,ptr) reduction ( + : FloodArea)
        for (i = 0; i < Parptr->xsz; i++) for (j = 0; j < Parptr->ysz; j++)
        {
            ptr = i + j * Parptr->xsz;
            if (Statesptr->latlong == ON) dA = Arrptr->dA[ptr]; // if latlong is on change dA to local cell area
            if (Statesptr->porosity == ON)
            {
                if (Arrptr->H[ptr] > 0.01) FloodArea += dA * Arrptr->paerial[ptr]; // If porosity used, scale flooded area by porosity (TJF)
            }
            else if (Statesptr->SGC == ON)
            {
                if (Arrptr->H[ptr] - Arrptr->SGCbfH[ptr] > Solverptr->DepthThresh) FloodArea += dA; // If sub-grid used remove channel depth
            }
            else
            {
                if (Arrptr->H[ptr] > 0.01) FloodArea += dA; // standard ara calculation
            }
        }
        Solverptr->FArea = FloodArea;

        fprintf(Fptr->mass_fp, "%-12.3f %-10.4f %-10.4f %-10li %12.4e %12.4e  %-11.3f %-10.3f %-11.3f %12.4e %12.4e %12.4e\n", Solverptr->t, Solverptr->Tstep, Solverptr->MinTstep, Solverptr->itCount, Solverptr->FArea, Solverptr->vol2, BCptr->Qin, Solverptr->Hds, BCptr->Qout, Solverptr->Qerror, Solverptr->Verror, Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss));
        fflush(Fptr->mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.

        Solverptr->vol1 = Solverptr->vol2;
        Parptr->MassTotal += Parptr->MassInt;
    }
    // Regular output

    if (Solverptr->t >= Parptr->SaveTotal)
    {

        double Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;
        time(&Solverptr->time_check);
        Comp_time = difftime(Solverptr->time_check, Solverptr->time_start) / 60;
        if (Comp_time != 0 && Statesptr->comp_out == ON) // only of t is not zero (can't divide by zero)
        {
            Model_Comp_Ratio = ((Solverptr->t / 60) / Comp_time);
            Model_time_left = (Solverptr->Sim_Time - Solverptr->t) / 60;
            Est_Time_Fin = (Model_time_left / Model_Comp_Ratio);
            Est_Time_Tot = Comp_time + Est_Time_Fin;
            printf("T(mins): M: %.1lf, C: %.1lf, M/C: %.2lf, ETot: %.1lf, EFin: %.1lf\n", (Solverptr->t / 60.0), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
        }
        write_regular_output(Fnameptr, Solverptr, Statesptr, Parptr, Arrptr);
        // update interval counter
        Parptr->SaveTotal += Parptr->SaveInt;
        Parptr->SaveNo += 1;
    }
    if(*verbose==ON)cout << "Time:    " << Solverptr->t << endl;
}

void Config::write_regular_output(Fnames* Fnameptr, Solver* Solverptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr)
{
    // write out water depth
    if (Statesptr->save_depth == ON)
    {
        // output binary of ascii rasters
        if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, Parptr->SaveNo, ".wdb", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
        else write_ascfile(Fnameptr->resrootname, Parptr->SaveNo, ".wd", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
    }
    return;
}
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
// BOUNDARY CONDITIONS
// Calculate Qx and Qy at edges of the domain in response to boundary
// conditions
void Config::BCs(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr)
{
    int i, j, BCi, numBCs, p0, p1, sign, dir, edge, pTQ;
    double h0, h1, z0, z1, hflow, fn, dh, Sf, alpha = 1e-6, g, q0;
    double Qlim, * qptr, * qoldptr;



    //double shl,shr,aux1,ul,ur,vl,vr,ubarra,vbarra,cbarra,a1,a2,a3,a1m,a2m,a3m,a1l,a1r,a3l,a3r,epsilon,dhu,dhv; CCS threw unreferenced local variable warning so I guess no longer needed?
    //double alfa1,alfa2,alfa3,e11,e12,e13,e21,e22,e23,e31,e32,e33,f1lp,f2lp,f3lp,f1rp,f2rp,f3rp;	CCS threw unreferenced local variable warning so I guess no longer needed?
    //double g1lp,g2lp,g3lp,g1rp,g2rp,g3rp; CCS threw unreferenced local variable warning so I guess no longer needed?

    g = Solverptr->g;

    // BCs default to zero flux
    for (j = 0; j < Parptr->ysz; j++) {
        Arrptr->Qx[j * (Parptr->xsz + 1)] = 0.0;  //zero West boundary
        Arrptr->Qx[Parptr->xsz + j * (Parptr->xsz + 1)] = 0.0; //zero East boundary
    }
    for (i = 0; i < Parptr->xsz; i++) {
        Arrptr->Qy[i] = 0.0;  //zero North boundary
        Arrptr->Qy[i + Parptr->ysz * (Parptr->xsz + 1)] = 0.0; //zero South boundary
    }
    numBCs = 2 * Parptr->xsz + 2 * Parptr->ysz;
    for (BCi = 0; BCi < numBCs; BCi++)
    {
        if (BCptr->BC_Ident[BCi] != 0 || Statesptr->Roe == ON) // do nothing if closed boundary and not the Roe solver
        {
            // First for each edge number work out where it is on the boundary,
            // the associated edge pixels and whether it's facing in the x or y
            // direction
            if (BCi <= Parptr->xsz - 1)
            {	// N(j=0) edge
                p0 = BCi;
                p1 = BCi + Parptr->xsz;
                sign = -1;
                qptr = Arrptr->Qy + BCi;
                qoldptr = Arrptr->Qyold + BCi;
                q0 = *qoldptr; // Set old value of Q for acceleration version
                dir = 1;
                pTQ = BCi;
                edge = 1;
            }
            else if (BCi >= Parptr->xsz && BCi <= Parptr->xsz + Parptr->ysz - 1)
            {  // E edge
                p0 = Parptr->xsz - 1 + (BCi - Parptr->xsz) * Parptr->xsz;
                p1 = Parptr->xsz - 2 + (BCi - Parptr->xsz) * Parptr->xsz;
                sign = 1;
                qptr = Arrptr->Qx + Parptr->xsz + (BCi - Parptr->xsz) * (Parptr->xsz + 1);
                qoldptr = Arrptr->Qxold + Parptr->xsz + (BCi - Parptr->xsz) * (Parptr->xsz + 1);
                dir = 2;
                pTQ = Parptr->xsz + (BCi - Parptr->xsz) * (Parptr->xsz + 1);
                q0 = *qoldptr; // Set old value of Q for acceleration version
                edge = 2;
            }
            else if (BCi >= Parptr->xsz + Parptr->ysz && BCi <= 2 * Parptr->xsz + Parptr->ysz - 1)
            {  // S(j=ysz-1) edge
                p0 = 2 * Parptr->xsz + Parptr->ysz - 1 - BCi + (Parptr->ysz - 1) * Parptr->xsz;
                p1 = 2 * Parptr->xsz + Parptr->ysz - 1 - BCi + (Parptr->ysz - 2) * Parptr->xsz;
                sign = 1;
                qptr = Arrptr->Qy + 2 * Parptr->xsz + Parptr->ysz - 1 - BCi + (Parptr->ysz) * (Parptr->xsz + 1);
                qoldptr = Arrptr->Qyold + 2 * Parptr->xsz + Parptr->ysz - 1 - BCi + (Parptr->ysz) * (Parptr->xsz + 1);
                dir = 1;
                pTQ = 2 * Parptr->xsz + Parptr->ysz - 1 - BCi + (Parptr->ysz) * (Parptr->xsz + 1);
                q0 = *qoldptr; // Set old value of Q for acceleration version
                edge = 3;
            }
            else
            {   // W edge
                p0 = 0 + (numBCs - 1 - BCi) * Parptr->xsz;
                p1 = 1 + (numBCs - 1 - BCi) * Parptr->xsz;
                sign = -1;
                qptr = Arrptr->Qx + (numBCs - 1 - BCi) * (Parptr->xsz + 1);
                qoldptr = Arrptr->Qxold + (numBCs - 1 - BCi) * (Parptr->xsz + 1);
                dir = 2;
                pTQ = (numBCs - 1 - BCi) * (Parptr->xsz + 1);
                q0 = *qoldptr; // Set old value of Q for acceleration version
                edge = 4;
            }

            // Now calculate flows
            if (BCptr->BC_Ident[BCi] == 1 && Arrptr->H[p0] > Solverptr->DepthThresh)
            {	// FREE boundary
                hflow = Arrptr->H[p0];
                h0 = Arrptr->H[p0];
                h1 = Arrptr->H[p1];
                z0 = Arrptr->DEM[p0];
                z1 = Arrptr->DEM[p1];
                if (Arrptr->Manningsn != NULL) fn = Arrptr->Manningsn[p0]; else fn = Parptr->FPn;

                if (!(Arrptr->ChanMask[p0] != -1))
                {
                    if (Statesptr->acceleration == ON || Statesptr->Roe == ON) // Use semi-implicit accceleration formulation
                    {
                        // if BCptr->BC_Val[BCi] is -1 there is no floodplain slope specifed... use local slope from elevation model (likely to be unstable)
                        if (BCptr->BC_Val[BCi] < -0.999) // ie -1 (done like this as double)
                        {
                            dh = z0 + h0 - z1 - h1;
                            //dh=-h0;  
                            Sf = -dh / Parptr->dx;
                        }
                        else // use floodplain slope
                        {
                            Sf = BCptr->BC_Val[BCi];
                        }
                        // multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
                        // fabs on Sf and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
                        *qoldptr = sign * (fabs(q0) + fabs(g * Solverptr->Tstep * hflow * Sf)) / (1 + g * Solverptr->Tstep * hflow * fn * fn * fabs(q0) / (pow(hflow, (10. / 3.))));
                        *qptr = *qoldptr * Parptr->dx; // potntially unnessasary?? check for repetition in UpdateQs (jcn)
                    }
                    else if (Statesptr->Roe == ON)
                    {
                        if (edge == 2) *qptr = 1.75 * sqrt(g) * pow(Arrptr->H[i - 1 + (Parptr->ysz) * (Parptr->xsz)], 1.5); // East
                        else if (edge == 4)*qptr = 1.75 * sqrt(g) * pow(Arrptr->H[i + 1 + (Parptr->ysz) * (Parptr->xsz)], 1.5); // West
                        else if (edge == 3) *qptr = 1.75 * sqrt(g) * pow(Arrptr->H[i + (Parptr->ysz - 1) * (Parptr->xsz)], 1.5);	// South			
                        else if (edge == 1) *qptr = 1.75 * sqrt(g) * pow(Arrptr->H[i + (Parptr->ysz + 1) * (Parptr->xsz)], 1.5);	// North
                    }
                    else // origional lisflood-fp boundary
                    {
                        // if BCptr->BC_Val[BCi] is -1 there is no floodplain slope specifed... use local slope from elevation model (origional lisflood)
                        if (BCptr->BC_Val[BCi] < -0.999) // ie -1 (done like this as double)
                        {
                            // origional lisflood-fp boundary
                            dh = -z0 - h0 + z1 + h1;
                            Sf = sqrt(fabs(dh) / Parptr->dx);
                        }
                        else
                        {
                            Sf = BCptr->BC_Val[BCi]; // sqrt of user specified slope calulated in input.cpp for qlim and adaptive
                            dh = BCptr->BC_Val[BCi] * BCptr->BC_Val[BCi] * Parptr->dx; // backcalulate dh for sqrt(user slope)... OK so this is not very efficient could be stored somewhere  
                        }

                        if (fabs(dh) < Solverptr->dhlin)
                        {
                            Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (fabs(dh) / Parptr->dx);
                            alpha = (pow(hflow, (5. / 3.)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
                        }
                        else alpha = pow(hflow, (5. / 3.)) / (2. * fn * Sf);

                        if (dh < 0) Sf = -Sf;
                        *qptr = sign * pow(hflow, (5. / 3.)) * Sf * Parptr->dy / fn;
                    }
                    if (Statesptr->adaptive_ts == ON)
                    {
                        Solverptr->Tstep = Tool::getmin(Solverptr->Tstep, (0.25 * Parptr->dy * Parptr->dy / alpha));
                        // TJF: added to record Tstep
                        if (dir == 1) Arrptr->TRecy[pTQ] = Solverptr->Tstep;
                        if (dir == 2) Arrptr->TRecx[pTQ] = Solverptr->Tstep;
                    }
                    else if (Statesptr->qlim == ON)
                    {
                        Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
                        if (fabs(*qptr) > Qlim)
                        {
                            if (*qptr > 0) *qptr = Qlim;
                            if (*qptr < 0) *qptr = -Qlim;
                            // TJF: added to record Qlim
                            if (dir == 1) Arrptr->LimQy[pTQ] = *qptr;
                            if (dir == 2) Arrptr->LimQx[pTQ] = *qptr;
                        }
                    }
                }
                else *qptr = 0.0;
                if (dir == 1 && sign == -1 && *qptr > 0) *qptr = 0.0;
                if (dir == 1 && sign == 1 && *qptr < 0) *qptr = 0.0;
                if (dir == 2 && sign == -1 && *qptr > 0) *qptr = 0.0;
                if (dir == 2 && sign == 1 && *qptr < 0) *qptr = 0.0;
            }

            // HFIX boundary
            if (BCptr->BC_Ident[BCi] == 2 && Arrptr->H[p0] > Solverptr->DepthThresh)
            {
                hflow = Tool::getmax(Arrptr->H[p0], BCptr->BC_Val[BCi] - Arrptr->DEM[p0]);
                h0 = BCptr->BC_Val[BCi];
                h1 = Arrptr->H[p0];
                z1 = Arrptr->DEM[p0];

                if (Arrptr->Manningsn != NULL) fn = Arrptr->Manningsn[p0]; else fn = Parptr->FPn;

                if (!(Arrptr->ChanMask[p0] != -1))
                {
                    if (Statesptr->acceleration == ON) // Use semi-implicit accceleration formulation
                    {
                        dh = h0 - z1 - h1;
                        // change slops direction depending on the edge
                        if (edge == 1 || edge == 4) Sf = -dh / Parptr->dx;
                        else Sf = dh / Parptr->dx;
                        // implement momentum equation
                        *qoldptr = (q0 - g * Solverptr->Tstep * hflow * Sf) / (1 + g * Solverptr->Tstep * hflow * fn * fn * fabs(q0) / (pow(hflow, (10. / 3.))));
                        *qptr = *qoldptr * Parptr->dx;
                    }
                    else // origional lisflood boundary
                    {
                        dh = -h0 + z1 + h1;
                        Sf = sqrt(fabs(dh) / Parptr->dx);

                        if (fabs(dh) < Solverptr->dhlin)
                        {
                            Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (fabs(dh) / Parptr->dx);
                            alpha = (pow(hflow, (5. / 3.)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
                        }
                        else alpha = pow(hflow, (5. / 3.)) / (2. * fn * Sf);

                        if (dh < 0) Sf = -Sf;
                        *qptr = sign * pow(hflow, (5. / 3.)) * Sf * Parptr->dy / fn;
                    }
                    if (Statesptr->adaptive_ts == ON && hflow > 0.0)
                    {
                        Solverptr->Tstep = Tool::getmin(Solverptr->Tstep, (0.25 * Parptr->dy * Parptr->dy / alpha));
                        // TJF: added to record Tstep
                        if (dir == 1) Arrptr->TRecy[pTQ] = Solverptr->Tstep;
                        if (dir == 2) Arrptr->TRecx[pTQ] = Solverptr->Tstep;
                    }
                    else if (Statesptr->qlim == ON)
                    {
                        Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
                        if (fabs(*qptr) > Qlim)
                        {
                            //printf("Q limited %lf->%lf at t=%i\n",fabs(*qptr),Qlim,ts);
                            if (*qptr > 0) *qptr = Qlim;
                            if (*qptr < 0) *qptr = -Qlim;
                            // TJF: added to record Qlim
                            if (dir == 1) Arrptr->LimQy[pTQ] = *qptr;
                            if (dir == 2) Arrptr->LimQx[pTQ] = *qptr;
                        }
                    }
                }
            }

            // HVAR boundary
            if (BCptr->BC_Ident[BCi] == 3 && Arrptr->H[p0] > Solverptr->DepthThresh)
            {
                //h0=BCVarlist[(int)BC_Val[BCi]][ts];
                h0 = Tool::InterpBC(BCptr->BCVarlist[(int)BCptr->BC_Val[BCi]], Solverptr->t);
                hflow = Tool::getmax(Arrptr->H[p0], h0 - Arrptr->DEM[p0]);
                h1 = Arrptr->H[p0];
                z1 = Arrptr->DEM[p0];

                if (Arrptr->Manningsn != NULL) fn = Arrptr->Manningsn[p0]; else fn = Parptr->FPn;
                if (!(Arrptr->ChanMask[p0] != -1)) {
                    dh = -h0 + z1 + h1;
                    Sf = sqrt(fabs(dh) / Parptr->dx);
                    if (fabs(dh) < Solverptr->dhlin)
                    {
                        Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (fabs(dh) / Parptr->dx);
                        alpha = (pow(hflow, (5. / 3.)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
                    }
                    else alpha = pow(hflow, (5. / 3.)) / (2. * fn * Sf);

                    if (dh < 0) Sf = -Sf;
                    *qptr = sign * pow(hflow, (5. / 3.)) * Sf * Parptr->dy / fn;
                    if (Statesptr->qlim == ON)
                    {
                        Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
                        if (fabs(*qptr) > Qlim)
                        {
                            if (*qptr > 0) *qptr = Qlim;
                            if (*qptr < 0) *qptr = -Qlim;
                            // TJF: added to record Qlim
                            if (dir == 1) Arrptr->LimQy[pTQ] = *qptr;
                            if (dir == 2) Arrptr->LimQx[pTQ] = *qptr;
                        }
                    }
                }
            }
        }


        // QFIX boundary
        if (BCptr->BC_Ident[BCi] == 4)
        {
            *qptr = -BCptr->BC_Val[BCi] * sign * Parptr->dx;
        }
        // QVAR boundary
        if (BCptr->BC_Ident[BCi] == 5)
        {
            *qptr = -sign * Tool::InterpBC(BCptr->BCVarlist[(int)BCptr->BC_Val[BCi]], Solverptr->t) * Parptr->dx;
            hflow = Tool::getmax(Arrptr->H[p0], h0 - Arrptr->DEM[p0]);
        }
    }
    return;
}
//-----------------------------------------------------------------------------
// GENERAL BINARY RASTER FILE WRITE ROUTINE
// MT new multi purpose version - eliminates a lot of repetition
// and implements new filenaming
void Config::write_binrasterfile(char* root, int SaveNumber, char* extension, double* data, double* dem, int outflag, States* Statesptr, Pars* Parptr)
/*
Purpose: writes data that would be in an ascii raster to a binary file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *dem		-	pointer to dem
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
    int i, j, tempi;
    double temp, no_data = -9999;
    FILE* fp;
    char fnam[800];
    //char tmp_sys_com[255];

    // check if there is a savenumber to add and create filename
    if (SaveNumber >= 0 && SaveNumber <= 9999)	sprintf(fnam, "%s-%.4d%s", root, SaveNumber, extension);
    else if (SaveNumber > 9999)			sprintf(fnam, "%s-%d%s", root, SaveNumber, extension);
    else 									sprintf(fnam, "%s%s", root, extension);

    // open file
    fp = fopen(fnam, "wb");

    // check file opened ok, if not give warning and exit function
    if (fp == NULL)
    {
        printf("Problems writing to file %s\n", fnam);
        return;
    }

    //Output alternative header or default header
    //if(Statesptr->alt_ascheader==ON)
    //{
      //for(i=0;i<6;i++) 
      //{
       // fprintf(fp,"%s",Parptr->ascheader[i]);
       //   size_t fwrite ( const void * ptr, size_t size, size_t count, FILE * stream );
      //    fwrite (ptr, double , sizeof(buffer) , fp );
      //}
    //} 
    if (outflag == 0 || outflag == 3)
    {
        fwrite(&Parptr->xsz, sizeof(int), 1, fp);
        fwrite(&Parptr->ysz, sizeof(int), 1, fp);
        fwrite(&Parptr->blx, sizeof(double), 1, fp);
        fwrite(&Parptr->bly, sizeof(double), 1, fp);
        fwrite(&Parptr->dx, sizeof(double), 1, fp);
        fwrite(&no_data, sizeof(double), 1, fp);
    }

    // output data switched by type
    switch (outflag)
    {
    case 0:	// normal cell output
        fwrite(data, sizeof(double), Parptr->xsz * Parptr->ysz, fp);

        break;
    case 1:	// flux ouput - ie cell boundaries not cells
      // Edited to output the correct number of rows as Qx includes fluxes across boundaries (TJF)
      // Origin offset by dx*0.5 to output fluxes at boundaries when read into GIS (TJF)
        tempi = Parptr->xsz + 1;
        fwrite(&tempi, sizeof(int), 1, fp);

        fwrite(&Parptr->ysz, sizeof(int), 1, fp);
        temp = Parptr->blx - (Parptr->dx / 2.0);
        fwrite(&temp, sizeof(double), 1, fp);
        fwrite(&Parptr->bly, sizeof(double), 1, fp);
        fwrite(&Parptr->dx, sizeof(double), 1, fp);
        fwrite(&no_data, sizeof(double), 1, fp);

        fwrite(data, sizeof(double), (Parptr->xsz + 1) * (Parptr->ysz + 1), fp);

        break;
    case 2:	// flux ouput - ie cell boundaries not cells
      // Edited to output the correct number of rows as Qy includes fluxes across boundaries (TJF)
      // Origin offset by dy*0.5 to output fluxes at boundaries when read into GIS (TJF)
        fwrite(&Parptr->xsz, sizeof(int), 1, fp);
        tempi = Parptr->ysz + 1;
        fwrite(&tempi, sizeof(int), 1, fp);
        fwrite(&Parptr->blx, sizeof(double), 1, fp);
        temp = Parptr->bly - (Parptr->dx / 2.0);
        fwrite(&temp, sizeof(double), 1, fp);
        fwrite(&Parptr->dx, sizeof(double), 1, fp);
        fwrite(&no_data, sizeof(double), 1, fp);

        fwrite(data, sizeof(double), (Parptr->xsz + 1) * (Parptr->ysz + 1), fp);

        break;
    case 3:	// special output for water elevation (add DEM to data)
        for (j = 0; j < Parptr->ysz; j++)
        {
            for (i = 0; i < Parptr->xsz; i++)
            {

                if (data[i + j * Parptr->xsz] <= 0.01) // set to null for very shallow depth
                {
                    // PB: Fixed bug in generation of elev files so these are only created when H<=0.01 (not H<=0).  
                    // Without this you get significant lateral curvature caused by thin films of water ponding in drying elements.
                    temp = NULLVAL;
                }
                else
                {
                    temp = dem[i + j * Parptr->xsz] + data[i + j * Parptr->xsz];
                }
                fwrite(&temp, sizeof(double), 1, fp);
            }
        }
        break;
    }

    // close file
    fclose(fp);

    return;
}
//-----------------------------------------------------------------------------
// GENERAL ASCII RASTER FILE WRITE ROUTINE
// MT new multi purpose version - eliminates a lot of repetition
// and implements new filenaming
void Config::write_ascfile(char* root, int SaveNumber, char* extension, double* data, double* dem, int outflag, States* Statesptr, Pars* Parptr)
/*
Purpose: writes data to a standard ascii file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *dem		-	pointer to dem
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
    int i, j;
    double temp;
    FILE* fp;
    char fnam[800];
    char tmp_sys_com[255];

    // check if there is a savenumber to add and create filename
    if (SaveNumber >= 0 && SaveNumber <= 9999)	sprintf(fnam, "%s-%.4d%s", root, SaveNumber, extension);
    else if (SaveNumber > 9999)				sprintf(fnam, "%s-%d%s", root, SaveNumber, extension);
    else 										sprintf(fnam, "%s%s", root, extension);

    // open file
    fp = fopen(fnam, "wb");

    // check file opened ok, if not give warning and exit function
    if (fp == NULL)
    {
        printf("Problems writing to file %s\n", fnam);
        return;
    }

    //Output alternative header or default header
    if (Statesptr->alt_ascheader == ON)
    {
        for (i = 0; i < 6; i++)
        {
            fprintf(fp, "%s", Parptr->ascheader[i]);
        }
    }
    else if (outflag == 0 || outflag == 3)
    {
        fprintf(fp, "ncols         %i\n", Parptr->xsz);
        fprintf(fp, "nrows         %i\n", Parptr->ysz);
        fprintf(fp, "xllcorner     %lf\n", Parptr->blx);
        fprintf(fp, "yllcorner     %lf\n", Parptr->bly);
        fprintf(fp, "cellsize      %lf\n", Parptr->dx);
        fprintf(fp, "NODATA_value  %lf\n", NULLVAL);
    }

    // output data switched by type
    switch (outflag)
    {
    case 0:	// normal cell output
        for (j = 0; j < Parptr->ysz; j++)
        {
            for (i = 0; i < Parptr->xsz; i++)
            {
                double T = data[i + j * (Parptr->xsz )];
                fprintf(fp, "%.3f\t", data[i + j * Parptr->xsz]);
            }
            // output end of line
            fprintf(fp, "\n");
        }
        break;
    case 1:	// flux ouput - ie cell boundaries not cells
      // Edited to output the correct number of rows as Qx includes fluxes across boundaries (TJF)
      // Origin offset by dx*0.5 to output fluxes at boundaries when read into GIS (TJF)
        fprintf(fp, "ncols         %i\n", Parptr->xsz + 1);
        fprintf(fp, "nrows         %i\n", Parptr->ysz);
        fprintf(fp, "xllcorner     %lf\n", Parptr->blx - (Parptr->dx / 2.0));
        fprintf(fp, "yllcorner     %lf\n", Parptr->bly);
        fprintf(fp, "cellsize      %lf\n", Parptr->dx);
        fprintf(fp, "NODATA_value  %lf\n", NULLVAL);

        for (j = 0; j < Parptr->ysz; j++)
        {
            for (i = 0; i < Parptr->xsz + 1; i++)
            {
                fprintf(fp, "%.3f\t", data[i + j * (Parptr->xsz + 1)]);
            }
            // output end of line
            fprintf(fp, "\n");
        }
        break;
    case 2:	// flux ouput - ie cell boundaries not cells
      // Edited to output the correct number of rows as Qy includes fluxes across boundaries (TJF)
      // Origin offset by dy*0.5 to output fluxes at boundaries when read into GIS (TJF)
        fprintf(fp, "ncols         %i\n", Parptr->xsz);
        fprintf(fp, "nrows         %i\n", Parptr->ysz + 1);
        fprintf(fp, "xllcorner     %lf\n", Parptr->blx);
        fprintf(fp, "yllcorner     %lf\n", Parptr->bly - (Parptr->dx / 2.0));
        fprintf(fp, "cellsize      %lf\n", Parptr->dx);
        fprintf(fp, "NODATA_value  %lf\n", NULLVAL);

        for (j = 0; j < Parptr->ysz + 1; j++)
        {
            for (i = 0; i < Parptr->xsz; i++)
            {
                fprintf(fp, "%.3f\t", data[i + j * (Parptr->xsz + 1)]);
            }
            // output end of line
            fprintf(fp, "\n");
        }
        break;
    case 3:	// special output for water elevation (add DEM to data)
        for (j = 0; j < Parptr->ysz; j++)
        {
            for (i = 0; i < Parptr->xsz; i++)
            {

                if (data[i + j * Parptr->xsz] <= 0.01) // set to null for very shallow depth
                {
                    // PB: Fixed bug in generation of elev files so these are only created when H<=0.01 (not H<=0).  
                    // Without this you get significant lateral curvature caused by thin films of water ponding in drying elements.
                    temp = NULLVAL;
                }
                else
                {
                    temp = dem[i + j * Parptr->xsz] + data[i + j * Parptr->xsz];
                }
                fprintf(fp, "%11.3f", temp);
            }
            // output end of line
            fprintf(fp, "\n");
        }
        break;
    }

    // close file
    fclose(fp);

    // check if we need to zip the file up
    if (Statesptr->call_gzip == ON)
    {
        sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", fnam);
        system(tmp_sys_com);
    }

    return;
}
//-----------------------------------------------------------------------------
// FILE OUTPUT
void Config::fileoutput(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr)
{
    // output binary or ascii rasters
    if (Statesptr->binary_out == ON) // output binary of ascii rasters
    {
        // Write time of initial flood inundation
        write_binrasterfile(Fnameptr->resrootname, -1, ".inittmb", Arrptr->initHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write total inundation time
        write_binrasterfile(Fnameptr->resrootname, -1, ".totaltmb", Arrptr->totalHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write maximum depth
        write_binrasterfile(Fnameptr->resrootname, -1, ".maxb", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write time of maximum depth
        write_binrasterfile(Fnameptr->resrootname, -1, ".maxtmb", Arrptr->maxHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        if (Statesptr->voutput == ON)
        {
            write_binrasterfile(Fnameptr->resrootname, -1, ".maxVxb", Arrptr->maxVx, Arrptr->DEM, 1, Statesptr, Parptr);
            write_binrasterfile(Fnameptr->resrootname, -1, ".maxVyb", Arrptr->maxVy, Arrptr->DEM, 2, Statesptr, Parptr);
        }
        if (Statesptr->hazard == ON)
        {
            // Write maximum V Vd and Hazard
            write_binrasterfile(Fnameptr->resrootname, -1, ".maxVcb", Arrptr->maxVc, Arrptr->DEM, 0, Statesptr, Parptr);
            write_binrasterfile(Fnameptr->resrootname, -1, ".maxVcdb", Arrptr->maxVcH, Arrptr->DEM, 0, Statesptr, Parptr);
            write_binrasterfile(Fnameptr->resrootname, -1, ".maxHazb", Arrptr->maxHaz, Arrptr->DEM, 0, Statesptr, Parptr);
        }
    }
    else
    {
        // Write time of initial flood inundation
        write_ascfile(Fnameptr->resrootname, -1, ".inittm", Arrptr->initHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write total inundation time
        write_ascfile(Fnameptr->resrootname, -1, ".totaltm", Arrptr->totalHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write maximum depth
        write_ascfile(Fnameptr->resrootname, -1, ".max", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
        // Write time of maximum depth
        write_ascfile(Fnameptr->resrootname, -1, ".maxtm", Arrptr->maxHtm, Arrptr->DEM, 0, Statesptr, Parptr);
        if (Statesptr->voutput == ON)
        {
            write_ascfile(Fnameptr->resrootname, -1, ".maxVx", Arrptr->maxVx, Arrptr->DEM, 1, Statesptr, Parptr);
            write_ascfile(Fnameptr->resrootname, -1, ".maxVy", Arrptr->maxVy, Arrptr->DEM, 2, Statesptr, Parptr);
        }
        if (Statesptr->hazard == ON)
        {
            // Write maximum V Vd and Hazard
            write_ascfile(Fnameptr->resrootname, -1, ".maxVc", Arrptr->maxVc, Arrptr->DEM, 0, Statesptr, Parptr);
            write_ascfile(Fnameptr->resrootname, -1, ".maxVcd", Arrptr->maxVcH, Arrptr->DEM, 0, Statesptr, Parptr);
            write_ascfile(Fnameptr->resrootname, -1, ".maxHaz", Arrptr->maxHaz, Arrptr->DEM, 0, Statesptr, Parptr);
        }
    }
}

void Config::LoadDEM(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    char dum[800];
    double no_data_value = -9999;
    int i, j;

    fp = fopen(Fnameptr->demfilename, "rb");
    if (fp == NULL) { fprintf(stderr, "No DEM file found. Aborting.\n"); exit(0); }

    if (*verbose == ON) printf("\nLoading DEM:\t%s\n", Fnameptr->demfilename);

    fscanf(fp, "%s %i", dum, &Parptr->xsz);
    fscanf(fp, "%s %i", dum, &Parptr->ysz);
    fscanf(fp, "%s %lf", dum, &Parptr->blx);
    fscanf(fp, "%s %lf", dum, &Parptr->bly);
    fscanf(fp, "%s %lf", dum, &Parptr->dx);

    Parptr->dx_sqrt = sqrt((double)Parptr->dx); // sqrt now for later use in flooplain calcs - small speed increase
    Parptr->dy = Parptr->dx; Parptr->dA = Parptr->dx * Parptr->dy;
    Parptr->tlx = Parptr->blx; Parptr->tly = Parptr->bly + Parptr->ysz * Parptr->dy;
    fscanf(fp, "%s %lf", dum, &no_data_value);

    if (*verbose == ON)
    {
        printf("%ix%i\nBL corner\t(%lf,%lf)\nNODATA_value\t%lf\n",
            Parptr->xsz, Parptr->ysz, Parptr->blx, Parptr->bly, no_data_value);
    }

    // allocate memory for arrays, Note the () at the end ensures all elements are initialised to zero
    Arrptr->H = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->temp_H = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->maxH = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->totalHtm = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->Qx = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
    Arrptr->Qy = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
    if (Statesptr->voutput == ON)
    {
        Arrptr->Vx = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
        Arrptr->Vy = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
        Arrptr->maxVx = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
        Arrptr->maxVy = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
    }
    if (Statesptr->hazard == ON)
    {
        Arrptr->maxVc = new double[Parptr->xsz * Parptr->ysz]();
        Arrptr->maxVcH = new double[Parptr->xsz * Parptr->ysz]();
        Arrptr->maxHaz = new double[Parptr->xsz * Parptr->ysz]();
    }
    Arrptr->Qxold = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
    Arrptr->Qyold = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();

    // allocate memory for velocity arrays U and V
    // currently only used in acceleration version - initialised under all conditions as may want them
    // for other lisflood versions (TJF)
    Arrptr->U = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();
    Arrptr->V = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)]();

    // allocate memory for none zero arrays
    Arrptr->maxHtm = new double[Parptr->xsz * Parptr->ysz];
    Arrptr->initHtm = new double[Parptr->xsz * Parptr->ysz];
    //Arrptr->TRecx = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)];
    //Arrptr->TRecy = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)];
    Arrptr->LimQx = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)];
    Arrptr->LimQy = new double[(Parptr->xsz + 1) * (Parptr->ysz + 1)];

    Arrptr->ChanMask = new int[Parptr->xsz * Parptr->ysz];
    Arrptr->SegMask = new int[Parptr->xsz * Parptr->ysz];

    Arrptr->DEM = new double[Parptr->xsz * Parptr->ysz];

    // allocate memory for flow direction array for routing very shallow flows from rainfall component CCS 13/03/2012
    if (Statesptr->routing == ON)
    {
        Arrptr->FlowDir = new int[Parptr->xsz * Parptr->ysz]();
        for (i = 0; i < Parptr->xsz * Parptr->ysz; i++) Arrptr->FlowDir[i] = (int)NULLVAL;
    }

    // allocate memory for lat long arrays
    Arrptr->dx = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->dy = new double[Parptr->xsz * Parptr->ysz]();
    Arrptr->dA = new double[Parptr->xsz * Parptr->ysz]();



    // set initial values of elements for some arrays to NULLVAL
    for (i = 0; i < Parptr->xsz * Parptr->ysz; i++) Arrptr->maxHtm[i] = NULLVAL;
    for (i = 0; i < Parptr->xsz * Parptr->ysz; i++) Arrptr->initHtm[i] = NULLVAL;
    //for (i = 0; i < (Parptr->xsz + 1) * (Parptr->ysz + 1); i++) Arrptr->TRecx[i] = NULLVAL;
    //for (i = 0; i < (Parptr->xsz + 1) * (Parptr->ysz + 1); i++) Arrptr->TRecy[i] = NULLVAL;
    for (i = 0; i < (Parptr->xsz + 1) * (Parptr->ysz + 1); i++) Arrptr->LimQx[i] = NULLVAL;
    for (i = 0; i < (Parptr->xsz + 1) * (Parptr->ysz + 1); i++) Arrptr->LimQy[i] = NULLVAL;

    // set initial values of elements for mask arrays to NULLVAL
    for (i = 0; i < Parptr->xsz * Parptr->ysz; i++) Arrptr->ChanMask[i] = -1;
    for (i = 0; i < Parptr->xsz * Parptr->ysz; i++) Arrptr->SegMask[i] = -1;

    for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
    {
        fscanf(fp, "%lf", Arrptr->DEM + i + j * Parptr->xsz);
        if ((int)Arrptr->DEM[i + j * Parptr->xsz] == no_data_value) Arrptr->DEM[i + j * Parptr->xsz] = 1e10;
    }
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");

    return;
}
//-----------------------------------------------------------------------------
// LOAD TIME EVAPORATION FROM .evap FILE
void Config::LoadEvap(Fnames* Fnameptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j, ndata;
    char buff[255], units[80];

    fp = fopen(Fnameptr->evapfilename, "r");
    if (fp == NULL) return;
    if (*verbose == ON) printf("\nLoading time varying evaporation:\t%s\n", Fnameptr->evapfilename);

    //while(!feof(fp)) // removed this because it used to write over the evaportaion data if you had a carage return at the end of the file (JCN)
    // this will need to be change if you ever what to import multiple evaporation time series (can't think why you would do spatially varying this way though)
    //{
    j = 0; // skip 1st comment line
    do {
        buff[j] = fgetc(fp);
        if (feof(fp)) break;
    } while (buff[j++] != '\n');
    //if(feof(fp) || buff[0]=='\n') break; // removed because of removal of while loop see comment above

    fscanf(fp, "%i%s", &ndata, units);

    Arrptr->evap = new double[ndata * 2 + 2]();
    Arrptr->evap[ndata * 2 + 1] = -1;

    for (i = 0; i < ndata; i++) fscanf(fp, "%lf%lf", Arrptr->evap + i * 2, Arrptr->evap + i * 2 + 1);

    // convert time to seconds
    if (!strcmp(units, "hours")) for (i = 0; i < ndata; i++) *(Arrptr->evap + i * 2 + 1) *= 3600;
    else if (!strcmp(units, "days")) for (i = 0; i < ndata; i++) *(Arrptr->evap + i * 2 + 1) *= (3600 * 24);

    // convert evaporation rate from mm/day to m/second
    for (i = 0; i < ndata; i++) *(Arrptr->evap + i * 2) /= (1000 * 24 * 3600);
    //}
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");
    return;
}
void Config::LoadRain(Fnames* Fnameptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j, ndata;
    char buff[255], units[80];

    fp = fopen(Fnameptr->rainfilename, "r");
    if (fp == NULL) return;
    if (*verbose == ON) printf("\nLoading time varying rainfall:\t%s\n", Fnameptr->rainfilename);

    //while(!feof(fp))// removed this because it used to write over the evaportaion data if you had a carage return at the end of the file (JCN)
    // this will need to be change if you ever what to import multiple evaporation time series (can't think why you would do spatially varying this way though)
    //{
    j = 0; // skip 1st comment line
    do {
        buff[j] = fgetc(fp);
        if (feof(fp)) break;
    } while (buff[j++] != '\n');

    //if(feof(fp) || buff[0]=='\n') break;

    fscanf(fp, "%i%s", &ndata, units);

    Arrptr->rain = new double[ndata * 2 + 2];
    Arrptr->rain[ndata * 2 + 1] = -1;

    for (i = 0; i < ndata; i++) fscanf(fp, "%lf%lf", Arrptr->rain + i * 2, Arrptr->rain + i * 2 + 1);

    // convert time to seconds
    if (!strcmp(units, "hours")) for (i = 0; i < ndata; i++) *(Arrptr->rain + i * 2 + 1) *= 3600;
    else if (!strcmp(units, "days")) for (i = 0; i < ndata; i++) *(Arrptr->rain + i * 2 + 1) *= (3600 * 24);

    // convert rainfall rate from mm/hr to m/second
    for (i = 0; i < ndata; i++) *(Arrptr->rain + i * 2) /= (1000 * 3600);

    //}
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");
    return;
}
void Config::FPInfiltration(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr)
{
    int i, j;
    double cell_inf, h0;

    // Calculate Infiltration
    for (j = 0; j < Parptr->ysz; j++)
    {
        for (i = 0; i < Parptr->xsz; i++)
        {
            if (Arrptr->ChanMask[i + j * Parptr->xsz] == -1) { //only on non-channel cells
                h0 = Arrptr->H[i + j * Parptr->xsz];
                if (h0 > Solverptr->DepthThresh)
                {
                    cell_inf = Parptr->InfilRate * Solverptr->Tstep; //rate for depth, not area
                    h0 -= cell_inf;

                    //check for -ve depths
                    if (h0 < 0)
                    {
                        cell_inf += h0;
                        h0 = 0;
                    }
                    Arrptr->H[i + j * Parptr->xsz] = h0;
                    //for mass-balance
                    Parptr->InfilTotalLoss += cell_inf * Parptr->dA;
                }
            }
        }
    }

    return;
}
//-----------------------------------------------------------------------------
// FLOODPLAIN EVAPORATION, MDW
void Config::Evaporation(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr)
{
    int i, j;
    double cell_evap, h0, evap_rate;

    evap_rate = Tool::InterpBC(Arrptr->evap, Solverptr->t);//constant rate across whole floodplain
    cell_evap = evap_rate * Solverptr->Tstep; //rate for depth, not area

    // Calculate Evaporation
    for (j = 0; j < Parptr->ysz; j++)
    {
        for (i = 0; i < Parptr->xsz; i++)
        {
            if (Arrptr->ChanMask[i + j * Parptr->xsz] == -1)
            { //only on non-channel cells
                h0 = Arrptr->H[i + j * Parptr->xsz];
                if (h0 > Solverptr->DepthThresh)
                {
                    h0 -= cell_evap;

                    //check for -ve depths 
                    if (h0 < 0)
                    {
                        cell_evap += h0;
                        h0 = 0;
                    }
                    Arrptr->H[i + j * Parptr->xsz] = h0;
                    //for mass-balance
                    Parptr->EvapTotalLoss += cell_evap * Parptr->dA;
                }
            }
        }
    }

    return;
}
//-----------------------------------------------------------------------------
// PRECIPITATION INPUT, TJF
// Adds rainfall (only used when routing scheme is disabled; if routing is enabled, rainfall
// is added in Routing function)
void Config::Rainfall(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr)
{
    int i, j, p0;
    double cell_rain, h0, rain_rate, loc_rainfall_total = 0.0;

    rain_rate = Tool::InterpBC(Arrptr->rain, Solverptr->t);//constant rate across whole floodplain

    // Add rainfall depth to water heights
#pragma omp parallel for private(i, p0, h0, cell_rain) reduction (+:loc_rainfall_total)  // Parallelised CCS May 2013
    for (j = 0; j < Parptr->ysz; j++)
    {
        for (i = 0; i < Parptr->xsz; i++)
        {
            p0 = i + j * Parptr->xsz; // location of cell
            cell_rain = rain_rate * Solverptr->Tstep; //rate for depth, not area.  Has to be inside loop due to negative rain legacy support.

            if (Arrptr->DEM[p0] != 1e10)
            {
                h0 = Arrptr->H[p0];
                h0 += cell_rain;

                //check for -ve depths (legacy support for pre-evap code when rainfall could be negative)
                if (h0 < 0)
                {
                    cell_rain += h0;
                    h0 = 0;
                }
                Arrptr->H[p0] = h0;
                loc_rainfall_total += cell_rain * Parptr->dA; // mass balance for local cell (cumulative)		
            }
        }
    }
    Parptr->RainTotalLoss += loc_rainfall_total; // Update domain mass balance
    return;
}
void Config::set_verbose() {
    *verbose = ON;
    cout << "   开启LISFLOOD-FP中的verbose" << endl;
    cout << "   " << endl;
}

double Config::get_save_time()
{
    return Parptr->SaveInt;
}

double Config::get_current_time()
{
    return Solverptr->t;
}

double Config::get_sim_time()
{
    return Solverptr->Sim_Time;
}
Water_API int Config::get_PSNum() {
    return BCptr->numPS;
}
Water_API int* Config::get_BCPSNum()
{
    return &(BCptr->numPS);
}
char* Config::get_PSName()
{
    return BCptr->PS_Name;
}
Config_pipe::Config_pipe(char* config_file, char* sheet)
{
    inpfile = new char[80];
    outfile = new char[80];
    rptfile = new char[80];
    readParamFile_SWMM(config_file, sheet);
}

void Config_pipe::init()
{
}

void Config_pipe::update_time()
{
}

void Config_pipe::finalize()
{
}

void Config_pipe::readParamFile_SWMM(char* config_file, char* sheet)
{
    excel = new Excel();
    excel->read_excel(config_file, sheet);
    if (excel->sheet) {//如果表存在
        size_t maxRows = excel->sheet->GetTotalRows();//获取工作表总行数
        size_t maxCols = excel->sheet->GetTotalCols();//总列数
        for (size_t r = 0; r < maxRows; ++r)
        {
            for (size_t c = 0; c < maxCols; ++c)
            {
                BasicExcelCell* cell = excel->sheet->Cell(r, c);
                if (cell->Type() == BasicExcelCell::STRING)//找到标识符（字符串类型）
                {
                    //读数据
                    string temp = cell->GetString();
                    if(temp=="inpfile")strcpy(inpfile, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "outputfile")strcpy(outfile, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "rptfile")strcpy(rptfile, excel->sheet->Cell(r, c + 1)->GetString());
                    if (temp == "river_pipe_index") {
                        rivercouple = 1;
                        strcpy(river_pipe_index_name, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                    if (temp == "land_pipe_index") {
                        landcouple = 1;
                        strcpy(land_pipe_index_name, excel->sheet->Cell(r, c + 1)->GetString());
                    }
                }
            }
        }
    }
    excel->excel_close();
    this->ResultsFlag = &SaveResultsFlag;
    this->Routing = &DoRouting;
    this->Runoff = &DoRunoff;
    this->Exception_Count = &ExceptionCount;
    swmm_open(inpfile, rptfile, outfile);
    this->OpenFlag =new int(true);
}
//----------------------------------------------------------------------------
// SUM Qs INTO A CELL AND UPDATE DEPTH ACCORDINGLY
void Config::UpdateH(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr)
{
    int i, j, pi;
    double* qxptr0, * qyptr0, * qyptr1, * hptr,*dptr;
    int* mptr,*rptr0,*rptr1;
    double dV, himp, qtmp;
    double dAPorTemp;
    // 加入河道影响
    for (j = 0; j < Parptr->ysz; j++) {
        double* temp_hptr,* temp_hptr1;
        temp_hptr = Arrptr->H + j * Parptr->xsz;
        temp_hptr1 = Arrptr->temp_H + j * Parptr->xsz;
        for (i = 0; i < Parptr->xsz; i++) {
            if (*temp_hptr1 != 0) {
                *temp_hptr = *temp_hptr1;
                * temp_hptr1 = 0;
            }
            temp_hptr++;
            temp_hptr1++;
        }
    }
    // Insert point sources ((MT) moved before dV check, in case flow out of cell results in negative H reset and the inflow would prevent this)
  // H point sources moved back to after UpdateH (JCN)
    BCptr->Qpoint = 0.0;
    for (pi = 0; pi < BCptr->numPS; pi++)
    {
        if (BCptr->PS_Ident[pi] == 4) // QFIX
        {
            Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz] += BCptr->PS_Val[pi] * Parptr->dx * Solverptr->Tstep / Parptr->dA;
            BCptr->Qpoint += BCptr->PS_Val[pi] * Parptr->dx;
        }
        if (BCptr->PS_Ident[pi] == 5) // QVAR
        {
            qtmp = Tool::InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]], Solverptr->t);
            Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz] += qtmp * Parptr->dx * Solverptr->Tstep / Parptr->dA;
            BCptr->Qpoint += qtmp * Parptr->dx;
        }
    }

    // Calculate dV (+ve => outflow) and update NewH
#pragma omp parallel for private( i,qxptr0,qyptr0,qyptr1,rptr0,rptr1,hptr,mptr,dV,dAPorTemp)
    for (j = 0; j < Parptr->ysz; j++)
    {
        qxptr0 = Arrptr->Qx + j * (Parptr->xsz + 1);
        qyptr0 = Arrptr->Qy + j * (Parptr->xsz + 1);
        qyptr1 = Arrptr->Qy + (j + 1) * (Parptr->xsz + 1);
        hptr = Arrptr->H + j * Parptr->xsz;
        rptr0 = Arrptr->RiverFlag + j * Parptr->xsz;
        mptr = Arrptr->ChanMask + j * Parptr->xsz;
        for (i = 0; i < Parptr->xsz; i++)
        {
            if (*mptr == -1)
            {
                if (Statesptr->porosity == ON)
                {
                    dV = Solverptr->Tstep * (*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
                    dAPorTemp = PorArea(i, j, Parptr, Arrptr);
                    if (dAPorTemp == 0.0) (*hptr) += 0.0;
                    else (*hptr) += dV / dAPorTemp;

                    if (*hptr < 0) *hptr = 0;
                }
                else
                {
                    dV = Solverptr->Tstep * (*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
                    double temp = dV / Parptr->dA;
                    if (Statesptr->river2d_couple == 1)
                    {
                        if (*rptr0 == 1) {
                            printf("");
                        }
                        else
                        {
                            (*hptr) += temp;
                        }
                    }
                    else
                    {
                        (*hptr) += temp;
                    }
                    if (*hptr < 0) *hptr = 0;
                }
            }
            qxptr0++;
            qyptr0++;
            qyptr1++;
            rptr0++;
            hptr++;
            mptr++;
        }
    }

    // Point source HVAR and HFIX
    for (pi = 0; pi < BCptr->numPS; pi++)
    {
        if (BCptr->PS_Ident[pi] == 2) // HFIX
        {
            himp = BCptr->PS_Val[pi] - Arrptr->DEM[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz];
            if (himp < 0.0) himp = 0.0;
            BCptr->Qpoint += (himp - Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz]) * Parptr->dA / Solverptr->Tstep;
            Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz] = himp;
        }
        if (BCptr->PS_Ident[pi] == 3) // HVAR
        {
            himp = Tool::InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]], Solverptr->t) - Arrptr->DEM[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz];
            if (himp < 0.0) himp = 0.0;
            BCptr->Qpoint += (himp - Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz]) * Parptr->dA / Solverptr->Tstep;
            Arrptr->H[BCptr->xpi[pi] + BCptr->ypi[pi] * Parptr->xsz] = himp;
        }
    }

    return;
}
//-----------------------------------------------------------------------------------
// Calculates velocity and hazard
void Config::UpdateV(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr)
{
    int i, j, p0, pxy0, px1, pyl;
    double Vc, Xv, Yv, Haz;

#pragma omp parallel for private( i,p0,pxy0,px1,pyl,Xv,Yv,Vc,Haz)
    // calc at Ps
    for (j = 0; j < Parptr->ysz; j++)
    {
        for (i = 0; i < Parptr->xsz; i++)
        {
            p0 = i + j * Parptr->xsz;
            pxy0 = i + j * (Parptr->xsz + 1);
            px1 = i + 1 + j * (Parptr->xsz + 1);
            pyl = i + (j + 1) * (Parptr->xsz + 1);

            Xv = Tool::getmax(fabs(Arrptr->Vx[pxy0]), fabs(Arrptr->Vx[px1]));
            Yv = Tool::getmax(fabs(Arrptr->Vy[pxy0]), fabs(Arrptr->Vy[pyl]));

            Vc = sqrt(Xv * Xv + Yv * Yv);
            Haz = Arrptr->H[p0] * (Vc + 1.5); // Changed to equation from DEFRA 2006 (ALD)
            Arrptr->maxHaz[p0] = Tool::getmax(Arrptr->maxHaz[p0], Haz);
            if (Vc > Arrptr->maxVc[p0])
            {
                Arrptr->maxVc[p0] = Vc;
                Arrptr->maxVcH[p0] = Arrptr->H[p0];
            }
        }
    }
    return;
}

// Calculate net channel and floodplain flow in and out of domain boundary
void Config::BoundaryFlux(States* Statesptr, Pars* Parptr, Solver* Solverptr, BoundCs* BCptr, Arrays* Arrptr)
{
    int i;
    BCptr->Qin = 0.0; BCptr->Qout = 0.0;
    double qina = 0.0, qinb = 0.0, qinc = 0.0, qind = 0.0, qouta = 0.0, qoutb = 0.0, qoutc = 0.0, qoutd = 0.0;


    // ***** Note ..... Qy is positive north to south

#pragma omp parallel // private (Arrptr)
    {
#pragma omp sections private (i)
        {
            // North Boundary
            // loop through boundary and add positive Qy fluxes to Qin or subtract negative Qy fluxes from Qout
#pragma omp section 
            for (i = 0; i < Parptr->xsz; i++) if (Arrptr->ChanMask[i] == -1)
            {
                if (Arrptr->Qy[i] > 0) qina += Arrptr->Qy[i];
                else qouta -= Arrptr->Qy[i];
            }
            // South boundary
            // loop through boundary and add positive Qy fluxes to Qout or subtract negative Qy fluxes from Qin
#pragma omp section
            for (i = 0; i < Parptr->xsz; i++) if (Arrptr->ChanMask[i + (Parptr->ysz - 1) * Parptr->xsz] == -1)
            {
                if (Arrptr->Qy[i + Parptr->ysz * (Parptr->xsz + 1)] > 0) qoutb += Arrptr->Qy[i + Parptr->ysz * (Parptr->xsz + 1)];
                else qinb -= Arrptr->Qy[i + Parptr->ysz * (Parptr->xsz + 1)];
            }
            // ***** Note ..... Qx is positive west to east

            //	West boundary
            // loop through boundary and add positive Qx fluxes to Qin or subtract negative Qx fluxes from Qout
#pragma omp section
            for (i = 0; i < Parptr->ysz; i++) if (Arrptr->ChanMask[i * Parptr->xsz] == -1)
            {
                if (Arrptr->Qx[i * (Parptr->xsz + 1)] > 0) qinc += Arrptr->Qx[i * (Parptr->xsz + 1)];
                else qoutc -= Arrptr->Qx[i * (Parptr->xsz + 1)];
            }

            // East boundary
            // loop through boundary and add positive Qx fluxes to Qout or subtract negative Qx fluxes from Qin
#pragma omp section
            for (i = 0; i < Parptr->ysz; i++) if (Arrptr->ChanMask[Parptr->xsz - 1 + i * Parptr->xsz] == -1)
            {
                if (Arrptr->Qx[Parptr->xsz + i * (Parptr->xsz + 1)] > 0) qoutd += Arrptr->Qx[Parptr->xsz + i * (Parptr->xsz + 1)];
                else qind -= Arrptr->Qx[Parptr->xsz + i * (Parptr->xsz + 1)];
            }
        }
    }

    BCptr->Qout = qouta + qoutb + qoutc + qoutd;
    BCptr->Qin = qina + qinb + qinc + qind;

    //// Channel flows
    //if (Statesptr->ChannelPresent == ON)
    //{
    //    BCptr->Qout += BCptr->QChanOut;
    //    for (chseg = 0; chseg < (int)ChannelSegmentsVecPtr->size(); chseg++) // CCS
    //    {
    //        csp = ChannelSegments + chseg;
    //        for (i = 0; i < csp->chsz; i++)
    //        {
    //            if (csp->Q_Ident[i] == 4) BCptr->Qin += csp->Q_Val[i]; // add up fixed flows QFIXs
    //            else if (csp->Q_Ident[i] == 5) BCptr->Qin += InterpBC(csp->QVarlist[(int)csp->Q_Val[i]], Solverptr->t);  // add up QVARs
    //        }
    //    }
    //}


    // Point sources
    BCptr->Qin += BCptr->Qpoint;

    // calculate volume in and volume out
    BCptr->VolInMT += BCptr->Qin * Solverptr->Tstep;
    BCptr->VolOutMT += BCptr->Qout * Solverptr->Tstep;


    return;
}
//-----------------------------------------------------------------------------------
// CALCULATE VOLUME IN CELL WHEN POROSITY SCALES THE AREA AVAILABLE FOR STORAGE
double Config::PorArea(int i, int j, Pars* Parptr, Arrays* Arrptr)
{
    double h0, por0;
    int p0, pH0;
    double dAPor;
    int zsz, maxelev, zlev;

    zsz = Parptr->zsz;
    maxelev = (int)Parptr->maxelev;
    zlev = (int)Parptr->zlev;

    p0 = i + j * Parptr->xsz;

    h0 = Arrptr->H[p0];


    // Calculate area based on the porosity value
    if (Parptr->Por_Ident == 1 || Parptr->Por_Ident == 3)
    {
        por0 = Arrptr->paerial[p0];
        dAPor = Parptr->dA * por0;
    }
    else if (Parptr->Por_Ident == 2 || Parptr->Por_Ident == 4)
    {
        pH0 = (int)(h0 / zlev);
        if (pH0 > (maxelev / zlev)) pH0 = (maxelev / zlev);
        por0 = Arrptr->paerial[i + j * Parptr->xsz + pH0 * Parptr->xsz * Parptr->ysz];
        dAPor = Parptr->dA * por0;

    }


    return(dAPor);
}
//-----------------------------------------------------------------------------------
// CHECK FOR DRYING ELEMENTS
// If dV is going to make the water depth -ve, scale all flows out of cell
// so that the water volume in the cell goes to 0 in
void Config::DryCheck(Pars* Parptr, Solver* Solverptr, Arrays* Arrptr)
{
    int i, j;
    double WDweight, dV;
    double* hptr, * q1, * q2, * q3, * q4, cv;

    // Check for drying elements
  //#pragma omp parallel for private( i,hptr,q1,q2,q3,q4,cv,WDweight,dV)
    for (j = 0; j < Parptr->ysz; j++)
    {
        hptr = Arrptr->H + j * Parptr->xsz;
        for (i = 0; i < Parptr->xsz; i++, hptr++) if (*hptr > Solverptr->DepthThresh)
        {
            q1 = Arrptr->Qx + j * (Parptr->xsz + 1) + i;
            q2 = Arrptr->Qx + j * (Parptr->xsz + 1) + i + 1;
            q3 = Arrptr->Qy + j * (Parptr->xsz + 1) + i;
            q4 = Arrptr->Qy + (j + 1) * (Parptr->xsz + 1) + i;

            dV = Solverptr->Tstep * (*q1 - *q2 + *q3 - *q4);
            cv = *hptr * Parptr->dA;

            if (cv + dV < 0)
            {
                WDweight = -cv / dV; // 0.5 for improved drying stability
                if (*q1 < 0) *q1 *= WDweight;
                if (*q2 > 0) *q2 *= WDweight;
                if (*q3 < 0) *q3 *= WDweight;
                if (*q4 > 0) *q4 *= WDweight;
            }
        }
    }

    return;
}
double Config_pipe::get_current_time()
{
    double s_time = get_current_time_swmm(&s_time) * MSECperDAY / 1000;
    return s_time;
}
double Config_pipe::get_Elapsed_Time()
{
    double s_time = get_current_time_swmm(&s_time);
    return s_time;
}

double Config_pipe::get_Total_Duration()
{
    double t = get_TotalDuration();
    return t;
}

double Config_pipe::get_New_Routing_Time()
{
    double t = get_NewRoutingTime() ;
    return  t;
}
void Config_pipe::set_inflow(double* a)
{
    set_inflow_swmm(a);
}

void Config_pipe::get_overflow(double* b)
{
    get_overflow_swmm(b);
}

char* Config_pipe::get_pipe_name()
{
    return get_pipe_name_swmm(inpfile);
}

int Config_pipe::get_pipe_number()
{
    return get_pipe_number_swmm();
}

double Config_pipe::get_pipe_elevation(int j)
{
    return get_pipe_elevation_swmm(j);
}

void Config_pipe::setOutletDepth(int j, double z)
{
    setOutletDepth_swmm(j, z);
}

bool Config_pipe::isOutlet(int j)
{
    return isOutlet_swmm(j);
}

double Config_pipe::get_Link_Flow(int j)
{
    return get_Link_Flow_swmm(j);
}

int Config_pipe::get_upstream_Link(char* ps_name)
{
    return get_upstream_Link_swmm(ps_name);
}

double Config_pipe::get_position_node_x(int j)
{
    return get_position_node_x_swmm(j);
}

double Config_pipe::get_position_node_y(int j)
{
    return get_position_node_y_swmm(j);
}

bool Config_pipe::isRiverCouple()
{
    if (rivercouple == 1) return true;
    return false;
}

bool Config_pipe::isLandCouple()
{
    return landcouple;
}

char* Config_pipe::get_river_pipe_index_name()
{
    return river_pipe_index_name;
}

char* Config_pipe::get_land_pipe_index_name()
{
    return land_pipe_index_name;
}

BasicExcelWorksheet* Excel::read_excel(const char* excel_name, const char* sheet_name)
{
    BasicExcel* e = new BasicExcel();
    xls = e;
    // Load a workbook with one sheet, display its contents and
    // save into another file.
    if (!e->Load(excel_name))
    {
        printf("文件载入错误");
        return NULL;
    }
    sheet = e->GetWorksheet(sheet_name);
    return sheet;
}

void Excel::excel_close()
{
    xls->~BasicExcel();
}
