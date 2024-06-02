#include "Land.h"
Region::Region(const char* config, const char* sheet)
{
    parameter = new Config(config, sheet);

    if (parameter)
    {
        cout << "Creating surfure module..." << endl;
        cout << "   " << endl;
        // Qlim solver
        parameter->Statesptr->adaptive_ts = OFF;
        parameter->Statesptr->qlim = ON;
        cout << "Using Qlim formulation for floodplain flow\n" << endl;
        cout << "   " << endl;
    }
}
Region::Region(Config* lisf)
{
    parameter = lisf;
    if (parameter)
    {
        cout << "   Creating surfure module..." << endl;
        cout << "   " << endl;
        // Qlim solver
        parameter->Statesptr->adaptive_ts = OFF;
        parameter->Statesptr->qlim = ON;
        cout << "   Using Qlim formulation for floodplain flow" << endl;
        cout << "   " << endl;
    }
    parameter->init();
}
void Region::init()
{

    if (parameter->Statesptr->startfile == ON) {
        LoadStart(
            parameter->Fnameptr,
            parameter->Statesptr,
            parameter->Parptr,
            parameter->Arrptr,
            parameter->verbose);
    }  
    //LoadBCs(
    //    parameter->Fnameptr,
    //    parameter->Statesptr,
    //    parameter->Parptr, 
    //    parameter->BCptr, 
    //    parameter->Arrptr, 
    //    parameter->verbose);
    //LoadBCVar(
    //    parameter->Fnameptr,
    //    parameter->Statesptr, 
    //    parameter->Parptr, 
    //    parameter->BCptr, 
    //    parameter->Arrptr, 
    //    parameter->verbose);

    LoadManningsn(
        parameter->Fnameptr, 
        parameter->Parptr, 
        parameter->Arrptr,
        parameter->verbose);

    if (parameter->Statesptr->land_couple==ON)
    {
        Load_land_pipe_index(parameter->Fnameptr,
            parameter->Statesptr);
    }
}
void Region::update()
{
    if (parameter->Solverptr->t > 0.0) parameter->Solverptr->Tstep = parameter->Solverptr->InitTstep;
    if (parameter->Statesptr->BuildingFlag==0)FloodplainQ(parameter->Statesptr, parameter->Parptr, parameter->Solverptr, parameter->Arrptr);//这个地方不太好，应该加个建筑物flag栅格，然后内部判断是建筑物还是地块,类似于2维河道和地块的判断，


}
void Region::finalize()
{
}

void Region::LoadStart(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j, gr;
    char dum[800];
    double no_data_value = -9999;

    fp = fopen(Fnameptr->startfilename, "r");
    if (fp == NULL)
    {
        if (*verbose == ON) printf("\nWARNING: Unable to load startfile:\t%s\t\n", Fnameptr->startfilename);
        return;
    }

    if (*verbose == ON) printf("Loading initial depths:\t%s\t", Fnameptr->startfilename);


    for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
    fscanf(fp, "%s %lf", dum, &no_data_value);

    for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
    {
        fscanf(fp, "%lf", Arrptr->H + i + j * Parptr->xsz);
        // if no_data set depth to zero
        if ((int)Arrptr->H[i + j * Parptr->xsz] == no_data_value) Arrptr->H[i + j * Parptr->xsz] = 0.0;
        else if (Statesptr->startelev == ON) // convert water surface elevation to depth is this is being used
        {
            Arrptr->H[i + j * Parptr->xsz] = Tool::getmax(Arrptr->H[i + j * Parptr->xsz] - Arrptr->DEM[i + j * Parptr->xsz], 0.0);
        }
    }
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");

    return;
}
void Region::Load_land_pipe_index(Fnames* Fnameptr, States* Statesptr)
{
    ifstream file(Fnameptr->land_pipe_index_name);
    if (!file.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return;
    }
    int num = parameter->get_PSNum();
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        // 使用空格分隔键和值
        if (iss >> key >> value) {
            // 在这里你可以对key和value进行处理，比如存储到map中
            for (int i = 0; i < num; i++)
            {
                if (strcmp((parameter->get_PSName() + i * 80), key.c_str()) == 0)
                {
                    land_pipe_index[key] = i;
                }
            }
        }
    }
    file.close();
}
void Region::LoadManningsn(Fnames* Fnameptr, Pars* Parptr, Arrays* Arrptr, int* verbose)
{
    FILE* fp;
    int i, j;
    char dum[800];
    double no_data_value = -9999;

    fp = fopen(Fnameptr->nfilename, "r");
    if (fp == NULL) return;

    if (*verbose == ON) printf("Loading floodplain Manning's n data:\t%s\t", Fnameptr->nfilename);

    for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
    fscanf(fp, "%s %lf", dum, &no_data_value);

    Arrptr->Manningsn = new double[Parptr->xsz * Parptr->ysz];
    for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
    {
        fscanf(fp, "%lf", Arrptr->Manningsn + i + j * Parptr->xsz);
        if ((int)Arrptr->Manningsn[i + j * Parptr->xsz] == no_data_value) Arrptr->Manningsn[i + j * Parptr->xsz] = Parptr->FPn;
    }
    fclose(fp);

    if (*verbose == ON) printf("Done.\n\n");
    return;
}

//-----------------------------------------------------------------------------
// LOAD TIME VARYING BOUNDARY CONDITIONS FROM .bdy FILE

//-----------------------------------------------------------------------------
// CALCULATES FLOODPLAIN, OVERBANK AND WEIR FLOWS
void Region::FloodplainQ(States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr)
{
  int i, j;
  double h0, h1, ThreadTS, TmpTstep;
  double* hptr0, * hptr1, * qptr, * TSptr;
  int* rptr0, *rptr1;
  TmpTstep = Solverptr->Tstep;

  // Calculate Qx
#pragma omp parallel for private( i, h0, h1, hptr0,qptr,rptr0,rptr1,TSptr,ThreadTS) 
  for (j = 0; j < Parptr->ysz; j++)
  {
    hptr0 = Arrptr->H + j * Parptr->xsz;
    qptr = Arrptr->Qx + j * (Parptr->xsz + 1) + 1;
    rptr0 = Arrptr->RiverFlag + j * Parptr->xsz;
    rptr1 = Arrptr->RiverFlag + j * Parptr->xsz+1;
    // initialise thread time step for openMP
    ThreadTS = TmpTstep;
    TSptr = &ThreadTS;
    for (i = 0; i < Parptr->xsz - 1; i++)
    {
      h0 = *hptr0;
      h1 = *(hptr0 + 1);
      *qptr = 0.0;
      if (h0 > Solverptr->DepthThresh || h1 > Solverptr->DepthThresh) {
          if(Statesptr->river2d_couple==1){
              *qptr = CalcFPQx(i, j, Statesptr, Parptr, Solverptr, Arrptr, TSptr);
              if((*rptr0 == 1 && *rptr1 == 1)) {//|| ((*rptr0 == 0 && *rptr1 == 1)&& *qptr>0)|| ((*rptr0 == 1 && *rptr1 == 0) && *qptr < 0) 建筑物也就加个bptr类似这样写，注释里的这段话意思是如果（（当前点是地块&&右边是河流）&&自地块流向河流||（当前点是河流&&右边是地块）&&自地块流向河流）），哪位好心人想写下建筑物与地块耦合逻辑也可以在这里写，下面y轴类似的
                  *qptr = 0;
              }
              else
              {
                  
              }
          }
          else
          {
              *qptr = CalcFPQx(i, j, Statesptr, Parptr, Solverptr, Arrptr, TSptr);
          }
      }
      qptr++;
      hptr0++;
      rptr0++;
      rptr1++;
    }
#pragma omp critical
      {
        Solverptr->Tstep = Tool::getmin(Solverptr->Tstep, ThreadTS);
      }
    
  }
  // Calculate Qy
//#pragma omp section
#pragma omp parallel for private( i, h0, h1, hptr0,hptr1,qptr,rptr0,rptr1,TSptr,ThreadTS)
  for (j = 0; j < Parptr->ysz - 1; j++)
  {
    hptr0 = Arrptr->H + j * Parptr->xsz;
    hptr1 = Arrptr->H + (j + 1) * Parptr->xsz;
    qptr = Arrptr->Qy + (j + 1) * (Parptr->xsz + 1);
    rptr0 = Arrptr->RiverFlag + j * Parptr->xsz;
    rptr1 = Arrptr->RiverFlag + (j + 1) * Parptr->xsz;
    // initialise thread time step for openMP
    ThreadTS = TmpTstep;
    TSptr = &ThreadTS;
    for (i = 0; i < Parptr->xsz; i++)
    {
      h0 = *hptr0;
      h1 = *hptr1;
      *qptr = 0.0;
      if (h0 > Solverptr->DepthThresh || h1 > Solverptr->DepthThresh)
      {
          if (Statesptr->river2d_couple == 1) {
              *qptr = CalcFPQy(i, j, Statesptr, Parptr, Solverptr, Arrptr, TSptr);
              if ((*rptr0 == 1 && *rptr1 == 1)) {// || ((*rptr0 == 0 && *rptr1 == 1) && *qptr > 0) || ((*rptr0 == 1 && *rptr1 == 0) && *qptr < 0)
                  *qptr = 0;
              }
              else
              {
                  
              }
          }
          else
          {
              *qptr = CalcFPQy(i, j, Statesptr, Parptr, Solverptr, Arrptr, TSptr);
          }
      }
      hptr0++;
      hptr1++;
      qptr++;
      rptr0++;
      rptr1++;
    }

#pragma omp critical
      {
        Solverptr->Tstep = Tool::getmin(Solverptr->Tstep, ThreadTS);
      }
  }
  return;
}
double Region::CalcFPQx(int i, int j, States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr, double* TSptr)
{
    double z0, z1, h0, h1, Sf, hflow, fn, dh, Qlim, Q, alpha;
    int p0, p1, pTQ;

    p0 = i + j * Parptr->xsz;
    p1 = i + 1 + j * Parptr->xsz;
    pTQ = i + j * (Parptr->xsz + 1) + 1;

    z0 = Arrptr->DEM[p0];
    z1 = Arrptr->DEM[p1];
    h0 = Arrptr->H[p0];
    h1 = Arrptr->H[p1];

    if (Arrptr->Manningsn != NULL) fn = 0.5 * (Arrptr->Manningsn[p0] + Arrptr->Manningsn[p1]);
    else fn = Parptr->FPn;

    if (z0 + h0 > z1 + h1 && h0 > Solverptr->DepthThresh) // Flow from 0->1
    {
        dh = z0 + h0 - z1 - h1;
        Sf = sqrt(dh / Parptr->dx);
        hflow = Tool::getmax(z0 + h0, z1 + h1) - Tool::getmax(z0, z1);
        hflow = Tool::getmax(hflow, 0);
        hflow = Tool::getmin(hflow, Solverptr->MaxHflow);
        // added to record Hflow
        //Arrptr->Hflowx[pTQ] = hflow;

        if (Tool:: MaskTest(Arrptr->ChanMask[p0], Arrptr->ChanMask[p1]) && hflow > Solverptr->DepthThresh) {
            if (dh < Solverptr->dhlin)
            {
                Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (dh / Parptr->dx);
                alpha = (pow(hflow, (5.0 / 3.0)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
            }
            else alpha = pow(hflow, (5.0 / 3.0)) / (2. * fn * Sf);

            Q = (pow(hflow, (5.0 / 3.0)) * Sf * Parptr->dy / fn);
            if (Statesptr->adaptive_ts == ON)
            {
                *TSptr = Tool::getmin(*TSptr, (0.25 * Parptr->dy * Parptr->dy / alpha));
                // MT: added to record Tstep
                Arrptr->TRecx[pTQ] = *TSptr;
            }
            else
            {
                // flow limiter
                Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
                if (fabs(Q) > Qlim)
                {
                    if (Q > 0) Q = Qlim;
                    if (Q < 0) Q = -Qlim;
                    // MT added to record Qlim
                    Arrptr->LimQx[pTQ] = Q;
                }
            }
        }
        else Q = 0.0;
    }
    else if (z0 + h0<z1 + h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
    {
        dh = z1 + h1 - z0 - h0;
        Sf = sqrt(dh / Parptr->dx);
        hflow = Tool::getmax(z0 + h0, z1 + h1) - Tool::getmax(z0, z1);
        hflow = Tool::getmax(hflow, 0);
        hflow = Tool::getmin(hflow, Solverptr->MaxHflow);
        // added to record Hflow
        //Arrptr->Hflowx[pTQ] = hflow;

        if (Tool::MaskTest(Arrptr->ChanMask[p0], Arrptr->ChanMask[p1]) && hflow > Solverptr->DepthThresh)
        {
            if (dh < Solverptr->dhlin)
            {
                Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (dh / Parptr->dx);
                alpha = (pow(hflow, (5.0 / 3.0)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
            }
            else alpha = pow(hflow, (5.0 / 3.0)) / (2. * fn * Sf);

            Q = (-pow(hflow, (5.0 / 3.0)) * Sf * Parptr->dy / fn);

            // flow limiter
            Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
            if (fabs(Q) > Qlim)
            {
                if (Q > 0) Q = Qlim;
                if (Q < 0) Q = -Qlim;
                // MT added to record Qlim
                Arrptr->LimQx[pTQ] = Q;
            }
        }
        else Q = 0.0;

    }
    else Q = 0.0;
    // option to save V's
    if (Statesptr->voutput == ON)
    {
        if (Q != 0)
        {
            Arrptr->Vx[pTQ] = Q / Parptr->dx / hflow;
            Arrptr->maxVx[pTQ] = Tool::getmax(Arrptr->maxVx[pTQ], fabs(Arrptr->Vx[pTQ]));
        }
        else Arrptr->Vx[pTQ] = 0.0;
    }
    return(Q);
}
double Region::CalcFPQy(int i, int j, States* Statesptr, Pars* Parptr, Solver* Solverptr, Arrays* Arrptr, double* TSptr)
{
    double z0, z1, h0, h1, Sf, hflow, fn, dh, Qlim, Q, alpha;
    int p0, p1, pTQ;

    p0 = i + j * Parptr->xsz;
    p1 = i + (j + 1) * Parptr->xsz;
    pTQ = i + (j + 1) * (Parptr->xsz + 1);

    z0 = Arrptr->DEM[p0];
    z1 = Arrptr->DEM[p1];
    h0 = Arrptr->H[p0];
    h1 = Arrptr->H[p1];

    if (Arrptr->Manningsn != NULL) fn = 0.5 * (Arrptr->Manningsn[p0] + Arrptr->Manningsn[p1]);
    else fn = Parptr->FPn;


    if (z0 + h0 > z1 + h1 && h0 > Solverptr->DepthThresh)
    {
        dh = z0 + h0 - z1 - h1;
        Sf = sqrt(dh / Parptr->dx);
        hflow = Tool::getmax(z0 + h0, z1 + h1) - Tool::getmax(z0, z1);
        hflow = Tool::getmax(hflow, 0);
        hflow = Tool::getmin(hflow, Solverptr->MaxHflow);
        // added to record Hflow
        //Arrptr->Hflowy[pTQ] = hflow;

        if (Tool::MaskTest(Arrptr->ChanMask[p0], Arrptr->ChanMask[p1]) && hflow > Solverptr->DepthThresh)
        {
            if (dh < Solverptr->dhlin)
            {
                Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (dh / Parptr->dx);
                alpha = (pow(hflow, (5.0 / 3.0)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
            }
            else alpha = pow(hflow, (5.0 / 3.0)) / (2. * fn * Sf);

            Q = (pow(hflow, (5.0 / 3.0)) * Sf * Parptr->dy / fn);

            // flow limiter
            Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
            if (fabs(Q) > Qlim)
            {
                if (Q > 0) Q = Qlim;
                if (Q < 0) Q = -Qlim;
                // MT added to record Qlim
                Arrptr->LimQy[pTQ] = Q;
            }
            
        }
        else Q = 0.0;

    }
    else if (z0 + h0<z1 + h1 && h1>Solverptr->DepthThresh)
    {
        dh = z1 + h1 - z0 - h0;
        Sf = sqrt(dh / Parptr->dx);
        hflow = Tool::getmax(z0 + h0, z1 + h1) - Tool::getmax(z0, z1);
        hflow = Tool::getmax(hflow, 0);
        hflow = Tool::getmin(hflow, Solverptr->MaxHflow);
        // added to record Hflow
        // Arrptr->Hflowy[pTQ] = hflow;

        if (Tool::MaskTest(Arrptr->ChanMask[p0], Arrptr->ChanMask[p1]) && hflow > Solverptr->DepthThresh)
        {
            if (dh < Solverptr->dhlin)
            {
                Sf = sqrt(Parptr->dx / Solverptr->dhlin) * (dh / Parptr->dx);
                alpha = (pow(hflow, (5.0 / 3.0)) * Parptr->dx_sqrt) / (fn * sqrt(Solverptr->dhlin));
            }
            else alpha = pow(hflow, (5.0 / 3.0)) / (2. * fn * Sf);

            Q = (-pow(hflow, (5.0 / 3.0)) * Sf * Parptr->dy / fn);

            // flow limiter
            Qlim = Solverptr->Qlimfact * Parptr->dA * fabs(dh) / (8 * Solverptr->Tstep);
            if (fabs(Q) > Qlim)
            {
                if (Q > 0) Q = Qlim;
                if (Q < 0) Q = -Qlim;
                // MT added to record Qlim
                Arrptr->LimQy[pTQ] = Q;
            }
            
        }
        else Q = 0.0;
    }
    else Q = 0.0;
    // option to save V's
    if (Statesptr->voutput == ON)
    {
        if (Q != 0)
        {
            Arrptr->Vy[pTQ] = Q / Parptr->dx / hflow;
            Arrptr->maxVy[pTQ] = Tool::getmax(Arrptr->maxVy[pTQ], fabs(Arrptr->Vy[pTQ]));
        }
        else Arrptr->Vy[pTQ] = 0.0;
    }
    return(Q);
}



void Region::set_inflow_lisflood(int j, double value)
{
    double* data = get_BCVar(j);   // GET BCVAR
    int k = 0;
    while (data[k] != -1)//change Q
    {
        if (k % 2 == 0) {
            data[k] = value;// overflow  into LISFLOOD unit m^2/s
        }
        k++;
    }
}
int Region::get_gird_x(int j)
{
    return parameter->BCptr->xpi[j];
}
int Region::get_gird_y(int j)
{
    return parameter->BCptr->ypi[j];
}
void Region::set_waterdepth(int row, int col,double t)
{
       parameter->Arrptr->H[row + (parameter->Parptr->ysz - col) * (parameter->Parptr->xsz)] = t;
}
double Region::get_elevation(int row, int col)
{
    return parameter->Arrptr->H[row + (parameter->Parptr->ysz - col) * (parameter->Parptr->xsz)];;
}
double Region::get_waterdepth(int row, int col)
{
    return parameter->Arrptr->H[row + col * (parameter->Parptr->xsz)];;
}
double* Region::get_BCVar(int j)
{
    int nbdy = parameter->BCptr->PS_Val[j];
    return parameter->BCptr->BCVarlist[nbdy];
}
double Region::get_dx()
{
    return parameter->Parptr->dx;
}
map<string, int> Region::get_land_pipe_index()
{
    return land_pipe_index;
}

