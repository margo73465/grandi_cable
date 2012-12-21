/*********************************************************************
grandi_cable_3c6.c

  Program integrates the Grandi et al. 2011 atrial cell model in a 1D cable.
  Based on matlab code downloaded from website listed in article

  TKM removed Q10 powers, assume T=310.

  Version 1.0 reads in initial condition file 
  Version 1.1 has adiabatic elimination of CSQN plus analytical calculation 
    of CaSR (as in LRd). It also allows adaptive time-stepping and uses
	dt_max=0.01 ms (version 1.0 used 1.0 b/c stiff CSQN).
  Version 1.2 is cleaned up: removed temp stuff, cell geometry parameters, 
    Mgi, I_ClCFTR, Ko dependence of gKr and gKi
 
  Version 2.0 has a remapping of the variable y (and ydot) to remove
    un-used state variables
  Version 2.1 has separated out the variables y and ydot to individual variables 
  Version 2.2 uses analytical formulations of gating variables
  Version 2.3 has hard-coded many parameters
  Version 2.4 has small changes to speed things up (added RToF; changed calcs of INaK, 
    ICa, Incx, I_pca, and J_serca)
 
  Version 3.0 has look-up table for gating variables
  Version 3.1 has look-up table for other V-dependent exponentials
  Version 3.2 has calculations of APD, dV/dt_max, and CV
  Version 3.3 solves PDE 5 times per dt_max and has different stim amp
  Version 3.4 reads in new IC file
  Version 3.5 has dynamic variations in [K+]i (based on conservation of Nai and Ki)
  Version 3.6 uses openMP and only calculates INa,Late in case of CAF

 
  Trine Krogh-Madsen, Nov 2012

**********************************************************************/

#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

const double AF = 0.;
const double ISO = 0.;
const double RA = 0.;

const double pi = 3.141592653589793;


/* Stimulation parameters */
const double dt_max = 0.01;
const double T = 1000.;   //BCL
const double t_end = 500000.-.01; 
const double stim_dur = 1.0;
const double stim_amp = -48.0;
const int N_BCL = 100000;  // BCL/dt
const int N_stim = 100;    // stim_dur/dt

/* Cable parameters */
const double dx = 0.0125; //cm
const int DIM = 100;
const int rec1 = 3;
const int rec2 = 49;


int main(int argc, char *argv[]) {
  
  if(argc != 7) printf("Incorrect number of arguments.\n");
  
  printf("Int main okay. \n");

  //Timer
  time_t start, stop;
  time(&start);

  double D;
  D = atof(argv[1]);

  /* files */  
  FILE *outputI;
  FILE *outputI2;
  FILE *outputAPD;
  FILE *outputCV;
  FILE *outputVdmax;
  FILE *inputIC;
  outputI2 = fopen(argv[2],"w");
  outputI = fopen(argv[3],"w");

  printf("Opening files okay. \n");

  // Physical constants
  const double Frdy = 96485.;   // [C/mol]  
  const double FoRT = Frdy/8314./310.;
  const double RToF = 1./FoRT;
  const double Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10;//
  const double eps = 0.000005;

	
  // Cell geometry
  double Vcell = pi*10.25*10.25*100.*1.0e-15;    // [L]
  double Vmyo = 0.65*Vcell; 
  double Vsr = 0.035*Vcell; 
  double Vsl = 0.02*Vcell; 
  double Vjunc = 0.0539*.01*Vcell; 
  double J_ca_juncsl = 1./1.2134e12; // [L/msec] = 8.2413e-13
  double J_ca_slmyo = 1./2.68510e11; // [L/msec] = 3.2743e-12
  double J_na_juncsl = 3./1.6382e14; // [L/msec] = 6.1043e-13
  double J_na_slmyo = 3./1.8308e12;  // [L/msec] = 5.4621e-11
	
  // Fractional currents in compartments
  const double Fjunc = 0.11;   
  const double Fsl = 1.-Fjunc;
  const double Fjunc_CaL = 0.9; 
  const double Fsl_CaL = 1.-Fjunc_CaL;
	
  // Fixed ion concentrations     
  const double Ko = 5.4;   // Extracellular K   [mM]
  const double Nao = 140.;  // Extracellular Na  [mM]
  const double Cao = 1.8;  // Extracellular Ca  [mM]

  double ena_junc, ena_sl, ek, eca_junc, eca_sl, ecl ;            // [mV]
	
  double GNa = 23.*(1.-0.1*AF);  // [mS/uF]
  double GNaB = 0.597e-3;    // [mS/uF] 
  double IbarNaK = 1.26;     // [uA/uF]
  double KmNaip = 11.*(1.-0.25*ISO);         // [mM]11
  double KmKo = 1.5;         // [mM]1.5
  double GNaL = 0.0025*AF;

  // K current parameters
  double pNaK = 0.01833;      
  double gkp = 0.002;
	
  // Cl current parameters
  double GClCa = 0.0548;   // [mS/uF]
  double GClB = 9.0e-3;        // [mS/uF]
  double KdClCa = 100.0e-3;    // [mM]
	
  // I_Ca parameters
  double pNa = (1.+0.5*ISO)*(1.-0.5*AF)*0.75e-8;       // [cm/sec]
  double pCa = (1.+0.5*ISO)*(1.-0.5*AF)*2.7e-4;       // [cm/sec]
  double pK = (1.+0.5*ISO)*(1.-0.5*AF)*1.35e-7;        // [cm/sec]
	
  // Ca transport parameters
  double IbarNCX = (1.+0.4*AF)*3.15;      // [uA/uF]5.5 before - 9 in rabbit
  double KmCai = 3.59e-3;    // [mM]
  double KmCao = 1.3;        // [mM]
  double KmNaiCubed = 12.29*12.29*12.29;      // [mM]
  double KmNaoCubed = 87.5*87.5*87.5;  // [mM^3]
  double ksat = 0.27;        // [none]  
  double nu = 0.35;          // [none]
  double Kdact = 0.384e-3;   // [mM] 0.256 rabbit
  double IbarSLCaP = 0.0471; // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
  double KmPCap = 5.2282e-06;  //pow(0.5e-3,1.6);     // [mM] 
  double GCaB = 6.0643e-4;    // [uA/uF] 3
	
  // SR flux parameters
  double Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
  double Kmf = (2.5-1.25*ISO)*0.246e-3;          // [mM] default
  double Kmr = 1.7;               // [mM]L cytosol
  double koCa = 10.+20.*AF+10.*ISO*(1.-AF);               // [mM^-2 1/ms]   //default 10   modified 20
  double kom = 0.06;              // [1/ms]     
  double kim = 0.005;             // [1/ms]

  double temp, temp2, temp3, b, c;
  double I_Na_junc, I_Na_sl, I_Na, I_NaL_junc, I_NaL_sl, I_NaL, I_nabk_junc, I_nabk_sl, I_nabk;
  double I_nak_junc, I_nak_sl, I_nak;
  double I_kr, I_ks, I_kp, I_kur, I_ki;
  double I_ClCa_junc, I_ClCa_sl, I_ClCa, I_Clbk;
  double ibarca_j, ibarca_sl, ibark, ibarna_j, ibarna_sl, I_Ca_junc, I_Ca_sl, I_Ca;
  double I_CaK, I_CaNa_junc, I_CaNa_sl, I_CaNa, I_Catot;
  double Ka_junc, Ka_sl, s1_junc, s1_sl, s2_junc, s2_sl, s3_junc, s3_sl, I_ncx_junc, I_ncx_sl, I_ncx;
  double I_pca_junc, I_pca_sl, I_pca;
  double I_cabk_junc, I_cabk_sl, I_cabk;
  double I_Na_tot, I_Cl_tot, I_Ca_tot, I_tot, I_stim, Istim;
  double kCaSR, RI, J_SRCarel, J_serca, J_SRleak, J_CaB_cytosol;
  double J_CaB_junction, J_CaB_sl;
  double I_Na_tot_junc, I_Na_tot_sl, I_Na_tot_sl2, I_Na_tot_junc2;
  double I_K_tot, I_Ca_tot_junc, I_Ca_tot_sl, I_to;
  float r;

  double t = 0.0;
  double dt = dt_max, Vstart, dV[DIM], Vn[DIM], VP[DIM];
  double V_min = -100.0;
  int ilow, tl = 1600;
  double linext, Vt, Vtab[tl][34], etab[34];   
  double APDstart=0., dVmax=0., CVstart;
  double V_max = -100.0, V_low = 100.0, V_90; 
	
  int z, i = 0, ti, ii, k = 1, cc;
  int n_BCL = N_BCL, N = -1;
  int n_stim = 1;
  int stim_on = 0;
	
  double Vo[DIM], gate_m[DIM], gate_h[DIM], gate_j[DIM], gate_d[DIM], gate_f[DIM], gate_gj[DIM], gate_gsl[DIM];
  double gate_ml[DIM], gate_hl[DIM], gate_xto[DIM], gate_yto[DIM], gate_xr[DIM], gate_xs[DIM];
  double gate_xkur[DIM], gate_ykur[DIM], RYRr[DIM], RYRo[DIM], RYRi[DIM];
  double buf_TnC_low[DIM], buf_TnC_high[DIM], buf_Mg[DIM], buf_CaM[DIM], buf_myoCa[DIM], buf_myoMg[DIM], buf_SR[DIM], buf_SLLj[DIM], buf_SLLsl[DIM], buf_SLHj[DIM], buf_SLHsl[DIM]; 
  double buf_Naj[DIM], buf_Nasl[DIM], Naj[DIM], Nasl[DIM], Nai[DIM], Caj[DIM], Casl[DIM], Cai[DIM], CaSR[DIM];
  double dRYRr, dRYRo, dRYRi, dbuf_Naj, dbuf_Nasl; 
  double dbuf_TnC_low, dbuf_TnC_high, dbuf_Mg, dbuf_CaM, dbuf_myoCa, dbuf_myoMg, dbuf_SR; 
  double dbuf_SLLj, dbuf_SLLsl, dbuf_SLHj, dbuf_SLHsl;
  double dNaj, dNasl, dNai, dCaj, dCasl, dCai, Ki[DIM];
  double INa1,INa2;
		
  ecl = RToF*log(15./150.);

  // make look-up table; each column contains V-dep of an expression, the order of which is: 
  // V, mss, taum, hss, tauh, jss, tauj, mlss, tauml, hlinf, xrss, tauxr, rkr, xsss, tauxs, rkp, xtoss, tauxto, ytoss, tauyto, xkurss, tauxkur, ykurss, tauykur, kiss, dss, taud, fss, tauf, 6 diff exps
  for (i=0; i<tl; i++) {
    Vt = V_min+0.1*i;
    Vtab[i][0] = V_min+0.1*i;
    Vtab[i][1] = 1./((1.+exp(-(56.86+Vt)/9.03))*(1.+exp(-(56.86+Vt)/9.03)));
    Vtab[i][2] = 0.1292*exp(-((Vt+45.79)/15.54)*((Vt+45.79)/15.54))+0.06487*exp(-((Vt-4.823)/51.12)*((Vt-4.823)/51.12));
    Vtab[i][3] = 1./((1.+exp((Vt+71.55)/7.43))*(1.+exp((Vt+71.55)/7.43)));
    if (Vt>=-40.) {
      Vtab[i][4] = 1./(0.77/(0.13*(1.+exp(-(Vt+10.66)/11.1))));
      Vtab[i][6] = 1./((0.6*exp(0.057*Vt))/(1.+exp(-0.1*(Vt+32.))));
    }
    else {
      Vtab[i][4] = 1./( 0.057*exp(-(Vt+80.)/6.8) + 2.7*exp(0.079*Vt)+3.1e5*exp(0.3485*Vt) );
      Vtab[i][6] = 1./( (((-2.5428*1.0e4*exp(0.2444*Vt)-6.948e-6*exp(-0.04391*Vt))*(Vt+37.78))/(1.+exp(0.311*(Vt+79.23)))) + ((0.02424*exp(-0.01052*Vt))/(1.+exp(-0.1378*(Vt+40.14)))) );
    }	
    Vtab[i][5] = 1./((1.+exp((Vt+71.55)/7.43))*(1.+exp((Vt+71.55)/7.43)));		
    Vtab[i][7] = (0.32*(Vt+47.13)/(1.-exp(-0.1*(Vt+47.13)))) / ( 0.32*(Vt+47.13)/(1.-exp(-0.1*(Vt+47.13))) + 0.08*exp(-Vt/11.) );
    Vtab[i][8] = 1./ ( 0.32*(Vt+47.13)/(1.-exp(-0.1*(Vt+47.13))) + 0.08*exp(-Vt/11.) );
    Vtab[i][9] = 1./(1.+exp((Vt+91.)/6.1));
    Vtab[i][10] = 1./(1.+exp(-(Vt+10.)/5.));
    Vtab[i][11] = 550./(1.+exp((-22.-Vt)/9.))*6./(1.+exp((Vt+11.)/9.))+230./(1.+exp((Vt+40.)/20.));
    Vtab[i][12] = 1./(1.+exp((Vt+74.)/24.));
    Vtab[i][13] = 1./(1.+exp(-(Vt+40.*ISO+3.8)/14.25));
    Vtab[i][14] = 990.1/(1.+exp(-(Vt+40.*ISO+2.436)/14.12));
    Vtab[i][15] = 1./(1.+exp(7.488-Vt/5.98));
    Vtab[i][16] = 1./(1.+exp(-(Vt+1.0)/11.0));
    Vtab[i][17] = 3.5*exp(-((Vt/30.0)*(Vt/30.0)))+1.5;
    Vtab[i][18] = 1./(1.+exp((Vt+40.5)/11.5));
    Vtab[i][19] = 25.635*exp(-(((Vt+52.45)/15.8827)*((Vt+52.45)/15.8827)))+24.14;
    Vtab[i][20] = 1./(1.+exp((Vt+6.)/-8.6));
    Vtab[i][21] = 9./(1.+exp((Vt+5.)/12.0))+0.5;
    Vtab[i][22] = 1./(1.+exp((Vt+7.5)/10.));
    Vtab[i][23] = 590./(1.+exp((Vt+60.)/10.0))+3050.;
    Vtab[i][24] = exp((nu-1.)*Vt*FoRT);
    Vtab[i][25] = 1./(1.+exp(-(Vt+3.*ISO+9.)/6.)); 
    if ( (ISO<.5 && Vt>-9.05 && Vt<-8.95) || (ISO<.5 && Vt>-12.05 && Vt<-11.95) ) Vtab[i][26] = 2.38095;
    else Vtab[i][26] = Vtab[i][25]*(1.-exp(-(Vt+3.*ISO+9.)/6.))/(0.035*(Vt+3.*ISO+9.)); 
    Vtab[i][27] = 1./(1.+exp((Vt+3.*ISO+30.)/7.))+0.2/(1.+exp((50.-Vt-3.*ISO)/20.)); 
    Vtab[i][28] = 1./(0.0197*exp( -(0.0337*(Vt+3.*ISO+25.))*(0.0337*(Vt+3.*ISO+25.)) )+0.02);
    Vtab[i][29] = exp(-0.1*Vt*FoRT);
    Vtab[i][30] = exp(-Vt*FoRT);
    Vtab[i][31] = exp(2.*Vt*FoRT);
    Vtab[i][32] = exp(Vt*FoRT);
    Vtab[i][33] = exp(nu*Vt*FoRT);
  }
  i=0;
	
	
  /* initial conditions */
  for (z=0; z<DIM; z++)  {	
    inputIC = fopen("grandi_IC_cell_3c3_def_1Hz_2000AP.dat","r");
    fscanf(inputIC,"%f",&r);
    Vo[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_m[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_h[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_j[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_d[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_f[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_gj[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_gsl[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_ml[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_hl[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_xto[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_yto[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_xr[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_xs[z] = r;
    fscanf(inputIC,"%f",&r);
    RYRr[z] = r;
    fscanf(inputIC,"%f",&r);
    RYRo[z] = r;
    fscanf(inputIC,"%f",&r);
    RYRi[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_Naj[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_Nasl[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_TnC_low[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_TnC_high[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_Mg[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_CaM[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_myoCa[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_myoMg[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_SR[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_SLLj[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_SLLsl[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_SLHj[z] = r;
    fscanf(inputIC,"%f",&r);
    buf_SLHsl[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_xkur[z] = r;
    fscanf(inputIC,"%f",&r);
    gate_ykur[z] = r;
    fscanf(inputIC,"%f",&r);
    Naj[z] = r;
    fscanf(inputIC,"%f",&r);
    Nasl[z] = r;
    fscanf(inputIC,"%f",&r);
    Nai[z] = r;
    fscanf(inputIC,"%f",&r);
    Caj[z] = r;
    fscanf(inputIC,"%f",&r);
    Casl[z] = r;
    fscanf(inputIC,"%f",&r);
    Cai[z] = r;
    fscanf(inputIC,"%f",&r);
    CaSR[z] = r;
    fclose(inputIC);
    Ki[z] = 120.;
    dV[z] = 0.;
    Vn[z] = Vo[z];
    VP[z] = Vo[z];
  }	

	
  while (t < t_end) {

    printf("Entering time loop okay. \n");

    /* compute stimulus current */
    Istim = 0.0;
    if (n_BCL == N_BCL) {
      stim_on = 1;
      n_BCL = 1;
      N++;
    }
    if (n_stim <= N_stim && stim_on == 1) {
      Istim = stim_amp;
      n_stim ++;
    }
    if (n_stim == N_stim) {
      stim_on = 0;
      n_stim = 1;
    }

	  
    /* ODE all inner grid points */
    #pragma omp parallel for private(I_stim,Vstart,k,dt,ii,ilow,linext,etab,ti, \
      ena_junc,ena_sl,eca_junc,eca_sl,ek, \
      I_Na,I_Na_junc,I_Na_sl,I_NaL,I_NaL_junc,I_NaL_sl,I_nabk_junc,I_nabk_sl,I_nabk, \
      temp,temp2,temp3,I_nak_junc,I_nak_sl,I_nak,I_kr,I_ks,I_kp,I_to,I_kur,I_ki, \
      I_ClCa_junc,I_ClCa_sl,I_ClCa,I_Clbk,ibarca_j,ibarca_sl,ibark,ibarna_j,ibarna_sl, \
      I_Ca_junc,I_Ca_sl,I_Ca,I_CaK,I_CaNa_junc,I_CaNa_sl,I_CaNa,I_Catot,Ka_junc,Ka_sl, \
      s1_junc,s1_sl,s2_junc,s2_sl,s3_junc,s3_sl,I_ncx_junc,I_ncx_sl,I_ncx, \
      I_pca_junc,I_pca_sl,I_pca,I_cabk_junc,I_cabk_sl,I_cabk,kCaSR,RI,J_SRCarel,J_serca,J_SRleak, \
      J_CaB_cytosol,J_CaB_junction,J_CaB_sl,b,c,I_Na_tot_junc,I_Na_tot_sl,I_Na_tot_sl2,I_Na_tot_junc2, \
      I_K_tot,I_Ca_tot_junc,I_Ca_tot_sl,I_Na_tot,I_Cl_tot,I_Ca_tot,I_tot, \
      dRYRr,dRYRo,dRYRi,dbuf_Naj,dbuf_Nasl,dbuf_TnC_low,dbuf_TnC_high,dbuf_Mg,dbuf_CaM, \
      dbuf_myoCa,dbuf_myoMg,dbuf_SR,dbuf_SLLj,dbuf_SLLsl,dbuf_SLHj,dbuf_SLHsl,dNaj,dNasl,dNai, \
      dCaj,dCasl,dCai)
    for (z=1; z<DIM-1; z++) {  
			
      /* compute stimulus current */	
      I_stim = 0.0;
      if (z<8) I_stim = Istim;
		  
      /* compute number of adaptive time steps */	
      Vstart = Vo[z];	  	  
      k = fabs(dV[z]/dt_max);
      if (dV[z]/dt_max > 0.1) {
        k = k+3;
        if (k > 10) k = 10;
      }		
      else {
        if (k > 10) k = 10;
        if (k < 1) k = 1;
      }
      dt = dt_max/k;  
	  
	
      for (ii=0; ii<k; ii++) {  //solve ODE's k times
  
	  ena_junc = RToF*log(Nao/Naj[z]);     
	  ena_sl = RToF*log(Nao/Nasl[z]);     
	  eca_junc = RToF/2.*log(Cao/Caj[z]);   
	  eca_sl = RToF/2.*log(Cao/Casl[z]);  
	  ek = RToF*log(Ko/Ki[z]);	


	  // extrapolate rate constants and other V-deps from table
	  ilow = fabs((Vo[z]-V_min)/0.1);
	  linext = (Vo[z]-Vtab[ilow][0])/0.1;
	  for (ti=1; ti<34; ti++) {
	    etab[ti] = (Vtab[ilow+1][ti]-Vtab[ilow][ti])*linext+Vtab[ilow][ti];
	  }	
		
	  // I_Na
	  I_Na_junc = Fjunc*GNa*gate_m[z]*gate_m[z]*gate_m[z]*gate_h[z]*gate_j[z]*(Vo[z]-ena_junc);
	  I_Na_sl = Fsl*GNa*gate_m[z]*gate_m[z]*gate_m[z]*gate_h[z]*gate_j[z]*(Vo[z]-ena_sl);
	  I_Na = I_Na_junc+I_Na_sl;
	  gate_m[z] = etab[1]-(etab[1]-gate_m[z])*exp(-dt/etab[2]);
	  gate_h[z] = etab[3]-(etab[3]-gate_h[z])*exp(-dt/etab[4]);
	  gate_j[z] = etab[5]-(etab[5]-gate_j[z])*exp(-dt/etab[6]);
          if (z==rec1) INa1 = I_Na;
          if (z==rec2) INa2 = I_Na;
		
	  // Late I_Na
	  if (AF<eps) {
	    I_NaL_junc = 0.;
	    I_NaL_sl = 0.;
	    I_NaL = 0.;
	  }
	  else{
       	    I_NaL_junc = Fjunc*GNaL*gate_ml[z]*gate_ml[z]*gate_ml[z]*gate_hl[z]*(Vo[z]-ena_junc);
	    I_NaL_sl = Fsl*GNaL*gate_ml[z]*gate_ml[z]*gate_ml[z]*gate_hl[z]*(Vo[z]-ena_sl);
	    I_NaL = I_NaL_junc+I_NaL_sl;
	    gate_ml[z] = etab[7]-(etab[7]-gate_ml[z])*exp(-dt/etab[8]);
	    gate_hl[z] = etab[9]-(etab[9]-gate_hl[z])*exp(-dt/600.);
	  }
		  
	  // I_nabk: Na Background Current
	  I_nabk_junc = Fjunc*GNaB*(Vo[z]-ena_junc);
	  I_nabk_sl = Fsl*GNaB*(Vo[z]-ena_sl);
	  I_nabk = I_nabk_junc+I_nabk_sl;
	  
	  // I_nak: Na/K Pump Current
	  temp = (exp(Nao/67.3)-1.)/7.;
	  temp2 = IbarNaK*Ko/(1.+0.1245*etab[29]+0.0365*temp*etab[30])/(Ko+KmKo);
	  temp3 = KmNaip/Naj[z];
	  I_nak_junc = Fjunc*temp2/(1.+temp3*temp3*temp3*temp3);
	  temp3 = KmNaip/Nasl[z];
	  I_nak_sl = Fsl*temp2/(1.+temp3*temp3*temp3*temp3);
	  I_nak = I_nak_junc+I_nak_sl;
	  
	  // I_kr: Rapidly Activating K Current
	  I_kr = 0.035*gate_xr[z]*etab[12]*(Vo[z]-ek);
	  gate_xr[z] = etab[10]-(etab[10]-gate_xr[z])*exp(-dt/etab[11]);
		
	  // I_ks: Slowly Activating K Current
	  temp2 = RToF*log((Ko+pNaK*Nao)/(Ki[z]+pNaK*Nai[z]));
	  temp = (1.+AF+2.*ISO)*0.0035;
	  I_ks = temp*gate_xs[z]*gate_xs[z]*(Vo[z]-temp2);
	  gate_xs[z] = etab[13]-(etab[13]-gate_xs[z])*exp(-dt/etab[14]);

	  //I_kp: Plateau K current
	  I_kp = gkp*etab[15]*(Vo[z]-ek);
	  
	  // I_to: Transient Outward K Current (fast component)
	  temp = (1.0-0.7*AF)*0.165; //nS/pF maleckar; //human atrium
	  I_to = temp*gate_xto[z]*gate_yto[z]*(Vo[z]-ek);
	  gate_xto[z] = etab[16]-(etab[16]-gate_xto[z])*exp(-dt/etab[17]);
	  gate_yto[z] = etab[18]-(etab[18]-gate_yto[z])*exp(-dt/etab[19]);

	  // I_kur: Ultra rapid delayed rectifier Outward K Current
	  temp = (1.0-0.5*AF)*(1.+2.*ISO)*0.045*(1.+0.2*RA); //nS/pF maleckar 0.045
	  I_kur = temp*gate_xkur[z]*gate_ykur[z]*(Vo[z]-ek);
	  gate_xkur[z] = etab[20]-(etab[20]-gate_xkur[z])*exp(-dt/etab[21]);
	  gate_ykur[z] = etab[22]-(etab[22]-gate_ykur[z])*exp(-dt/etab[23]);
	  
	  // I_ki: Time-Independent K Current
	  temp = (1.02/(1.+exp(0.2385*(Vo[z]-ek-59.215)))) / ( 1.02/(1.+exp(0.2385*(Vo[z]-ek-59.215))) + (0.49124*exp(0.08032*(Vo[z]+5.476-ek))+exp(0.06175*(Vo[z]-ek-594.31)))/(1.+exp(-0.5143*(Vo[z]-ek+4.753))) );
	  I_ki = (1.+AF)*0.0525*temp*(Vo[z]-ek);
	  
	  // I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
	  I_ClCa_junc = Fjunc*GClCa/(1.+KdClCa/Caj[z])*(Vo[z]-ecl);
	  I_ClCa_sl = Fsl*GClCa/(1.+KdClCa/Casl[z])*(Vo[z]-ecl);
	  I_ClCa = I_ClCa_junc+I_ClCa_sl;
	  I_Clbk = GClB*(Vo[z]-ecl);
	  	  
	  // I_Ca: L-type Calcium Current
	  temp = Vo[z]*Frdy*FoRT;
	  ibarca_j = pCa*4.*temp * (0.341*Caj[z]*etab[31]-0.341*Cao) /(etab[31]-1.);
	  ibarca_sl = pCa*4.*temp * (0.341*Casl[z]*etab[31]-0.341*Cao) /(etab[31]-1.);
	  ibark = pK*temp*(0.75*Ki[z]*etab[32]-0.75*Ko) /(etab[32]-1.);
	  ibarna_j = pNa*temp *(0.75*Naj[z]*etab[32]-0.75*Nao) /(etab[32]-1.);
	  ibarna_sl = pNa*temp *(0.75*Nasl[z]*etab[32]-0.75*Nao) /(etab[32]-1.);
	  I_Ca_junc = (Fjunc_CaL*ibarca_j*gate_d[z]*gate_f[z]*(1.-gate_gj[z]))*0.45;
	  I_Ca_sl = (Fsl_CaL*ibarca_sl*gate_d[z]*gate_f[z]*((1.-gate_gsl[z])))*0.45;
	  I_Ca = I_Ca_junc+I_Ca_sl;
	  I_CaK = (ibark*gate_d[z]*gate_f[z]*(Fjunc_CaL*(1.-gate_gj[z])+Fsl_CaL*(1.-gate_gsl[z])))*0.45;
	  I_CaNa_junc = (Fjunc_CaL*ibarna_j*gate_d[z]*gate_f[z]*(1.-gate_gj[z]))*0.45;
	  I_CaNa_sl = (Fsl_CaL*ibarna_sl*gate_d[z]*gate_f[z]*(1.-gate_gsl[z]))*0.45;
	  I_CaNa = I_CaNa_junc+I_CaNa_sl;
	  I_Catot = I_Ca+I_CaK+I_CaNa;
	  gate_d[z] = etab[25]-(etab[25]-gate_d[z])*exp(-dt/etab[26]);
	  gate_f[z] = etab[27]-(etab[27]-gate_f[z])*exp(-dt/etab[28]);
	  temp = 1./(1.7*Caj[z]+11.9e-3);
	  temp2 = 1.7*Caj[z]*temp;
	  gate_gj[z] = temp2-(temp2-gate_gj[z])*exp(-dt/temp);
	  temp = 1./(1.7*Casl[z]+11.9e-3);
	  temp2 = 1.7*Casl[z]*temp;
	  gate_gsl[z] = temp2-(temp2-gate_gsl[z])*exp(-dt/temp);
		
	  // I_ncx: Na/Ca Exchanger flux
	  Ka_junc = 1./(1.+(Kdact/Caj[z])*(Kdact/Caj[z]));
	  Ka_sl = 1./(1.+(Kdact/Casl[z])*(Kdact/Casl[z]));
	  temp3 = Naj[z]*Naj[z]*Naj[z];	
	  temp2 = Nasl[z]*Nasl[z]*Nasl[z];	
	  s1_junc = etab[33]*temp3*Cao;
	  s1_sl = etab[33]*temp2*Cao;
	  temp = Nao*Nao*Nao;	
	  s2_junc = etab[24]*temp*Caj[z];
	  s3_junc = KmCai*temp*(1.+temp3/KmNaiCubed) + KmNaoCubed*Caj[z]*(1.+Caj[z]/KmCai) + KmCao*temp3 + temp3*Cao + temp*Caj[z];
	  s2_sl = etab[24]*temp*Casl[z];
	  s3_sl = KmCai*temp*(1.+temp2/KmNaiCubed) + KmNaoCubed*Casl[z]*(1.+Casl[z]/KmCai) + KmCao*temp2 + temp2*Cao + temp*Casl[z];	  
	  I_ncx_junc = Fjunc*IbarNCX*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1.+ksat*etab[24]);
	  I_ncx_sl = Fsl*IbarNCX*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1.+ksat*etab[24]);
	  I_ncx = I_ncx_junc+I_ncx_sl;
	  
	  // I_pca: Sarcolemmal Ca Pump Current
	  temp2 = pow(Caj[z],1.6);	
	  temp3 = pow(Casl[z],1.6);	
	  I_pca_junc = Fjunc*IbarSLCaP*temp2/(KmPCap+temp2);
	  I_pca_sl = Fsl*IbarSLCaP*temp3/(KmPCap+temp3);
	  I_pca = I_pca_junc+I_pca_sl;
	  
	  // I_cabk: Ca Background Current
	  I_cabk_junc = Fjunc*GCaB*(Vo[z]-eca_junc);
	  I_cabk_sl = Fsl*GCaB*(Vo[z]-eca_sl);
	  I_cabk = I_cabk_junc+I_cabk_sl;
	  
	  // SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
	  kCaSR = 15.-14./(1.+pow((0.45/CaSR[z]),2.5));
	  temp = koCa/kCaSR;
	  temp2 = 0.5*kCaSR;
	  RI = 1.-RYRr[z]-RYRo[z]-RYRi[z];
	  dRYRr = (kim*RI-temp2*Caj[z]*RYRr[z])-(temp*Caj[z]*Caj[z]*RYRr[z]-kom*RYRo[z]);   // R
	  dRYRo = (temp*Caj[z]*Caj[z]*RYRr[z]-kom*RYRo[z])-(temp2*Caj[z]*RYRo[z]-kim*RYRi[z]);// O
	  dRYRi = (temp2*Caj[z]*RYRo[z]-kim*RYRi[z])-(kom*RYRi[z]-temp*Caj[z]*Caj[z]*RI);   // I
	  J_SRCarel = 25.*RYRo[z]*(CaSR[z]-Caj[z]);          // [mM/ms]	 
	  temp = pow((Cai[z]/Kmf),1.787);
	  temp2 = pow((CaSR[z]/Kmr),1.787);
	  J_serca = Vmax_SRCaP*(temp-temp2)/(1.+temp+temp2);
	  J_SRleak = (1.0+0.25*AF)*5.348e-6*(CaSR[z]-Caj[z]);           //   [mM/ms]
	  
	  // Sodium and Calcium Buffering
	  dbuf_Naj = 0.1e-3*Naj[z]*(7.561-buf_Naj[z])-1.0e-3*buf_Naj[z];        // NaBj      [mM/ms]
	  dbuf_Nasl = 0.1e-3*Nasl[z]*(1.65-buf_Nasl[z])-1.0e-3*buf_Nasl[z];       // NaBsl     [mM/ms]
	  
	  // Cytosolic Ca Buffers
	  dbuf_TnC_low = 32.7*Cai[z]*(70.0e-3-buf_TnC_low[z])-(1.+0.5*ISO)*19.6e-3*buf_TnC_low[z];            // TnCL      [mM/ms]
	  dbuf_TnC_high = 2.37*Cai[z]*(140.0e-3-buf_TnC_high[z]-buf_Mg[z])-0.032e-3*buf_TnC_high[z]; // TnCHc     [mM/ms]
	  dbuf_Mg = 3.0e-3*(140.0e-3-buf_TnC_high[z]-buf_Mg[z])-3.33e-3*buf_Mg[z];   // TnCHm     [mM/ms]
	  dbuf_CaM = 34.*Cai[z]*(24.0e-3-buf_CaM[z])-238.0e-3*buf_CaM[z];                 // CaM       [mM/ms]
	  dbuf_myoCa = 13.8*Cai[z]*(140.0e-3-buf_myoCa[z]-buf_myoMg[z])-0.46e-3*buf_myoCa[z];    // Myosin_ca [mM/ms]
	  dbuf_myoMg = 0.0157*(140.0e-3-buf_myoCa[z]-buf_myoMg[z])-0.057e-3*buf_myoMg[z];      // Myosin_mg [mM/ms]
	  dbuf_SR = 100.*Cai[z]*(19.*.9e-3-buf_SR[z])-60.0e-3*buf_SR[z];                    // SRB       [mM/ms]
	  J_CaB_cytosol = dbuf_TnC_low+dbuf_TnC_high+dbuf_Mg+dbuf_CaM+dbuf_myoCa+dbuf_myoMg+dbuf_SR;
	  
	  // Junctional and SL Ca Buffers
	  dbuf_SLLj = 100.*Caj[z]*(4.6e-4*Vmyo/Vjunc-buf_SLLj[z])-1300.0e-3*buf_SLLj[z];       // SLLj      [mM/ms]
	  dbuf_SLLsl = 100.*Casl[z]*(37.4e-3*Vmyo/Vsl-buf_SLLsl[z])-1300.0e-3*buf_SLLsl[z];      // SLLsl     [mM/ms]
	  dbuf_SLHj = 100.*Caj[z]*(1.65e-4*Vmyo/Vjunc-buf_SLHj[z])-30.0e-3*buf_SLHj[z];      // SLHj      [mM/ms]
	  dbuf_SLHsl = 100.*Casl[z]*(13.4e-3*Vmyo/Vsl-buf_SLHsl[z])-30.0e-3*buf_SLHsl[z];     // SLHsl     [mM/ms]
	  J_CaB_junction = dbuf_SLLj+dbuf_SLHj;
	  J_CaB_sl = dbuf_SLLsl+dbuf_SLHsl;
	  
	  // Ion concentrations
	  // SR Ca Concentrations
	  temp2 = 2.6*CaSR[z]/(.65+CaSR[z]);
	  temp = dt*(J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel));
	  b = 2.6-temp2-temp-CaSR[z]+.65;
	  c = .65*(temp2+temp+CaSR[z]);
	  CaSR[z] = (sqrt(b*b+4.*c)-b)/2.;
		
	  // Sodium Concentrations
	  I_Na_tot_junc = I_Na_junc+I_nabk_junc+3.*I_ncx_junc+3.*I_nak_junc+I_CaNa_junc+I_NaL_junc;   // [uA/uF]
	  I_Na_tot_sl = I_Na_sl+I_nabk_sl+3.*I_ncx_sl+3.*I_nak_sl+I_CaNa_sl+I_NaL_sl;   // [uA/uF]
	  I_Na_tot_sl2 = 3.*I_ncx_sl+3.*I_nak_sl+I_CaNa_sl;   // [uA/uF]
	  I_Na_tot_junc2 = 3.*I_ncx_junc+3.*I_nak_junc+I_CaNa_junc;   // [uA/uF]	  
	  dNaj = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(Nasl[z]-Naj[z])-dbuf_Naj;
	  dNasl = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(Naj[z]-Nasl[z])+J_na_slmyo/Vsl*(Nai[z]-Nasl[z])-dbuf_Nasl;
	  dNai = J_na_slmyo/Vmyo*(Nasl[z]-Nai[z]);             // [mM/msec] 
	  
	  // Potassium Concentration
	  I_K_tot = I_to+I_kr+I_ks+I_ki-2.*I_nak+I_CaK+I_kp+I_kur;     // [uA/uF] //SVP: added IKur
	  Ki[z] = 120.+9.2-Nai[z]; 	
		  
	  // Calcium Concentrations
	  I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2.*I_ncx_junc;                   // [uA/uF]
	  I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2.*I_ncx_sl;            // [uA/uF]
	  dCaj = -I_Ca_tot_junc*Cmem/(Vjunc*2.*Frdy)+J_ca_juncsl/Vjunc*(Casl[z]-Caj[z])-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  // Ca_j
	  dCasl = -I_Ca_tot_sl*Cmem/(Vsl*2.*Frdy)+J_ca_juncsl/Vsl*(Caj[z]-Casl[z]) + J_ca_slmyo/Vsl*(Cai[z]-Casl[z])-J_CaB_sl;   // Ca_sl
	  dCai = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(Casl[z]-Cai[z]);

	  // Membrane Potential
	  I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          // [uA/uF]
	  I_Cl_tot = I_ClCa+I_Clbk;                        // [uA/uF]
	  I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
	  I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
	  Vo[z] = Vo[z] - (I_tot+I_stim)*dt;
	  			
	  // Integration of remaining vars
	  RYRr[z] = RYRr[z]+dt*dRYRr;
	  RYRo[z] = RYRo[z]+dt*dRYRo;
	  RYRi[z] = RYRi[z]+dt*dRYRi;
	  buf_Naj[z] = buf_Naj[z]+dt*dbuf_Naj;
	  buf_Nasl[z] = buf_Nasl[z]+dt*dbuf_Nasl;
	  buf_TnC_low[z] = buf_TnC_low[z]+dt*dbuf_TnC_low;
	  buf_TnC_high[z] = buf_TnC_high[z]+dt*dbuf_TnC_high;
	  buf_Mg[z] = buf_Mg[z]+dt*dbuf_Mg;
	  buf_CaM[z] = buf_CaM[z]+dt*dbuf_CaM;
	  buf_myoCa[z] = buf_myoCa[z]+dt*dbuf_myoCa;
	  buf_myoMg[z] = buf_myoMg[z]+dt*dbuf_myoMg;
	  buf_SR[z] = buf_SR[z]+dt*dbuf_SR;
	  buf_SLLj[z] = buf_SLLj[z]+dt*dbuf_SLLj;
	  buf_SLLsl[z] = buf_SLLsl[z]+dt*dbuf_SLLsl;
	  buf_SLHj[z] = buf_SLHj[z]+dt*dbuf_SLHj;
	  buf_SLHsl[z] = buf_SLHsl[z]+dt*dbuf_SLHsl;
	  Naj[z] = Naj[z]+dt*dNaj;
	  Nasl[z] = Nasl[z]+dt*dNasl;
	  Nai[z] = Nai[z]+dt*dNai;
	  Caj[z] = Caj[z]+dt*dCaj;
	  Casl[z] = Casl[z]+dt*dCasl;
	  Cai[z] = Cai[z]+dt*dCai;
		
	} //end ODE loop		
																		   
	dV[z] = Vo[z]-Vstart;  
		
		
	/* ODE boundaries */
	Vo[0] = Vo[2];
	Vo[DIM-1] = Vo[DIM-3];
		
    } //end z-loop
    printf("End of z-loop.\n");
	  
	/* PDE all inner grid points */ 
	for (cc=0; cc<5; cc++) {
	  for (z=1; z<DIM-1; z++)  {
		Vn[z] = dt_max/(dx*dx)*D*(Vo[z+1]+Vo[z-1]-2.0*Vo[z]) + Vo[z];
	  }
	  
	  /* boundaries */
	  Vn[0] = Vn[2];
	  Vn[DIM-1] = Vn[DIM-3];
	  
	  /* update */
	  for (z=0; z<DIM; z++) {
		Vo[z] = Vn[z];
	  }
	}
	  
	printf("PDE calculated.\n");


	// APD
	if (Vn[49]>=-60. && VP[49]<-60.) APDstart = t;
	if (Vn[49]<=-60. && VP[49]>-60.) {
	  outputAPD = fopen(argv[4],"a");
	  fprintf(outputAPD,"%.2f\t%.2f\n",t,t-APDstart);
	  fclose(outputAPD);
	  printf("APD recorded.\n");
	  }

	/*APD_90
	if (n_BCL == N_BCL-1){
	  V_low = Vn[49];
	  APDstart = t+2*dt_max;
	  V_max = -100.;
	}
	if (Vn[49] > V_max){
	  V_max = Vn[49];
	  V_90 = V_max - (V_max - V_low)*0.90;
	}
	if (Vn[49] <= V_90 && Vstart > V_90 && t-APDstart>0){
	  outputAPD = fopen("grandi_APD90_cable_3c6_1Hz_def_D0c5_dx0c0125_dtmax0c01.dat","a");
	  fprintf(outputAPD, "%.2f \t %.2f \n",t,t-APDstart);
	  fclose(outputAPD);
	  }*/
		  
	// Vdot_max 
	if (dV[49]/dt_max>dVmax) dVmax = dV[49]/dt_max;
	if (Vn[49]<=0. && VP[49]>0.) {
	  outputVdmax = fopen(argv[5],"a");
       	  fprintf(outputVdmax,"%.2f\t%.2f\n",t,dVmax); 
	  fclose(outputVdmax);
	  dVmax = 0.;
	  printf("Vdmax recorded.\n");
        }

	// CV
	if (Vn[19]>=-40. && VP[19]<-40.) CVstart = t;
	if (Vn[79]>=-40. && VP[79]<-40.) {
	  outputCV = fopen(argv[6],"a");
	  fprintf(outputCV,"%.2f\t%.2f\n",t,60.*dx/(t-CVstart)*1000.);
	  fclose(outputCV);
	  printf("CV recorded.\n");
	}
	  
	/* update */
	for (z=0; z<DIM; z++) {
	  VP[z] = Vn[z];
	}
																   
	/* write to file */ 	  
	if (i==10)  {
	  fprintf(outputI,"%.1f\t%.2f\t%f\t%.2f\t%f\n",t,Vo[rec1],Cai[rec1],Vo[rec2],Cai[rec2]);
	  fprintf(outputI2,"%.1f\t%f\t%f\t%f\t%f\n",t,Ki[rec1],Ki[rec2],INa1,INa2);
	  i=0;
	}
	  
    t = N*T + n_BCL*dt_max;
    n_BCL++;
    //if (t<2.*T) i++;
    if (t>t_end-T) i++;
	  
  }

  fclose(outputI);
  fclose(outputI2);

  time(&stop);
  printf("Finished in %.0f seconds \n", difftime(stop,start));
  return 0;

}
