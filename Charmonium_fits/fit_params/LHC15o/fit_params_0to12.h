//-----------------------------------------------------------------
//
// systematic uncertainties display parameters
//
//-----------------------------------------------------------------

    double display_NJPsi_range = 50000;
    double display_Sigma_range = 0.006;
    double display_Mass_range = 0.005;

//-----------------------------------------------------------------
//
// crystal ball function + variable width gaussian
//
//-----------------------------------------------------------------

//Tail: pp13TeV. Fit Range: 2.4 4.7
double pp13_fit0_par0_CB2_VWG[] = {1.,0.,10.}; // initial value, lower limit, upper limit
double pp13_fit0_par1_CB2_VWG[] = {1.,-10.,10.};
double pp13_fit0_par2_CB2_VWG[] = {1.,-10.,10.};
double pp13_fit0_par3_CB2_VWG[] = {1.,-10.,10.};
double pp13_fit0_par4_CB2_VWG[] = {1.,-10.,10.};

double pp13_fit0_par0_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit0_par1_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit0_par2_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit0_par3_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit0_par4_variance_CB2_VWG[] = {0.,0.,0.};

double pp13_fit0_mean_jpsi_CB2_VWG[] = {3.10, 3.07, 3.20};
double pp13_fit0_width_jpsi_CB2_VWG[] = {0.07, 0.04, 0.10};
double pp13_fit0_norm_psi2S_CB2_VWG[] = {0.1, 0.0, 0.5};


//Tail: pp13TeV. Fit Range: 2.2 4.4
double pp13_fit1_par0_CB2_VWG[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double pp13_fit1_par1_CB2_VWG[] = {1,-10.,10.};
double pp13_fit1_par2_CB2_VWG[] = {0.1,-10.,10.};
double pp13_fit1_par3_CB2_VWG[] = {1,-10.,10.};
double pp13_fit1_par4_CB2_VWG[] = {0.1,-10.,10.};

double pp13_fit1_par0_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit1_par1_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit1_par2_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit1_par3_variance_CB2_VWG[] = {0.,0.,0.};
double pp13_fit1_par4_variance_CB2_VWG[] = {0.,0.,0.};

double pp13_fit1_mean_jpsi_CB2_VWG[] = {3.10, 3.07, 3.20};
double pp13_fit1_width_jpsi_CB2_VWG[] = {0.07, 0.04, 0.10};
double pp13_fit1_norm_psi2S_CB2_VWG[] = {0.1, 0, 0.5};

// Tail: geant3. Fit Range: 2.2 4.5
double geant3_fit0_par0_CB2_VWG[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_CB2_VWG[] = {1,-10.,10.};
double geant3_fit0_par2_CB2_VWG[] = {1,-10.,10.};
double geant3_fit0_par3_CB2_VWG[] = {1,-10.,10.};
double geant3_fit0_par4_CB2_VWG[] = {1,-10.,10.};

double geant3_fit0_par0_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit0_par1_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit0_par2_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit0_par3_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit0_par4_variance_CB2_VWG[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_CB2_VWG[] = {3.10, 3.07, 3.20};
double geant3_fit0_width_jpsi_CB2_VWG[] = {0.07, 0.04, 0.10};
double geant3_fit0_norm_psi2S_CB2_VWG[] = {0.1, 0, 0.5};


// Tail: geant3. Fit Range: 2.2 4.4
double geant3_fit1_par0_CB2_VWG[] = {1.,0.,10.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_CB2_VWG[] = {1,-10.,10.};
double geant3_fit1_par2_CB2_VWG[] = {1,-10.,10.};
double geant3_fit1_par3_CB2_VWG[] = {1,-10.,10.};
double geant3_fit1_par4_CB2_VWG[] = {0.1,-5.,5.};

double geant3_fit1_par0_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit1_par1_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit1_par2_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit1_par3_variance_CB2_VWG[] = {0.,0.,0.};
double geant3_fit1_par4_variance_CB2_VWG[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_CB2_VWG[] = {3.10, 3.07, 3.20};
double geant3_fit1_width_jpsi_CB2_VWG[] = {0.07, 0.04, 0.10};
double geant3_fit1_norm_psi2S_CB2_VWG[] = {0.1, 0, 0.5};

//-----------------------------------------------------------------
//
// crystal ball function + pol2/pol3
//
//-----------------------------------------------------------------

double pp13_fit0_par0_CB2_pol[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double pp13_fit0_par1_CB2_pol[] = {1.,-100.,100.};
double pp13_fit0_par2_CB2_pol[] = {1.,-100.,100.};
double pp13_fit0_par3_CB2_pol[] = {1.,-100.,100.};
double pp13_fit0_par4_CB2_pol[] = {1.,-100.,100.};
double pp13_fit0_par5_CB2_pol[] = {1.,-100.,100.};
double pp13_fit0_par6_CB2_pol[] = {0.1,-50.,50.};

double pp13_fit0_par0_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par1_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par2_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par3_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par4_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par5_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit0_par6_variance_CB2_pol[] = {0.,0.,0.};

double pp13_fit0_mean_jpsi_CB2_pol[] = {3.10, 3.07, 3.20};
double pp13_fit0_width_jpsi_CB2_pol[] = {0.07, 0.04, 0.10}; //0.09, 0.05, 0.15
double pp13_fit0_norm_psi2S_CB2_pol[] = {0.1, 0, 0.5};

//Tail: pp13TeV. Fit Range: 2.2 4.4
double pp13_fit1_par0_CB2_pol[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double pp13_fit1_par1_CB2_pol[] = {0.1,-100.,100.};
double pp13_fit1_par2_CB2_pol[] = {0.1,-100.,100.};
double pp13_fit1_par3_CB2_pol[] = {0.1,-100.,100.};
double pp13_fit1_par4_CB2_pol[] = {0.1,-100.,100.};
double pp13_fit1_par5_CB2_pol[] = {0.1,-100.,100.};
double pp13_fit1_par6_CB2_pol[] = {0.1,-5.,5.};

double pp13_fit1_par0_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par1_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par2_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par3_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par4_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par5_variance_CB2_pol[] = {0.,0.,0.};
double pp13_fit1_par6_variance_CB2_pol[] = {0.,0.,0.};

double pp13_fit1_mean_jpsi_CB2_pol[] = {3.10, 3.07, 3.20};
double pp13_fit1_width_jpsi_CB2_pol[] = {0.07, 0.04, 0.10};
double pp13_fit1_norm_psi2S_CB2_pol[] = {0.1, 0, 0.5};


//Tail: geant3. Fit Range: 2.2 4.5
double geant3_fit0_par0_CB2_pol[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_CB2_pol[] = {1.,-10.,10.}; // 1.,-10000.,10000.
double geant3_fit0_par2_CB2_pol[] = {1.,-10.,10.};
double geant3_fit0_par3_CB2_pol[] = {1.,-10.,10.}; // 1.,-10000.,10000.
double geant3_fit0_par4_CB2_pol[] = {1.,-10.,10.};
double geant3_fit0_par5_CB2_pol[] = {1.,-10.,10.};
double geant3_fit0_par6_CB2_pol[] = {0.1,-10.,10.};

double geant3_fit0_par0_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par1_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par2_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par3_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par4_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par5_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit0_par6_variance_CB2_pol[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_CB2_pol[] = {3.10, 3.07, 3.20};
double geant3_fit0_width_jpsi_CB2_pol[] = {0.07, 0.04, 0.10};
double geant3_fit0_norm_psi2S_CB2_pol[] = {0.1, 0, 0.5};

//Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit1_par0_CB2_pol[] = {1, 0.,10.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_CB2_pol[] = {1,-10.,10.};
double geant3_fit1_par2_CB2_pol[] = {1,-10.,10.};
double geant3_fit1_par3_CB2_pol[] = {1,-10.,10.};
double geant3_fit1_par4_CB2_pol[] = {1,-10.,10.};
double geant3_fit1_par5_CB2_pol[] = {1,-10.,10.};
double geant3_fit1_par6_CB2_pol[] = {0.1,-10.,10.};

double geant3_fit1_par0_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par1_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par2_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par3_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par4_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par5_variance_CB2_pol[] = {0.,0.,0.};
double geant3_fit1_par6_variance_CB2_pol[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_CB2_pol[] = {3.10, 3.07, 3.15};
double geant3_fit1_width_jpsi_CB2_pol[] = {0.07, 0.04, 0.09};
double geant3_fit1_norm_psi2S_CB2_pol[] = {0.1, 0, 0.5};

//-----------------------------------------------------------------
//
// CB2 function + double exponential
//
//-----------------------------------------------------------------

double pp13_fit0_par0_CB2_exp[] = {1.,0.,10000.}; // initial value, lower limit, upper limit
double pp13_fit0_par1_CB2_exp[] = {1.,0.,50.};
double pp13_fit0_par2_CB2_exp[] = {1.,0.,1000000.};
double pp13_fit0_par3_CB2_exp[] = {1.,0.,50.};

double pp13_fit0_par0_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit0_par1_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit0_par2_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit0_par3_variance_CB2_exp[] = {0.,0.,0.};

double pp13_fit0_mean_jpsi_CB2_exp[] = {3.10, 3.07, 3.20};
double pp13_fit0_width_jpsi_CB2_exp[] = {0.07, 0.04, 0.10};
double pp13_fit0_norm_psi2S_CB2_exp[] = {0.1, 0.0, 0.5};

double pp13_fit1_par0_CB2_exp[] = {1.,0.,10000.}; // initial value, lower limit, upper limit
double pp13_fit1_par1_CB2_exp[] = {1.,0.,50.};
double pp13_fit1_par2_CB2_exp[] = {1.,0.,10000.};
double pp13_fit1_par3_CB2_exp[] = {1.,0.,50.};

double pp13_fit1_par0_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit1_par1_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit1_par2_variance_CB2_exp[] = {0.,0.,0.};
double pp13_fit1_par3_variance_CB2_exp[] = {0.,0.,0.};

double pp13_fit1_mean_jpsi_CB2_exp[] = {3.10, 3.07, 3.20};
double pp13_fit1_width_jpsi_CB2_exp[] = {0.07, 0.04, 0.10};
double pp13_fit1_norm_psi2S_CB2_exp[] = {0.1, 0.0, 0.5};

//Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit0_par0_CB2_exp[] = {1,0.,10000.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_CB2_exp[] = {1,0.,50.};
double geant3_fit0_par2_CB2_exp[] = {1,0.,10000.};
double geant3_fit0_par3_CB2_exp[] = {1,0.,50.};

double geant3_fit0_par0_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit0_par1_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit0_par2_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit0_par3_variance_CB2_exp[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_CB2_exp[] = {3.10, 3.07, 3.15};
double geant3_fit0_width_jpsi_CB2_exp[] = {0.07, 0.04, 0.09};
double geant3_fit0_norm_psi2S_CB2_exp[] = {0.1, 0, 0.5};

//Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit1_par0_CB2_exp[] = {1, 0.,10000.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_CB2_exp[] = {1, 0.,50.};
double geant3_fit1_par2_CB2_exp[] = {1, 0.,10000.};
double geant3_fit1_par3_CB2_exp[] = {1, 0.,50.};
//double geant3_fit1_par4_CB2_exp[] = {0.1,-100.,100.};
//double geant3_fit1_par5_CB2_e[] = {0.1,-100.,100.};
//double geant3_fit1_par6_CB2_pol[] = {0.1,-50.,50.};

double geant3_fit1_par0_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit1_par1_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit1_par2_variance_CB2_exp[] = {0.,0.,0.};
double geant3_fit1_par3_variance_CB2_exp[] = {0.,0.,0.};
//double geant3_fit1_par4_variance_CB2_pol[] = {0.,0.,0.};
//double geant3_fit1_par5_variance_CB2_pol[] = {0.,0.,0.};
//double geant3_fit1_par6_variance_CB2_pol[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_CB2_exp[] = {3.10, 3.07, 3.15};
double geant3_fit1_width_jpsi_CB2_exp[] = {0.07, 0.04, 0.09};
double geant3_fit1_norm_psi2S_CB2_exp[] = {0.1, 0, 0.5};


//-----------------------------------------------------------------
//
// NA60 function + variable width gaussian
//
//-----------------------------------------------------------------

double geant3_fit0_par0_NA60_VWG[] = {1, 0, 10.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_NA60_VWG[] = {1, -10.,10.}; //2.5,0.,8.
double geant3_fit0_par2_NA60_VWG[] = {1, -10.,10.}; // 1.,-5.,5.
double geant3_fit0_par3_NA60_VWG[] = {1, -10.,10.}; // 1.,-5.,5.
double geant3_fit0_par4_NA60_VWG[] = {1, -10.,10.}; // 1.,-5.,5.

double geant3_fit0_par0_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit0_par1_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit0_par2_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit0_par3_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit0_par4_variance_NA60_VWG[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_NA60_VWG[] = {3.10, 3.07, 3.20};
double geant3_fit0_width_jpsi_NA60_VWG[] = {0.07, 0.04, 0.10};
double geant3_fit0_norm_psi2S_NA60_VWG[] = {0.1, 0, 0.5};


// Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit1_par0_NA60_VWG[] = {1.,0.,10.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_NA60_VWG[] = {0.1,-10.,10.};
double geant3_fit1_par2_NA60_VWG[] = {0.1,-10.,10.};
double geant3_fit1_par3_NA60_VWG[] = {0.1,-10.,10.};
double geant3_fit1_par4_NA60_VWG[] = {0.1,-10.,10.};

double geant3_fit1_par0_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit1_par1_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit1_par2_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit1_par3_variance_NA60_VWG[] = {0.,0.,0.};
double geant3_fit1_par4_variance_NA60_VWG[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_NA60_VWG[] = {3.10, 3.07, 3.20};
double geant3_fit1_width_jpsi_NA60_VWG[] = {0.07, 0.04, 0.10};
double geant3_fit1_norm_psi2S_NA60_VWG[] = {0.1, 0, 0.5};


//-----------------------------------------------------------------
//
// NA60 function + pol2/pol3
//
//-----------------------------------------------------------------


// Tail: geant3. Fit Range: 2.2 4.5
double geant3_fit0_par0_NA60_pol[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_NA60_pol[] = {1.,-10.,10.};
double geant3_fit0_par2_NA60_pol[] = {1.,-10.,10.};
double geant3_fit0_par3_NA60_pol[] = {1.,-10.,10.};
double geant3_fit0_par4_NA60_pol[] = {1.,-10.,10.};
double geant3_fit0_par5_NA60_pol[] = {1.,-10.,10.};
double geant3_fit0_par6_NA60_pol[] = {0.1,-5.,5.};

double geant3_fit0_par0_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par1_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par2_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par3_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par4_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par5_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit0_par6_variance_NA60_pol[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_NA60_pol[] = {3.10, 3.07, 3.20};
double geant3_fit0_width_jpsi_NA60_pol[] = {0.07, 0.04, 0.10};
double geant3_fit0_norm_psi2S_NA60_pol[] = {0.1, 0, 0.5};


// Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit1_par0_NA60_pol[] = {1.,0.,5.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_NA60_pol[] = {1,-10.,10.};
double geant3_fit1_par2_NA60_pol[] = {1,-10.,10.};
double geant3_fit1_par3_NA60_pol[] = {0.1,-10.,10.};
double geant3_fit1_par4_NA60_pol[] = {0.1,-10.,10.};
double geant3_fit1_par5_NA60_pol[] = {0.1,-10.,10.};
double geant3_fit1_par6_NA60_pol[] = {1.,-5.,5.};

double geant3_fit1_par0_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par1_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par2_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par3_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par4_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par5_variance_NA60_pol[] = {0.,0.,0.};
double geant3_fit1_par6_variance_NA60_pol[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_NA60_pol[] = {3.10, 3.07, 3.20};
double geant3_fit1_width_jpsi_NA60_pol[] = {0.07, 0.04, 0.10};
double geant3_fit1_norm_psi2S_NA60_pol[] = {0.1, 0, 0.5};

//-----------------------------------------------------------------
//
// NA60 function + exp
//
//-----------------------------------------------------------------

// Tail: geant3. Fit Range: 2.2 4.5
double geant3_fit0_par0_NA60_exp[] = {1.,0.,10000.}; // initial value, lower limit, upper limit
double geant3_fit0_par1_NA60_exp[] = {1.,0.,50.};
double geant3_fit0_par2_NA60_exp[] = {1.,0.,10000.};
double geant3_fit0_par3_NA60_exp[] = {1.,0.,50.};

double geant3_fit0_par0_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit0_par1_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit0_par2_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit0_par3_variance_NA60_exp[] = {0.,0.,0.};

double geant3_fit0_mean_jpsi_NA60_exp[] = {3.10, 3.07, 3.20};
double geant3_fit0_width_jpsi_NA60_exp[] = {0.07, 0.04, 0.10};
double geant3_fit0_norm_psi2S_NA60_exp[] = {0.1, 0, 0.5};

// Tail: geant3. Fit Range: 2.4 4.7
double geant3_fit1_par0_NA60_exp[] = {1.,0,100000.}; // initial value, lower limit, upper limit
double geant3_fit1_par1_NA60_exp[] = {1,0.,50.};
double geant3_fit1_par2_NA60_exp[] = {1,0.,100000.};
double geant3_fit1_par3_NA60_exp[] = {1,0.,100.};

double geant3_fit1_par0_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit1_par1_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit1_par2_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit1_par3_variance_NA60_exp[] = {0.,0.,0.};
double geant3_fit1_par6_variance_NA60_exp[] = {0.,0.,0.};

double geant3_fit1_mean_jpsi_NA60_exp[] = {3.10, 3.07, 3.20};
double geant3_fit1_width_jpsi_NA60_exp[] = {0.07, 0.04, 0.10};
double geant3_fit1_norm_psi2S_NA60_exp[] = {0.1, 0, 0.5};
