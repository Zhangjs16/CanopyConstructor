/*  Canopy Constructor algorithm v.1.0, published on 29/06/2020:

 1. Background: Inspired by the individual-based forest growth simulator TROLL (https://github.com/TROLL-code/TROLL) as well as canopy filling algorithms (Taubert et al. 2015, PNAS), and programmed to be compatible with TROLL v.2.5.
        - To allow for easy co-development and consistency, the overall style is in keeping with the TROLL programming layout (including a provision of all functions in one single file).
        - Some functions are taken directly from TROLL. However, unlike in TROLL, most tree-based functions have been programmed as non-member functions of the Tree class, since they need to be applicable to both the initialised trees and tests with shifted crowns/trees. This could be potentially transformed into a more object-oriented style by incorporating a class for potential trees in the future
 2. Basic functioning: the script constructs a virtual forest from a number of data sources that are provided at initialization. The core functions are ConstructCanopy() and FitCanopy().
        - ConstructCanopy() function converts empirical data into virtual tree distributions and creates an initial canopy with randomly tree dimensions. It then calls FitCanopy() which optimizes the fit.
        - There are two main options: a) If allometric scaling rules, ground data (geolocated tree trunks) and ALS information (canopy height model) are provided, then the algorithm draws tree crowns from the allometric rules and swaps them until a best fit is achieved. This can be inversed to obtain allometric rules. b) If a crown packing matrix is provided, i.e. a description of within-canopy crown packing in addition to ground data, then the algorithm draws trees until the canopy is well-packed.
 3. Extended functioning (in development): the script also considers leaf physiological information and imports the TROLL photosynthetic module. This introduces a biological viability condition when swapping trees, which should further narrow down the potential configurations. Furthermore, this can be used to initialise TROLL v.2.5 from a reconstructed forest
 4. User input and error handling:
        - the Canopy Constructor automatically outputs all provided input sheets, as they were read into the variables, and generates new ones from default values for those that were not provided. This can be used to ensure that the input was read correctly, or to generate default input sheets for a first run.
        - users are flexible in what parameters they provide for the global and detailed input sheets, but parameter names have to be exact (as in R functions, for example). Parameters that were not provided (or not correctly), will be automatically replaced by default values, with a WARNING attached.
        - for the data files (CHM, inventory, climate, crown packing), this is not the case. Files will be checked for consistency with the parameter files, and if they are not consistent, the Canopy Constructor algorithm is exited. The only exception to this are parameters in the inventory file beyond the minimum requirements (x, y coordinates and dbh), e.g. height, CR etc. These do not have to be provided, or do not have to be correct and will be overwritten by allometric functions etc.
        - errors will always be returned if input files are provided, but not readable
 */

/* overall TODO: optimize tree locating in step 2, physiology module, species identities, repeated ALS acquisitions, event() function that allows to simulate events on reconstructed canopy (selective logging, etc.) */
/* technical TODO: throw error for input parameters outside allowed range, error checks for climate input file, R-interface via Rcpp, MCMC functionality, more object-oriented and fewer global variables */

#define exchangepositions_step2
#define accept_equalfits

# include <cstdio>
# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <limits>
# include <ctime>
# include <cmath>
# include <algorithm>
# include <vector>
# include <sstream>

/* gsl libraries for stochasticity */

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>

using namespace std;

/* user control */
int seed_rng;
int ConstructorStep;

/* Global constants (e.g. PI and various derivatives...) */
# define PI 3.141592654
# define twoPi 6.2831853071
# define Pis2 1.570796327
# define iPi 0.3183099

/* buffer for inputs */
char buffer[512], inputfile_inventory[256], inputfile_climate[256], inputfile_global[256], inputfile_detail[256],inputfile_chm[256], *bufi_inventory(0), *bufi_climate(0), *bufi_global(0), *bufi_detail(0), *bufi_chm(0);
char output[256], *bufo(0);
char inputfile_cp[256], *bufi_cp(0);

/* file output streams */
fstream output_input_global;
fstream output_input_detailed;
fstream output_input_CHM;
fstream output_input_inventory;
fstream output_input_inventorycutoff;
fstream output_input_cp;
fstream output_input_cpcutoff;

fstream output_chm;
fstream output_sumstat;
fstream output_agbquarterha;
fstream output_sdd;

fstream output_vertical;
fstream output_verticalfromtop;
fstream output_vertical_ALS;
fstream output_verticalfromtop_ALS;
fstream output_twodim;
fstream output_field;
fstream output_trees;
fstream output_LAI3D;
fstream output_volume;
fstream output_troll;

/********************************/
/* Parameters of the simulation */
/********************************/

/* global parameters */
int sites,              /* number of sites */
    cols,               /* nb of columns */
    rows,               /* nb of lines */
    height_max,         /* maximum height */
    nb_parametersets;   /* number of sets of detailed parameters to be varied */

int mincol_absolute, minrow_absolute;

int nb_parametersdetailed;      /* nb of parameters in one set for abc routine */
int parameterset;               /* current parameter set*/
int paramIDcurrent;             /* saves the current parameter set ID, can be passed on to the output */

float a_sdd;                    /* parameter for stem diameter distribution, either power law exponent (if <= 0) or shape parameter of Weibull distribution */
float b_sdd;                    /* parameter for stem diameter distribution, either inverse of rate of exponential decline for large diameters (> dbh_transition from power law) or scale parameter of Weibull */
float dbh_transition;           /* diameter where power law transitions into exponential function */
float fraction_drawspowerlaw;   /* if powerlaw distribution is assumed, the fraction of draws within the exponential tale */
float nbavg_sdd;                   /* the average number of trees per ha*/
int nb_sdd;                   /* the total tree number for the plot */
float crown_gap_fraction;   /* this is the fraction of gaps in a crown */
int icrown_gap_fraction;    /* when we loop over a crown, every n steps, a gap is preserved. icrown_gap_fraction is this n */

int trees_fromground;       /* nb of trees read in from ground census */
float dbh_cutoff_fromground;    /* the lowest dbh measured on the ground */
int height_cutoff_fromground;   /* the lowest height measured on the ground */

float mean_height_emp, sd_height_emp;

float dbh_limittop, dbh_limitbottom, dbhlog_binsize;
int dbhlog_binnb;
float dbh_limitstoppingrule, stoppingrule;

int nbsteps_dissimilarity, nbsteps_carbonbalance, nbsteps_mae, nbsteps_combined;

float crowndisplacement_factor;

int flag_ApplyMedianfilter, flag_Prefitting, flag_PreventCrownpiercing, flag_OutputReduced, flag_OutputTROLL, flag_LAIgradient, flag_powerlaw;

vector<float> parameters_detailed;  // detailed parameters for abc inference, columns are the parameters, rows are the parameter sets (paramID, if several sets are provided), access: parameter + parameterset * nb_parametersdetailed

/* random generator */
gsl_rng *gslrand;

/* LAI3D field */
float **LAI3D(0);   /* leaf density (per volume unit) */
int **Voxel3D(0);   /* a voxel field where the number of contributions of different trees are counted */

float **correlation_structure;
int trees_sampled;
int trees_fitted;
int trees_fitted10;
/* derived fields */
vector<int> chm_field;
float **transmittance_simulated(0);
int **ALS_sampling(0);
int **ALS_echos(0);

/* and three more vector fields */
vector<int> chm_field_random;

int sum_abserror;
int sum_abserrormin, sum_abserrormax;
vector<int> hist_CHM_sim;
float dissimilarity_min, dissimilarity_max;

/* now the output summary stats */
int sum_abserror_random, sum_abserror_dist, sum_abserror_spatial, sum_abserror_final;
float dissimilarity_random, dissimilarity_dist, dissimilarity_spatial, dissimilarity_final;
float rate_spatial, rate_dist, rate_final;

float carbonstarv;
float carbonstarv_min, carbonstarv_max;
int sum_abserror_physiology;
float dissimilarity_physiology;
float carbonstarv_random,carbonstarv_spatial, carbonstarv_dist, carbonstarv_physiology, carbonstarv_final;
float rate_physiology;

float dbh_limitcarbon;   /* only trees above this threshold will actually be checked for carbon starvation */
int nbtrees_carbonstarv;
vector<float> maxCR_perheight;
float heightminimum_threshold;

int nbtrees_total;
int nbtrees_total_threshold;
int nbingrowth_random, nbingrowth_final;
int nbmaxovl_random, nbmaxovl_final;
int nbmaxovl10_random, nbmaxovl10_final;
int CA_exterior, CA_full;
int CA_exterior_random, CA_exterior_final;
int CA_full_random, CA_full_final;

vector<int> trees_dbhsorted;
vector<int> trees_dbhlogbins;
vector<int> trees_sitesnspecies_empirical;
vector<float> trees_traits_empirical;
vector<float> trees_draw;
vector<float> CP_twodimemp_params;
vector<float> CP_twodimemp;

/* these are two vectors that can store values associated with trees (up to int limits) and sort them */
int **trees_sortedheight(0);
int *chm_empirical(0);
int *hist_CHM_emp(0);
float *hist_CHM_empcompcumulative(0);

int iteration_final;

int nbattempt;
int nbsuccess;

/* LookUp_Crown_site */
int LookUp_Crown_site[2601];               /* lookup table to fill crown cylinder sequentially from inside to outside, allowing for smooth crown growth */
int LookUp_neighbour_site[10201];            /* lookup table to search for neighbours in the tree's surroundings that might be affected by additional LAI */

/* allometries for height, crown radius and crown depth, as well as leaf and wood traits */
float a_height, b_height, sigma_height;
float a_CR, b_CR, sigma_CR;
float a_CD, b_CD, sigma_CD;
float mean_N, mean_P, mean_LMA, mean_wsg;
float sigma_N, sigma_P, sigma_LMA, sigma_wsg;
float max_LAI, sigma_LAI;
float corr_CR_height, corr_N_P, corr_N_LMA, corr_P_LMA,
cov_N_P, cov_N_LMA, cov_P_LMA;

int covariance_status;      /* if one of N, P, or LMA has 0 variation, the Cholesky decomposition fails, we then use no correlation at all */
gsl_matrix *mcov_N_P_LMA;   /* covariance matrices for crown radius - tree height covariance and leaf_properties, respectively */
gsl_vector *mu_N_P_LMA, *variation_N_P_LMA; /* the mean values of the distributions as well as the result vector for the multivariate draw */

float hetscedast_height;
float hetscedast_CR;
float *sigma_heightbin;
float *sigma_CRbin;

float maxoverlap_percentage = 1.0;
float shape_crown;
float timestep = 1.0;   /* transposing concept from TROLL model, here we calculate all quantities over a year */

/*********************************************/
/* Environmental variables of the simulation */
/*********************************************/

int iterperyear = 12;   /* number of values for climate variables, imported from TROLL model */

/* Climate input data  */
/* these climate input data are given in the input file, its structure depends on the timestep and scenario used for the simulation */

float daily_light[24];    /* normalized (ie between 0 and 1) daily light variation (used if DAILYLIGHT defined) */
float daily_vpd[24];      /* normalized (ie between 0 and 1) daily vpd variation (used if DAILYLIGHT defined) */
float daily_T[24];        /* normalized (ie between 0 and 1) daily T variation (used if DAILYLIGHTdefined) */

float *Temperature(0);                      /* in degree Celsius */
float *DailyMaxTemperature(0);              /* in degree Celsius */
float *NightTemperature(0);                 /* in degree Celsius */
float *Rainfall(0);                         /* in mm */
float *WindSpeed(0);                        /* in m/s */
float *MaxIrradiance(0);                    /* in W/m2 */
float *MeanIrradiance(0);                   /* in W/m2 */
float *SaturatedVapourPressure(0);          /* in kPa */
float *VapourPressure(0);                   /* in kPa */
float *VapourPressureDeficit(0);            /* in kPa */
float *DailyVapourPressureDeficit(0);       /* in kPa */
float *DailyMaxVapourPressureDeficit(0);    /* in kPa */

float Wmax_avg,     /* avg Wmax across the year */
Tmax_avg,           /* avg Tmax across the year */
VPDmax_avg,         /* avg VPDmax across the year */
Tnight_avg,         /* avg tnight across the year */
temp_avg;           /* avg temp across the year */

float klight,       /* original light absorption rate or extinction cefficient used in Beer-Lambert law (leaf angle distribution) */
kpar,               /* effective light absorption rate or extinction cefficient used in Beer-Lambert law to compute light within the canopy */
phi,                /* apparent quantum yield (in micromol C/micromol photon). phi is the quantum yield multiplied by leaf absorbance (in the literature, the quantum yield is often given per absorbed light flux, so one should multiply incident PPFD by leaf absorbance, see Poorter et al American Journal of Botany (assumed to be 0.91 for tropical tree species). Even though holding phi constant across species is widely assumed, in reality, phi varies across species and environmental conditions (see eg Domingues et al 2014 Plant Ecology & Diversity) */
theta,              /* parameter of the Farquhar model set to 0.7 in this version. For some authors, it should be species-dependent or temperature dependent, but these options are not implemented here */
g1,                 /* g1 parameter of Medlyn et al's model of stomatal conductance */
alpha,              /* apparent quantum yield to electron transport in mol e-/mol photons, equal to the true quantum yield multiplied by the leaf absorbance */
Cair,               /* atmosphericCO2 concentration, if we aim at making CO2 vary (scenarios), CO2 will have to have the same status as other climatic variables  */
iCair;              /* inverse of Cair */

int nbTbins;                    /*nb of bins for the temperature lookup tables */
float iTaccuracy;                /* inverse of accuracy of a temperature bin (e.g. Taccuracy us 0.1 or 0.5 °C, so iTaccuracy is 10.0 or 2.0, resp) */
float *LookUp_KmT(0);                   /* lookup table for fast comput of Farquhar as a function of T */
/* !! leaf temperature must be comprised between 0°C and 60°C
 (T_leaf is stored every 0.5°C, so 120 values in total */
float *LookUp_GammaT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_VcmaxT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_JmaxT(0);                   /* lookup table for fast comput of Farquhar */
float *LookUp_Rday(0);                 /* lookup table for fast comput of Farquhar */
float *LookUp_flux_absorption(0);      /* lookup table for faster computation of PPFD / new in v.2.4: absorption flux */
float *LookUp_flux(0);                 /* lookup table for faster computation of PPFD / new in v.2.4: averaging instead of top value*/
float *LookUp_VPD(0);                  /* lookup table for faster computation of VPD / new in v.2.4:averaging instead of top value */
float *LookUp_T(0);                    /* lookup table for faster computation of T / new in v.2.4: averaging instead of top value */
float *LookUp_Rstem(0);                /* lookup table for faster computation of Rstem */
float *LookUp_Rnight(0);                /* lookup table for faster computation of Rstem */

/**********************************/
/* Various mathematical functions */
/**********************************/

float flor(float f) {
    if(f>0.) return f;
    else return 0.;
}
float florif(int i) {
    if(i>0) return float(i);
    else return 0.;
}
float maxf(float f1, float f2) {
    if(f1>f2) return f1;
    else return f2;
}
float minf(float f1, float f2) {
    if(f1<f2) return f1;
    else return f2;
}
int min(int i1, int i2) {
    if(i1<i2) return i1;
    else return i2;
}
int max(int i1, int i2) {
    if(i1>i2) return i1;
    else return i2;
}
int sgn(float f) {
    if(f>0.0) return 1;
    else return -1;
}

/*#########################*/
/*####### FUNCTIONS #######*/
/*#########################*/

/* Initialisation */

/* Initialisation of parameter files and output streams */
/* Error handling: if no or incomplete files are provided, default values are assumed. This allows simplified tests of the Canopy Constructor with reduced sets of parameters */
void InitialiseGlobal(int flag_global, int &error_global);
void InitialiseDetail(int flag_detail, int &error_detail);
void InitialiseOutput(int flag_output, int &error_output);

/* Initialisation of data files */
/* Error handling: if files are not provided, default options are assumed. If files are provided, but do not correspond to expected format, then an error is thrown and the Canopy Constructor algorithm exited. This is to ensure internal consistency. */
void InitialiseCHM(int flag_chm, int &error_chm);
void InitialiseInventory(int flag_inventory, int &error_inventory);
void InitialiseClimate(int flag_climate, int &error_climate);
void InitialiseCP(int flag_cp, int &error_cp);

/* main routines */
void ConstructCanopy();
void ConstructInitialCanopy();
void FitCanopy(int iteration_minimum, int nbsteps, int dissimilarity_flag, int mae_flag, int carbonstarv_flag, vector<int> &sites_random);

void SupplementTraits(float dbh, int species_label, float &height, float &CR, float &CD, float &Nmass, float &Pmass, float &LMA, float &wsg);

void FindBetterCrown(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstav_flag);
void FindBetterPosition(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstarv_flag);
#ifdef exchangepositions_step2
void FindBetterPositions(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstarv_flag);
#endif


int TreeIngrowth(int site_crown, float height_tree, float CR_tree);
int TreeIngrowthFromabove(int site_crown, float height_tree, float CR_tree);

/* functions for updating fields */
void GetCanopyHeight(vector<int> &chm_patch, int site_crown, int extent);
void SetCanopyHeight(vector<int> chm_patch, int site_crown, int extent);


void UpdateCanopyAroundTree(int site_crown, int site_crown_updated, float height_tree, float height_tree_updated, float CR_tree, float CR_tree_updated,float CD_tree, float CD_tree_updated, int extent_tree);
void UpdateCanopyNewPosition(int site_crown, int site_crown_updated, float height_tree, float CR_tree, float CD_tree, int extent_tree);
int GetTreeExtent(int site_crown, int site_crown_updated, float CR_tree, float CR_tree_updated);
int EvaluateCanopy(int mae_flag, int dissimilarity_flag, int carbonstarv_flag, float dissimilarityprevious, float dissimilaritynew, int abserror_previous, int abserror_new, float meanheight_previous, float meanheight_new, float CA_exterior_previous, float CA_exterior_new, float carbonstarv_previous, float carbonstarv_new);

void CanopyDistanceFull(float &dissimilarity, int &sum_abserror);

void CanopyDevMetrics(int site_crown, int extent, int remove);
float CalcDissimilarity(int hbin_size);
float CanopyMeanheight();

/* two functions that calculate the LAI contribution of a tree that has not yet been accepted, the remove flag at 0 adds LAI, and removes LAI at 1 */
void CalcLAItrial(bool remove, int site_tree, float height, float CR, float CD, float LAI, float dbh);
void FillVoxel(int increment, int site_crown, float height, float CR, float CD);
float GetVoxeloverlap(int site_tree, float height, float CR, float CD, bool multiple);
int CalcCrownAreaExterior(int site_crown, float CR);
int CalcCV(int site_tree, float height, float CR, float CD);
int CalcCV_generalized(float height, float CR, float CD);

/* three functions to a) get a list of trees affected, b) to update the leaves in the crown and c) to reverse the update */
/* all three functions are overloaded for three different cases: crown swapping, crown shifting/updating, tree shifting */
int GetNbTreesThreshold(vector<int> &list_trees_potential, vector<int> &list_trees_examined, vector<float> &heights_list_trees_potential, int site_tree, float dbh, float height, float CR, int col_crowncenter, int row_crowncenter, bool addtree, int &functioncalls);
void GetTreeList(int site_tree1, int site_tree2, float CR1_previous, float CR1_new, float CR2_previous, float CR2_new, float CD1_previous, float CD1_new, float CD2_previous, float CD2_new, float height1_previous, float height1_new, float height2_previous, float height2_new, int crown_displacement1_previous, int crown_displacement1_new, int crown_displacement2_previous, int crown_displacement2_new, vector<int> &list_trees);
void GetTreeList(int site_tree, float CR_previous, float CR_new, float CD_previous, float CD_new, float height_previous, float height_new, float crown_displacement_previous, float crown_displacement_new, vector<int> &list_trees);
void GetTreeList(int site_original, int site_new, vector<int> &list_trees);

void UpdateLeaves();

void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site1, float height1_new, float CR1_new, float CD1_new, int crown_displacement1_new, float Nmass1_new, float Pmass1_new, float LMA1_new, int site2, float height2_new, float CR2_new, float CD2_new, int crown_displacement2_new, float Nmass2_new, float Pmass2_new, float LMA2_new);
void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site_tree, float height_new, float CR_new, float CD_new, int crown_displacement_new);
void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site_original, int site_new);

void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site_tree1, float height1_new, float CR1_new, float CD1_new, int crown_displacement1_new, int site_tree2, float height2_new, float CR2_new, float CD2_new, int crown_displacement2_new);
void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site, float height_new, float CR_new, float CD_new, int crown_displacement_new);
void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site_original, int site_new);


void GetCVcumulatedemp(vector<int> &CV_cumulatedemp, int &CV_totalemp);
void SaveCHM(vector<int> &chm_field_safe);
void CalcSumstatCanopy(int &sum_abserror_safe, float &dissimilarity_safe, float &carbonstarv_safe, int &nbingrowth_safe, int &nbmaxovl_safe, int &nbmaxovl10_safe, int &CA_exterior_safe, int &CA_full_safe);
void UpdateArrays();

int CalcCrownDisplacement(int site_tree, int height);
void UpdateCorrelation();
void UpdateCHM();
void UpdateCHMradius(int site_tree, float crown_radius);
void AllocMem();
void FreeMem();

int PlaceTreeWithincanopy_avg(float height_test, float CR_test, float CD_test, vector<int> &sites_free);
void CalcCV_layers(vector<int> &CVtree_layers, int site_crown, float height, float CR, float CD);

float GetAvgHeight(int site_crown, float height, float CR, float CD);
void GetCrownHeights(vector<int> &crown_heights, float height, float CR, float CD);
float UpdateLAIabove_eff(int site_tree,float height, float CR, float CD, float dbh);
void CheckAllometricScaling();
void OverallCanopyStatistics();
void CreateFullCPmatrix(vector<float> &CP_twodim_fullextent, vector<float> &CP_twodim, vector<int> &voxtotal_incan_twodim, int nbsamples_minimum, int height_canopy_max);

float GPPleaf(float PPFD, float VPD, float T,float Vcmax, float Jmax);
float Rdayleaf(float T, float Rdark);
float dailyRdayleaf(float T, float Rdark);
float dailyGPPleaf(float PPFD, float VPD, float T, float Vcmax, float Jmax);
float treeLAImax(float Vcmax, float Jmax, float Rdark);

/* Postprocessing */
float GetCrownExposure(bool canopy, int site_crown, float height, float CR, float CD, float radiusfactor);
float GetCanopyStatus(int site_crown, float height, float CR, float CD);

/* Output functions */
void OutputSumstatEmp();
void OutputSumstatFinal();
void OutputAGBtrees_finescale();
void OutputAGBtrees_distributions();
void OutputInputCHMinventory();
void OutputTrees();
void OutputTroll();

void CalcHistCA(vector<int> &hist_CA);
void CalcFieldCA(vector<int> &field_CA);
void CalcFieldAGB(vector<float> &field_AGB);


/**********************/
/* Treeshape calculations */
/**********************/

void GetPPFDabove(int height, int site, float noinput, float (&PPFD)[2]);
void AddupVolume(int height, int site, int voxel_volume, int &crown_volume);
void AddupVolumeLayers(int height, int site, int voxel_volume, vector<int> &crown_volume_layers);
void UpdateVoxelValue(int height, int site, int increment, int &increment_sum);
void GetVoxelValue(int height, int site, int noinput, int (&voxels)[2]);
void GetVoxelValue_countmultiple(int height, int site, int noinput, int (&voxels)[2]);
void GetCanopyEnvironment(int height, int site, float dens, float (&canopy_environment_cumulated)[4]);
void UpdateLAI3D(int height, int site, float dens, float &LA_cumulated);
void ModifyNoModification_float(float input, float &input_unmodified, float CD, float height, int layer_fromtop);
void ModifyNoModification_int(int input, int &input_unmodified, float CD, float height, int layer_fromtop);
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop);
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens);
void GetDensityUniform(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens);
int GetCrownIntarea(float radius);
float GetCrownAreaFilled(float crown_radius, float fraction_filled_target);
float GetRadiusSlope(float CR, float crown_extent, float crownshell_extent_toplayer);
float GetRadiusCylinder(float CR, float crown_extent, float crownshell_extent_toplayer );
template <typename I, typename O, typename M, typename F>
void LoopLayerUpdateCrownStatistic_template(int row_center, int col_center, float height, float CR, float CD, float fraction_filled_target, int shell_fromtop, int layers_filled, float GetRadiusLayer(float, float, float), I CrownStatistic_input, O &CrownStatistic_output, M ModifyCrownStatistic_input, F UpdateCrownStatistic_output);
template <typename I, typename O, typename F>
void CircleAreaUpdateCrownStatistic_template(int row_center, int col_center, int pos_start, int pos_end, float fraction_filled_target, float &fraction_filled_actual, int height_layer, I CrownStatistic_input, O &CrownStatistic_output, F UpdateCrownStatistic);


/* Tree class */
class Tree {
    
public:
    int   t_site,            /* location */
            t_age;           /* analogy to TROLL where this indicates whether tree is alive or not */
    
    int     t_CrownVolume,
    t_CrownDisplacement;
    
    int   t_species_label;
    
    float t_dbh,
            t_Tree_Height,
            t_Crown_Radius,
            t_Crown_Depth,
            t_LAI,
            t_LAImax,
            t_LAIabove,
            t_Pmass,
            t_Nmass,
            t_LMA,
            t_wsg,
            t_dens_top,
            t_dens_belowtop,
            t_dens,
            t_overlap,
            t_overlap_multiple;
    
    float   t_dev_height,
            t_dev_CR,
            t_dev_CD;
    
    float   t_GPP,
            t_NPP;
    
    Tree(){
        t_age = 0;
        t_dbh = t_Tree_Height = t_Crown_Radius = t_Crown_Depth = t_LAI = t_LAIabove = t_dens = 0.0;
        t_dens_top = t_dens_belowtop = 0.0;
    };  /* constructor */
    
    virtual ~Tree() {
    };	/* destructor */
    
    void SetZero(int site);
    void InitialiseTree(int site, float dbh, float height, float CR, float CD, float Nmass, float Pmass, float LMA, float wsg, int species_label);
    void CrownVolume();
    void UpdateDensity();
};
                        
Tree *T=NULL;


/*############################################
 ############################################
 ############     MAIN PROGRAM    ###########
 ############################################
 ############################################*/

int main(int argc,char *argv[]) {

    int flag_global = 0, flag_detail = 0, flag_output = 0, flag_inventory = 0, flag_chm = 0, flag_climate = 0,  flag_cp = 0;
    
    for(int argn=1;argn<argc;argn++){ /* Arguments of the input and output files */
        if(*argv[argn] == '-'){
            switch(*(argv[argn]+1)){
                case 'o':
                       bufo = argv[argn]+2;
                       flag_output = 1;
                       break;
                case 'g':
                    bufi_global = argv[argn]+2;
                    flag_global = 1;
                    break;
                case 'd':
                    bufi_detail = argv[argn]+2;
                    flag_detail = 1;
                    break;
                case 'f':
                    bufi_inventory = argv[argn]+2;
                    flag_inventory = 1;
                    break;
                case 'c':
                    bufi_chm = argv[argn]+2;
                    flag_chm = 1;
                    break;
                case 'm':
                    bufi_climate = argv[argn]+2;
                    flag_climate = 1;
                    break;
                case 'v':
                    bufi_cp = argv[argn]+2;
                    flag_cp = 1;
                    break;
            }
        }
    }
    
    cout << "\n<--------------------------------------------------------------------------------------------------------------->";
    cout << "\n<------------------------------------- Starting Canopy Constructor v.1.0 --------------------------------------->";
    cout << "\n<--------------------------------------------------------------------------------------------------------------->" << endl;
    
    cout << "\n#######################################";
    cout << "\n##### Initialising output streams #####";
    cout << "\n#######################################" << endl;
    int error_output;
    InitialiseOutput(flag_output, error_output);
  
    cout << "\n##########################################";
    cout << "\n##### Initialising global parameters #####";
    cout << "\n##########################################" << endl;
    int error_global;
    InitialiseGlobal(flag_global, error_global);

    cout << "\n############################################";
    cout << "\n##### Initialising detailed parameters #####";
    cout << "\n############################################" << endl;
    int error_detail;
    InitialiseDetail(flag_detail, error_detail);
    
    cout << "\n#############################";
    cout << "\n##### Allocating memory #####";
    cout << "\n#############################" << endl;
    AllocMem();

    cout << "\n##########################################";
    cout << "\n##### Initialising empirical canopy  #####";
    cout << "\n##########################################" << endl;
    int error_chm;
    InitialiseCHM(flag_chm, error_chm);
    
    cout << "\n############################$$###########";
    cout << "\n##### Initialising field inventory  #####";
    cout << "\n############################$$###########" << endl;
    int error_inventory;
    InitialiseInventory(flag_inventory, error_inventory);
    
    cout << "\n###########################################";
    cout << "\n##### Initialising climate parameters #####";
    cout << "\n###########################################" << endl;
    int error_climate;
    InitialiseClimate(flag_climate, error_climate);

    cout << "\n############################################";
    cout << "\n##### Initialising crown packing stats #####";
    cout << "\n############################################" << endl;
    int error_cp;
    InitialiseCP(flag_cp, error_cp);
    
    int initialisation_success;
    
    if(error_output < 2 && error_global < 2 && error_detail < 2 && error_chm < 2 && error_inventory < 2 && error_climate < 2 && error_cp < 2) {
        initialisation_success = 1;
        
        cout << "\n<--------------------------------------------------------------------------------------------------------------->";
        cout << "\n<----------------------------------------- INITIALISATION SUCCESSFUL ------------------------------------------->";
        cout << "\n<--------------------------------------------------------------------------------------------------------------->" << endl;
        
        string ConstructorStep_explanation;
        if(error_chm == 1){
            ConstructorStep = 0;
            flag_Prefitting = 0;
            ConstructorStep_explanation = "Creating random canopy, without fitting to CHM.";
        } else {
            if(error_cp == 1){
                ConstructorStep = 1;
                ConstructorStep_explanation = "Fitting canopy to field inventory and CHM.";
            } else {
                ConstructorStep = 2;
                ConstructorStep_explanation = "Inferring field inventory from crown packing. Fitting canopy to CHM.";
            }
        }
        cout << "Constructor step: " << ConstructorStep << ". " << ConstructorStep_explanation << endl;
        
        cout << "\n####################################";
        cout << "\n##### Starting new simulations #####";
        cout << "\n####################################" << endl;

        /* Initialising random number generator (gsl) */
        const gsl_rng_type *Trandgsl;
        gsl_rng_env_setup();
        Trandgsl = gsl_rng_default;
        gslrand = gsl_rng_alloc (Trandgsl);

        gsl_rng_set(gslrand, seed_rng);

        OutputSumstatEmp();

        double start_time;
        double stop_time = clock();
        double duration = 0.0;

        for(parameterset = 0; parameterset < nb_parametersets; parameterset++){
            start_time = stop_time;
            
            cout << "\n#############################";
            cout << "\n##### New parameter set #####";
            cout << "\n#############################" << endl;
            cout << "\nParameter set: " << parameterset << endl;

            /* reset the Tree class */
            for(int site = 0; site < sites; site++){
                T[site].SetZero(site);
            }
            
            /* reset the chm */
            for(int s = 0; s < sites; s++){
                chm_field[s] = 0;
            }
            
            /* reset the LAI3D field */
            for(int h=height_max;h>=0;h--){
                for(int site=0;site<sites;site++){
                    LAI3D[h][site] = 0.0;
                }
            }
            
            /* reset the Voxel field */
            for(int h=height_max;h>=0;h--){
                for(int site=0;site<sites;site++){
                    Voxel3D[h][site] = 0;
                }
            }
            
            /* retrieve parameter set ID */
            paramIDcurrent = int(lround(parameters_detailed[0 + parameterset * nb_parametersdetailed]));
            
            /* retrieve stem diameter distribution details */
            a_sdd = parameters_detailed[1 + parameterset * nb_parametersdetailed];
            b_sdd = parameters_detailed[2 + parameterset * nb_parametersdetailed];
            
            /* determine which distribution is used to draw diameters */
            if(a_sdd <= 0){
                flag_powerlaw = 1;
                /* for clarity's sake, we ranem a_sdd and b_sdd */
                float slope_pl = a_sdd;
                float invrate_exp = b_sdd;
                
                dbh_transition = slope_pl/invrate_exp;
                
                /* we now calculate the fraction of trees that are drawn from the exponential tale of the power law distribution */
                /* we simply calculate the area under the power law distribution, compared to the area under the exponential distribution */
                float density_transition = exp(invrate_exp * dbh_transition); // or exp(a_sdd)
                float alpha_pl = density_transition/pow(dbh_transition,slope_pl); // we normalize the power law through the continuity condition at dbh_transition
                float density_cumulated_exponential = (exp(invrate_exp * dbh_limittop) - exp(invrate_exp * dbh_transition))/invrate_exp;
                float density_cumulated_powerlaw = (pow(dbh_transition,slope_pl + 1.0) - pow(dbh_limitbottom, slope_pl + 1.0)) * alpha_pl/ (slope_pl + 1.0);
                fraction_drawspowerlaw = density_cumulated_powerlaw/(density_cumulated_exponential + density_cumulated_powerlaw);
                
                //cout << "Slope_pl: " << slope_pl << " invrate_exp: " << invrate_exp << " dbh_transition: " << dbh_transition << " density_transition: " << density_transition << " alpha_pl: " << alpha_pl << " fraction_drawspowerlaw: " << fraction_drawspowerlaw << endl;
            }
            
            nbavg_sdd = parameters_detailed[3 + parameterset * nb_parametersdetailed];
            nb_sdd = int(nbavg_sdd * float(sites)/float(10000));
            
            /* retrieve parameters */
            a_height = parameters_detailed[4 + parameterset * nb_parametersdetailed];
            b_height = parameters_detailed[5 + parameterset * nb_parametersdetailed];
            sigma_height = parameters_detailed[6 + parameterset * nb_parametersdetailed];
            float hetscedast_height_factor = parameters_detailed[7 + parameterset * nb_parametersdetailed];
            a_CR = parameters_detailed[8 + parameterset * nb_parametersdetailed];
            b_CR = parameters_detailed[9 + parameterset * nb_parametersdetailed];
            sigma_CR = parameters_detailed[10 + parameterset * nb_parametersdetailed];
            float hetscedast_CR_factor = parameters_detailed[11 + parameterset * nb_parametersdetailed];
            a_CD = parameters_detailed[12 + parameterset * nb_parametersdetailed];
            b_CD = parameters_detailed[13 + parameterset * nb_parametersdetailed];
            sigma_CD = parameters_detailed[14 + parameterset * nb_parametersdetailed];
            float hetscedast_CD_factor = parameters_detailed[15 + parameterset * nb_parametersdetailed]; // not used at the moment
            mean_LMA = parameters_detailed[16 + parameterset * nb_parametersdetailed];
            sigma_LMA = parameters_detailed[17 + parameterset * nb_parametersdetailed];
            mean_N = parameters_detailed[18 + parameterset * nb_parametersdetailed];
            sigma_N = parameters_detailed[19 + parameterset * nb_parametersdetailed];
            mean_P = parameters_detailed[20 + parameterset * nb_parametersdetailed];
            sigma_P = parameters_detailed[21 + parameterset * nb_parametersdetailed];
            mean_wsg = parameters_detailed[22 + parameterset * nb_parametersdetailed];
            sigma_wsg = parameters_detailed[23 + parameterset * nb_parametersdetailed];
            corr_N_P = parameters_detailed[24 + parameterset * nb_parametersdetailed];
            corr_N_LMA = parameters_detailed[25 + parameterset * nb_parametersdetailed];
            corr_P_LMA = parameters_detailed[26 + parameterset * nb_parametersdetailed];
            shape_crown = parameters_detailed[27 + parameterset * nb_parametersdetailed];
            crown_gap_fraction = parameters_detailed[28 + parameterset * nb_parametersdetailed];

            hetscedast_height = (hetscedast_height_factor - 1.0) * sigma_height;
            hetscedast_CR = (hetscedast_CR_factor - 1.0) * sigma_CR;

            // convert correlations to covariances
            cov_N_P = corr_N_P * sigma_N * sigma_P;
            cov_N_LMA = corr_N_LMA * sigma_N * sigma_LMA;
            cov_P_LMA = corr_P_LMA * sigma_P * sigma_LMA;
            
            if(cov_N_P == 0.0 || cov_N_LMA == 0.0 || cov_P_LMA == 0.0){
                cerr << "\nCovariance matrix N,P,LMA could not be decomposed. Using uncorrelated variation of trait values instead" << endl;
                covariance_status = 0;
            } else {
                cout << "Correlation status. corr_N_P: " << corr_N_P << " cov_N_LMA: " << corr_N_LMA << " corr_P_LMA: " << corr_P_LMA << endl;
                covariance_status = 1;
                // Initialise covariance matrix for N, P, LMA
                
                mcov_N_P_LMA = gsl_matrix_alloc(3,3);
                gsl_matrix_set(mcov_N_P_LMA,0,0,sigma_N*sigma_N);
                gsl_matrix_set(mcov_N_P_LMA,0,1,cov_N_P);
                gsl_matrix_set(mcov_N_P_LMA,0,2,cov_N_LMA);
                gsl_matrix_set(mcov_N_P_LMA,1,0,cov_N_P);
                gsl_matrix_set(mcov_N_P_LMA,1,1,sigma_P*sigma_P);
                gsl_matrix_set(mcov_N_P_LMA,1,2,cov_P_LMA);
                gsl_matrix_set(mcov_N_P_LMA,2,0,cov_N_LMA);
                gsl_matrix_set(mcov_N_P_LMA,2,1,cov_P_LMA);
                gsl_matrix_set(mcov_N_P_LMA,2,2,sigma_LMA*sigma_LMA);
                
                
                // Cholesky decomposition for multivariate draw
                
                cout << "\nCovariance matrix N,P,LMA: " << endl;
                for(int mrow=0; mrow<3;mrow++){
                    cout << gsl_matrix_get(mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 2) << endl;
                }
                
                gsl_linalg_cholesky_decomp1(mcov_N_P_LMA);
                
                cout << "\nCovariance matrix N,P,LMA (after Cholesky decomposition) " << endl;
                for(int mrow=0; mrow<3;mrow++){
                    cout << gsl_matrix_get(mcov_N_P_LMA, mrow, 0) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 1) << "\t" << gsl_matrix_get(mcov_N_P_LMA, mrow, 2) << endl;
                }
                
                // allocating mean and result vectors for multivariate draw
                mu_N_P_LMA = gsl_vector_alloc(3);
                for(int j=0;j<3;j++){
                    gsl_vector_set(mu_N_P_LMA,j,0.0);
                }
                variation_N_P_LMA = gsl_vector_alloc(3);
            }

            
            /* the input parameters for height and CR heteroscedasticity are given as linear decrease or increase in variation on log10 coordinates */
            /* i.e. a slope of -0.05 would mean a decrease of 0.05 for increase of 1 in log(dbh) */
            /* to convert into parameters per bin, we simply divide by 1.0/dbhlog_binsize, or equally: multiply by dbhlog_binsize */
            hetscedast_height *= dbhlog_binsize;
            hetscedast_CR *= dbhlog_binsize;
            
            cout << endl;
            cout << "Hetscedas_height: " << hetscedast_height << endl;
            cout << "Hetscedas_CR: " << hetscedast_CR << endl;

            /* additional computations for crown gap fraction */
            crown_gap_fraction = maxf(crown_gap_fraction,0.000001);     // crown_gap_fraction is prevented from becoming zero in order to avoid division by zero. Given that crown area is currently limited to 1963, the lowest crown_gap_fraction that could potentially have an effect would be 1/1963, which is ~ 0.0005
            crown_gap_fraction = minf(crown_gap_fraction,0.499999);     // we also limit to a max gap fraction of 0.5
            icrown_gap_fraction = int(0.01+1.0/crown_gap_fraction);     // we add 0.01 before division to make sure that we do not encounter rounding problems

            for(int bin=0;bin<dbhlog_binnb;bin++){
                sigma_heightbin[bin] = maxf(sigma_height + hetscedast_height * float(bin), 0.0);
                //cout << "Hetscedas_height: " << sigma_heightbin[bin] << endl;
            }
            
            for(int bin=0;bin<dbhlog_binnb;bin++){
                sigma_CRbin[bin] = maxf(sigma_CR + hetscedast_CR * float(bin),0.0);
                //cout << "Hetscedas_CR: " << sigma_CRbin[bin] << endl;
            }
            //int bin_cutoff_fromground = int(log10(dbh_cutoff_fromground*100.0)/dbhlog_binsize);
            //float sigma_heightbin_cutoff_fromground = sigma_heightbin[bin_cutoff_fromground];
            height_cutoff_fromground = 0;
            cout << "DBH_cutoff at: " << dbh_cutoff_fromground << " | Height_cutoff at: " << height_cutoff_fromground << endl;
            
            cout << "Setting crown packing matrix" << endl;
            int height_canopy_max = height_max + 1;
            for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
                for(int height_within_canopy = 0; height_within_canopy < height_canopy_max; height_within_canopy++){
                    CP_twodimemp[height_canopy + height_within_canopy * height_canopy_max] = CP_twodimemp_params[parameterset + nb_parametersets * (height_canopy + height_within_canopy * height_canopy_max)];
                
                }
            }
            
            /*********************/
            /* CONSTRUCT CANOPY */
            /*********************/
            
            ConstructCanopy();

            /******************/
            /* POSTPROCESSING */
            /******************/
            
            CheckAllometricScaling();
            
            /* update the geometries */
            for(int s=0; s < sites; s++){
                if(T[s].t_age > 0){
                    int site_crown = s + T[s].t_CrownDisplacement;
                    
                    T[s].CrownVolume();
                    T[s].t_overlap = GetVoxeloverlap(site_crown, T[s].t_Tree_Height, T[s].t_Crown_Radius, T[s].t_Crown_Depth, 0);
                    T[s].t_overlap_multiple = GetVoxeloverlap(site_crown, T[s].t_Tree_Height, T[s].t_Crown_Radius, T[s].t_Crown_Depth, 1);
                }
            }
            
            OutputSumstatFinal();
            OutputAGBtrees_finescale();
            OutputAGBtrees_distributions();
            
            stop_time = clock();
            duration += stop_time-start_time;
            
            float durf = float(duration/double(CLOCKS_PER_SEC));        /* output of the effective CPU time */
            
            float percentage = 100.0 * float(parameterset)/float(nb_parametersets);
            if(int(percentage) > 0 && int(percentage)%5 == 0) cout << "Percentage completion: " << percentage << " Time elapsed: " << durf << " seconds." << endl;
            
        }

        float durf = float(duration/double(CLOCKS_PER_SEC));        /* output of the effective CPU time */

        cout << "End of simulation. Time: " << durf << " seconds" << endl;
        
        cout << "\n<--------------------------------------------------------------------------------------------------------------->";
        cout << "\n<---------------------------------------------- POSTPROCESSING ------------------------------------------------->";
        cout << "\n<--------------------------------------------------------------------------------------------------------------->" << endl;
        
        output_input_global.close();
        output_input_detailed.close();
        
        output_agbquarterha.close();
        output_sumstat.close();
        output_chm.close();
        output_sdd.close();

        if(flag_OutputReduced == 0){
            OverallCanopyStatistics();
            OutputInputCHMinventory();
            OutputTrees();
            
            output_input_CHM.close();
            output_input_inventory.close();
            output_input_inventorycutoff.close();
            output_input_cp.close();
            output_input_cpcutoff.close();
            output_field.close();
            output_LAI3D.close();
            output_vertical.close();
            output_verticalfromtop.close();
            output_vertical_ALS.close();
            output_verticalfromtop_ALS.close();
            output_troll.close();
            output_trees.close();
            output_twodim.close();
            output_volume.close();
        }
        
        if(flag_OutputTROLL == 1){
            OutputTroll();
        }
        
        cout << "Output streams closed." << endl;
        
    } else {
        initialisation_success = 0;
        cout << "\n<-------------------------------->";
        cout << "\n<---- INITIALISATION FAILURE ---->";
        cout << "\n<-------------------------------->" << endl;
    }
    cout << "\n##########################";
    cout << "\n##### FREEING MEMORY #####";
    cout << "\n##########################" << endl;
    FreeMem();
    
    cout << "Memory freed." << endl;
    
    if(initialisation_success) exit(EXIT_SUCCESS);
    else exit(EXIT_FAILURE);
}

/* function definitions */
void Tree::SetZero(int site) {
    
    /* This initialises the trees with basic values */
    t_age = 0;
    t_site = site;
    t_dbh = 0.0;
    t_Tree_Height = 0.0;
    t_Crown_Radius = 0.0;
    t_Crown_Depth = 0.0;
    t_CrownVolume = 0;
    t_wsg = 0.0;
    t_Pmass = 0.0;
    t_Nmass = 0.0;
    t_LMA = 0.0;
    t_dens_top = 0.0;
    t_dens_belowtop = 0.0;
    t_dens = 0.0;
    t_LAIabove = 0.0;
    t_LAImax = 0.0;
    t_LAI = 0.0;
    t_overlap = 0.0;
    t_overlap_multiple = 0.0;
    t_dev_height = 0.0;
    t_dev_CR = 0.0;
    t_dev_CD = 0.0;
    
    t_GPP = 0.0;
    t_NPP = 0.0;
    
}

void Tree::InitialiseTree(int site, float dbh, float height, float CR, float CD, float Nmass, float Pmass, float LMA, float wsg, int species_label) {
    
    /* This initialises the trees with basic values */
    t_age = 1;
    t_site = site;
    t_dbh = dbh;
    t_Tree_Height = minf(height,float(height_max));
    t_Crown_Radius = CR;
    t_Crown_Depth = CD;
    t_Nmass = Nmass;
    t_Pmass = Pmass;
    t_LMA = LMA;
    t_wsg = wsg;
    
    t_species_label = species_label;
  
    float height_mean = b_height * dbh / (dbh + a_height); //* exp(deviation_height);
    float CD_mean = a_CD + b_CD * height;
    float CR_mean = exp(a_CR + b_CR * log(dbh));
    
    t_dev_height = log(height/height_mean);
    t_dev_CD = log(CD/CD_mean);
    t_dev_CR = log(CR/CR_mean);

    if(Nmass > 0.0 && Pmass > 0.0 && LMA > 0.0){
        /* calculate the Farquhar parameters */
        
        float SLA=10000.0/LMA;
        float Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
        float Jmaxm=pow(10.0, minf((-1.50+0.41*log10(Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(Pmass*1000.0)+0.32*log10(SLA)))); // this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
        float Vcmax=Vcmaxm*LMA;
        float Jmax=Jmaxm*LMA;
        float Narea = LMA * Nmass;
        float Parea = LMA * Pmass;
        
        float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2)); //t_Rdark corresponds to leaf maintenance respiration. Originally from Table 6 in Atkin et al 2015 New phytologist v.2.0, updated to reflect published corretion
        
        t_LAImax = treeLAImax(Vcmax, Jmax, Rdark);
        t_LAIabove = 0.0;
        t_LAI = 0.0;        // TODO: should we set the LAI already here?
    } else {
        t_LAImax = 0.0;
        t_LAI = 0.0;
    }
    
    if(t_Crown_Depth > 0.0 && t_LAI > 0.0){
        UpdateDensity();
    } else {
        t_dens = t_dens_top = t_dens_belowtop = 0.0;
    }
    
    t_GPP = 0.0;
    t_NPP = 0.0;
    t_CrownDisplacement = 0;
}

void Tree::UpdateDensity(){
    if(flag_LAIgradient == 1) GetDensitiesGradient(t_LAI, t_Crown_Depth, t_dens_top, t_dens_belowtop, t_dens);
    else GetDensityUniform(t_LAI, t_Crown_Depth, t_dens_top, t_dens_belowtop, t_dens);
}

void Tree::CrownVolume(){
    t_CrownVolume = CalcCV(t_site + t_CrownDisplacement, t_Tree_Height, t_Crown_Radius, t_Crown_Depth);
}

/*#############################################
 ###   Farquhar von Caemmerer Berry model  ###
 #############################################*/

/* This function returns the leaf-level carbon assimilation rate in micromoles C/m^2/s according to Farquhar-von Caemmerer-Berry model */

float GPPleaf(float PPFD, float VPD, float T, float Vcmax, float Jmax) {
    /* v.2.3.0: theta defined as a global variable */
    /* Parameters for Farquhar model, with temperature dependencies */
    int convT= int(iTaccuracy*T); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    
    //if(convT>500 || isnan(convT) || convT <0) cout << t_site << " | convT: " << convT << " | T: " << T << " | PPFD: " << PPFD << " | VPD: " << VPD << endl;
    float KmT    =LookUp_KmT[convT];
    float GammaT =LookUp_GammaT[convT];
    
    // float g1 = -3.97 * t_wsg + 6.53 (Lin et al. 2015)
    
    float t_fci =g1/(g1+sqrt(VPD));
    float VcmaxT = Vcmax*LookUp_VcmaxT[convT];
    float JmaxT  = Jmax*LookUp_JmaxT[convT];
    
    /* FvCB model */
    
    float I = alpha * PPFD;
    float J = (I+JmaxT-sqrt((JmaxT+I)*(JmaxT+I)-4.0*theta*JmaxT*I))*0.5/theta;
    float A = minf(VcmaxT/(t_fci+KmT),0.25*J/(t_fci+2.0*GammaT))*(t_fci-GammaT);
    
    return A;
}



/* NEW in v. 2.4.0: Separate function for calculation of Rday */
/* separation necessitates extra function calls, but allows for decoupling of photosynthesis and leaf respiration in special cases (i.e. no photosynthesis at very low ppfd, but still respiration*/
/* in future: maybe higher respiration for very low ppfd?, no factor 0.4 */

float Rdayleaf(float T, float Rdark) {
    int convT= int(iTaccuracy*T);
    //if(T < 0 || isnan(T) || T > 50) cout << t_site << " species: " << t_s->s_name << " convT: " << T << endl;
    float Rleaf = Rdark*LookUp_Rday[convT];
    return Rleaf;
}

float dailyRdayleaf(float T, float Rdark) {
    float Rdayleaf_daily = 0.0;
    
    for(int i=0; i<24; i++){
        //cout << t_site << " i: " << i << " Rdayleaf_daily: " << Rdayleaf_daily << endl;
        Rdayleaf_daily += Rdayleaf(T*daily_T[i], Rdark);
    }
    
    Rdayleaf_daily *= 0.0417;
    return Rdayleaf_daily;
}

float dailyGPPleaf(float PPFD, float VPD, float T, float Vcmax, float Jmax) {
    float dailyA=0.0;
    
    for(int i=0; i<24; i++) {
        float ppfd_halfhour = PPFD * daily_light[i];
        float vpd_halfhour = VPD * daily_vpd[i];
        float t_halfhour = T * daily_T[i];
        if(ppfd_halfhour > 0.1) dailyA += GPPleaf(ppfd_halfhour,vpd_halfhour,t_halfhour, Vcmax, Jmax);
    }
    dailyA*=0.0417;                                 // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    return dailyA;
}

float treeLAImax(float Vcmax, float Jmax, float Rdark){
    
    /* the maximum LAI is assumed to lie between 0 and 10, and designates the LAI where adding additional leaves results in net losses in terms of carbon balance */
    /* we narrow down the range where the maximum LAI can lie sequentially, by halfing the possible range */
    /* letting i run from 0 to n, we thus have a precision of 1/(2^n) */
    /* importantly, we here use actually absorbed PPFD as well instead of incident PPFD, since we calculate not for a perfectly illuminated leaf, but one that is distributed realistically in the lower canopy layers */
    
    float LAI_lowerbound = 0.0;
    float LAI_upperbound = 10.0;
    float LAImax_temp = 0.5 * (LAI_lowerbound + LAI_upperbound);
    
    for(int i=0;i<10;i++){
        float absorb_prev = LAImax_temp;
        float absorb_delta = 0.5;                                           /* we calculate everything, assuming a medium leaf density in the lower canopy, TODO: if the crown has a non-cylindric shape, this should probably be taken as a third of LAImax, or if LAI_gradient is activated, as 25% of LAImax */
        int intabsorb = int(absorb_prev*20.0)+ 400*int(absorb_delta*20.0);  /* Attention: technically, both absorb_prev and absorb_delta should be controlled (e.g. whether they are < 9.95 and 19.95, but both conditions are automatically fulfilled here, since LAImax considered here is 10.0 and absorb_delta is fixed */
        
        /* get PPFD, VPD, and temperature at each discretisation step */
        float PPFD_LAI = Wmax_avg * LookUp_flux_absorption[intabsorb];
        float VPD_LAI = VPDmax_avg * LookUp_VPD[intabsorb];
        float Tmp_LAI = Tmax_avg - LookUp_T[intabsorb];
        
        /* calculate the GPP */
        
        float GPP_LAI = dailyGPPleaf(PPFD_LAI, VPD_LAI, Tmp_LAI, Vcmax, Jmax);
        float Rday_LAI = dailyRdayleaf(Tmp_LAI, Rdark);
        /* assuming that one third of the leaves are mature and that the rest of the leaves have half the respiration/assimilation rates, we derive a factor 0.66 */
        float effLA = 0.66 * 189.3 * timestep;

        GPP_LAI *= effLA;
        Rday_LAI *= effLA * 0.4;
        /* get the night respiration */
        int convTnight = int(iTaccuracy*Tnight_avg);

        float Rnight_LAI = Rdark * effLA * LookUp_Rnight[convTnight];
        
        /* we add up the two components of leaf respiration, and multiply the result by 1.5 (fine root respiration), since in TROLL, this cannot be separated from leaf respiration */
        
        /* calculate effective npp */
        float Rleaf_LAI = 1.5 * (Rday_LAI + Rnight_LAI);    // cf. explanation for tree respiration
        float NPP_leaf = 0.7 * (GPP_LAI - Rleaf_LAI);   // cf. explanation for t_NPP calculation
 
        /* update boundaries */
        if(NPP_leaf > 0.0) LAI_lowerbound = LAImax_temp;
        else LAI_upperbound = LAImax_temp;
        
        /* now update LAImax_temp */
        LAImax_temp = 0.5 * (LAI_lowerbound + LAI_upperbound);
    }
    
    return(LAImax_temp);
}


void GettreeNPPGPP(float &treeGPP, float &treeNPP, float Vcmax, float Jmax, float Rdark, int site_tree, int crown_displacement, float wsg, float height, float CR, float CD, float dbh, float LAI) {
    int site_crown = site_tree + crown_displacement;
    
    float leafarea = 0.0;
    float treeRday = 0.0;
    
    if(LAI > 0.0){
        int crown_base = int(height-CD),
        crown_top = int(height),
        row_crowncenter = site_crown/cols,
        col_crowncenter = site_crown%cols;

        float CR_mean = exp(a_CR + b_CR * log(dbh));
        float multiplier_CR = CR/CR_mean;
        float fraction_filled_general = 1.0 - crown_gap_fraction;
        float fraction_filled_target = minf(fraction_filled_general/(multiplier_CR * multiplier_CR),1.0);

        int max_shells = min(crown_top - crown_base + 1, 4);
        int layers_filled = 0;
        
        for(int shell_fromtop = 0; shell_fromtop < max_shells; shell_fromtop++){
            float canopy_environment_cumulated[4] = {0.0,0.0,0.0,0.0};

            LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, LAI, canopy_environment_cumulated, LAI2dens, GetCanopyEnvironment);
    
            float leafarea_layer = canopy_environment_cumulated[0];
            float ileafarea_layer;
            if(leafarea_layer > 0.0) ileafarea_layer = 1.0/leafarea_layer;
            else ileafarea_layer = 0.0;
 
            float PPFD_layer = canopy_environment_cumulated[1] * ileafarea_layer;
            float VPD_layer = canopy_environment_cumulated[2] * ileafarea_layer;
            float Tmp_layer = canopy_environment_cumulated[3] * ileafarea_layer;
            treeGPP += leafarea_layer * dailyGPPleaf(PPFD_layer, VPD_layer, Tmp_layer, Vcmax, Jmax);
            treeRday += leafarea_layer * dailyRdayleaf(Tmp_layer, Rdark);
            leafarea += leafarea_layer;
        }
    }
    /* assuming that one third of the leaves are mature and that the rest of the leaves have half the respiration/assimilation rates, we derive a factor 0.66 */
    float effLA = 0.66 * 189.3 * timestep;
    
    treeGPP *= effLA;
    treeRday *= effLA * 0.4;
    
    /* tree wood respiration */
    float sapwood_area =  0.0001 * 2.0 * leafarea / (0.066 + 0.017 * height - 0.18 + 1.6 * wsg); /* from Fyllas et al. 2014, based on inversion of pipe model, multiplication with 0.0001 to convert cm2 to m2 */
    float sapwood_minimum;
    if(dbh < 0.01) sapwood_minimum = dbh * dbh * 0.25 * PI;
    else sapwood_minimum = 0.005 * (dbh - 0.005) * PI;
    sapwood_area = maxf(sapwood_area, sapwood_minimum);
    
    /* get the night respiration */
    int convTnight = int(iTaccuracy*Tnight_avg);
    float treeRnight = Rdark * effLA * LookUp_Rnight[convTnight] * leafarea;
    
    int convT= int(iTaccuracy*temp_avg); // temperature data at a resolution of Taccuracy=0.1°C -- stored in lookup tables ranging from 0°C to 50°C ---
    float treeRstem = sapwood_area * (height-CD) * LookUp_Rstem[convT];
    
    /* we add up the two components of leaf respiration, and multiply the result by 1.5 (fine root respiration), since in TROLL, this cannot be separated from leaf respiration */
    
    /* calculate effective npp */
    float treeRleaf = 1.5 * (treeRday + treeRnight);
    treeRstem *= 1.5;
    treeNPP = 0.7 * (treeGPP - treeRleaf - treeRstem);
    
    //cout << "NPP_GPP conventional: " << site_tree << "\tsite_crown: " <<  site_crown << "\theight: " << height << "\tCR: " << CR << "\tCD: " << CD << "\tLeafarea: " << leafarea << " Rday: " << treeRday << " Rnight: " << treeRnight << " Rstem: " << treeRstem << endl;
}

float GetCanopyStatus(int site_crown, float height, float CR, float CD){
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int crown_area_int = GetCrownIntarea(CR);
    
    /* loop across the crown */
    
    int abovecrown_count = 0;
    int abovecrownheight_fullsum = 0;
    
    for(int i = 0; i < crown_area_int; i++){
        int site_relative = LookUp_Crown_site[i];
        
        /* choose orientation of tree depending on site, non-random in order to save calculation time */
        //            int row, col;
        //            if(row_crowncenter%2==0) row = row_crowncenter + site_relative/51 - 25;
        //            else row = row_crowncenter - site_relative/51 + 25;
        //            if(col_crowncenter%2==0) col = col_crowncenter + site_relative%51 - 25;
        //            else  col = col_crowncenter - site_relative%51 + 25;
        
        int row = row_crowncenter + site_relative/51 - 25;
        int col = col_crowncenter + site_relative%51 - 25;
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
            int site=col+cols*row;
            abovecrownheight_fullsum += int(chm_field[site]);
            abovecrown_count++;
        }
    }
    
    float heightcrownabove_avg;
    if(abovecrown_count > 0) heightcrownabove_avg = float(abovecrownheight_fullsum)/float(abovecrown_count);
    else heightcrownabove_avg = 0.0;
 
    float heightcrown_avg = GetAvgHeight(site_crown, height, CR, CD);
    float heightdiff = heightcrownabove_avg - heightcrown_avg;
    
    return(heightdiff);
}

float GetCrownExposure(bool canopy, int site_crown, float height, float CR, float CD, float radiusfactor){
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    int crown_area_int = GetCrownIntarea(CR);

    /* loop across the crown */

    int abovecrown_count = 0;
    int abovecrownheight_fullsum = 0;

    for(int i = 0; i < crown_area_int; i++){
        int site_relative = LookUp_Crown_site[i];
        
        /* choose orientation of tree depending on site, non-random in order to save calculation time */
        //            int row, col;
        //            if(row_crowncenter%2==0) row = row_crowncenter + site_relative/51 - 25;
        //            else row = row_crowncenter - site_relative/51 + 25;
        //            if(col_crowncenter%2==0) col = col_crowncenter + site_relative%51 - 25;
        //            else  col = col_crowncenter - site_relative%51 + 25;
        
        int row = row_crowncenter + site_relative/51 - 25;
        int col = col_crowncenter + site_relative%51 - 25;
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
            int site=col+cols*row;
            int height_canopy = chm_field[site];
            if(height_canopy > 0){
                abovecrownheight_fullsum += height_canopy;
                abovecrown_count++;
            }
        }
    }

    /* now calculate the surrounding crown height */
    /* deduct the crown itself */

    int canopy_count = -abovecrown_count;
    int canopyheight_fullsum = -abovecrownheight_fullsum;

    float radius_surrounding = maxf(radiusfactor * CR, CR);

    int radius_surrounding_int = int(radius_surrounding) + 1;
    for(int col=max(0,col_crowncenter-radius_surrounding_int);col<min(cols,col_crowncenter+radius_surrounding_int+1);col++) {
        for(int row=max(0,row_crowncenter-radius_surrounding_int);row<min(rows,row_crowncenter+radius_surrounding_int+1);row++) {
            if((col-col_crowncenter)*(col-col_crowncenter)+(row-row_crowncenter)*(row-row_crowncenter)<radius_surrounding_int*radius_surrounding_int){
                int height_canopy = chm_field[col + row * cols];
                /* condition to exclude canopy spots that have not yet been filled by trees */
                if(height_canopy > 0){
                    canopyheight_fullsum += height_canopy;
                    canopy_count++;
                }
            }
        }
    }

    float heightcanopy_avg = 0.0;
    if(canopy_count > 0) heightcanopy_avg = float(canopyheight_fullsum)/float(canopy_count);

    float heightcrown_avg;
    if(canopy == 1){
        heightcrown_avg = 0.0;
        if(abovecrown_count > 0) heightcrown_avg = float(abovecrownheight_fullsum)/float(abovecrown_count);
    } else{
        heightcrown_avg = GetAvgHeight(site_crown, height, CR, CD);
    }

    float heightdiff = heightcrown_avg - heightcanopy_avg;
    
    return(heightdiff);
}

void FillVoxel(int increment, int site_crown, float height, float CR, float CD){
    int crown_base = int(height - CD),
    crown_top = int(height);
    int row_crowncenter = (site_crown)/cols;
    int col_crowncenter = (site_crown)%cols;
    
    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int increment_sum = 0;

    int layers_filled = min(crown_top - crown_base,3);

    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, increment, increment_sum, ModifyNoModification_int, UpdateVoxelValue);
}

int CalcCrownAreaExterior(int site_crown, float CR){

    int row_crowncenter = site_crown/cols;
    int col_crowncenter = site_crown%cols;

    int crown_area_int = GetCrownIntarea(CR);
    
    int crown_area_exterior = 0;
    
    for(int i = 0; i < crown_area_int; i++){
        int site_relative = LookUp_Crown_site[i];
    
        int row = row_crowncenter + site_relative/51 - 25;
        int col = col_crowncenter + site_relative%51 - 25;
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
        } else {
            crown_area_exterior++;
        }
    }
    
    return(crown_area_exterior);
}

/* the following functions are not member functions of the Tree class object, but operate at the tree level and could be converted to member functions. They're therefore defined next to the other tree level functions. Whether conversion to Tree members adds any benefit in terms of performance should be tested. */
/* the main function is a template to loop through one crown shell (i.e. a 1m layer of the tree) and do something within the tree's canopy, such as allocating leaves or computing the flux. It is supported by a second template that simulates the actual loop. The crown shell can be bent via a function to simulate the umbrella like crown shapes */
/* variables to be provided to the template */
/*  a) crown properties, including the crown position (row, col), the tree height, the crown radius, the crown depth, the fraction of filled voxels (inverse of gap fraction) as well as the layer/shell counting from the tree top */
/*  b) a function "GetRadiusLayer" that creates the "umbrella" shape. It defines the change in crown radius between the top of the crown and the second layer above the crown base. Here by default a linear slope */
/*  c) crown statistics to be updated, one as input, one as output. These statistics are updated by the function "ModifyCrownStatistic" and then applied across the crown with the function "UpdateCrownStatistic". While they always have to be provided, they do not always have to be used. */
/*  EXAMPLE 1: for adding volume to a voxel field, the input statistic would be +1, for removing volume, -1, the ModifyCrownStatistic function empty, and the update function would simply add the input value to the Voxel3D field */
/*  EXAMPLE 2: for adding leaves to the LAI field, the input statistic would be the tree LAI, or for removing, -LAI, ModifyCrownStatistic would convert the LAI to density within a specific layer, and the update function would simply add the resulting density values to the LAI3D field */
/*  Input and output variables can be separate types (e.g. int and float) and of different length (e.g. input can be a single variable, output can be a vector. This is needed, for example, to compute PPFD, VPD, Tmp and leafarea in Fluxh) */
/*  In the current implementation, crowns below 3m in crown depth are simply treated as cylinders, this could be changed in future implementations */

template <typename I, typename O, typename M, typename F>
void LoopLayerUpdateCrownStatistic_template(int row_center, int col_center, float height, float CR, float CD, float fraction_filled_target, int shell_fromtop, int layers_filled, float GetRadiusLayer(float, float, float), I CrownStatistic_input, O &CrownStatistic_output, M ModifyCrownStatistic_input, F UpdateCrownStatistic_output){

    int crown_top = int(height);
    
    /* we start out with 0 actually filled voxels. As a result, the first voxel will always be filled */
    float fraction_filled_actual = 0.0;

    if(CD <= 3.0){

        I CrownStatistic_input_modified;

        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_modified, CD, height, shell_fromtop);

        int crown_intarea_previous = 0;
        int crown_intarea = GetCrownIntarea(CR);
        layers_filled = max(layers_filled,0);

        int heightfill_top = crown_top - shell_fromtop;
        int heightfill_base = max(heightfill_top - layers_filled,0);
        for(int h = heightfill_top; h >= heightfill_base; h--){
            CircleAreaUpdateCrownStatistic_template(row_center, col_center, crown_intarea_previous, crown_intarea, fraction_filled_target, fraction_filled_actual, h, CrownStatistic_input_modified,CrownStatistic_output, UpdateCrownStatistic_output);
        }

    } else {
        
        /* This function computes the extent of the crown at every height layer, given a specific function */
        /* it separates out the innermost sector (a slowly increasing cylinder), and the surrounding parts of the crown */

        /* first the metrics with respect to the internal crown structure (i.e. z coordinate with respect to crown base) */
        float crownshell_base = height - CD + 2.0;              /* lower reference point for the crown slope function is two layers up from the crown base */
        float crownshell_extent = height - crownshell_base;          /* this is the extent from the "base layer" to the top */
        float crownshell_extent_toplayer = floor(crownshell_extent);      /* this is the extent to the lower limit of the toplayer */

        /* then we translate the crown coordinates into discretised variables with respect to the absolute location in the voxel field, as needed for location in the voxel field, with layers defined from top to bottom */
        int height_innermost = crown_top - shell_fromtop;
        int height_toplayer = int(crownshell_base + crownshell_extent_toplayer) - shell_fromtop;
        int height_baselayer = int(crownshell_base + 1.0) - shell_fromtop;
        
        /* now calculate the two modifications of the input statistic */
        I CrownStatistic_input_innermost;
        I CrownStatistic_input_outer;
        
        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_innermost, CD, height, shell_fromtop);
        ModifyCrownStatistic_input(CrownStatistic_input, CrownStatistic_input_outer, CD, crownshell_base, shell_fromtop);
        
        /* now do calculations */
        /* first the inner crown shell section that grows dynamically */
        int crown_intarea_previous = 0;
        
        if(!(layers_filled > 0) || height_innermost != height_toplayer){
            float radius_innermost = GetRadiusLayer(CR, crownshell_extent, crownshell_extent_toplayer);
            int crown_intarea_innermost = GetCrownIntarea(radius_innermost);
            CircleAreaUpdateCrownStatistic_template(row_center, col_center, crown_intarea_previous, crown_intarea_innermost, fraction_filled_target, fraction_filled_actual, height_innermost, CrownStatistic_input_innermost,CrownStatistic_output, UpdateCrownStatistic_output);
            if(!(layers_filled > 0)) crown_intarea_previous = crown_intarea_innermost;
        }

        /* now loop through the outer crown shell cylinders */
        for(int h_outer = height_toplayer; h_outer >= height_baselayer; h_outer--){
            /* calculating the radius of the current layer depending on the respective slopes, to be replaced by function  */
            //float radius_height = CR - crown_slope * (h_outer - height_baselayer);    /* for the lowest layer, i.e. h == height_baselayer, radius = t_Crown_Radius */
            int extent_layerouter = h_outer - height_baselayer;
            float radius_height = GetRadiusLayer(CR, crownshell_extent, extent_layerouter);
            int crown_intarea = GetCrownIntarea(radius_height);
            CircleAreaUpdateCrownStatistic_template(row_center, col_center, crown_intarea_previous, crown_intarea, fraction_filled_target, fraction_filled_actual, h_outer, CrownStatistic_input_outer, CrownStatistic_output,UpdateCrownStatistic_output);
            if(!(layers_filled > 0)) crown_intarea_previous = crown_intarea;
        }
        
        /* plus a function to potentially fill up the crown below the layers considered */
         
        if(layers_filled > 0){
            /* crown area is now uniformly at maximum extent */
            int crown_intarea = GetCrownIntarea(CR);
             
            int heightfill_top = height_baselayer - 1;
            int heightfill_base = max(height_baselayer - layers_filled,0);
            for(int h = heightfill_top; h >= heightfill_base; h--){
                CircleAreaUpdateCrownStatistic_template(row_center, col_center, crown_intarea_previous, crown_intarea, fraction_filled_target, fraction_filled_actual, h, CrownStatistic_input_outer, CrownStatistic_output,UpdateCrownStatistic_output);
            }
        }
    }
}

/* This is the template that simulates the actual circling through a crown shell, between a given start and stop position (i.e. between a starting crown area and a stopping crown area) */
/* If the crown shells are bent, i.e. extend across several canopy layers, the starting position within a lower layer is the stopping position of the layer just above */

template <typename I, typename O, typename F>
void CircleAreaUpdateCrownStatistic_template(int row_center, int col_center, int pos_start, int pos_end, float fraction_filled_target, float &fraction_filled_actual, int height_layer, I CrownStatistic_input, O &CrownStatistic_output, F UpdateCrownStatistic){
 
    for(int i = pos_start; i < pos_end; i++){
        
        if(fraction_filled_actual > fraction_filled_target){
            fraction_filled_actual = (fraction_filled_actual * float(i))/(float(i) + 1.0);
        } else {
            fraction_filled_actual = (fraction_filled_actual * float(i) + 1.0)/(float(i) + 1.0);
            
            int site_relative = LookUp_Crown_site[i];
            int row, col;
            
            row = row_center + site_relative/51 - 25;
            col = col_center + site_relative%51 - 25;
            
            if(row >= 0 && row < rows && col >= 0 && col < cols){
                int site=col+cols*row;
                UpdateCrownStatistic(height_layer, site, CrownStatistic_input, CrownStatistic_output);
            }
        }
    }
}

/* linear decrease of crown radius */
float GetRadiusSlope(float CR, float crown_extent, float crownshell_extent_toplayer){
    float crown_slope = CR * (1.0 - shape_crown) / crown_extent;
    float radius = CR - crown_slope * float(crownshell_extent_toplayer);
    return(radius);
}

/* not currently used, but simply returns the input radius */
float GetRadiusCylinder(float CR, float crown_extent, float crownshell_extent_toplayer){
    return(CR);
}

/* converts floating point crown area into integer value, imposing lower and upper limits */
int GetCrownIntarea(float crown_radius){
    /* crown area */
    float crown_area = PI * crown_radius * crown_radius;
    int crown_intarea = int(crown_area);                                   // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea,1);                                 // minimum area of crown (1)
    crown_intarea = min(crown_intarea,1963);                              // maximum area of crown (radius 25)
    
    return(crown_intarea);
}

float GetCrownAreaFilled(float crown_radius, float fraction_filled_target){
    /* for now calculated explicitly */
    /* ideally, we would replace the loop by an equation that expresses the underlying logic */
    
    float crown_area = PI * crown_radius * crown_radius;
    int crown_intarea = int(crown_area);     // floor of crown_area to bound area accumulation
    crown_intarea = max(crown_intarea,1);                                 // minimum area of crown (1)
    crown_intarea = min(crown_intarea,1963);                              // maximum area of crown (radius 25)

    int crown_intarea_gaps = 0;
    float fraction_filled_actual = 0.0;
    
    for(int i = 0; i < crown_intarea; i++){
        if(fraction_filled_actual > fraction_filled_target){
            fraction_filled_actual = (fraction_filled_actual * float(i))/(float(i) + 1.0);
            crown_intarea_gaps++;
        } else {
            fraction_filled_actual = (fraction_filled_actual * float(i) + 1.0)/(float(i) + 1.0);
        }
    }
    
    /* now determine crown_area_filled, depending on whether the next voxel is filled or not filled */
    float crown_area_filled;
    if(fraction_filled_actual > fraction_filled_target){
        crown_area_filled= float(crown_intarea - crown_intarea_gaps);
    } else {
        crown_area_filled = crown_area - float(crown_intarea_gaps);
    }
    
    return(crown_area_filled);
}

/* deduces within-crown densities from LAI with a gradient from 50% in top layer to 25% in belowtop and 25% in all shells underneath (1 layer for umbrella-like shape) */
void GetDensitiesGradient(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens){
    if(CD < 2.0){
        dens_top = dens_belowtop = dens = LAI/ CD;
    } else if (CD < 3.0){
        dens_top = 0.5 * LAI;
        dens_belowtop = dens = 0.5 * LAI / (CD - 1.0);
    } else {
        dens_top = 0.5 * LAI;
        dens_belowtop = 0.25 * LAI;
        dens = 0.25 * LAI;
    }
}

/* deduces within-crown density from LAI, assuming uniform leaf distribution */
void GetDensityUniform(float LAI, float CD, float &dens_top, float &dens_belowtop, float &dens){
    float crownshells_limit = minf(CD, 3.0);
    dens_top = dens_belowtop = dens = LAI / crownshells_limit;
}

/* modifier for GPP calculation where we need the leaves per layer to weight our results */
/* LAI is input, dens_layer output */

/* a dummy function when no modification is needed */
void ModifyNoModification_float(float input, float &input_unmodified, float CD, float height, int layer_fromtop){
    input_unmodified = input;
}

void ModifyNoModification_int(int input, int &input_unmodified, float CD, float height, int layer_fromtop){
    input_unmodified = input;
}

/* a modifying function that converts LAI to the density of a specific layer, using the GetDensity functions above */
void LAI2dens(float LAI, float &dens_layer, float CD, float height, int layer_fromtop){
    
    int crown_top = int(height);
    int crown_base = int(height - CD);
    
    float dens_top, dens_belowtop, dens;
    if(flag_LAIgradient == 1) GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens);
    else GetDensityUniform(LAI, CD, dens_top,dens_belowtop, dens);
    dens_top = dens_belowtop = dens;
    
    if(CD < 3.0 && crown_top == crown_base){
        dens_layer = dens_top * CD;
    } else if(CD < 3.0 && (crown_top - layer_fromtop == crown_base)){
        float fraction_belowbase = float(crown_base+1) - (height - CD);
        dens_layer = dens * fraction_belowbase;
    } else {
        float fraction_layer = height - floor(height);    /* this is the fraction that each layer apart from the topmost layer will extend into the voxel above */
        float fraction_layer_fromabove = 1.0 - fraction_layer;                              /* the inverse of the fraction above */
        
        if(layer_fromtop == 0) dens_layer = dens_top * fraction_layer;
        else if(layer_fromtop == 1) dens_layer = dens_top * fraction_layer_fromabove + dens_belowtop * fraction_layer;
        else if(layer_fromtop == 2) dens_layer = dens_belowtop * fraction_layer_fromabove + dens * fraction_layer;
        else dens_layer = dens * fraction_layer_fromabove;
    }
}

/* a modifying function that converts LAI to the density of a specific layer, using the GetDensity functions above */
/* the difference to LAI2dens is that densities are cumulated across layers so that they can be directly allocated to the LAI3D field without summing them up afterwards */
void LAI2dens_cumulated(float LAI, float &dens_layer, float CD, float height, int layer_fromtop){
    
    int crown_top = int(height);
    int crown_base = int(height - CD);
    
    float dens_top, dens_belowtop, dens;
    if(flag_LAIgradient == 1) GetDensitiesGradient(LAI, CD, dens_top, dens_belowtop, dens);
    else GetDensityUniform(LAI, CD, dens_top, dens_belowtop, dens);
    dens_top = dens_belowtop = dens;
    
    if(CD < 3.0 && crown_top == crown_base){
        dens_layer = LAI;         /* full LAI allocation */
    } else if(CD < 3.0 && (crown_top - layer_fromtop == crown_base)){
        dens_layer = LAI;
    } else {
        float fraction_layer = height - floor(height);    /* this is the fraction that each layer apart from the topmost layer will extend into the voxel above */
        if(layer_fromtop == 0) dens_layer = dens_top * fraction_layer;
        else if(layer_fromtop == 1) dens_layer = dens_top + dens_belowtop * fraction_layer;
        else if(layer_fromtop == 2) dens_layer = dens_top + dens_belowtop + dens * fraction_layer;
        else dens_layer = LAI;
    }
}


/* the actual Update function for CalcLAI() */
void UpdateLAI3D(int height, int site, float dens, float &LA_cumulated){
    float lai3d = LAI3D[height][site];
    lai3d += dens;
    if(lai3d < 0.000001){
        lai3d = 0.0;
    }
    LAI3D[height][site] = lai3d;
    LA_cumulated += dens;
}

/* the PPFD retrieval function for leafarea_max() */
void GetPPFDabove(int height, int site, float noinput, float (&ppfd_CA)[2]){
    /* first get voxel field densities */
    float absorb_prev = LAI3D[height+1][site];

    /* translate into absorptance value */
    absorb_prev = minf(absorb_prev,19.95);
    int intabsorb = int(absorb_prev*20.0);
    float flux = LookUp_flux[intabsorb];

    /* obtain PPFD for the voxel, and also record the circled area */
    ppfd_CA[0] += Wmax_avg * flux;
    ppfd_CA[1] += 1.0; /* add area */
}

void AddHeight(int height, int site, int noinput, int (&height_avg)[2]){
    /* obtain PPFD for the voxel, and also record the circled area */
    height_avg[0] += height;
    height_avg[1]++; /* add area */
}
    
void AddNegativeDiff(int height, int site, int noinput, int (&negdiff)[2]){
    /* obtain PPFD for the voxel, and also record the circled area */
    int height_emp = chm_empirical[site];
    
    int diff = height_emp - height;
    if(diff < 0) negdiff[0] += diff;
    negdiff[1] += 1;
}

void AddupVolume(int height, int site, int voxel_volume, int &crown_volume){
    crown_volume += voxel_volume;
}

void AddupVolumeLayers(int height, int site, int voxel_volume, vector<int> &crown_volume_layers){
    crown_volume_layers[height]++;
}

void UpdateVoxelValue(int height, int site, int increment, int &increment_sum){
    Voxel3D[height][site] += increment;
    increment_sum += increment;
}

void GetVoxelValue(int height, int site, int noinput, int (&voxels)[2]){
    if(Voxel3D[height][site] > 1){
        voxels[0]++;
    }
    voxels[1]++;
}

void GetVoxelValue_countmultiple(int height, int site, int noinput, int (&voxels)[2]){
    int trees_overlapping = Voxel3D[height][site];
    if(trees_overlapping > 1){
        voxels[0] += (trees_overlapping - 1);
    }
    voxels[1]++;
}


/* the PPFD, VPD, Tmp and leafarea_layer retrieval function for Fluxh() */
void GetCanopyEnvironment(int height, int site, float dens, float (&canopy_environment_cumulated)[4]){
    /* this function adds to the environmental variables provided in canopy_environment_cumulated */
    /* first get voxel field densities */
    float absorb_prev = LAI3D[height+1][site];
    float absorb_curr = LAI3D[height][site];
    float absorb_delta = absorb_curr - absorb_prev;
    if(absorb_delta < 0.001) absorb_delta = 0.0;    /* eliminate rounding errors */

    /* translate into absorptance value */
    absorb_delta = minf(absorb_delta,9.95);
    absorb_prev = minf(absorb_prev,19.95);
    int intabsorb = int(absorb_prev*20.0) + 400*int(absorb_delta*20.0);
    
    /* obtain PPFD, VPD and T for the voxel */
    float PPFD_voxel = Wmax_avg * LookUp_flux_absorption[intabsorb];
    float VPD_voxel = VPDmax_avg * LookUp_VPD[intabsorb];
    float T_voxel = Tmax_avg - LookUp_T[intabsorb];
    
    //cout << site << " height: " << height << " dens: " << dens << " absorb_prev: " << absorb_prev << " absorb_curr: " << absorb_curr << " absorb_delta : " << absorb_delta << " intabsorb: " << intabsorb << endl;
    
    /* add the three variables up, weighted by leaf density inside voxel*/
    canopy_environment_cumulated[0] += dens;
    canopy_environment_cumulated[1] += PPFD_voxel * dens;
    canopy_environment_cumulated[2] += VPD_voxel * dens;
    canopy_environment_cumulated[3] += T_voxel * dens;
}

/* this is a function specifically written to calculate packing densities */
void AddCrownVolumeLayer(int row_center, int col_center, float height, float CR, float CD, int crownvolume[70]){

    int crown_top = int(height);
    int crown_base = int(height-CD);
    
    if(CD <= 3.0){
        /* for the smallest crowns it is simply cylinder rings being filled up */
        int crown_intarea = GetCrownIntarea(CR);
        for(int h = crown_top; h >= crown_base; h--){
            crownvolume[h] += crown_intarea;
        }
        
    } else {
        
        /* For the rest of the crown, we go through different crown shells. We separate out the innermost sector (a slowly increasing cylinder), and the surrounding parts of the crown */

        /* first the metrics with respect to the internal crown structure (i.e. z coordinate with respect to crown base) */
        float crownshell_base = height - CD + 2.0;              /* lower reference point for the crown slope function is two layers up from the crown base */
        float crownshell_extent = height - crownshell_base;          /* this is the extent from the "base layer" to the top */
        float crownshell_extent_toplayer = floor(crownshell_extent);      /* this is the extent to the lower limit of the toplayer */

        /* then we translate the crown coordinates into discretised variables with respect to the absolute location in the voxel field, as needed for location in the voxel field, with layers defined from top to bottom */
        int shell_fromtop = 0;
        int height_innermost = crown_top - shell_fromtop;
        int height_toplayer = int(crownshell_base + crownshell_extent_toplayer) - shell_fromtop;
        int height_baselayer = int(crownshell_base + 1.0) - shell_fromtop;
        
        /* now do calculations */
        /* first the inner crown shell section that grows dynamically */
        float radius_innermost = GetRadiusSlope(CR, crownshell_extent, crownshell_extent_toplayer);
        int crown_intarea_innermost = GetCrownIntarea(radius_innermost);
        for(int h = height_innermost; h >= crown_base; h--){
            crownvolume[h] += crown_intarea_innermost;
        }

        /* now loop through the outer crown shell cylinders */
        for(int h_outer = height_toplayer; h_outer >= crown_base; h_outer--){
            /* calculating the radius of the current layer depending on the respective slopes, to be replaced by function  */
            //float radius_height = CR - crown_slope * (h_outer - height_baselayer);    /* for the lowest layer, i.e. h == height_baselayer, radius = t_Crown_Radius */
            int extent_layerouter = max(h_outer - height_baselayer,0);                  /* we also fill up underneath the baselayer */
            float radius_height = GetRadiusSlope(CR, crownshell_extent, extent_layerouter);
            int crown_intarea = GetCrownIntarea(radius_height);
            
            crownvolume[h_outer] += (crown_intarea - crown_intarea_innermost);
            
        }
    }
}

int CalcCV(int site_crown, float height, float CR, float CD){
    
    int crown_base = int(height - CD),
    crown_top = int(height),
    row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int voxel_volume = 1;
    int crown_volume = 0;
    int layers_filled = min(crown_top - crown_base,3);

    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, voxel_volume, crown_volume, ModifyNoModification_int, AddupVolume);
    
    return(crown_volume);
}

float GetVoxeloverlap(int site_crown, float height, float CR, float CD, bool multiple){
    int crown_base = int(height - CD),
    crown_top = int(height),
    row_crowncenter = (site_crown)/cols,
    col_crowncenter = (site_crown)%cols;
    
    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int layers_filled = min(crown_top - crown_base,3);
    int noinput = 0;
    int voxels[2] = {0,0};
    
    if(multiple == 0){
        LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, noinput, voxels, ModifyNoModification_int, GetVoxelValue);
    } else {
        LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, noinput, voxels, ModifyNoModification_int, GetVoxelValue_countmultiple);
    }
    
    int nbvoxel_ovl = voxels[0];
    int nbvoxel_total = voxels[1];
    
    float overlap = float(nbvoxel_ovl)/float(nbvoxel_total);
    return(overlap);
}

float GetAvgHeight(int site_crown, float height, float CR, float CD){
    float height_tree_avg;
    
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;

    /* array, first element saves cumulated height, second the cumulated crown area, ratio yields avg crown height */
    int height_avg[2] = {0,0};

    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int layers_filled = 0;
    int noinput = 0;
    
    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, noinput, height_avg, ModifyNoModification_int, AddHeight);
    
    height_tree_avg = float(height_avg[0]);
    if(height_avg[1] > 0) height_tree_avg *= 1.0/float(height_avg[1]);
    
    return(height_tree_avg);
}

float GetAvgNegdiff(int site_crown, float height, float CR, float CD){
    int negdiff[2] = {0,0};
    
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;

    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int layers_filled = 0;
    int noinput = 0;
    
    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, noinput, negdiff, ModifyNoModification_int, AddNegativeDiff);
    
    float avg_negdiff = float(negdiff[0]);
    if(negdiff[1] > 0.0) avg_negdiff *= 1.0/float(negdiff[1]);
    return(avg_negdiff);
}

void GetCrownHeights(vector<int> &crown_heights, float height, float CR, float CD){

    for(int i = 0; i < crown_heights.size(); i++){
        crown_heights[i] = 0;
    }
    
    if(CD < 3.0){
        int crown_area_int = GetCrownIntarea(CR);
        
        int crown_area_nogaps_int = 0;
        int height_int = int(height);
        
        /* calculate the crown area without gaps */
        for(int i = 0; i < crown_area_int; i++){
            crown_area_nogaps_int++;
            crown_heights[i] = height_int;
        }
    
    } else {
        /* the height variable to sum up heights across the crown */
        float crownshell_base = height - CD + 1.5;
        float crown_extent = height - crownshell_base;
        float crown_slope = CR * (1.0 - shape_crown) / (crown_extent - 0.5);
        /* area for the innermost (and innermost) cylinder: this is the part of the crown that keeps pushing out and creates new layers */
        int crown_extent_toplayer = int(crown_extent - 0.5);                 /* the extent gives the number of full layers both upwards and downwards not including the central layer (i.e. half of it, ~ 0.5) */
        float radius_innermost = CR - crown_slope * float(crown_extent_toplayer);
        
        /* calculate the innermost/uppermost cylinder */
        int crown_area_previous = 0;
        
        /* crown area */
        int crown_area_innermost_int = GetCrownIntarea(radius_innermost);
        
        /* height at which we allocate */
        int height_innermost = int(height);
        
        for(int i = crown_area_previous; i < crown_area_innermost_int; i++){
            crown_heights[i] = height_innermost;
                
        }
        
        /* update previous crown area */
        crown_area_previous = crown_area_innermost_int;
        
        /* now loop through the outer crown cylinders */
        int height_toplayer = int(crownshell_base + 0.5 + float(crown_extent_toplayer));
        int height_baselayer = int(crownshell_base + 1.5);
        
        for(int h = height_toplayer; h >= height_baselayer; h--){
            /* calculating the radius of the current layer depending on the respective slopes */
            
            float radius_height = CR - crown_slope * (h - height_baselayer);    /* for the lowest layer, i.e. h == height_baselayer, radius = t_Crown_Radius */
            
            /* crown area */
            int crown_area_int = GetCrownIntarea(radius_height);
            
            for(int i = crown_area_previous; i < crown_area_int; i++){
                crown_heights[i] += h;
            }
            
            /* update previous crown area */
            crown_area_previous = crown_area_int;
        }
    }
};
    
float UpdateLAIabove_eff(int site_crown, float height, float CR, float CD, float dbh){

    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;

    float ppfd_CA[2] = {0.0,0.0};
    
    float CR_mean = exp(a_CR + b_CR * log(dbh));
    float multiplier_CR = CR/CR_mean;
    float fraction_filled_general = 1.0 - crown_gap_fraction;
    float fraction_filled_target = minf(fraction_filled_general/(multiplier_CR * multiplier_CR),1.0);
    int shell_fromtop = 0;
    int layers_filled = 0;
    float noinput = 0;

    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, noinput, ppfd_CA, ModifyNoModification_float, GetPPFDabove);
    
    float ppfd_experienced = ppfd_CA[0];
    float crown_area_looped = ppfd_CA[1];

    if(crown_area_looped > 0.0){
        ppfd_experienced *= 1.0/crown_area_looped;
    } else {
        cout << "Warning. Crown with no values at site " << site_crown << endl;
        ppfd_experienced = Wmax_avg;
    }
    
    //if(ppfd_experienced < 0) cout << "ppfd: " << ppfd_experienced << " or crown_area_looped: " << crown_area_looped << endl;
    float la_experienced_eff = -log(ppfd_experienced/Wmax_avg)/kpar;
    if(isnan(la_experienced_eff)){
        
        cout << "ppfd: " << ppfd_experienced << " or crown_area_looped: " << crown_area_looped << endl;
    }
    
    if(la_experienced_eff < 0.0001) la_experienced_eff = 0.0;
    
    return(la_experienced_eff);
}

void CalcLAItrial(bool remove, int site_crown, float height, float CR, float CD, float LAI, float dbh) {
    
    /* only allocate if there are leaves to be allocated */
    if(LAI > 0.0){
        if(remove == 1) LAI = -LAI;
        
        int crown_base = int(height-CD),
        crown_top = int(height),
        row_crowncenter = site_crown/cols,
        col_crowncenter = site_crown%cols;
        
        float leafarea_cumulated = 0.0;   /* currently, an output variable is required by LoopLayerUpdateCrownStatistic_template, we here use leafarea_cumulated as control variable */
        float CR_mean = exp(a_CR + b_CR * log(dbh));
        float multiplier_CR = CR/CR_mean;
        float fraction_filled_general = 1.0 - crown_gap_fraction;
        float fraction_filled_target = minf(fraction_filled_general/(multiplier_CR * multiplier_CR),1.0);

        int max_shells = min(crown_top - crown_base + 1, 4);
        int layers_filled = 0;

        for(int shell_fromtop = 0; shell_fromtop < max_shells; shell_fromtop++){
            if(shell_fromtop == max_shells - 1) layers_filled = 100;
            LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, LAI, leafarea_cumulated, LAI2dens_cumulated, UpdateLAI3D);
        }
    }
}

int TreeIngrowth(int site_crown, float height_tree, float CR_tree){
    /* we calculate whether there is a larger tree underneath or a smaller tree above, which indicates that trees grow "through" each other, i.e. a smaller tree pierces the canopy of a larger tree */
    /* while this might happen naturally (trees leaning on another tree), it should be rare. A crown fitting algorithm could, however, preferentially select such tree combinations, since small trees above large trees create better fitting canopies and thus falsely render heterogeneous large crowns through an assemblage of small trees above larger trees */
    int treeingrowth = 0;
    
    int extent = 25;
    extent = extent - ceill(CR_tree);
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int min_row = max(0, row_crowncenter - extent);
    int min_col = max(0, col_crowncenter - extent);
    int max_row = min(rows, row_crowncenter + extent);
    int max_col = min(cols, col_crowncenter + extent);
    
    for(int r = min_row; r < max_row; r++){
        for(int c = min_col; c < max_col; c++){
            int s = c + r * cols;
            if(T[s].t_age > 0){
                /* calculate the distance */
                int site_crowncompare = s + T[s].t_CrownDisplacement;
                float height_treecompare = T[s].t_Tree_Height;
                float CR_treecompare = T[s].t_Crown_Radius;
                int row_crowncompare = site_crowncompare/cols;
                int col_crowncompare = site_crowncompare%cols;
                int col_dist = col_crowncompare - col_crowncenter;
                int row_dist = row_crowncompare - row_crowncenter;
                
                int dist_crowns_squared = col_dist * col_dist + row_dist * row_dist;
                
                if(dist_crowns_squared <= extent){
                    float dist_crowns = sqrt(dist_crowns_squared);
                    
                    /* now separate the two cases depending on a height comparison */
                    if(height_treecompare < height_tree){
                        /* we need to test whether the other tree (smaller) has a crown area large enough to encompass the tree's crown area */
                        /* to this effect, we extend the lign linking both trees till the edge of the tree's crown and test whether the combined distance is smaller than the crown radius of the tree underneath */
                        /* in case that CR_tree > CR_treecompare, this is trivially false */
                        float dist_crownsedge = dist_crowns + CR_tree;
                        if(dist_crownsedge <= CR_treecompare){
                            treeingrowth = 1;
                        }
                    } else {
                        /* we need to test whether the other tree (taller) has a crown area small enough to be encompassed by the tree's crown area */
                        float dist_crownsedge = dist_crowns + CR_treecompare;
                        if(dist_crownsedge <= CR_tree){
                            treeingrowth = 1;
                        }
                    }
                }
            }
        }
    }
    return(treeingrowth);
}

int TreeIngrowthFromabove(int site_crown, float height_tree, float CR_tree){
    /* we calculate whether there is a larger tree underneath or a smaller tree above, which indicates that trees grow "through" each other, i.e. a smaller tree pierces the canopy of a larger tree */
    /* while this might happen naturally (trees leaning on another tree), it should be rare. A crown fitting algorithm could, however, preferentially select such tree combinations, since small trees above large trees create better fitting canopies and thus falsely render heterogeneous large crowns through an assemblage of small trees above larger trees */

    int treeingrowth = 0;
    
    int extent = 25;
    extent = extent - ceill(CR_tree);
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int min_row = max(0, row_crowncenter - extent);
    int min_col = max(0, col_crowncenter - extent);
    int max_row = min(rows, row_crowncenter + extent);
    int max_col = min(cols, col_crowncenter + extent);
    
    for(int r = min_row; r < max_row; r++){
        for(int c = min_col; c < max_col; c++){
            int s = c + r * cols;
            if(T[s].t_age > 0){
                /* calculate the distance */
                int site_crowncompare = s + T[s].t_CrownDisplacement;
                int row_crowncompare = site_crowncompare/cols;
                int col_crowncompare = site_crowncompare%cols;
                int col_dist = col_crowncompare - col_crowncenter;
                int row_dist = row_crowncompare - row_crowncenter;
                
                float dist_crowns = sqrt(col_dist * col_dist + row_dist * row_dist);
                
                float height_treecompare = T[s].t_Tree_Height;
                float CR_treecompare = T[s].t_Crown_Radius;
                /* now separate the two cases depending on a height comparison */
                if(height_treecompare < height_tree){
                    /* we need to test whether the other tree (smaller) has a crown area large enough to encompass the tree's crown area */
                    /* to this effect, we extend the lign linking both trees till the edge of the tree's crown and test whether the combined distance is smaller than the crown radius of the tree underneath */
                    /* in case that CR_tree > CR_treecompare, this is trivially false */
                    float dist_crownsedge = dist_crowns + CR_tree;
                    if(dist_crownsedge <= CR_treecompare){
                        treeingrowth = 1;
                        //cout << "Site_crown: " << site_crown << " Col: " << col_crowncenter << " Row: " << row_crowncenter << " Height: " << height_tree << " CR: " << CR_tree << " site_crown2: " << site_crowncompare << " Col2: " << col_crowncompare << " Row2: " << row_crowncompare << " Height2: " << height_treecompare << " CR2: " << CR_treecompare << endl;
                    }
                }
            }
        }
    }
    
    return(treeingrowth);
}

void GetCanopyHeight(vector<int> &chm_patch, int site_crown, int extent){
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int min_row = max(0, row_crowncenter - extent);
    int min_col = max(0, col_crowncenter - extent);
    int max_row = min(rows, row_crowncenter + extent);
    int max_col = min(cols, col_crowncenter + extent);
    
    for(int r = min_row; r < max_row; r++){
        for(int c = min_col; c < max_col; c++){
            int site = c + r * cols;
            int height_sim = chm_field[site];
            chm_patch.push_back(height_sim);
        }
    }
}

void SetCanopyHeight(vector<int> chm_patch, int site_crown, int extent){
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int min_row = max(0, row_crowncenter - extent);
    int min_col = max(0, col_crowncenter - extent);
    int max_row = min(rows, row_crowncenter + extent);
    int max_col = min(cols, col_crowncenter + extent);
    
    int i = 0;
    
    for(int r = min_row; r < max_row; r++){
        for(int c = min_col; c < max_col; c++){
        
            int site = c + r * cols;
            int height_sim = chm_patch[i];
            chm_field[site] = height_sim;
            i++;
        }
    }
}

void FindBetterCrown(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstarv_flag){
    /* calculate statistics */
    nbattempt++;
    
    /* now get the dbh bin */
    float dbh = T[site_tree].t_dbh;
    float log10dbh = log10(dbh*100.0);
    int bin_tree = int(log10dbh/dbhlog_binsize);
    int bindex = 0;
    for(int bin = 0; bin < bin_tree; bin++){
        bindex += trees_dbhlogbins[bin];
    }
    
    /* now we need to randomly pick a crown from the same dbh bin to exchange crowns only for similarly sized trees */
    /* to do so we get the size of the bin */
    int nbtrees_bin = trees_dbhlogbins[bin_tree];
    int treeexchange = int(gsl_rng_uniform_int(gslrand, nbtrees_bin));
    int site_treeexchange = trees_dbhsorted[bindex + treeexchange];
    
    /* we then get the tree's properties */
    int crowndisplacement_tree = T[site_tree].t_CrownDisplacement;
    float dev_height_tree = T[site_tree].t_dev_height;
    float dev_CR_tree = T[site_tree].t_dev_CR;
    float dev_CD_tree = T[site_tree].t_dev_CD;
    float height_tree = T[site_tree].t_Tree_Height;
    float CR_tree = T[site_tree].t_Crown_Radius;
    float CD_tree = T[site_tree].t_Crown_Depth;
    float Nmass_tree = T[site_tree].t_Nmass;
    float Pmass_tree = T[site_tree].t_Pmass;
    float LMA_tree = T[site_tree].t_LMA;
    int site_crown = site_tree + crowndisplacement_tree;
    
    /* Now we differentiate between two cases */
    if(site_tree == site_treeexchange){
        /* 1) in case we draw the tree itself, we will simply consider moving the tree crown a bit and drawing new dimensions */
        /* get new crown displacement */
        int crowndisplacement_updated = CalcCrownDisplacement(site_tree, height_tree);
        int site_crown_updated = site_tree + crowndisplacement_updated;
        
        /* now calculate the distance between old and new tree tree centers, add 1m to be outside of it */
        int col_crown = site_crown%cols;
        int row_crown = site_crown/cols;
        int col_crown_updated = site_crown_updated%cols;
        int row_crown_updated = site_crown_updated/cols;
        
        float crown_shift = sqrt((row_crown_updated - row_crown) * (row_crown_updated - row_crown) + (col_crown_updated - col_crown) * (col_crown_updated - col_crown));
        crown_shift += 1.0;
        
        /* now update the dimension */
        /* this function computes a new crown radius and height and then assesses the canopy fit */
        
        float height_tree_updated = height_tree;
        float CR_tree_updated = CR_tree;
        float CD_tree_updated = CD_tree;
        float dev_height_tree_updated = dev_height_tree;
        float dev_CR_tree_updated = dev_CR_tree;
        float dev_CD_tree_updated = dev_CD_tree;
        
        float sigma_heightbin_tree = sigma_heightbin[bin_tree];
        float sigma_CRbin_tree = sigma_CRbin[bin_tree];
        
        /* draw new dimensions */
        if(nbtrees_bin <= 10){
           float sumheight_bin_prev = correlation_structure[bin_tree][1];
           float sumCR_bin_prev = correlation_structure[bin_tree][3];
           
           float mean_height_prev = sumheight_bin_prev/float(nbtrees_bin);
           float mean_CR_prev = sumCR_bin_prev/float(nbtrees_bin);
           
           int success_correlation = 0;
           int tries = 0;
           int max_tries = 10;
           
           while(success_correlation == 0 && tries < max_tries){
               tries++;

               float dev_height_tree_test = float(gsl_ran_gaussian(gslrand, sigma_heightbin_tree));
               float dev_CR_tree_test = float(gsl_ran_gaussian(gslrand, sigma_CRbin_tree));

               float dev_CD_tree_test = float(gsl_ran_gaussian(gslrand, sigma_CD));
               
               if((mean_height_prev < 0 && dev_height_tree_test < 0) || (mean_height_prev > 0 && dev_height_tree_test > 0) || (mean_CR_prev < 0 && dev_CR_tree_test < 0) || (mean_CR_prev > 0 && dev_CR_tree_test > 0)){
               } else {
                   success_correlation = 1;
                   dev_height_tree_updated = dev_height_tree_test;
                   dev_CR_tree_updated = dev_CR_tree_test;
                   dev_CD_tree_updated = dev_CD_tree_test;
                   height_tree_updated  = b_height * dbh / (dbh + a_height) * exp(dev_height_tree_updated);
                   height_tree_updated  = minf(height_tree_updated , height_max-1);
                   height_tree_updated = maxf(height_tree_updated , 1.0);   // minimum 1m height!
                   
                   CR_tree_updated = exp(a_CR + b_CR * log(dbh) + dev_CR_tree_updated);
                   CR_tree_updated = minf(CR_tree_updated, 25.0);
                   
                   CD_tree_updated = (a_CD + b_CD * height_tree_updated) * exp(dev_CD_tree_updated);
                   CD_tree_updated = minf(CD_tree_updated, height_tree_updated * 0.5);
                   //cout << bin_tree << " Update!" << endl;
               }
           }
        }
        
        /* the radius for a CHM update (needs to be the largest possible radius, including possible crown displacement) */
        float radius_update =  maxf(CR_tree, CR_tree_updated) + crown_shift;
        int extent_tree = int(radius_update) + 2;

        //        vector<int> chm_field_save(chm_field);
        vector<int> chm_patch_tree;
        chm_patch_tree.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
        vector<int> hist_CHM_sim_save(hist_CHM_sim);
        
        /* safe the statistics */
        float dissimilarityprevious_test = CalcDissimilarity(1);
        int abserror_previous_test = sum_abserror;
        float carbonstarv_previous_test = carbonstarv;
        int nbtrees_carbonstarv_previous_test = nbtrees_carbonstarv;
        
        float meanheight_previous = CanopyMeanheight();
        float CA_exterior_previous = CA_exterior;
        
        UpdateCanopyAroundTree(site_crown, site_crown_updated, height_tree, height_tree_updated, CR_tree, CR_tree_updated, CD_tree, CD_tree_updated, extent_tree);
        
        /* calculate the potential new canopy statistics */
        float dissimilaritynew_test = CalcDissimilarity(1);
        int abserror_new_test = sum_abserror;
        
        /* update the LAI field */
        vector<int> list_trees;
        vector<float> list_trees_properties;
        
        /* since it is computationally intensive, we only compute this when it is needed */
        if(carbonstarv_flag == 1){
            
            /* get all trees affected */
            GetTreeList(site_tree, CR_tree, CR_tree_updated, CD_tree, CD_tree_updated, height_tree, height_tree_updated, crowndisplacement_tree, crowndisplacement_updated, list_trees);
            
            list_trees_properties.reserve(4 * list_trees.size());
            /* update them and save results in list_trees_properties */
            UpdateLeaves_fromList(list_trees, list_trees_properties, site_tree, height_tree_updated, CR_tree_updated, CD_tree_updated, crowndisplacement_updated);
        }
        
        float carbonstarv_new_test = carbonstarv;
        
        /* now we test whether the new configuration is better than the previous one */
        float meanheight_new = CanopyMeanheight();
        float CA_exterior_new = CA_exterior;
        int improved = EvaluateCanopy(mae_flag, dissimilarity_flag, carbonstarv_flag, dissimilarityprevious_test, dissimilaritynew_test, abserror_previous_test, abserror_new_test, meanheight_previous, meanheight_new, CA_exterior_previous, CA_exterior_new, carbonstarv_previous_test, carbonstarv_new_test);
        
        if(flag_PreventCrownpiercing == 1 && improved == 1){
            int ingrowth_tree = TreeIngrowth(site_crown_updated, height_tree, CR_tree);
            if(ingrowth_tree == 1) improved = 0;
        }
        
        if(improved == 1){
            /* update the success rate */
            nbsuccess++;
            T[site_tree].t_CrownDisplacement = crowndisplacement_updated;
            T[site_tree].t_dev_height = dev_height_tree_updated;
            T[site_tree].t_dev_CR = dev_CR_tree_updated;
            T[site_tree].t_dev_CD = dev_CD_tree_updated;
            T[site_tree].t_Tree_Height = height_tree_updated;
            T[site_tree].t_Crown_Radius = CR_tree_updated;
            T[site_tree].t_Crown_Depth = CD_tree_updated;
            
            correlation_structure[bin_tree][1] -= (dev_height_tree + dev_height_tree_updated);
            correlation_structure[bin_tree][2] -= (dev_height_tree * dev_height_tree + dev_height_tree_updated * dev_height_tree_updated);
            correlation_structure[bin_tree][3] -= (dev_CR_tree + dev_CR_tree_updated);
            correlation_structure[bin_tree][4] -= (dev_CR_tree * dev_CR_tree  + dev_CR_tree_updated * dev_CR_tree_updated);
            correlation_structure[bin_tree][5] -= (dev_CR_tree  * dev_height_tree + dev_CR_tree_updated * dev_height_tree_updated);

            if(carbonstarv_flag == 1){
                for(int i=0; i < list_trees.size(); i++){
                    int site_tree_update = list_trees[i];
                    T[site_tree_update].t_LAIabove = list_trees_properties[0 + i*4];
                    T[site_tree_update].t_LAI = list_trees_properties[1 + i*4];
                    T[site_tree_update].t_GPP = list_trees_properties[2 + i*4];
                    T[site_tree_update].t_NPP = list_trees_properties[3 + i*4];
                }
            }
        } else {
            /* otherwise, undo the changes */
            FillVoxel(-1, site_crown_updated, height_tree_updated, CR_tree_updated, CD_tree_updated);
            FillVoxel(1, site_crown, height_tree, CR_tree, CD_tree);
            
            //chm_field = chm_field_save;
            SetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
            hist_CHM_sim = hist_CHM_sim_save;
            sum_abserror = abserror_previous_test;
            CA_exterior = CA_exterior_previous;
            if(carbonstarv_flag == 1){
                ReverseLeaves(list_trees, list_trees_properties, site_tree, height_tree_updated, CR_tree_updated, CD_tree_updated, crowndisplacement_updated);
                carbonstarv = carbonstarv_previous_test;
                nbtrees_carbonstarv = nbtrees_carbonstarv_previous_test;
            }
            //UpdateCHM();
        }
        
    } else {
        
        /* 2) in the other case, we will exchange the deviations of the two trees, move the tree crowns at the same time, and evaluate the new goodness of fit */
        /* first, we get the properties of the second tree */
        int crowndisplacement_treeexchange = T[site_treeexchange].t_CrownDisplacement;
        float dev_height_treeexchange = T[site_treeexchange].t_dev_height;
        float dev_CR_treeexchange = T[site_treeexchange].t_dev_CR;
        float dev_CD_treeexchange = T[site_treeexchange].t_dev_CD;
        float height_treeexchange = T[site_treeexchange].t_Tree_Height;
        float CR_treeexchange = T[site_treeexchange].t_Crown_Radius;
        float CD_treeexchange = T[site_treeexchange].t_Crown_Depth;
        int site_crownexchange = site_treeexchange + crowndisplacement_treeexchange;
        float Nmass_treeexchange = T[site_treeexchange].t_Nmass;
        float Pmass_treeexchange = T[site_treeexchange].t_Pmass;
        float LMA_treeexchange = T[site_treeexchange].t_LMA;
        /* now compute the factors for intraspecific variation */
        float factor_height_tree = exp(dev_height_tree);
        float factor_CR_tree = exp(dev_CR_tree);
        float factor_CD_tree = exp(dev_CD_tree);
        
        float factor_height_treeexchange = exp(dev_height_treeexchange);
        float factor_CR_treeexchange = exp(dev_CR_treeexchange);
        float factor_CD_treeexchange = exp(dev_CD_treeexchange);
        
        /* exchange the factors and update the trees */
        /* original tree */
        float height_tree_updated = height_tree * factor_height_treeexchange/factor_height_tree;
        height_tree_updated = minf(height_tree_updated, height_max-1);
        height_tree_updated = maxf(height_tree_updated, 1.0);   // minimum 1m height!
        
        float CR_tree_updated = CR_tree * factor_CR_treeexchange/factor_CR_tree;
        CR_tree_updated = minf(CR_tree_updated, 25.0);
        
        // we update the crown depth directly, since the tree's height has changed
        float CD_tree_updated = (a_CD + b_CD * height_tree_updated) * factor_CD_treeexchange;
        CD_tree_updated = minf(CD_tree_updated, height_tree_updated * 0.5);
        
        /* tree to be exchanged */
        float height_treeexchange_updated = height_treeexchange * factor_height_tree/factor_height_treeexchange;
        height_treeexchange_updated = minf(height_treeexchange_updated, height_max-1);
        height_treeexchange_updated = maxf(height_treeexchange_updated, 1.0);   // minimum 1m height!
        
        float CR_treeexchange_updated = CR_treeexchange * factor_CR_tree/factor_CR_treeexchange;
        CR_treeexchange_updated = minf(CR_treeexchange_updated, 25.0);
        
        // we update the crown depth directly, since the tree's height has changed
        float CD_treeexchange_updated = (a_CD + b_CD * height_treeexchange_updated) * factor_CD_tree;
        CD_treeexchange_updated = minf(CD_treeexchange_updated, height_treeexchange_updated * 0.5);
        
        /* get new crown displacement */
        int crowndisplacement_tree_updated = CalcCrownDisplacement(site_tree, height_tree_updated);
        int crowndisplacement_treeexchange_updated = CalcCrownDisplacement(site_treeexchange, height_treeexchange_updated);
        
        int site_crown_updated = site_tree + crowndisplacement_tree_updated;
        int site_crownexchange_updated = site_treeexchange + crowndisplacement_treeexchange_updated;
        
        /* now we update the statistics, remove and add */
        /* this has to be done in succession, in case the two tree crowns overlap */
        int extent_tree = GetTreeExtent(site_crown, site_crown_updated, CR_tree, CR_tree_updated);
        int extent_treeexchange = GetTreeExtent(site_crownexchange, site_crownexchange_updated, CR_treeexchange, CR_treeexchange_updated);

        vector<int> chm_patch_tree;
        chm_patch_tree.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
        vector<int> chm_patch_treeexchange;
        chm_patch_treeexchange.reserve(4 * extent_treeexchange * extent_treeexchange);
        GetCanopyHeight(chm_patch_treeexchange, site_crownexchange, extent_treeexchange);
        vector<int> hist_CHM_sim_save(hist_CHM_sim);
        
        /* then update the fields */
        float dissimilarityprevious_test = CalcDissimilarity(1);
        int abserror_previous_test = sum_abserror;
        float carbonstarv_previous_test = carbonstarv;
        int nbtrees_carbonstarv_previous_test = nbtrees_carbonstarv;
        
        float meanheight_previous = CanopyMeanheight();
        float CA_exterior_previous = CA_exterior;

        UpdateCanopyAroundTree(site_crown, site_crown_updated, height_tree, height_tree_updated, CR_tree, CR_tree_updated, CD_tree, CD_tree_updated, extent_tree);
        UpdateCanopyAroundTree(site_crownexchange, site_crownexchange_updated, height_treeexchange, height_treeexchange_updated, CR_treeexchange, CR_treeexchange_updated, CD_treeexchange, CD_treeexchange_updated, extent_treeexchange);

        //UpdateCHM();
        
        /* calculate the new theoretical canopy statistics */
        float dissimilaritynew_test = CalcDissimilarity(1);
        int abserror_new_test = sum_abserror;
        
        /* update the LAI field */
        vector<int> list_trees;
        vector<float> list_trees_properties;

        //if(carbonstarv_flag == 1) UpdateLeaves();
        
        /* since it is computationally intensive, we only compute this when it is needed */
        if(carbonstarv_flag == 1){
            GetTreeList(site_tree, site_treeexchange, CR_tree, CR_tree_updated, CR_treeexchange, CR_treeexchange_updated, CD_tree, CD_tree_updated, CD_treeexchange, CD_treeexchange_updated, height_tree, height_tree_updated, height_treeexchange, height_treeexchange_updated, crowndisplacement_tree, crowndisplacement_tree_updated, crowndisplacement_treeexchange, crowndisplacement_treeexchange_updated, list_trees);
            
            list_trees_properties.reserve(4 * list_trees.size());
            
            UpdateLeaves_fromList(list_trees, list_trees_properties, site_tree, height_tree_updated, CR_tree_updated, CD_tree_updated, crowndisplacement_tree_updated, Nmass_tree, Pmass_tree, LMA_tree, site_treeexchange, height_treeexchange_updated, CR_treeexchange_updated, CD_treeexchange_updated, crowndisplacement_treeexchange_updated, Nmass_treeexchange, Pmass_treeexchange, LMA_treeexchange);
        }
        
        float carbonstarv_new_test = carbonstarv;
        
        float meanheight_new = CanopyMeanheight();
        float CA_exterior_new = CA_exterior;

        /* now we test whether the new configuration is better than the previous one */
        int bettercrowns = EvaluateCanopy(mae_flag, dissimilarity_flag, carbonstarv_flag, dissimilarityprevious_test, dissimilaritynew_test, abserror_previous_test, abserror_new_test, meanheight_previous, meanheight_new, CA_exterior_previous, CA_exterior_new, carbonstarv_previous_test, carbonstarv_new_test);
 
        if(flag_PreventCrownpiercing == 1 && bettercrowns == 1){
            int ingrowth_tree = TreeIngrowth(site_crown_updated, height_tree_updated, CR_tree_updated);
            int ingrowth_treeexchange = TreeIngrowth(site_crownexchange_updated, height_treeexchange_updated, CR_treeexchange_updated);
            if(ingrowth_tree == 1 || ingrowth_treeexchange == 1) bettercrowns = 0;
        }

        if(bettercrowns == 1){
            //if(abserror_new_test < abserror_previous_test){
            /* update the success rate */
            nbsuccess++;
            
            /* if so, use the updated dimensions */
            T[site_tree].t_Tree_Height = height_tree_updated;
            T[site_tree].t_Crown_Radius = CR_tree_updated;
            T[site_tree].t_Crown_Depth = CD_tree_updated;
            T[site_tree].t_CrownDisplacement = crowndisplacement_tree_updated;
            
            /* and exchange the variance values */
            T[site_tree].t_dev_height = dev_height_treeexchange;
            T[site_tree].t_dev_CR = dev_CR_treeexchange;
            T[site_tree].t_dev_CD = dev_CD_treeexchange;
            
            /* updated dimensions for second tree */
            T[site_treeexchange].t_Tree_Height = height_treeexchange_updated;
            T[site_treeexchange].t_Crown_Radius = CR_treeexchange_updated;
            T[site_treeexchange].t_Crown_Depth = CD_treeexchange_updated;
            T[site_treeexchange].t_CrownDisplacement = crowndisplacement_treeexchange_updated;
            
            /* and exchange of the variance values for the second tree */
            T[site_treeexchange].t_dev_height = dev_height_tree;
            T[site_treeexchange].t_dev_CR = dev_CR_tree;
            T[site_treeexchange].t_dev_CD = dev_CD_tree;

            if(carbonstarv_flag == 1){
                for(int i=0; i < list_trees.size(); i++){
                    int site_tree_update = list_trees[i];
                    T[site_tree_update].t_LAIabove = list_trees_properties[0 + i*4];
                    T[site_tree_update].t_LAI = list_trees_properties[1 + i*4];
                    T[site_tree_update].t_GPP = list_trees_properties[2 + i*4];
                    T[site_tree_update].t_NPP = list_trees_properties[3 + i*4];
                }
            }
        } else {
            /* otherwise, undo the changes */
            FillVoxel(-1, site_crown_updated, height_tree_updated, CR_tree_updated, CD_tree_updated);
            FillVoxel(1, site_crown, height_tree, CR_tree, CD_tree);
            FillVoxel(-1, site_crownexchange_updated, height_treeexchange_updated, CR_treeexchange_updated, CD_treeexchange_updated);
            FillVoxel(1, site_crownexchange, height_treeexchange, CR_treeexchange, CD_treeexchange);
            
            SetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
            SetCanopyHeight(chm_patch_treeexchange, site_crownexchange, extent_treeexchange);
            hist_CHM_sim = hist_CHM_sim_save;
            sum_abserror = abserror_previous_test;
            CA_exterior = CA_exterior_previous;

            if(carbonstarv_flag == 1){
                ReverseLeaves(list_trees, list_trees_properties, site_tree, height_tree_updated, CR_tree_updated, CD_tree_updated, crowndisplacement_tree_updated, site_treeexchange, height_treeexchange_updated, CR_treeexchange_updated, CD_treeexchange_updated, crowndisplacement_treeexchange_updated);
                carbonstarv = carbonstarv_previous_test;
                nbtrees_carbonstarv = nbtrees_carbonstarv_previous_test;
            }
        }
    }
}

int GetTreeExtent(int site_crown, int site_crown_updated, float CR_tree, float CR_tree_updated){
    /* first get the distance between the previous and the updated crown site */
    int col_crown = site_crown%cols;
    int row_crown = site_crown/cols;
    int col_crown_updated = site_crown_updated%cols;
    int row_crown_updated = site_crown_updated/cols;
    
    float crown_shift = sqrt((row_crown_updated - row_crown) * (row_crown_updated - row_crown) + (col_crown_updated - col_crown) * (col_crown_updated - col_crown));
    crown_shift += 1.0;
    
    /* the radius for a CHM update (needs to be the largest possible radius, including possible crown displacement) */
    float radius_update = maxf(CR_tree, CR_tree_updated) + crown_shift;
    
    int extent_tree = int(radius_update) + 2;
    return(extent_tree);
}

void UpdateCanopyAroundTree(int site_crown, int site_crown_updated, float height_tree, float height_tree_updated, float CR_tree, float CR_tree_updated,float CD_tree, float CD_tree_updated, int extent_tree){
    int CA_exterior_tree_previous = CalcCrownAreaExterior(site_crown, CR_tree);
    int CA_exterior_tree_new = CalcCrownAreaExterior(site_crown_updated, CR_tree_updated);
    CA_exterior += (CA_exterior_tree_new - CA_exterior_tree_previous);
    
    /* update the voxel field */
    FillVoxel(-1, site_crown, height_tree, CR_tree, CD_tree);
    FillVoxel(1, site_crown_updated, height_tree_updated, CR_tree_updated, CD_tree_updated);
    
    CanopyDevMetrics(site_crown, extent_tree, 1);
    UpdateCHMradius(site_crown, (CR_tree + 1.0));
    UpdateCHMradius(site_crown_updated, (CR_tree_updated + 1.0));
    CanopyDevMetrics(site_crown, extent_tree, 0);
}

#ifdef exchangepositions_step2
void UpdateCanopyNewPosition(int site_crown, int site_crown_updated, float height_tree, float CR_tree, float CD_tree, int extent_tree){
    int CA_exterior_tree_previous = CalcCrownAreaExterior(site_crown, CR_tree);
    int CA_exterior_tree_new = CalcCrownAreaExterior(site_crown_updated, CR_tree);
    CA_exterior += (CA_exterior_tree_new - CA_exterior_tree_previous);
    
    /* update the voxel field */
    FillVoxel(-1, site_crown, height_tree, CR_tree, CD_tree);
    FillVoxel(1, site_crown_updated, height_tree, CR_tree, CD_tree);
    
//    float dissimilarity_test = 0.0;
//    int sum_abserror_test = 0;
//    CanopyDistanceFull(dissimilarity_test, sum_abserror_test);
//
//    if(sum_abserror != sum_abserror_test) cout << nbattempt << " col " << site_crown%cols << " row: " << site_crown/cols << " updated col: " << site_crown_updated%cols << " updated row: " << site_crown_updated/cols << " distance: " << sqrt(pow(site_crown_updated%cols - site_crown%cols,2) + pow(site_crown_updated/cols - site_crown/cols,2)) << " radius: " << CR_tree << " extent: " << extent_tree << " sumabserror: " << sum_abserror << " recalc: " << sum_abserror_test << endl;
    
    CanopyDevMetrics(site_crown, extent_tree, 1);
    UpdateCHMradius(site_crown, (CR_tree + 1.0));
    CanopyDevMetrics(site_crown, extent_tree, 0);
    
//    dissimilarity_test = 0.0;
//    sum_abserror_test = 0;
//    CanopyDistanceFull(dissimilarity_test, sum_abserror_test);
//
//     if(sum_abserror != sum_abserror_test) cout << nbattempt << " col " << site_crown%cols << " row: " << site_crown/cols << " updated col: " << site_crown_updated%cols << " updated row: " << site_crown_updated/cols << " distance: " << sqrt(pow(site_crown_updated%cols - site_crown%cols,2) + pow(site_crown_updated/cols - site_crown/cols,2)) << " radius: " << CR_tree << " extent: " << extent_tree << " sumabserror: " << sum_abserror << " recalc: " << sum_abserror_test << endl;
    
    CanopyDevMetrics(site_crown_updated, extent_tree, 1);
    UpdateCHMradius(site_crown_updated, (CR_tree + 1.0));
    CanopyDevMetrics(site_crown_updated, extent_tree, 0);
    
//    dissimilarity_test = 0.0;
//    sum_abserror_test = 0;
//    CanopyDistanceFull(dissimilarity_test, sum_abserror_test);
//
//    if(sum_abserror != sum_abserror_test) cout << nbattempt << " col " << site_crown%cols << " row: " << site_crown/cols << " updated col: " << site_crown_updated%cols << " updated row: " << site_crown_updated/cols << " distance: " << sqrt(pow(site_crown_updated%cols - site_crown%cols,2) + pow(site_crown_updated/cols - site_crown/cols,2)) << " radius: " << CR_tree << " extent: " << extent_tree << " sumabserror: " << sum_abserror << " recalc: " << sum_abserror_test << endl;
}
#endif

int EvaluateCanopy(int mae_flag, int dissimilarity_flag, int carbonstarv_flag, float dissimilarityprevious, float dissimilaritynew, int abserror_previous, int abserror_new, float meanheight_previous, float meanheight_new, float CA_exterior_previous, float CA_exterior_new, float carbonstarv_previous, float carbonstarv_new){
    int success;
    
    /* four cases: each metric on its own, or all together, calibrated internally, otherwise nothing happens */
    if(dissimilarity_flag == 0 && mae_flag == 1 && carbonstarv_flag == 0){
#ifdef accept_equalfits
        if(abserror_new <= abserror_previous) success = 1;
        else success = 0;
#else
        if(abserror_new < abserror_previous) success = 1;
        else success = 0;
#endif
    } else if(dissimilarity_flag == 1 && mae_flag == 0 && carbonstarv_flag == 0){
#ifdef accept_equalfits
        if(dissimilaritynew <= dissimilarityprevious) success = 1;
        else success = 0;
#else
        if(dissimilaritynew < dissimilarityprevious) success = 1;
        else success = 0;
#endif
    } else if(dissimilarity_flag == 0 && mae_flag == 0 && carbonstarv_flag == 1){
#ifdef accept_equalfits
        if(carbonstarv_new <= carbonstarv_previous) success = 1;
        else success = 0;
#else
        if(carbonstarv_new < carbonstarv_previous) success = 1;
        else success = 0;
#endif
    } else {
        /* here we combine both metrics */

        float dissimilarity_previous_std;
        float dissimilarity_new_std;
        if(dissimilarity_flag == 1){
            dissimilarity_previous_std = (dissimilarityprevious - dissimilarity_min)/(dissimilarity_max - dissimilarity_min);
            dissimilarity_new_std = (dissimilaritynew - dissimilarity_min)/(dissimilarity_max - dissimilarity_min);
        } else {
            dissimilarity_previous_std = dissimilarity_new_std = 0.0;
        }
        
        float abserror_previous_std;
        float abserror_new_std;
        if(mae_flag == 1){
            abserror_previous_std = float(abserror_previous - sum_abserrormin)/float(sum_abserrormax - sum_abserrormin);
            abserror_new_std = float(abserror_new - sum_abserrormin)/float(sum_abserrormax - sum_abserrormin);
        } else {
            abserror_previous_std = 0.0;
            abserror_new_std = 0.0;
        }
          
        float carbonstarv_previous_std;
        float carbonstarv_new_std;
        if(carbonstarv_flag == 1){
            carbonstarv_previous_std = (carbonstarv_previous - carbonstarv_min)/(carbonstarv_max - carbonstarv_min);
            carbonstarv_new_std = (carbonstarv_new - carbonstarv_min)/(carbonstarv_max - carbonstarv_min);
        } else {
            carbonstarv_previous_std = 0.0;
            carbonstarv_new_std = 0.0;
        }

        float metric_combined_previous = sqrt(abserror_previous_std * abserror_previous_std + dissimilarity_previous_std * dissimilarity_previous_std + carbonstarv_previous_std * carbonstarv_previous_std);
        float metric_combined_new = sqrt(abserror_new_std * abserror_new_std + dissimilarity_new_std * dissimilarity_new_std + carbonstarv_new_std * carbonstarv_new_std);
        
#ifdef accept_equalfits
        if(metric_combined_new <= metric_combined_previous) success = 1;
        else success = 0;
#else
        if(metric_combined_new < metric_combined_previous) success = 1;
        else success = 0;
#endif

    }
    
    //if(((meanheight_previous - mean_height_emp) < 0 && meanheight_new < meanheight_previous) || ((meanheight_previous - mean_height_emp) > 0 && meanheight_new > meanheight_previous)) success = 0;
    
    if(((CA_exterior_previous - CA_exterior_random) < 0 && CA_exterior_new < CA_exterior_previous) || ((CA_exterior_previous - CA_exterior_random) > 0 && CA_exterior_new > CA_exterior_previous)) success = 0;
    
    return(success);
}
    
void FindBetterPosition(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstarv_flag){
    nbattempt++;
    
    /* this is a function to move a tree to another site, depending on whether that improves the canopy */
    float height_tree = T[site_tree].t_Tree_Height;
    float CR_tree = T[site_tree].t_Crown_Radius;
    float CD_tree = T[site_tree].t_Crown_Depth;
    float LMA_tree = T[site_tree].t_LMA;
    float Nmass_tree = T[site_tree].t_Nmass;
    float Pmass_tree = T[site_tree].t_Pmass;
    
    /* choose random site in the vicinity (exponential increase with tree height to ensure good fit), minimum 10m */
    float distance_max = height_tree*height_tree/10.0;
    //float distance_max = exp(height_tree/10.0);
    distance_max = maxf(distance_max, 5.0);
    distance_max = minf(distance_max, 200.0);
    
    float distance = gsl_rng_uniform(gslrand) * distance_max;
    float angle_move = float(twoPi*gsl_rng_uniform(gslrand));                                                /* angle in which tree is moved */

    int dist_cols = int(distance*cos(angle_move));
    int dist_rows = int(distance*sin(angle_move));
               
    int col_tree = site_tree%cols;
    int row_tree = site_tree/cols;
    int col_shifted = dist_cols + col_tree;
    int row_shifted = dist_rows + row_tree;
    
    /* we ensure that new site falls within simulated area via periodic boundary conditions */
    if(col_shifted < 0) col_shifted = cols - (abs(col_shifted)%cols);
    else col_shifted = col_shifted%cols;
    if(row_shifted < 0) row_shifted = rows - (abs(row_shifted)%rows);
    else row_shifted = row_shifted%rows;
    
    int site_tree_shifted = col_shifted + row_shifted * cols;
    
    if(T[site_tree_shifted].t_age == 0){
        /* the radius for a CHM update (needs to be the largest possible radius, including possible crown displacement) */
        float radius_update =  CR_tree;
        int extent_tree = int(radius_update) + 2;
        
        //vector<int> chm_field_save(chm_field);
        vector<int> chm_patch_tree;
        chm_patch_tree.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_tree, site_tree, extent_tree);
        vector<int> chm_patch_treeshifted;
        chm_patch_treeshifted.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_treeshifted, site_tree_shifted, extent_tree);
        vector<int> hist_CHM_sim_save(hist_CHM_sim);
        
        float dissimilarityprevious_test = CalcDissimilarity(1);
        int abserror_previous_test = sum_abserror;
        float carbonstarv_previous_test = carbonstarv;
        int nbtrees_carbonstarv_previous_test = nbtrees_carbonstarv;

        float meanheight_previous = CanopyMeanheight();
        float CA_exterior_previous = CA_exterior;
        
        UpdateCanopyNewPosition(site_tree, site_tree_shifted, height_tree, CR_tree, CD_tree, extent_tree);
        //UpdateCanopyAroundTree(site_tree, site_tree_shifted, height_tree, height_tree, CR_tree, CR_tree, CD_tree, CD_tree, extent_tree);
        
        float dissimilaritynew_test = CalcDissimilarity(1);
        int abserror_new_test = sum_abserror;
        
        /* update the LAI field */
        vector<int> list_trees;
        vector<float> list_trees_properties;
        
        /* since it is computationally intensive, we only compute this when it is needed */
        if(carbonstarv_flag == 1){
            GetTreeList(site_tree, site_tree_shifted, list_trees);
            
            list_trees_properties.reserve(4 * list_trees.size());
            UpdateLeaves_fromList(list_trees, list_trees_properties, site_tree, site_tree_shifted);
        }
        
        float carbonstarv_new_test = carbonstarv;
        
        float meanheight_new = CanopyMeanheight();
        float CA_exterior_new = CA_exterior;
        
        int move = EvaluateCanopy(mae_flag, dissimilarity_flag, carbonstarv_flag, dissimilarityprevious_test, dissimilaritynew_test, abserror_previous_test, abserror_new_test, meanheight_previous, meanheight_new, CA_exterior_previous, CA_exterior_new, carbonstarv_previous_test, carbonstarv_new_test);
        
        if(flag_PreventCrownpiercing == 1 && move == 1){
            int ingrowth_tree = TreeIngrowth(site_tree_shifted, height_tree, CR_tree);
            if(ingrowth_tree == 1) move = 0;
        }

        if(move == 1){
            nbsuccess++;
            
            float dbh = T[site_tree].t_dbh;
            float wsg = T[site_tree].t_wsg;
            int species_label = T[site_tree].t_species_label;
            
            T[site_tree_shifted].InitialiseTree(site_tree_shifted, dbh, height_tree, CR_tree, CD_tree, Nmass_tree, Pmass_tree, LMA_tree, wsg, species_label);
            if(carbonstarv_flag == 1){
                for(int i=0; i < list_trees.size(); i++){
                    int site_tree_update = list_trees[i];
                    if(site_tree_update == site_tree) site_tree_update = site_tree_shifted;
                    T[site_tree_update].t_LAIabove = list_trees_properties[0 + i*4];
                    T[site_tree_update].t_LAI = list_trees_properties[1 + i*4];
                    T[site_tree_update].t_GPP = list_trees_properties[2 + i*4];
                    T[site_tree_update].t_NPP = list_trees_properties[3 + i*4];
                }
            }
            
            T[site_tree].SetZero(site_tree);
        } else {
            FillVoxel(-1, site_tree_shifted, height_tree, CR_tree, CD_tree);
            FillVoxel(1, site_tree, height_tree, CR_tree, CD_tree);
            //chm_field = chm_field_save;
    
            SetCanopyHeight(chm_patch_tree, site_tree, extent_tree);
            SetCanopyHeight(chm_patch_treeshifted, site_tree_shifted, extent_tree);

            hist_CHM_sim = hist_CHM_sim_save;
            sum_abserror = abserror_previous_test;
            CA_exterior = CA_exterior_previous;
            if(carbonstarv_flag == 1){
                ReverseLeaves(list_trees, list_trees_properties, site_tree, site_tree_shifted);
                carbonstarv = carbonstarv_previous_test;
                nbtrees_carbonstarv = nbtrees_carbonstarv_previous_test;
            }
        }
    }
}
#ifdef exchangepositions_step2
void FindBetterPositions(int site_tree, int mae_flag, int dissimilarity_flag, int carbonstarv_flag){
    nbattempt++;
    
    /* this is a function to move a tree to another site, depending on whether that improves the canopy */
    float height_tree = T[site_tree].t_Tree_Height;
    float CR_tree = T[site_tree].t_Crown_Radius;
    float CD_tree = T[site_tree].t_Crown_Depth;
    float LMA_tree = T[site_tree].t_LMA;
    float Nmass_tree = T[site_tree].t_Nmass;
    float Pmass_tree = T[site_tree].t_Pmass;
    
    float distance_max = height_tree*height_tree/10.0;
    //float distance_max = exp(height_tree/10.0);
    distance_max = maxf(distance_max, 5.0);
    distance_max = minf(distance_max, 200.0);
    
    float distance = gsl_rng_uniform(gslrand) * distance_max;
    float angle_move = float(twoPi*gsl_rng_uniform(gslrand));                                                /* angle in which tree is moved */

    int dist_cols = int(distance*cos(angle_move));
    int dist_rows = int(distance*sin(angle_move));
               
    int col_tree = site_tree%cols;
    int row_tree = site_tree/cols;
    int col_shifted = dist_cols + col_tree;
    int row_shifted = dist_rows + row_tree;
    
    /* we ensure that new site falls within simulated area via periodic boundary conditions */
    if(col_shifted < 0) col_shifted = cols - (abs(col_shifted)%cols);
    else col_shifted = col_shifted%cols;
    if(row_shifted < 0) row_shifted = rows - (abs(row_shifted)%rows);
    else row_shifted = row_shifted%rows;
    
    int site_tree_shifted = col_shifted + row_shifted * cols;
    
    if(T[site_tree_shifted].t_age == 0){
        /* the radius for a CHM update (needs to be the largest possible radius, including possible crown displacement) */
        float radius_update =  CR_tree;
        int extent_tree = int(radius_update) + 2;
        
        //vector<int> chm_field_save(chm_field);
        vector<int> chm_patch_tree;
        chm_patch_tree.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_tree, site_tree, extent_tree);
        vector<int> chm_patch_treeshifted;
        chm_patch_treeshifted.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_treeshifted, site_tree_shifted, extent_tree);
        vector<int> hist_CHM_sim_save(hist_CHM_sim);
        
        float dissimilarityprevious_test = CalcDissimilarity(1);
        int abserror_previous_test = sum_abserror;
        float carbonstarv_previous_test = carbonstarv;
        int nbtrees_carbonstarv_previous_test = nbtrees_carbonstarv;

        float meanheight_previous = CanopyMeanheight();
        float CA_exterior_previous = CA_exterior;
        
        UpdateCanopyNewPosition(site_tree, site_tree_shifted, height_tree, CR_tree, CD_tree, extent_tree);
        //UpdateCanopyAroundTree(site_tree, site_tree_shifted, height_tree, height_tree, CR_tree, CR_tree, CD_tree, CD_tree, extent_tree);
        
        float dissimilaritynew_test = CalcDissimilarity(1);
        int abserror_new_test = sum_abserror;
        
        /* update the LAI field */
        vector<int> list_trees;
        vector<float> list_trees_properties;
        
        /* since it is computationally intensive, we only compute this when it is needed */
        if(carbonstarv_flag == 1){
            GetTreeList(site_tree, site_tree_shifted, list_trees);
            
            list_trees_properties.reserve(4 * list_trees.size());
            UpdateLeaves_fromList(list_trees, list_trees_properties, site_tree, site_tree_shifted);
        }
        
        float carbonstarv_new_test = carbonstarv;
        
        float meanheight_new = CanopyMeanheight();
        float CA_exterior_new = CA_exterior;
        
        int move = EvaluateCanopy(mae_flag, dissimilarity_flag, carbonstarv_flag, dissimilarityprevious_test, dissimilaritynew_test, abserror_previous_test, abserror_new_test, meanheight_previous, meanheight_new, CA_exterior_previous, CA_exterior_new, carbonstarv_previous_test, carbonstarv_new_test);
        
        if(flag_PreventCrownpiercing == 1 && move == 1){
            int ingrowth_tree = TreeIngrowth(site_tree_shifted, height_tree, CR_tree);
            if(ingrowth_tree == 1) move = 0;
        }

        if(move == 1){
            nbsuccess++;
            
            float dbh = T[site_tree].t_dbh;
            float wsg = T[site_tree].t_wsg;
            int species_label = T[site_tree].t_species_label;
            
            T[site_tree_shifted].InitialiseTree(site_tree_shifted, dbh, height_tree, CR_tree, CD_tree, Nmass_tree, Pmass_tree, LMA_tree, wsg, species_label);
            if(carbonstarv_flag == 1){
                for(int i=0; i < list_trees.size(); i++){
                    int site_tree_update = list_trees[i];
                    if(site_tree_update == site_tree) site_tree_update = site_tree_shifted;
                    T[site_tree_update].t_LAIabove = list_trees_properties[0 + i*4];
                    T[site_tree_update].t_LAI = list_trees_properties[1 + i*4];
                    T[site_tree_update].t_GPP = list_trees_properties[2 + i*4];
                    T[site_tree_update].t_NPP = list_trees_properties[3 + i*4];
                }
            }
            
            T[site_tree].SetZero(site_tree);
        } else {
            FillVoxel(-1, site_tree_shifted, height_tree, CR_tree, CD_tree);
            FillVoxel(1, site_tree, height_tree, CR_tree, CD_tree);
            //chm_field = chm_field_save;
    
            SetCanopyHeight(chm_patch_tree, site_tree, extent_tree);
            SetCanopyHeight(chm_patch_treeshifted, site_tree_shifted, extent_tree);

            hist_CHM_sim = hist_CHM_sim_save;
            sum_abserror = abserror_previous_test;
            CA_exterior = CA_exterior_previous;
            if(carbonstarv_flag == 1){
                ReverseLeaves(list_trees, list_trees_properties, site_tree, site_tree_shifted);
                carbonstarv = carbonstarv_previous_test;
                nbtrees_carbonstarv = nbtrees_carbonstarv_previous_test;
            }
        }
    } else {
        /* this is a function to move a tree to another site, depending on whether that improves the canopy */
        int site_treeexchange = site_tree_shifted;

        float height_treeexchange = T[site_treeexchange].t_Tree_Height;
        float CR_treeexchange = T[site_treeexchange].t_Crown_Radius;
        float CD_treeexchange = T[site_treeexchange].t_Crown_Depth;
        float LMA_treeexchange = T[site_treeexchange].t_LMA;
        float Nmass_treeexchange = T[site_treeexchange].t_Nmass;
        float Pmass_treeexchange = T[site_treeexchange].t_Pmass;

        int crowndisplacement_tree = T[site_tree].t_CrownDisplacement;
        int crowndisplacement_treeexchange = T[site_treeexchange].t_CrownDisplacement;

        /* now we test whether moving the tree creates a better canopy */
        int site_crown = site_tree + crowndisplacement_tree;
        int site_crown_updated = site_tree + crowndisplacement_treeexchange;
        int site_crownexchange = site_treeexchange + crowndisplacement_treeexchange;
        int site_crownexchange_updated = site_treeexchange + crowndisplacement_tree;

        /* now we update the statistics, remove and add */
        /* this has to be done in succession, in case the two tree crowns overlap */

        int extent_tree = GetTreeExtent(site_crown, site_crown_updated, CR_tree, CR_treeexchange);
        int extent_treeexchange = GetTreeExtent(site_crownexchange, site_crownexchange_updated, CR_treeexchange, CR_tree);

        vector<int> chm_patch_tree;
        chm_patch_tree.reserve(4 * extent_tree * extent_tree);
        GetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
        vector<int> chm_patch_treeexchange;
        chm_patch_treeexchange.reserve(4 * extent_treeexchange * extent_treeexchange);
        GetCanopyHeight(chm_patch_treeexchange, site_crownexchange, extent_treeexchange);
        vector<int> hist_CHM_sim_save(hist_CHM_sim);

        /* then update the fields */
        float dissimilarityprevious_test = CalcDissimilarity(1);
        int abserror_previous_test = sum_abserror;
        float carbonstarv_previous_test = carbonstarv;
        int nbtrees_carbonstarv_previous_test = nbtrees_carbonstarv;

        float meanheight_previous = CanopyMeanheight();
        float CA_exterior_previous = CA_exterior;

        UpdateCanopyAroundTree(site_crown, site_crown_updated, height_tree, height_treeexchange, CR_tree, CR_treeexchange, CD_tree, CD_treeexchange, extent_tree);
        UpdateCanopyAroundTree(site_crownexchange, site_crownexchange_updated, height_treeexchange, height_tree, CR_treeexchange, CR_tree, CD_treeexchange, CD_tree, extent_treeexchange);

        /* calculate the new theoretical canopy statistics */
        float dissimilaritynew_test = CalcDissimilarity(1);
        int abserror_new_test = sum_abserror;

        /* update the LAI field */
        vector<int> list_trees;
        vector<float> list_trees_properties;

        /* since it is computationally intensive, we only compute this when it is needed */
        if(carbonstarv_flag == 1){
            GetTreeList(site_tree, site_treeexchange, CR_tree, CR_treeexchange, CR_treeexchange, CR_tree, CD_tree, CD_treeexchange, CD_treeexchange, CD_tree, height_tree, height_treeexchange, height_treeexchange, height_tree, crowndisplacement_tree, crowndisplacement_treeexchange, crowndisplacement_treeexchange, crowndisplacement_tree, list_trees);

            list_trees_properties.reserve(4 * list_trees.size());

            UpdateLeaves_fromList(list_trees, list_trees_properties, site_tree, height_treeexchange, CR_treeexchange, CD_treeexchange, crowndisplacement_treeexchange, Nmass_treeexchange, Pmass_treeexchange, LMA_treeexchange, site_treeexchange, height_tree, CR_tree, CD_tree, crowndisplacement_tree, Nmass_tree, Pmass_tree, LMA_tree);
        }

        float carbonstarv_new_test = carbonstarv;

        float meanheight_new = CanopyMeanheight();
        float CA_exterior_new = CA_exterior;

        /* now we test whether the new configuration is better than the previous one */
        int betterpositions = EvaluateCanopy(mae_flag, dissimilarity_flag, carbonstarv_flag, dissimilarityprevious_test, dissimilaritynew_test, abserror_previous_test, abserror_new_test, meanheight_previous, meanheight_new, CA_exterior_previous, CA_exterior_new, carbonstarv_previous_test, carbonstarv_new_test);

        if(flag_PreventCrownpiercing == 1 && betterpositions == 1){
            int ingrowth_tree = TreeIngrowth(site_crown_updated, height_treeexchange, CR_treeexchange);
            int ingrowth_treeexchange = TreeIngrowth(site_crownexchange_updated, height_tree, CR_tree);
            if(ingrowth_tree == 1 || ingrowth_treeexchange == 1) betterpositions = 0;
        }

        if(betterpositions == 1){
            //if(abserror_new_test < abserror_previous_test){
            /* update the success rate */
            nbsuccess++;

            float dbh_tree = T[site_tree].t_dbh;
            float wsg_tree = T[site_tree].t_wsg;
            int species_label_tree = T[site_tree].t_species_label;

            float dbh_treeexchange = T[site_treeexchange].t_dbh;
            float wsg_treeexchange = T[site_treeexchange].t_wsg;
            int species_label_treeexchange = T[site_treeexchange].t_species_label;

            T[site_treeexchange].InitialiseTree(site_treeexchange, dbh_tree, height_tree, CR_tree, CD_tree, Nmass_tree, Pmass_tree, LMA_tree, wsg_tree, species_label_tree);
            T[site_tree].InitialiseTree(site_tree, dbh_treeexchange, height_treeexchange, CR_treeexchange, CD_treeexchange, Nmass_treeexchange, Pmass_treeexchange, LMA_treeexchange, wsg_treeexchange, species_label_treeexchange);

//            float dissimilarity_test = 0.0;
//            int  sum_abserror_test = 0;
//            //UpdateCHM();
//            CanopyDistanceFull(dissimilarity_test, sum_abserror_test);
//
//            if(abserror_new_test != sum_abserror_test) cout << site_tree << " dbh: " << dbh_tree << " site_exchange: " << dbh_treeexchange << " Dissimilarity new: " << dissimilaritynew_test << " full recalc: " << dissimilarity_test << " abserror_new: " << abserror_new_test << " full recalc: " << sum_abserror_test << endl;
//
//

            if(carbonstarv_flag == 1){
                for(int i=0; i < list_trees.size(); i++){
                    int site_tree_update = list_trees[i];
                    T[site_tree_update].t_LAIabove = list_trees_properties[0 + i*4];
                    T[site_tree_update].t_LAI = list_trees_properties[1 + i*4];
                    T[site_tree_update].t_GPP = list_trees_properties[2 + i*4];
                    T[site_tree_update].t_NPP = list_trees_properties[3 + i*4];
                }
            }
        } else {
            /* otherwise, undo the changes */
            FillVoxel(-1, site_crown_updated, height_treeexchange, CR_treeexchange, CD_treeexchange);
            FillVoxel(1, site_crown, height_tree, CR_tree, CD_tree);
            FillVoxel(-1, site_crownexchange_updated, height_tree, CR_tree, CD_tree);
            FillVoxel(1, site_crownexchange, height_treeexchange, CR_treeexchange, CD_treeexchange);

            SetCanopyHeight(chm_patch_tree, site_crown, extent_tree);
            SetCanopyHeight(chm_patch_treeexchange, site_crownexchange, extent_treeexchange);
            hist_CHM_sim = hist_CHM_sim_save;
            sum_abserror = abserror_previous_test;
            CA_exterior = CA_exterior_previous;

            if(carbonstarv_flag == 1){
                ReverseLeaves(list_trees, list_trees_properties, site_tree, height_treeexchange, CR_treeexchange, CD_treeexchange, crowndisplacement_treeexchange, site_treeexchange, height_tree, CR_tree, CD_tree, crowndisplacement_tree);
                carbonstarv = carbonstarv_previous_test;
                nbtrees_carbonstarv = nbtrees_carbonstarv_previous_test;
            }
        }
    }
}
#endif


template <typename N>
void SetParameter(string &parameter_name, string &parameter_value, N &parameter, N parameter_min, N parameter_max, N parameter_default, bool quiet) {
    /* idea for checking whether float from https://stackoverflow.com/questions/447206/c-isfloat-function */
    istringstream iss(parameter_value);
    N numeric;
    iss >> numeric;
    
    // Check the entire string was consumed and if either failbit or badbit is set
    bool isfloat = iss.eof() && !iss.fail();
    if(isfloat){
        if(numeric >= parameter_min && numeric <= parameter_max){
            parameter = numeric;
            if(!quiet) cout << parameter_name << ": " << parameter << endl;
        } else {
            parameter = parameter_default;
            if(!quiet) cout << "Warning. Value provided for '" << parameter_name << "' (" << numeric << ") is outside the allowed range (" << parameter_min << ", " << parameter_max << "). Set to default: " << parameter_default << endl;
        }
    } else {
        parameter = parameter_default;
        if(!quiet) cout << "Warning. Value provided for '" << parameter_name << "' (" << parameter_value << ") is not a " << typeid(numeric).name() << ". Set to default: " << parameter_default << endl;
    }
}

void AssignValueGlobal(string parameter_name, string parameter_value){
    /* we set parameters to values that have been read, or to their defaults, if outside of range or not the right type */
    bool quiet = 0;
    if(parameter_name == "cols"){
        SetParameter(parameter_name, parameter_value, cols, 0, 10000, 400, quiet);
    } else if(parameter_name == "rows"){
        SetParameter(parameter_name, parameter_value, rows, 0, 10000, 400, quiet);
    } else if(parameter_name == "height_max"){
        SetParameter(parameter_name, parameter_value, height_max, 0, 150, 80, quiet);
    } else if(parameter_name == "seed_rng"){
        SetParameter(parameter_name, parameter_value, seed_rng, 0, INT_MAX, 0, quiet);
        if(seed_rng < 0){
            int t = (int) time(NULL);
            seed_rng = 3*t+1;
            cout << "For seed_rng < 0, random initialization is assumed. We update seed to: " << seed_rng << endl;
        }
    } else if(parameter_name == "dbh_limitbottom"){
        SetParameter(parameter_name, parameter_value, dbh_limitbottom, float(0.01), float(35.0), float(0.01), quiet);
    } else if(parameter_name == "dbh_limittop"){
        SetParameter(parameter_name, parameter_value, dbh_limittop, float(0.01), float(35.0), float(1.5), quiet);
    } else if(parameter_name == "dbhlog_binsize"){
        SetParameter(parameter_name, parameter_value, dbhlog_binsize, float(0.0), float(1.0), float(0.025), quiet);
    } else if(parameter_name == "dbh_limitcarbon"){
        SetParameter(parameter_name, parameter_value, dbh_limitcarbon, float(0.01), float(35.0), float(0.1), quiet);
    } else if(parameter_name == "dbh_limitstoppingrule"){
        SetParameter(parameter_name, parameter_value, dbh_limitstoppingrule, float(0.01), float(35.0), float(0.1), quiet);
    } else if(parameter_name == "stoppingrule"){
        SetParameter(parameter_name, parameter_value, stoppingrule, float(0.0001), float(1.0), float(0.01), quiet);
    } else if(parameter_name == "nbsteps_dissimilarity"){
        SetParameter(parameter_name, parameter_value, nbsteps_dissimilarity, 0, 1000, 10, quiet);
    } else if(parameter_name == "nbsteps_mae"){
        SetParameter(parameter_name, parameter_value, nbsteps_mae, 0, 1000, 90, quiet);
    } else if(parameter_name == "nbsteps_carbonbalance"){
        SetParameter(parameter_name, parameter_value, nbsteps_carbonbalance, 0, 1000, 0, quiet);
    } else if(parameter_name == "nbsteps_combined"){
        SetParameter(parameter_name, parameter_value, nbsteps_combined, 0, 1000, 50, quiet);
    } else if(parameter_name == "crowndisplacement_factor"){
        SetParameter(parameter_name, parameter_value, crowndisplacement_factor, float(0.0), float(1.0), float(0.0), quiet);
    } else if(parameter_name == "flag_ApplyMedianfilter"){
        SetParameter(parameter_name, parameter_value, flag_ApplyMedianfilter, 0, 1, 1, quiet);
    } else if(parameter_name == "flag_Prefitting"){
        SetParameter(parameter_name, parameter_value, flag_Prefitting, 0, 1, 1, quiet);
    } else if(parameter_name == "flag_PreventCrownpiercing"){
        SetParameter(parameter_name, parameter_value, flag_PreventCrownpiercing, 0, 1, 0, quiet);
    } else if(parameter_name == "flag_LAIgradient"){
        SetParameter(parameter_name, parameter_value, flag_LAIgradient, 0, 1, 0, quiet);
    } else if(parameter_name == "flag_OutputReduced"){
        SetParameter(parameter_name, parameter_value, flag_OutputReduced, 0, 1, 0, quiet);
    } else if(parameter_name == "flag_OutputTROLL"){
        SetParameter(parameter_name, parameter_value, flag_OutputTROLL, 0, 1, 0, quiet);
    }
}

void AssignValueDetail(string parameter_name, string parameter_value, int parameterset){
    /* we set parameters to values that have been read, or to their defaults, if outside of range or not the right type */
    bool quiet = parameterset; // defaults to 0 for the first one, so output will be provided, otherwise quiet
    if(parameter_name == "paramID"){
        int paramIDtest;
        SetParameter(parameter_name, parameter_value, paramIDtest, 0, INT_MAX, parameterset, quiet);
        parameters_detailed[0 + parameterset * nb_parametersdetailed] = float(paramIDtest);
    } else if(parameter_name == "a_sdd"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[1 + parameterset * nb_parametersdetailed], -10000.0f, 10000.0f, -1.8f, quiet);
    } else if(parameter_name == "b_sdd"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[2 + parameterset * nb_parametersdetailed], -10000.0f, 10000.0f, -6.0f, quiet);
    } else if(parameter_name == "nbavg_sdd"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[3 + parameterset * nb_parametersdetailed], 0.0f, 10000.0f, 4000.0f, quiet);
    } else if(parameter_name == "a_height"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[4 + parameterset * nb_parametersdetailed], 0.0f, 10.0f, 0.3f, quiet);
    } else if(parameter_name == "b_height"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[5 + parameterset * nb_parametersdetailed], 0.0f, 1000.0f, 50.0f, quiet);
    } else if(parameter_name == "sigma_height"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[6 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.2f, quiet);
    } else if(parameter_name == "hetscedas_height"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[7 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 1.0f, quiet);
    } else if(parameter_name == "a_CR"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[8 + parameterset * nb_parametersdetailed], 0.0f, 10.0f, 2.3f, quiet);
    } else if(parameter_name == "b_CR"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[9 + parameterset * nb_parametersdetailed], 0.0f, 5.0f, 0.6f, quiet);
    } else if(parameter_name == "sigma_CR"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[10 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.25f, quiet);
    } else if(parameter_name == "hetscedas_CR"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[11 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 1.0f, quiet);
    } else if(parameter_name == "a_CD"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[12 + parameterset * nb_parametersdetailed], 0.0f, 0.0f, 0.0f, quiet);
    } else if(parameter_name == "b_CD"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[13 + parameterset * nb_parametersdetailed], -1.0f, 1.0f, 0.2f, quiet);
    } else if(parameter_name == "sigma_CD"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[14 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.0f, quiet);
    } else if(parameter_name == "hetscedas_CD"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[15 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 1.0f, quiet);
    } else if(parameter_name == "LMA"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[16 + parameterset * nb_parametersdetailed], 0.0f, 10000.0f, 95.0f, quiet);
    } else if(parameter_name == "sigma_LMA"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[17 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.35f, quiet);
    } else if(parameter_name == "Nmass"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[18 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.02f, quiet);
    } else if(parameter_name == "sigma_Nmass"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[19 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.25f, quiet);
    } else if(parameter_name == "Pmass"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[20 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.00055f, quiet);
    } else if(parameter_name == "sigma_Pmass"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[21 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.4f, quiet);
    } else if(parameter_name == "wsg"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[22 + parameterset * nb_parametersdetailed], 0.0f, 1.5f, 0.65f, quiet);
    } else if(parameter_name == "sigma_wsg"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[23 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.1f, quiet);
    } else if(parameter_name == "corr_Nmass_Pmass"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[24 + parameterset * nb_parametersdetailed], -1.0f, 1.0f, 0.65f, quiet);
    } else if(parameter_name == "corr_Nmass_LMA"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[25 + parameterset * nb_parametersdetailed], -1.0f, 1.0f, -0.43f, quiet);
    } else if(parameter_name == "corr_Pmass_LMA"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[26 + parameterset * nb_parametersdetailed], -1.0f, 1.0f, -0.39f, quiet);
    } else if(parameter_name == "shape_crown"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[27 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.8f, quiet);
    } else if(parameter_name == "gapfraction_crown"){
        SetParameter(parameter_name, parameter_value, parameters_detailed[28 + parameterset * nb_parametersdetailed], 0.0f, 1.0f, 0.0f, quiet);
    }
}

void InitialiseGlobal(int flag_global, int &error_global){
    /* Error handling: There are no fatal errors. If file cannot be read, default values are used. If parameters are lacking from the file, default values are used */
    /* reading input files */
    
    //vector<string> parameter_names{"cols","rows","height_max","seed_rng","dbh_limitbottom","dbh_limittop","dbhlog_binsize","dbh_limitcarbon","dbh_limitstoppingrule","stoppingrule","nbsteps_dissimilarity","nbsteps_mae","nbsteps_carbonbalance","nbsteps_combined","crowndisplacement_factor","flag_ApplyMedianfilter","flag_PreventCrownpiercing","flag_LAIgradient","flag_OutputReduced","flag_OutputTROLL"};
    // int nb_parameters = int(parameter_names.size());
    string parameter_names[21] = {"cols","rows","height_max","seed_rng","dbh_limitbottom","dbh_limittop","dbhlog_binsize","dbh_limitcarbon","dbh_limitstoppingrule","stoppingrule","nbsteps_dissimilarity","nbsteps_mae","nbsteps_carbonbalance","nbsteps_combined","crowndisplacement_factor","flag_ApplyMedianfilter","flag_Prefitting","flag_PreventCrownpiercing","flag_LAIgradient","flag_OutputReduced","flag_OutputTROLL"};
    int nb_parameters = 21;
    
    vector<string> parameter_values(nb_parameters,"");
    /* we update the values from the file */
    if (flag_global == 1) {
        if(flag_global) sprintf(inputfile_global,"%s",bufi_global);
        fstream InGlobal(inputfile_global, ios::in);
        /* first we initialise the potential parameter arrays with empty parameter values */
        
        if(InGlobal){
            error_global = 0;
            cout << "Reading in file: " << inputfile_global << endl;
            InGlobal.getline(buffer,256,'\n');
            string parameter_name, parameter_value;
            
            while (InGlobal >> parameter_name >> parameter_value){
                InGlobal.getline(buffer,256,'\n');
                for(int i = 0; i < nb_parameters; i++){
                    if(parameter_name == parameter_names[i]) parameter_values[i] = parameter_value;
                }
            }
        } else {
            error_global = 2;
            cout << "ERROR. General input file could not be read. Canopy Constructor will exit." << endl;
        }
        InGlobal.close();
    } else {
        error_global = 1;
        cout << "WARNING. No general input file provided. All global parameters will be set to default." << endl;
    }
    
    if(error_global < 2){
        /* now we assign values */
        for(int i = 0; i < nb_parameters; i++){
            AssignValueGlobal(parameter_names[i], parameter_values[i]);
        }

        cout << endl;
        
        /* now we write the parameters out */
        output_input_global << "Parameter\tValue" << endl;
        output_input_global << "cols\t" << cols << endl;
        output_input_global << "rows\t" << rows << endl;
        output_input_global << "height_max\t" << height_max << endl;
        output_input_global << "seed_rng\t" << seed_rng << endl;
        output_input_global << "dbh_limitbottom\t" << dbh_limitbottom << endl;
        output_input_global << "dbh_limittop\t" << dbh_limittop << endl;
        output_input_global << "dbhlog_binsize\t" << dbhlog_binsize << endl;
        output_input_global << "dbh_limitcarbon\t" << dbh_limitcarbon << endl;
        output_input_global << "dbh_limitstoppingrule\t" << dbh_limitstoppingrule << endl;
        output_input_global << "stoppingrule\t" << stoppingrule << endl;
        output_input_global << "nbsteps_dissimilarity\t" << nbsteps_dissimilarity << endl;
        output_input_global << "nbsteps_mae\t" << nbsteps_mae << endl;
        output_input_global << "nbsteps_carbonbalance\t" << nbsteps_carbonbalance << endl;
        output_input_global << "nbsteps_combined\t" << nbsteps_combined << endl;
        output_input_global << "crowndisplacement_factor\t" << crowndisplacement_factor << endl;
        output_input_global << "flag_ApplyMedianfilter\t" << flag_ApplyMedianfilter << endl;
        output_input_global << "flag_Prefitting\t" << flag_Prefitting << endl;
        output_input_global << "flag_PreventCrownpiercing\t" << flag_PreventCrownpiercing << endl;
        output_input_global << "flag_LAIgradient\t" << flag_LAIgradient << endl;
        output_input_global << "flag_OutputReduced\t" << flag_OutputReduced << endl;
        output_input_global << "flag_OutputTROLL\t" << flag_OutputTROLL << endl;
    
        /* update dependent parameters */
        sites = rows*cols;
        
        cout << "Simulation initialised with extension of: " << cols << " x " << rows << " | maximum height: " << height_max << endl;

        dbhlog_binnb = int(log10(dbh_limittop*100.0)/dbhlog_binsize + 1.0);
        cout << "Diameter binning for tree swapping" << endl;
        cout << "Logbin nb: " << dbhlog_binnb  << " Logbin size: " << dbhlog_binsize <<  " lower dbh_limit: " << dbh_limitbottom << " upper dbh_limit: " << dbh_limittop << " upper dbh_limit (log) " << log10(dbh_limittop*100.0) << endl;
        
        /* read in physiology information */
        klight = 0.5;
        kpar = klight * 0.9;    /* kpar is the klight factor times the absorptance of leaves */
        theta = 0.70; /* This could be added in the input file, just like kpar */
        phi = 0.093;
        alpha = 4.0 * phi;
        g1 = 3.77;
        Cair = 400;
        iCair = 1.0/Cair;
    }
}

void InitialiseCHM(int flag_chm, int &error_chm){
    /* Error handling: error == 0 means no error, error == 1 means nonfatal error, i.e. no file is provided and default options are used, error == 2 is fatal error (file was inconsistent with other parameter/data sets provided) */

    if(flag_chm == 1){
        if(flag_chm) sprintf(inputfile_chm,"%s",bufi_chm);
        fstream InCHM(inputfile_chm, ios::in);
        if(InCHM){
            error_chm = 0;
            cout << "Reading in file: " << inputfile_chm << endl;
            
            string line;
            
            getline(InCHM,line);
            istringstream firstlinestream(line);

            //vector<string> variable_names{"x","y","z"};
            //int nb_variables = int(variable_names.size());
            string variable_names[3] = {"x","y","z"};
            int nb_variables = 3;
            
            vector<int> variable_index;
            
            string variable_name;
            while(firstlinestream >> variable_name){
                int index = -1;
                for(int i = 0; i < nb_variables; i++){
                    if(variable_name == variable_names[i]) index = i;
                }
                if(index == -1){
                    cout << "Ignoring unknown variable: " << variable_name << ". Please provide dimensions as x, y, z." << endl;
                }
                variable_index.push_back(index);
            }
        
            mincol_absolute = INT_MAX;
            minrow_absolute = INT_MAX;
            
            int height_max_empirical = 0;
            int cols_empirical = 0, rows_empirical = 0;
            
            vector<int> CHMflat;
            CHMflat.reserve(sites * 3);

            int linenb = 0;
            int linenb_error = 0;
            
            /* first iteration just checks whether the file is good */
            while(getline(InCHM, line) && linenb_error == 0){
                istringstream linestream(line);
                linenb++;
                
                string variable_value;
                
                int col = -1, row = -1;
                float height_unrounded = -1.0;
                
                int v = 0;
                while(linestream >> variable_value){
                    bool quiet;
                    if(linenb < 5) quiet = 0;
                    else quiet = 1;
         
                    int index = variable_index[v];
                    
                    if(index >= 0){
                        string variable_name = variable_names[index];
                        if(variable_name == "x") SetParameter(variable_name, variable_value, col, 0, INT_MAX, -1, quiet);
                        else if(variable_name == "y") SetParameter(variable_name, variable_value, row, 0, INT_MAX, -1, quiet);
                        else if(variable_name == "z") SetParameter(variable_name, variable_value, height_unrounded, 0.0f, float(height_max), -1.0f, quiet);
                    }
                    v++;
                }
            
                if(col >= 0 && row >= 0 && height_unrounded >= 0.0){
                    if(col < mincol_absolute) mincol_absolute = col;
                    if(col > cols_empirical) cols_empirical = col;

                    if(row < minrow_absolute) minrow_absolute = row;
                    if(row > rows_empirical) rows_empirical = row;
                       
                    int height_rounded = int(lround(height_unrounded));
                    if(height_rounded > height_max_empirical) height_max_empirical = height_rounded;
                    
                    CHMflat.push_back(col);
                    CHMflat.push_back(row);
                    CHMflat.push_back(height_rounded);
                } else {
                    if(linenb_error == 0) linenb_error = linenb;
                }
            }
            
            InCHM.close();

            cols_empirical -= mincol_absolute;
            cols_empirical += 1;
            
            rows_empirical -= minrow_absolute;
            rows_empirical += 1;
            
            if(linenb == 0){
                error_chm = 1;
                cout << "WARNING. CHM input was empty. No fitting is performed." << endl;
            }
            
            if(linenb_error > 0){
                error_chm = 2;
                cout << "ERROR. Line: " << linenb_error << " of CHM file could not be read in. Canopy Constructor will exit." << endl;
            }
            if(height_max < height_max_empirical){
                error_chm = 2;
                cout << "ERROR. Empirical maximum height is: " << height_max_empirical << ", which is greater than initialised/allowed maximum height of: " << height_max << ". Canopy Constructor will exit." << endl;
            }
            if(cols_empirical != cols || rows_empirical != rows || CHMflat.size()/3 != sites){
                error_chm = 2;
                cout << "ERROR. empirical cols and rows values (" << cols_empirical << " | " << rows_empirical << ") and site number (" << CHMflat.size()/3 << ") do not correspond to initialised/allowed values (" << cols << " | " << rows << " and " << sites << ")  Canopy Constructor will exit." << endl;
            }
            
            if(error_chm == 0){
                cout << "Successfully read in empirical CHM" << endl;
                
                cout << "Minimum absolute coordinates (row/col): " << minrow_absolute << " | " << mincol_absolute << endl;
                
                for(int s=0;s<sites;s++){
                    int col = CHMflat[0 + s * 3];
                    int row = CHMflat[1 + s * 3];
                    int height = CHMflat[2 + s * 3];
                    
                    /* now subtract the minimum values to normalize to zero */
                    col -= mincol_absolute;
                    row -= minrow_absolute;
                    
                    int site = col + row * cols;
                    
                    chm_empirical[site] = height;
                }
                
                if(flag_ApplyMedianfilter){
                     /* we implement a naive version, relying on sorting the neighborhood */
                     vector<int> chm_empirical_filtered(sites,0);
                     
                     for(int r = 0; r < rows; r++){
                         for(int c = 0; c < cols; c++){
                             int site = c + r * cols;
                             vector<int> neighborhood;
                             neighborhood.reserve(5);
                             
                             neighborhood.push_back(chm_empirical[site]);
                             
                             /* use periodic boundary conditions */
                             int neighbor_left;
                             if(c > 0) neighbor_left = c - 1 + r * cols;
                             else neighbor_left = cols - 1 + r * cols;
                             
                             int neighbor_right;
                             if(c < cols - 1) neighbor_right = c + 1 + r * cols;
                             else neighbor_right = r * cols; /* + 0 */
                             
                             int neighbor_below;
                             if(r > 0) neighbor_below = c + (r - 1) * cols;
                             else neighbor_below = c + (rows - 1) * cols;
                             
                             int neighbor_above;
                             if(r < rows - 1) neighbor_above = c + (r + 1) * cols;
                             else neighbor_above = c;    /* + 0 * cols */
                             
                             neighborhood.push_back(chm_empirical[neighbor_left]);
                             neighborhood.push_back(chm_empirical[neighbor_right]);
                             neighborhood.push_back(chm_empirical[neighbor_below]);
                             neighborhood.push_back(chm_empirical[neighbor_above]);
                             
                             /* now sort and pick the middle entry */
                             /* here we do not need to find middle, because we work with 5 values */
                             sort(neighborhood.begin(), neighborhood.end());
                             int median_neighborhood = neighborhood[2];
                             
                             chm_empirical_filtered[site] = median_neighborhood;
                         }
                     }
                     
                     for(int site = 0; site < sites; site++){
                         chm_empirical[site] = chm_empirical_filtered[site];
                     }
                     
                     cout << "MEDIAN filter applied to empirical CHM (von Neumann neighborhood)" << endl;
                 }
                 
                 cout << "Calculate height distribution." << endl;
                 
                 for(int s=0;s<sites;s++){
                     int height_canopy = min(chm_empirical[s],height_max);
                     hist_CHM_emp[height_canopy]++;
                 }

                 cout << "Calculate complementary cumulative height distribution." << endl;
                 
                 int pixel_cumulative = 0;
                 for(int h=height_max; h >= 0; h--){
                     pixel_cumulative += hist_CHM_emp[h];
                     hist_CHM_empcompcumulative[h] += float(pixel_cumulative)/float(sites);
                 }
            }
        } else {
            error_chm = 2;
            cout << "ERROR. CHM input could not be read. Canopy Constructor will exit." << endl;
        }
    } else {
        error_chm = 1;
        cout << "WARNING. CHM input was not found. No fitting is performed." << endl;
    }
}
 

void InitialiseOutput(int flag_output, int &error_output){
    /* Error handling: There are no fatal errors. If output name could not be read, default values are used. If parameters are lacking from the file, default values are used */
    /*** Initialization of output streams ***/
    /****************************************/
    
    if(flag_output == 1){
        error_output = 0;
        char nnn[200];

        sprintf(nnn,"%s_input_global.txt",bufo);
        output_input_global.open(nnn, ios::out);

        sprintf(nnn,"%s_input_detailed.txt",bufo);
        output_input_detailed.open(nnn, ios::out);

        sprintf(nnn,"%s_input_CHM.txt",bufo);
        output_input_CHM.open(nnn, ios::out);

        sprintf(nnn,"%s_input_inventory.txt",bufo);
        output_input_inventory.open(nnn, ios::out);
        
        sprintf(nnn,"%s_input_inventorycutoff.txt",bufo);
        output_input_inventorycutoff.open(nnn, ios::out);

        sprintf(nnn,"%s_input_cp.txt",bufo);
        output_input_cp.open(nnn, ios::out);
        
        sprintf(nnn,"%s_input_cpcutoff.txt",bufo);
        output_input_cpcutoff.open(nnn, ios::out);

        sprintf(nnn,"%s_sumstat.txt",bufo);
        output_sumstat.open(nnn, ios::out);

        sprintf(nnn,"%s_agbquarterha.txt",bufo);
        output_agbquarterha.open(nnn, ios::out);

        sprintf(nnn,"%s_sdd.txt",bufo);
        output_sdd.open(nnn, ios::out);

        sprintf(nnn,"%s_chm.txt",bufo);
        output_chm.open(nnn, ios::out);

        if(!(flag_OutputReduced == 1)){
            sprintf(nnn,"%s_field.txt",bufo);
            output_field.open(nnn, ios::out);

            sprintf(nnn,"%s_LAI3D.txt",bufo);
            output_LAI3D.open(nnn, ios::out);

            sprintf(nnn,"%s_vertical.txt",bufo);
            output_vertical.open(nnn, ios::out);

            sprintf(nnn,"%s_verticalfromtop.txt",bufo);
            output_verticalfromtop.open(nnn, ios::out);

            sprintf(nnn,"%s_vertical_ALS.txt",bufo);
            output_vertical_ALS.open(nnn, ios::out);

            sprintf(nnn,"%s_verticalfromtop_ALS.txt",bufo);
            output_verticalfromtop_ALS.open(nnn, ios::out);

            sprintf(nnn,"%s_troll.txt",bufo);
            output_troll.open(nnn, ios::out);

            sprintf(nnn,"%s_trees.txt",bufo);
            output_trees.open(nnn, ios::out);

            sprintf(nnn,"%s_twodim.txt",bufo);
            output_twodim.open(nnn, ios::out);

            sprintf(nnn,"%s_volume.txt",bufo);
            output_volume.open(nnn, ios::out);
        }
    } else {
        error_output = 2;
        cout << "ERROR. No output to write to was given. This is a minimum for initialization with the Canopy Constructor. Please use option -o'name_outputstream'." << endl;
    }
}


void InitialiseDetail(int flag_detail, int &error_detail){
    /* Error handling: There are no fatal errors. If file cannot be read, default values are used. If parameters are lacking from the file, default values are used */
    
    //vector<string> parameter_names{"paramID", "a_sdd", "b_sdd", "nbavg_sdd", "a_height", "b_height", "sigma_height", "hetscedas_height", "a_CR", "b_CR", "sigma_CR", "hetscedas_CR", "a_CD", "b_CD", "sigma_CD", "hetscedas_CD", "LMA", "sigma_LMA", "Nmass", "sigma_Nmass", "Pmass", "sigma_Pmass", "wsg", "sigma_wsg", "corr_Nmass_Pmass", "corr_Nmass_LMA", "corr_Pmass_LMA", "shape_crown", "gapfraction_crown"}; // this type of vector initialization does not seem to work across all types of compilers/standards, so using in-built array instead
    // nb_parametersdetailed = int(parameter_names.size());
    
    string parameter_names[29] = {"paramID", "a_sdd", "b_sdd", "nbavg_sdd", "a_height", "b_height", "sigma_height", "hetscedas_height", "a_CR", "b_CR", "sigma_CR", "hetscedas_CR", "a_CD", "b_CD", "sigma_CD", "hetscedas_CD", "LMA", "sigma_LMA", "Nmass", "sigma_Nmass", "Pmass", "sigma_Pmass", "wsg", "sigma_wsg", "corr_Nmass_Pmass", "corr_Nmass_LMA", "corr_Pmass_LMA", "shape_crown", "gapfraction_crown"};
    
    nb_parametersdetailed = 29;
    int nb_parametersdetailed_fromfile = 0;

    if (flag_detail == 1) {
        if(flag_detail) sprintf(inputfile_detail,"%s",bufi_detail);
        fstream InDetail(inputfile_detail, ios::in);
        
        if(InDetail){
            cout << "Reading in file: " << inputfile_detail << endl;
            
            string line;
            getline(InDetail,line);
            istringstream firstlinestream(line);
            
            string parameter_name;
            
            vector<int> parameter_index;
            while(firstlinestream >> parameter_name){
                int index = -1;
                for(int i = 0; i < nb_parametersdetailed; i++){
                     if(parameter_name == parameter_names[i]){
                         index = i;
                         nb_parametersdetailed_fromfile++;
                     }
                }
                if(index == -1){
                    cout << "Ignoring unknown parameter: " << parameter_name << ". Parameters must be one of:";
                    for(int i = 0; i < nb_parametersdetailed; i++){
                        cout << "\t" << parameter_names[i];
                    }
                    cout << endl;
                }
                parameter_index.push_back(index);
            }
     
            //InDetail.getline(buffer,512,'\n');
            
            //for(int i=0; i < nb_parametersets; i++){
              
            nb_parametersets = 0;
            
            while(getline(InDetail, line)){
                istringstream linestream(line);
            
                /* first we initialise the potential parameter arrays with empty parameter values */
                vector<string> parameter_values(nb_parametersdetailed,"");
                for(int p = 0; p < nb_parametersdetailed; p++){
                    parameters_detailed.push_back(0.0);
                }

                /* we update the values from the file */
                int s = 0;
                string parameter_value;
                while(linestream >> parameter_value){
                    int index = parameter_index[s];
                    if(index >= 0){
                        parameter_values[index] = parameter_value;
                    }
                    s++;
                }
            
                /* now we assign values */
                for(int i = 0; i < nb_parametersdetailed; i++){
                    AssignValueDetail(parameter_names[i], parameter_values[i], nb_parametersets);
                }

                cout << endl;

                nb_parametersets++;
            }
            
            if(nb_parametersets == 0){
                error_detail = 1;
                cout << "WARNING. Detailed input file was empty. All detailed parameters will be set to default." << endl;
            } else {
                error_detail = 0;
                cout << "Read " << nb_parametersets << " parameter sets with " << nb_parametersdetailed_fromfile << " parameters. " << nb_parametersdetailed - nb_parametersdetailed_fromfile << " parameters have been set to default." << endl;
            }
            
            InDetail.close();
        } else {
            error_detail = 2;
            cout << "ERROR. Detailed input file could not be read. Canopy Constructor will exit." << endl;
        }
    } else {
        error_detail = 1;
        cout << "WARNING. No detailed input file provided. All detailed parameters will be set to default." << endl;
    }

    if(error_detail == 1){
        /* default parameterization */
        nb_parametersets = 1;
        int parameterset = 0;
        vector<string> parameter_values(nb_parametersdetailed,"");
        for(int p = 0; p < nb_parametersdetailed; p++){
            parameters_detailed.push_back(0.0);
        }
        for(int i = 0; i < nb_parametersdetailed; i++){
            AssignValueDetail(parameter_names[i], parameter_values[i], parameterset);
        }
        
        cout << endl;
        cout << "Default initialization, using " << nb_parametersets << " parameter sets with " << nb_parametersdetailed << " parameters." << endl;
    }
    
    if(error_detail < 2){
        /* now write the last parameter set to file */
        output_input_detailed << parameter_names[0];
        for(int i = 1; i < nb_parametersdetailed; i++){
            output_input_detailed << "\t" << parameter_names[i];
        }
        output_input_detailed << endl;
        
        output_input_detailed << parameters_detailed[0 + (nb_parametersets-1) * nb_parametersdetailed];
        for(int i = 1; i < nb_parametersdetailed; i++){
            output_input_detailed << "\t" << parameters_detailed[i + (nb_parametersets-1) * nb_parametersdetailed];
        }
    }
}

void AllocMem(){

    /*** Initialization of trees ***/
    /*******************************/
    
    if(NULL==(T=new Tree[sites])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Tree" << endl;
    }
    
    /* LookUp table that gives the sites within a crown in order of distance from the center */
    /* crowns can thus be assembled from inside out, in a radial fashion */
    
    int Crown_dist[2601];                                           // this saves the distances from the center of the crown
    int extent = 25;                                                // maximum extent of crowns (after test: either enlarge or allocate dynamically)
    int extent_full = 2*25 + 1;
    int index_crown = 0, xx, yy, site_rel, dist;
    Crown_dist[index_crown] = 0;                                    // this is the distance of the center of the crown from the center (0)
    LookUp_Crown_site[index_crown] = extent + extent_full * extent;      // this is the label of the site at the center of the crown ( x = extent, y = extent)
    
    /* loop over crown */
    for(int col = 0; col < extent_full; col++){
        for(int row = 0; row < extent_full; row++){
            xx = col - extent;                                      // distance from center (x = extent) in x direction
            yy = row - extent;                                      // distance from center (y = extent) in y direction
            if(!(xx == 0 & yy == 0)){
                site_rel = col + extent_full * row;
                dist = xx*xx + yy*yy;
                /* now order the arrays according to distance from center */
                /* index_crown saves last filled position in array */
                /* for every voxel we run through the array from position zero to last filled position, and check where to put the new value */
                for(int i=0; i < index_crown + 1; i++){
                    int temp = Crown_dist[i];
                    int site_rel_temp = LookUp_Crown_site[i];
                    if(dist <= temp){                               // if distance is smaller than at current array position, move everything up
                        Crown_dist[i] = dist;
                        LookUp_Crown_site[i] = site_rel;
                        dist = temp;
                        site_rel = site_rel_temp;
                    }
                }
                Crown_dist[index_crown+1] = dist;                   // the last value that has been pushed into the dist variable fills a new position
                index_crown = index_crown + 1;
            }
        }
    }
    
    /* LookUp table that gives the sites within a radius in order of distance from the center */
    /* neighbours can thus be assembled from inside out, in a radial fashion */
    
    int neighbour_dist[10201];                                           // this saves the distances from the center of the crown
    int nextent = 50;                                                // maximum extent of neighbour search
    int nextent_full = 2*50 + 1;
    int index_neighbour = 0, nxx, nyy, nsite_rel, ndist;
    neighbour_dist[index_neighbour] = 0;                                    // this is the distance of the center of the crown from the center (0)
    LookUp_neighbour_site[index_neighbour] = nextent + nextent_full * nextent;      // this is the label of the site at the center of the crown ( x = extent, y = extent)
    
    /* loop across neighbours */
    for(int col = 0; col < nextent_full; col++){
        for(int row = 0; row < nextent_full; row++){
            nxx = col - nextent;                                      // distance from center (x = extent) in x direction
            nyy = row - nextent;                                      // distance from center (y = extent) in y direction
            if(!(nxx == 0 & nyy == 0)){
                nsite_rel = col + nextent_full * row;
                ndist = nxx*nxx + nyy*nyy;
                /* now order the arrays according to distance from center */
                /* index_neighbour saves last filled position in array */
                /* for every voxel we run through the array from position zero to last filled position, and check where to put the new value */
                for(int i=0; i < index_neighbour + 1; i++){
                    int ntemp = neighbour_dist[i];
                    int nsite_rel_temp = LookUp_neighbour_site[i];
                    if(ndist <= ntemp){                               // if distance is smaller than at current array position, move everything up
                        neighbour_dist[i] = ndist;
                        LookUp_neighbour_site[i] = nsite_rel;
                        ndist = ntemp;
                        nsite_rel = nsite_rel_temp;
                    }
                }
                neighbour_dist[index_neighbour+1] = ndist;                   // the last value that has been pushed into the dist variable fills a new position
                index_neighbour = index_neighbour + 1;
            }
        }
    }
    
    /* now voxel fields and CHM */
    
    if (NULL==(LAI3D=new float*[height_max+1]))
        cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        if (NULL==(LAI3D[h]=new float[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        for(int site=0;site<sites;site++)
            LAI3D[h][site] = 0.0;
    
    if (NULL==(Voxel3D=new int*[height_max+1]))
        cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        if (NULL==(Voxel3D[h]=new int[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        for(int site=0;site<sites;site++)
            Voxel3D[h][site] = 0;

    if (NULL==(correlation_structure=new float*[dbhlog_binnb]))
        cerr<<"!!! Mem_Alloc\n";
    for(int bin=0;bin < dbhlog_binnb; bin++)
        if (NULL==(correlation_structure[bin]=new float[6]))
            cerr<<"!!! Mem_Alloc\n";
    for(int bin=0;bin < dbhlog_binnb; bin++)
        for(int i=0;i<6;i++)
            correlation_structure[bin][i] = 0.0;
    
    if (NULL==(transmittance_simulated=new float*[height_max+1]))
        cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        if (NULL==(transmittance_simulated[h]=new float[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        for(int site=0;site<sites;site++)
            transmittance_simulated[h][site] = 0.0;
    
    if (NULL==(ALS_sampling=new int*[height_max+1]))
        cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        if (NULL==(ALS_sampling[h]=new int[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        for(int site=0;site<sites;site++)
            ALS_sampling[h][site] = 0;
    
    if (NULL==(ALS_echos=new int*[height_max+1]))
        cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        if (NULL==(ALS_echos[h]=new int[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<(height_max+1);h++)
        for(int site=0;site<sites;site++)
            ALS_echos[h][site] = 0;
    
    if (NULL==(chm_empirical=new int[sites])) cerr<<"!!! Mem_Alloc\n";
    for(int s=0;s<sites;s++){
        chm_empirical[s] = 0;
    }

    chm_field.reserve(sites);
    for(int s=0;s<sites;s++){
        chm_field[s] = 0;
    }
    
    if (NULL==(trees_sortedheight=new int*[2])) cerr<<"!!! Mem_Alloc\n";
    for(int l=0;l<2;l++)
        if (NULL==(trees_sortedheight[l]=new int[sites]))
            cerr<<"!!! Mem_Alloc\n";
    for(int s=0;s<sites;s++){
        trees_sortedheight[0][s] = -1;
        trees_sortedheight[1][s] = -1;
    }
    
    if (NULL==(hist_CHM_emp=new int[height_max+1])) cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<height_max+1;h++){
        hist_CHM_emp[h] = 0;
    }

    if (NULL==(hist_CHM_empcompcumulative=new float[height_max+1])) cerr<<"!!! Mem_Alloc\n";
    for(int h=0;h<height_max+1;h++){
        hist_CHM_empcompcumulative[h] = 0.0;
    }

    hist_CHM_sim.reserve(height_max+1);
    for(int i = 0; i < height_max + 1; i++){
        hist_CHM_sim.push_back(0);
    }

    if (NULL==(sigma_heightbin=new float[dbhlog_binnb])) cerr<<"!!! Mem_Alloc\n";
    for(int bin=0;bin<dbhlog_binnb;bin++){
        sigma_heightbin[bin] = 0.0;
    }
    
    if (NULL==(sigma_CRbin=new float[dbhlog_binnb])) cerr<<"!!! Mem_Alloc\n";
    for(int bin=0;bin<dbhlog_binnb;bin++){
        sigma_CRbin[bin] = 0.0;
    }
    
    /* Climate variables */
    /* monthly averages */
    if(NULL==(Temperature=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Temperature" << endl;
    }

    if(NULL==(DailyMaxTemperature=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyMaxTemperature" << endl;
    }

    if(NULL==(NightTemperature=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc NightTemperature" << endl;
    }

    if(NULL==(Rainfall=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Rainfall" << endl;
    }

    if(NULL==(WindSpeed=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc WindSpeed" << endl;
    }

    if(NULL==(MaxIrradiance=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Irradiance" << endl;
    }

    if(NULL==(MeanIrradiance=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Irradiance" << endl;
    }

    if(NULL==(SaturatedVapourPressure=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc SaturatedVapourPressure" << endl;
    }

    if(NULL==(VapourPressure=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc VapourPressure" << endl;
    }

    if(NULL==(VapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc VapourPressureDeficit" << endl;
    }

    if(NULL==(DailyVapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyVapourPressureDeficit" << endl;
    }

    if(NULL==(DailyMaxVapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyMaxVapourPressureDeficit" << endl;
    }
   
    /* LookUp tables */
    nbTbins=500;
    float Taccuracy=0.1;
    iTaccuracy=1.0/Taccuracy;
    cout << "Built-in maximal temperature: " << float(nbTbins)*Taccuracy <<endl;
    if(NULL==(LookUp_KmT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_KmT" << endl;
    if(NULL==(LookUp_GammaT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_GammaT" << endl;
    if(NULL==(LookUp_VcmaxT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_VcmaxT" << endl;
    if(NULL==(LookUp_JmaxT=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_JmaxT" << endl;
    if(NULL==(LookUp_Rday=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rday" << endl;
    if(NULL==(LookUp_Rstem=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rstem" << endl;
    if(NULL==(LookUp_Rnight=new float[nbTbins])) cerr<<"!!! Mem_Alloc LookUp_Rnight" << endl;
    for(int i=0;i<nbTbins;i++) { // loop over "T" in GPPleaf()
        float temper=float(i)*Taccuracy;
        LookUp_KmT[i] = 404.0*exp(((temper-25.0)/(298*0.00831*(273+temper)))*59.36)*
        (1+210*1.0/248.0*exp(-(temper-25.0)/(298*0.00831*(273+temper))*35.94))*iCair;
        LookUp_GammaT[i]=37.0*exp(((temper-25.0)/(298*0.00831*(273+temper)))*23.4)*iCair;
        LookUp_VcmaxT[i]=exp(26.35-65.33/(0.00831*(temper+273.15)));
        LookUp_JmaxT[i]=exp(17.57-43.54/(0.00831*(temper+273.15)));
        LookUp_Rday[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper)));
        LookUp_Rstem[i]=39.6*378.7*timestep*exp(((temper-25.0)/10.0)*log(2.0));
        LookUp_Rnight[i]=exp((temper-25.0)*0.1*log(3.09-0.0215*(25.0+temper)));
        /* exp((temp-25)/10*log(2)) is the temperature dependency of Rstem, supposing a constant Q10=2, according to Ryan et al 1994 and Meir & Grace 2002
        exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight))) is the temperature dependencies used by Atkin 2015 (equ1) */
    }
   
    /* look up table for flux averaging/integration */
    /* division into absorption prior to current voxel (absorb_prev) and absorption in current voxel (absorb_delta) */
    /* prior absorption has a maximum of 20 m2/m3, while absorption within one voxel can maximally reach 10 m2/m3*/

    if(NULL==(LookUp_flux_absorption=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_flux" << endl;
    if(NULL==(LookUp_flux=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_flux" << endl;
    if(NULL==(LookUp_VPD=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_VPD" << endl;
    if(NULL==(LookUp_T=new float[80000])) cerr<<"!!! Mem_Alloc LookUp_VPD" << endl;
    for(int i=0;i<400;i++) { // loop over "absorb" in Fluxh()
        for(int j=0;j<200;j++){
            float absorb_prev=float(i)/20.0;
            float absorb_delta=float(j)/20.0;
            if(absorb_delta==0.0){
                // flux is now computed simply as absorption on a per m2 plant matter basis (needs to be modified in TROLL, weighted by contribution of each tree)
                // but since there is no plant matter in the case of absorb_delta == 0, PPFD should be zero here
                LookUp_flux_absorption[i+400*j] = 0.0;    // in this case
                LookUp_flux[i+400*j] = exp(-kpar * absorb_prev);
                // if the voxel does not contain any plant matter, values are constant across the voxel, e.g. just the top value calculated from absorb_prev
                LookUp_VPD[i+400*j] = 0.25 + sqrt(maxf(0.0 , 0.08035714*(7.0-absorb_prev)));
                LookUp_T[i+400*j] = 0.4285714 * (minf(7.0,absorb_prev));
            } else {
                // if the voxel does contain plant matter, then average values will be computed
                // for voxels of 1 unit length depth, this corresponds just to the integral over LAI, which can be decomposed into a constant absorb_prev and a linearly increasing absorb_delta
                // once LAI reaches the critical value of 7.0, VPD And T do not decrease anymore, hence the distinction between two cases
                // flux is now computed simply as absorption on a per m2 plant matter basis (needs to be modified in TROLL, weighted by contribution of each tree). This means taking the PPFD at the top of the voxel, and calculating the absorption inside the voxel, and then convert it into flux per square meter vegetation (absorb_delta given in m2/m2)
                LookUp_flux_absorption[i+400*j] = exp(-kpar * absorb_prev) * (1.0 - exp(-kpar * absorb_delta)) / absorb_delta;
                LookUp_flux[i+400*j] = exp(-kpar * absorb_prev) * (1.0 - exp(-kpar * absorb_delta)) / (kpar * absorb_delta);
                if(absorb_prev+absorb_delta >= 7){
                    // this could be calculated more exactly, but difference is negligible
                    LookUp_VPD[i+400*j] = 0.25;
                    LookUp_T[i+400*j] = 3.0;  // 0.4285714 * 7.0
                } else {
                    LookUp_VPD[i+400*j] = 0.25 + (0.188982/absorb_delta) * (pow((7.0 - absorb_prev),1.5) - pow((7.0 - absorb_prev - absorb_delta),1.5));
                    LookUp_T[i+400*j] = 0.4285714 * (absorb_prev + 0.5 * absorb_delta);
                }
            }
        }
    }
    
    int height_canopy_max = height_max + 1;
    /* crown packing across all parameter sets */
    CP_twodimemp_params.reserve(nb_parametersets * height_canopy_max * height_canopy_max);
    for(int parameterset = 0; parameterset < nb_parametersets; parameterset++){
        for(int height_canopy = 0; height_canopy < height_max + 1; height_canopy++){
            for(int h = 0; h < height_max + 1; h++){
                CP_twodimemp.push_back(0.0);
            }
        }
    }
    
    /* for one parameter set */
    CP_twodimemp.reserve(height_canopy_max * height_canopy_max);
    for(int height_canopy = 0; height_canopy < height_max + 1; height_canopy++){
       for(int h = 0; h < height_max + 1; h++){
           CP_twodimemp.push_back(0.0);
       }
    }
}

void CalcHistCA(vector<int> &hist_CA){
    /* clean up field */
    for(int h = 0; h < height_max + 1; h++){
        hist_CA[h] = 0;
    }
    
    for(int s = 0; s < sites; s++){
        if(T[s].t_age > 0 && T[s].t_Tree_Height > float(height_cutoff_fromground)){
            float height = T[s].t_Tree_Height;
            int int_height = int(height);
            float CR = T[s].t_Crown_Radius;
            float CD = T[s].t_Crown_Depth;
            
            if(CD < 3.0){
                float crown_area = PI * CR * CR;     // floor of crown_area to bound area accumulation
                int crown_area_int = int(crown_area);   // floor of crown_area to bound area accumulation
                crown_area_int = max(crown_area_int,1);                                 // minimum area of crown (1), important for radii below ~0.5
                crown_area_int = min(crown_area_int,1963);                              // maximum area of crown (radius 25)
                hist_CA[int_height] += crown_area_int;
                
            } else {
                /* the height variable to sum up heights across the crown */
                float crownshell_base = height - CD + 1.5;
                float crown_extent = height - crownshell_base;
                float crown_slope = CR * (1.0 - shape_crown) / (crown_extent - 0.5);
                /* area for the innermost (and innermost) cylinder: this is the part of the crown that keeps pushing out and creates new layers */
                int crown_extent_toplayer = int(crown_extent - 0.5);                 /* the extent gives the number of full layers both upwards and downwards not including the central layer (i.e. half of it, ~ 0.5) */
                float radius_innermost = CR - crown_slope * float(crown_extent_toplayer);
                
                /* calculate the innermost/uppermost cylinder */
                int crown_area_previous = 0;
                
                /* crown area */
                int crown_area_innermost_int = GetCrownIntarea(radius_innermost);
                
                /* height at which we allocate */
                int height_innermost = int(height);
                hist_CA[height_innermost] += crown_area_innermost_int;
                
                /* update previous crown area */
                crown_area_previous = crown_area_innermost_int;
                
                /* now loop through the outer crown cylinders */

                int height_toplayer = int(crownshell_base + 0.5 + float(crown_extent_toplayer));
                int height_baselayer = int(crownshell_base + 1.5);
                for(int h = height_toplayer; h >= height_baselayer; h--){
                    /* calculating the radius of the current layer depending on the respective slopes */
                    
                    float radius_height = CR - crown_slope * (h - height_baselayer);    /* for the lowest layer, i.e. h == height_baselayer, radius = t_Crown_Radius */
                    
                    /* crown area */
                    int crown_area_int = GetCrownIntarea(radius_height);
                    
                    hist_CA[h] += (crown_area_int - crown_area_previous);
                    
                    /* update previous crown area */
                    crown_area_previous = crown_area_int;
                }
            }
        }
    }
}

void CalcFieldCA(vector<int> &field_CA){
    /* clean up field */
    for(int s = 0; s < sites; s++){
        field_CA[s] = 0;
    }
    
    for(int s = 0; s < sites; s++){
        if(T[s].t_age > 0 && T[s].t_Tree_Height > float(height_cutoff_fromground)){
            int site_crowncenter = s + T[s].t_CrownDisplacement;
            
            int row_crowncenter = site_crowncenter/cols,
            col_crowncenter = site_crowncenter%cols;
            
            float crown_radius = T[s].t_Crown_Radius;
            float crown_area = PI * crown_radius * crown_radius;     // floor of crown_area to bound area accumulation
            int crown_area_int = int(crown_area);   // floor of crown_area to bound area accumulation
            crown_area_int = max(crown_area_int,1);                                 // minimum area of crown (1), important for radii below ~0.5
            crown_area_int = min(crown_area_int,1963);                              // maximum area of crown (radius 25)
         
            for(int i = 0; i < crown_area_int; i++){
                int site_relative = LookUp_Crown_site[i];
                
                /* choose orientation of tree depending on site, non-random in order to save calculation time */
                //        int row, col;
                //        if(row_crowncenter%2==0) row = row_crowncenter + site_relative/51 - 25;
                //        else row = row_crowncenter - site_relative/51 + 25;
                //        if(col_crowncenter%2==0) col = col_crowncenter + site_relative%51 - 25;
                //        else  col = col_crowncenter - site_relative%51 + 25;
                //
                int row = row_crowncenter + site_relative/51 - 25;
                int col = col_crowncenter + site_relative%51 - 25;
                
                if(row >= 0 && row < rows && col >= 0 && col < cols){
                    int site = col + row * cols;
                    field_CA[site]++;
                }
            }
        }
    }
}
    

void CalcFieldAGB(vector<float> &field_AGB){
    /* clean up field */
    for(int s = 0; s < sites; s++){
        field_AGB[s] = 0.0;
    }
    
    /* now we run through all the trees, calculate their biomass and then spread it out equally across the crown radius */
    for(int s = 0; s < sites; s++){
        if(T[s].t_age > 0){
            if(T[s].t_Tree_Height > float(height_cutoff_fromground)){
                float mass_cylinder = T[s].t_wsg*T[s].t_Tree_Height*T[s].t_dbh*T[s].t_dbh; // the mass of the stem if assumed to be cylindric. This is converted into real biomass using tapering factors, etc.
                float AGBtree = 0.0673*pow(mass_cylinder*10000, 0.976);
                
                int site_crowncenter = s + T[s].t_CrownDisplacement;
                
                int row_crowncenter = site_crowncenter/cols,
                col_crowncenter = site_crowncenter%cols;
                
                float crown_radius = T[s].t_Crown_Radius;
                int crown_area_int = GetCrownIntarea(crown_radius);

                vector<int> crown_heights(crown_area_int,0);
                GetCrownHeights(crown_heights, T[s].t_Tree_Height, T[s].t_Crown_Radius, T[s].t_Crown_Depth);
                
                int height_cumulated = 0;
                for(int i = 0; i < crown_area_int; i++){
                    height_cumulated += crown_heights[i];
                }
                
                for(int i = 0; i < crown_area_int; i++){
                    int site_relative = LookUp_Crown_site[i];
                    
                    int row = row_crowncenter + site_relative/51 - 25;
                    int col = col_crowncenter + site_relative%51 - 25;
                    
                    int height_site = crown_heights[i];
                    float height_perc = float(height_site)/float(height_cumulated);
                    float agb_site = height_perc * AGBtree;
                    
                    if(row >= 0 && row < rows && col >= 0 && col < cols){
                        int site = col + row * cols;
                        field_AGB[site] += agb_site;
                        
                    }
                }
            }
        }
    }
}
    
void CheckAllometricScaling(){
    /* linear least squares fitting of effective allometries */
    double a_CReff, b_CReff;
    double cov00, cov01, cov11, chisq;
    
    vector<double> dbh_census;
    vector<double> CR_census;
    int trees_alive = 0;
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0 && T[site].t_dbh > dbh_cutoff_fromground){
            trees_alive++;
            dbh_census.push_back(log(T[site].t_dbh));
            CR_census.push_back(log(T[site].t_Crown_Radius));
        }
    }
    
    cout << "\n#########################################################################";
    cout << "\n##### Checking preservation of crown radius allometry after fitting #####";
    cout << "\n#########################################################################" << endl;
    
    cout << "Linear least squares fit on log-log scales for crown radius allometry " << endl;
    gsl_fit_linear(&dbh_census.front(), 1, &CR_census.front(), 1, trees_alive, &a_CReff, &b_CReff, &cov00, &cov01, &cov11, &chisq);
    cout << "Initial parameters: " << a_CR << " | " << b_CR << " Final fit (inferred): " << a_CReff << " | " << b_CReff << endl;
}
    
// TODO: needs to be split up
void OverallCanopyStatistics(){
    cout << "\n#######################################################";
    cout << "\n###### Writing overall canopy statistics to file ######";
    cout << "\n#######################################################" << endl;

    /* exclude layers underneath cutoff */
    int height_lowerlimit = height_cutoff_fromground;

    /* we take the full return densities (e.g. around 17/6.8 for PP2012 and GP2012, and then 36/15 for PP2015 and 31/12 for GP2015) */
    /* alternative: taking the pulse density, e.g. first/last return density (~ 12) */
    // alternative: model second returns and third returns explicitly, for example with calibration: 40% second return, 10% third return (but these are field averages, ideally derived from Beer-Lambert
    float knir = klight;
    float transmittance_nir = 0.4;
    float mean_beam = 12.0;
    float sd_beam = 5.0;
    
    
    /* now calculate overall canopy statistics, no ALS */
    /* two options for statistics: one normalized with respect to the ground, the other one with respect to the top */
    int height_canopy_max = height_max + 1;
    
    vector<int> voxtotal(height_canopy_max, 0);
    vector<int> voxfilled(height_canopy_max, 0);
    vector<int> voxempty(height_canopy_max, 0);
    vector<int> voxtotal_incan(height_canopy_max, 0);
    vector<int> voxfilled_incan(height_canopy_max, 0);
    vector<int> voxempty_incan(height_canopy_max, 0);
    
    vector<int> voxtotal_incan_fromtop(height_canopy_max, 0);
    vector<int> voxfilled_incan_fromtop(height_canopy_max, 0);
    vector<int> voxempty_incan_fromtop(height_canopy_max, 0);
    
    vector<float> transmfull_geom(height_canopy_max, 0.0);
    vector<float> transmfilled_geom(height_canopy_max, 0.0);
    vector<float> transmfull_arithm(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm(height_canopy_max, 0.0);
    
    vector<float> transmfull_geom_incan(height_canopy_max, 0.0);
    vector<float> transmfilled_geom_incan(height_canopy_max, 0.0);
    vector<float> transmfull_arithm_incan(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm_incan(height_canopy_max, 0.0);
    
    vector<float> transmfull_geom_incan_fromtop(height_canopy_max, 0.0);
    vector<float> transmfilled_geom_incan_fromtop(height_canopy_max, 0.0);
    vector<float> transmfull_arithm_incan_fromtop(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm_incan_fromtop(height_canopy_max, 0.0);
    
    /* PAD is within the canopy by default */
    vector<float> PADsum_incan(height_canopy_max, 0.0);
    vector<float> PADsum_incan_fromtop(height_canopy_max, 0.0);
    
    /* now calculate overall LAI3D */
    /* two options for each statistic: one normalized with respect to the ground, the other one with respect to the top */
    
    for(int h = height_max-1; h > height_lowerlimit - 1; h--){
        float LAD_layer = 0.0;
        float transm_layer = 0.0;
        float transmfilled_layer = 0.0;
        int voxtotal_layer = 0;
        int voxfilled_layer = 0;
        int voxempty_layer = 0;
        
        for(int s = 0; s < sites; s++){
            float lad_voxel = LAI3D[h][s] - LAI3D[h+1][s];
            if(lad_voxel < 0.001) lad_voxel = 0.0;         // test this!
            float transm_voxel = exp(-kpar * lad_voxel);
            
            voxtotal_layer++;
            LAD_layer += lad_voxel;
            transm_layer += transm_voxel;
            
            if(lad_voxel > 0.0){
                voxfilled_layer++;
                transmfilled_layer += transm_voxel;
            } else {
                voxempty_layer++;
            }
        }

        voxtotal[h] = voxtotal_layer;
        voxfilled[h] = voxfilled_layer;
        voxempty[h] = voxempty_layer;
        
        if(voxtotal_layer > 0){
            transmfull_geom[h] = exp(-kpar * LAD_layer/float(voxtotal_layer));
            transmfull_arithm[h] = transm_layer/float(voxtotal_layer);
        } else {
            transmfull_geom[h] = -1.0;
            transmfull_arithm[h] = -1.0;
        }
        
        if(voxfilled_layer > 0){
            transmfilled_geom[h] = exp(-kpar * LAD_layer/float(voxfilled_layer));
            transmfilled_arithm[h] = transmfilled_layer/float(voxfilled_layer);
        } else {
            transmfilled_geom[h] = -1.0;
            transmfilled_arithm[h] = -1.0;
        }
        
    }
    
    /* now add field */
    vector<float> fieldLAI;
    fieldLAI.reserve(sites);
    
    vector<float> fieldtransm;
    fieldtransm.reserve(sites);
    
    for(int s = 0; s < sites; s++){
        int height_canopy = chm_field[s];
        
        float LAI_cum = 0.0;

        for(int h = height_canopy; h > height_lowerlimit - 1; h--){
            int dist_top = height_canopy - h;
            float lad_voxel = LAI3D[h][s] - LAI3D[h+1][s];
            if(lad_voxel < 0.001) lad_voxel = 0.0;         // test this!
            float transm_voxel = exp(-kpar * lad_voxel);
            
            LAI_cum += lad_voxel;
            //if(LAI > 5) cout << " Warning! LAI: " << LAI << endl;
            
            PADsum_incan[h] += lad_voxel;
            PADsum_incan_fromtop[dist_top] += lad_voxel;
            
            voxtotal_incan[h]++;
            voxtotal_incan_fromtop[dist_top]++;
            transmfull_arithm_incan[h] += transm_voxel;
            transmfull_arithm_incan_fromtop[dist_top] += transm_voxel;

            if(lad_voxel > 0.0){
                voxfilled_incan[h]++;
                voxfilled_incan_fromtop[dist_top]++;
                transmfilled_arithm_incan[h] += transm_voxel;
                transmfilled_arithm_incan_fromtop[dist_top] += transm_voxel;
            } else {
                voxempty_incan[h]++;
                voxempty_incan_fromtop[dist_top]++;
            }
            
        }
        
        fieldLAI[s] = LAI_cum;
        fieldtransm[s] = exp(-kpar * LAI_cum);
        // compare: field_LAI[s] = LAI3D[h][s];
        
    }
    
    for(int h = height_lowerlimit; h < height_max + 1; h++){
        
        int voxtotal_incan_layer = voxtotal_incan[h];
        int voxfilled_incan_layer = voxfilled_incan[h];
        
        if(voxtotal_incan_layer > 0){
            transmfull_geom_incan[h] = exp(-kpar * PADsum_incan[h]/float(voxtotal_incan_layer));
            transmfull_arithm_incan[h] = transmfull_arithm_incan[h]/float(voxtotal_incan_layer);
        } else {
            transmfull_geom_incan[h] = -1.0;
            transmfull_arithm_incan[h] = -1.0;
        }
        
        if(voxfilled_incan_layer > 0){
            transmfilled_geom_incan[h] = exp(-kpar * PADsum_incan[h]/float(voxfilled_incan_layer));
            transmfilled_arithm_incan[h] = transmfilled_arithm_incan[h]/float(voxfilled_incan_layer);
        } else {
            transmfilled_geom_incan[h] = -1.0;
            transmfilled_arithm_incan[h] = -1.0;
        }
    }
    
    for(int dist = 0; dist < height_max + 1; dist++){
        int voxtotal_incan_layer = voxtotal_incan_fromtop[dist];
        int voxfilled_incan_layer = voxfilled_incan_fromtop[dist];
        
        if(voxtotal_incan_layer > 0){
            transmfull_geom_incan_fromtop[dist] = exp(-kpar * PADsum_incan_fromtop[dist]/float(voxtotal_incan_layer));
            transmfull_arithm_incan_fromtop[dist] = transmfull_arithm_incan_fromtop[dist]/float(voxtotal_incan_layer);
        } else {
            transmfull_geom_incan_fromtop[dist] = -1.0;
            transmfull_arithm_incan_fromtop[dist] = -1.0;
        }
        
        if(voxfilled_incan_layer > 0){
            transmfilled_geom_incan_fromtop[dist] = exp(-kpar * PADsum_incan_fromtop[dist]/float(voxfilled_incan_layer));
            transmfilled_arithm_incan_fromtop[dist] = transmfilled_arithm_incan_fromtop[dist]/float(voxfilled_incan_layer);
        } else {
            transmfilled_geom_incan_fromtop[dist] = -1.0;
            transmfilled_arithm_incan_fromtop[dist] = -1.0;
        }
    }
    
    
    vector<int> CV(height_canopy_max, 0);
    vector<int> CVsqueezed(height_canopy_max, 0);
    vector<int> CV_fromtop(height_canopy_max, 0);
    vector<int> CVsqueezed_fromtop(height_canopy_max, 0);
    
    /* now add fields */
    vector<int> fieldCV;
    vector<int> fieldCV_squeezed;
    vector<int> fieldempty_incan;
    fieldCV_squeezed.reserve(sites);
    fieldCV.reserve(sites);
    fieldempty_incan.reserve(sites);
    
    /* calculate voxels per height layer */
    /*
    for(int h = 0; h < height_max+1; h++){
        for(int s = 0; s < sites; s++){
            int voxel_ovl = Voxel3D[h][s];
            voxtotal[h]++;
            
            if(voxel_ovl > 0){
                voxfilled[h]++;
            } else {
                voxempty[h]++;
            }
        }
    }*/
        
        
    /* max height */
    int height_canopy_maxrealized = 0;
    
    for(int s = 0; s < sites; s++){
        int height_canopy = chm_field[s];
        if(height_canopy > height_canopy_maxrealized) height_canopy_maxrealized = height_canopy;
        
        int voxcrown_squeezed_cum = 0;
        int voxcrown_cum = 0;
        int voxempty_cum = 0;
        
        for(int h = height_canopy; h > height_lowerlimit - 1; h--){

            int voxel_ovl = Voxel3D[h][s];
            int dist_top = height_canopy - h;
            
            if(voxel_ovl == 0){
                voxempty_cum++;
            } else {
                voxcrown_squeezed_cum++;
                voxcrown_cum += voxel_ovl;
                CVsqueezed[h]++;
                CV[h] += voxel_ovl;
                CVsqueezed_fromtop[dist_top]++;
                CV_fromtop[dist_top] += voxel_ovl;
            }
        }
        fieldCV_squeezed[s] = voxcrown_squeezed_cum;
        fieldCV[s] = voxcrown_cum;
        fieldempty_incan[s] = voxempty_cum;
    }
    
    /* now create statistics for top crown area */
    
    vector<int> CA(height_canopy_max, 0);
    CalcHistCA(CA);
    
    vector<int> fieldCA(sites, 0);;
    CalcFieldCA(fieldCA);
    
    vector<float> fieldAGB(sites, 0.0);
    CalcFieldAGB(fieldAGB);
    
    /* NOW SIMULATE THE LIDAR AND DO IT ALL OVER AGAIN */
    /* calculating the TROLL transmittance field from simulated LiDAR */
    /* first we draw from a distribution to calculate the sampling density, i.e. the number of beams per voxel column */
    /* then we loop over voxel column from top to bottom and calculate the number of hits given the density of the respective voxel */
    /* if beams are all extinct, NAs are returned (-1), if beams hit ground, all produce guaranteed returns */
    
    /* loop over the LAI3D field */
    for(int r = 0; r < rows; r++){
        for(int c = 0; c < cols; c++){
            int site = c + r * cols;
            int nbbeams = int(mean_beam + gsl_ran_gaussian(gslrand, sd_beam));
            nbbeams = max(nbbeams,1);
            
            /* we loop over the field from maximum height to 0. If we needed ground returns, there is a simple extension to h = -1 (the latter corresponding to ground returns) */

            ALS_sampling[height_max][site] = nbbeams;
            ALS_echos[height_max][site] = 0;
            
            for(int h = height_max - 1; h >= 0; h--){
                int hits;
                
                ALS_sampling[h][site] = nbbeams;
                
                if(nbbeams == 0){
                    /* if there is no beam reaching the voxel, transmittance is set to NA (i.e. -1.0) and hits to zero */
                    hits = 0;
                } else {
                    if(h > 0){
                        /* returns due to vegetation */
                        float LAI_above = LAI3D[h+1][site];
                        float LAI_current = LAI3D[h][site];
                        
                        float prob_hit;
                        if(LAI_above == 100.0 & LAI_current == 100.0){
                            /* stem returns */
                            hits = nbbeams;
                            nbbeams = 0;
                        } else {
                            /* leaf/twig returns */
                            float LAD = LAI_current - LAI_above;
                            if(LAD > 0.0) prob_hit = 1.0 - exp(-knir * LAD);
                            else prob_hit = 0.0;
                            hits = gsl_ran_binomial(gslrand, prob_hit, nbbeams);
                            //transmittance = exp(-knir * LAD);
                            if(hits == 0){
                            } else {
                                nbbeams -= hits;
                                
                                int hits_notextinct = gsl_ran_binomial(gslrand, transmittance_nir, hits);
                                nbbeams += hits_notextinct;
                            }
                        }
                        
                    } else {
                        /* ground returns */
                        hits = nbbeams;
                    }
                }
                ALS_echos[h][site] = hits;
                
            }
        }
    }
    
    /* now evaluate the simulated ALS field */
    vector<int> voxtotal_ALS(height_canopy_max, 0);
    vector<int> voxnona_ALS(height_canopy_max, 0);
    vector<int> voxfilled_ALS(height_canopy_max, 0);
    vector<int> voxempty_ALS(height_canopy_max, 0);
    vector<int> voxblocked_ALS(height_canopy_max, 0);
    vector<int> voxnotblocked_ALS(height_canopy_max, 0);
    vector<int> voxna_ALS(height_canopy_max, 0);
    vector<int> voxtotal_incan_ALS(height_canopy_max, 0);
    vector<int> voxnona_incan_ALS(height_canopy_max, 0);
    vector<int> voxfilled_incan_ALS(height_canopy_max, 0);
    vector<int> voxempty_incan_ALS(height_canopy_max, 0);
    vector<int> voxblocked_incan_ALS(height_canopy_max, 0);
    vector<int> voxnotblocked_incan_ALS(height_canopy_max, 0);
    vector<int> voxna_incan_ALS(height_canopy_max, 0);
  
    vector<int> voxtotal_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxnona_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxfilled_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxempty_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxblocked_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxnotblocked_incan_ALS_fromtop(height_canopy_max, 0);
    vector<int> voxna_incan_ALS_fromtop(height_canopy_max, 0);
    
    vector<int> echos_incan(height_canopy_max, 0);
    vector<int> transmissions(height_canopy_max, 0);
    vector<int> transmissionsfilled(height_canopy_max, 0);
    vector<int> pulses(height_canopy_max, 0);
    vector<int> pulsesfilled(height_canopy_max, 0);
    vector<int> transmissions_incan(height_canopy_max, 0);
    vector<int> transmissionsfilled_incan(height_canopy_max, 0);
    vector<int> pulses_incan(height_canopy_max, 0);
    vector<int> pulsesfilled_incan(height_canopy_max, 0);
    
    vector<int> echos_incan_fromtop(height_canopy_max, 0);
    vector<int> transmissions_incan_fromtop(height_canopy_max, 0);
    vector<int> transmissionsfilled_incan_fromtop(height_canopy_max, 0);
    vector<int> pulses_incan_fromtop(height_canopy_max, 0);
    vector<int> pulsesfilled_incan_fromtop(height_canopy_max, 0);
    
    /* Very crude estimate of PAD. Inferred via log of local transmittance, maximum value is 3.5. PAD is within the canopy by default */
    vector<float> PADsum_incan_ALS(height_canopy_max, 0.0);
    vector<float> PADsum_incan_ALS_fromtop(height_canopy_max, 0.0);
    
    /* calculate empties per height layer */
   
    for(int h = height_lowerlimit; h < height_max+1; h++){

        for(int s = 0; s < sites; s++){
            int nbsampling = ALS_sampling[h][s];
            int nbechos = ALS_echos[h][s];
            int nbtransmissions = nbsampling - nbechos;

            voxtotal_ALS[h]++;
            
            if(nbsampling > 0){
                voxnona_ALS[h]++;
                pulses[h] += nbsampling;
                transmissions[h] += nbtransmissions;
                
                if(nbechos > 0){
                    voxfilled_ALS[h]++;
                    pulsesfilled[h] += nbsampling;
                    transmissionsfilled[h] += nbtransmissions;
                } else {
                    voxempty_ALS[h]++;
                }
                
                if(nbtransmissions > 0){
                    voxnotblocked_ALS[h]++;
                } else {
                    voxblocked_ALS[h]++;
                }
                
            } else {
                voxna_ALS[h]++;
            }
        }
    }
    
    vector<int> fieldCHM_ALS;
    fieldCHM_ALS.reserve(sites);
    
    for(int s = 0; s < sites; s++){
        int h = height_max+1;
        int nbechos = 0;
        
        while(h > 0 && nbechos == 0){
            h--;
            nbechos = ALS_echos[h][s];
        }
        
        fieldCHM_ALS[s] = h;
    }
    
    for(int s = 0; s < sites; s++){
        //int height_canopy = chm_field[s];
        int height_canopy = fieldCHM_ALS[s];
        

        for(int h = height_canopy; h > height_lowerlimit - 1; h--){

            int dist_top = height_canopy - h;
            int nbsampling = ALS_sampling[h][s];
            int nbechos = ALS_echos[h][s];
            int nbtransmissions = nbsampling - nbechos;
            
            voxtotal_incan_ALS[h]++;
            voxtotal_incan_ALS_fromtop[dist_top]++;
            
            if(nbsampling > 0){
                voxnona_incan_ALS[h]++;
                voxnona_incan_ALS_fromtop[dist_top]++;
                
                pulses_incan[h] += nbsampling;
                echos_incan[h] += nbechos;
                transmissions_incan[h] += nbtransmissions;
                
                pulses_incan_fromtop[dist_top] += nbsampling;
                echos_incan_fromtop[dist_top] += nbechos;
                transmissions_incan_fromtop[dist_top] += nbtransmissions;
                
                
                if(nbechos > 0){
                    //cout << " transm_voxel < 1.0: " << transm_voxel << endl;
                    voxfilled_incan_ALS[h]++;
                    pulsesfilled_incan[h] += nbsampling;
                    transmissionsfilled_incan[h] += nbtransmissions;
                    
                    voxfilled_incan_ALS_fromtop[dist_top]++;
                    pulsesfilled_incan_fromtop[dist_top] += nbsampling;
                    transmissionsfilled_incan_fromtop[dist_top] += nbtransmissions;
                } else {
                    //cout << " transm_voxel >= 1.0: " << transm_voxel << endl;
                    voxempty_incan_ALS[h]++;
                    voxempty_incan_ALS_fromtop[dist_top]++;
                }
                
                if(nbtransmissions > 0){
                    voxnotblocked_incan_ALS[h]++;
                    voxnotblocked_incan_ALS_fromtop[dist_top]++;
                } else {
                    voxblocked_incan_ALS[h]++;
                    voxblocked_incan_ALS_fromtop[dist_top]++;
                }
                
                float transmittance_local = float(nbtransmissions)/float(nbsampling);
                //if(transmittance_local < 0.1737739) transmittance_local = 0.1737739;
                float PAD_local = log(transmittance_local)/(-knir * (1.0 - transmittance_nir));
                
                PADsum_incan_ALS[h] += PAD_local;
                PADsum_incan_ALS_fromtop[dist_top] += PAD_local;
                
            } else {
                voxna_incan_ALS[h]++;
                voxna_incan_ALS_fromtop[dist_top]++;
            }
        }
    }
    
    /* now compute the transmissions etc. */
    
    vector<float> transmfull_arithm_ALS(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm_ALS(height_canopy_max, 0.0);
    vector<float> transmfull_arithm_incan_ALS(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm_incan_ALS(height_canopy_max, 0.0);
    
    vector<float> transmfull_arithm_incan_ALS_fromtop(height_canopy_max, 0.0);
    vector<float> transmfilled_arithm_incan_ALS_fromtop(height_canopy_max, 0.0);

    for(int h = height_lowerlimit; h < height_max + 1; h++){

        if(pulses[h] > 0) transmfull_arithm_ALS[h] = float(transmissions[h])/float(pulses[h]);
        else transmfull_arithm_ALS[h] = -1;
        
        if(pulses_incan[h] > 0) transmfull_arithm_incan_ALS[h] = float(transmissions_incan[h])/float(pulses_incan[h]);
        else transmfull_arithm_incan_ALS[h] = -1;

        if(pulsesfilled[h] > 0) transmfilled_arithm_ALS[h] = float(transmissionsfilled[h])/float(pulsesfilled[h]);
        else transmfilled_arithm_ALS[h] = -1;
        
        if(pulsesfilled_incan[h] > 0) transmfilled_arithm_incan_ALS[h] = float(transmissionsfilled_incan[h])/float(pulsesfilled_incan[h]);
        else transmfilled_arithm_incan_ALS[h] = -1;
    }
    
    for(int dist = 0; dist < height_max + 1; dist++){
        if(pulses_incan_fromtop[dist] > 0) transmfull_arithm_incan_ALS_fromtop[dist] = float(transmissions_incan_fromtop[dist])/float(pulses_incan_fromtop[dist]);
        else transmfull_arithm_incan_ALS_fromtop[dist] = -1;
        
        if(pulsesfilled_incan_fromtop[dist] > 0) transmfilled_arithm_incan_ALS_fromtop[dist] = float(transmissionsfilled_incan_fromtop[dist])/float(pulsesfilled_incan_fromtop[dist]);
        else transmfilled_arithm_incan_ALS_fromtop[dist] = -1;
    }
    
        
    vector<int> fieldHeightsweighted(sites, 0);
        
    for(int s = 0; s < sites; s++){
        for(int h = 1; h < height_max + 1; h++){
            if(ALS_echos[h][s] > 0){
                fieldHeightsweighted[s] += h;
            }
        }
    }
        
        
    /*####################################################*/
    /*########### SPECIAL CANOPY METRICS #################*/
    /*####################################################*/
    
    /* Here we define a number of special canopy metrics that are supposed to capture more accurately space filling inside the canopy of tropical forests */
    /* This is essential for a good transition from constructing a canopy both from ground and above to from-above only */
    /* Furthermore, we would like to provide an accurate description of canopy structure, taking into account both vertical (depth in the canopy) and horizontal variation (by separating out different canopy heights) */
    /* The idea is to construct a two dimensional array (here as a one-dimensional vector) with canopy height on the x-axis and depth/within-canopy-height on the y-axis */
    /* If the y-axis is normalized to 1 (i.e. percentage height within the canopy), we both plot depth and canopy height together in a natural way */
    /* from this, we can derive relationships for the metric in question to pass on to canopy construction only from above */
    
    /* As before, we construct the arrays for several key metrics */
    /* In addition, we also add a crown number metric, measuring how many crowns we find in each layer */
    /* All of the metrics should be normalized by the total number of voxels in statistical software */
    /* Indexing on vector works like a 2dim array to be more efficient computationally. We separate into two dimensions (height_canopy as column, and within-canopy-height, typically h, as row) */
    
    int height_canopy_maxtwodim = height_canopy_max * height_canopy_max;
    
    vector<int> voxtotal_incan_twodim(height_canopy_maxtwodim, 0);
    vector<int> voxfilled_incan_twodim(height_canopy_maxtwodim, 0);
    vector<int> voxempty_incan_twodim(height_canopy_maxtwodim, 0);
    
    vector<float> PADsum_twodim(height_canopy_maxtwodim, 0.0);
    
    /* now compute the voxel field statistics */
    for(int s = 0; s < sites; s++){
        int height_canopy = chm_field[s];

        //for(int h = height_canopy; h > height_lowerlimit - 1; h--){
        for(int h = height_canopy; h >= 0; h--){
            float pad_voxel = LAI3D[h][s] - LAI3D[h+1][s];
            if(pad_voxel < 0.001) pad_voxel = 0.0;

            PADsum_twodim[height_canopy + h * height_canopy_max] += pad_voxel;
            voxtotal_incan_twodim[height_canopy + h * height_canopy_max]++;
            
            if(pad_voxel > 0.0){
                voxfilled_incan_twodim[height_canopy + h * height_canopy_max]++;

            } else {
                voxempty_incan_twodim[height_canopy + h * height_canopy_max]++;
            }
            
        }
    }

    vector<int> CV_twodim(height_canopy_maxtwodim, 0);
    vector<int> CVsqueezed_twodim(height_canopy_maxtwodim, 0);
    
    cout << "Outputting full crown packing matrix." << endl;
    
    /* same for the overlap field */
    for(int s = 0; s < sites; s++){
        int height_canopy = chm_field[s];

        for(int h = height_canopy; h >= 0; h--){
        //for(int h = height_canopy; h > height_lowerlimit - 1; h--){
            int voxel_ovl = Voxel3D[h][s];
            
            if(voxel_ovl > 0){
                CV_twodim[height_canopy + h * height_canopy_max] += voxel_ovl;
                CVsqueezed_twodim[height_canopy + h * height_canopy_max]++;
            }
        }
    }

    /* now compute the CV_twodim field, assuming a cutoff */
    vector<int> CV_twodimcutoff(height_canopy_maxtwodim, 0);
    
    /* first save the original voxel field */
    vector<int> Voxel3D_safe;
    Voxel3D_safe.reserve((height_max + 1)*sites);
    for(int h = 0; h < height_max + 1; h++){
        for(int s = 0; s < sites; s++){
            Voxel3D_safe.push_back(Voxel3D[h][s]);
            Voxel3D[h][s] = 0;
        }
    }
    
    float dbh_cutoff_fromground_output;
     if(ConstructorStep == 0){
         dbh_cutoff_fromground_output = 0.1;
         cout << "Outputting crown packing matrix with lower limit for dbh. Assuming default cutoff of 10cm." << endl;
     } else {
         dbh_cutoff_fromground_output = dbh_cutoff_fromground;
         cout << "Outputting crown packing matrix  with lower limit for dbh. Cutoff is: " << lround(dbh_cutoff_fromground_output * 100) << "cm." << endl;
     }
    
    for(int s = 0; s < sites; s++){
        if(T[s].t_dbh >= dbh_cutoff_fromground_output){
            int site_crown = s + T[s].t_CrownDisplacement;
            FillVoxel(1, site_crown, T[s].t_Tree_Height, T[s].t_Crown_Radius, T[s].t_Crown_Depth);
        }
    }
    
    /* same for the overlap field */
    for(int s = 0; s < sites; s++){
        int height_canopy = chm_field[s];

        for(int h = height_canopy; h >= 0; h--){
        //for(int h = height_canopy; h > height_lowerlimit - 1; h--){
            int voxel_ovl = Voxel3D[h][s];
            
            if(voxel_ovl > 0){
                CV_twodimcutoff[height_canopy + h * height_canopy_max] += voxel_ovl;
            }
        }
    }
    
    /* refill the voxel field */
    for(int h = 0; h < height_max + 1; h++){
        for(int s = 0; s < sites; s++){
            Voxel3D[h][s] = Voxel3D_safe[s + h * sites];
        }
    }
    
    Voxel3D_safe.clear();
    
    /*##############################################*/
    /*############## WRITE OUTPUTS #################*/
    /*##############################################*/
    
    output_vertical << "z" << "\t" << "voxtotal" << "\t" << "voxfilled" << "\t" << "voxempty" << "\t" << "voxtotal_incan" << "\t" << "voxfilled_incan" << "\t" << "voxempty_incan" << "\t" << "CA" << "\t" << "CV" << "\t" << "CVsqueezed" << /* "\t" << "CVtreelevel" << "\t" << "nbtrees" << */ "\t" << "transmfull_arithm" << "\t" << "transmfilled_arithm"  << "\t" << "transmfull_geom" << "\t" << "transmfilled_geom"  <<  "\t" << "PADsum_incan" << "\t" << "transmfull_arithm_incan"  << "\t" << "transmfilled_arithm_incan" << "\t" << "transmfull_geom_incan"  << "\t" << "transmfilled_geom_incan" <<  endl;
        
    for(int h = height_lowerlimit; h < height_max+1; h++){
        output_vertical << h << "\t" << voxtotal[h] << "\t" << voxfilled[h] << "\t" << voxempty[h] << "\t" << voxtotal_incan[h] << "\t" << voxfilled_incan[h] << "\t" << voxempty_incan[h] << "\t" << CA[h] << "\t" << CV[h] << "\t" << CVsqueezed[h] << "\t" << transmfull_arithm[h] << "\t" << transmfilled_arithm[h] << "\t"  << transmfull_geom[h] << "\t" << transmfilled_geom[h] << "\t" << PADsum_incan[h] << "\t" << transmfull_arithm_incan[h] << "\t" << transmfilled_arithm_incan[h] << "\t" << transmfull_geom_incan[h] << "\t" << transmfilled_geom_incan[h] << endl;
    }
    
    output_verticalfromtop << "d_cnpy" << "\t" << "voxtotal_incan" << "\t" << "voxfilled_incan" << "\t" << "voxempty_incan" << "\t" << "CV" << "\t" << "CV" << /* "\t" << "CVtreelevel" << "\t" << "nbtrees"  << */ "\t" << "PADsum_incan" << "\t" << "transmfull_arithm_incan" << "\t" << "transmfilled_arithm_incan" << "\t" << "transmfull_geom_incan" << "\t" << "transmfilled_geom_incan" << endl;
    
    for(int dist = 0; dist < height_max+1; dist++){
        output_verticalfromtop << dist << "\t" << voxtotal_incan_fromtop[dist] << "\t" << voxfilled_incan_fromtop[dist] << "\t" << voxempty_incan_fromtop[dist] << "\t" << CV_fromtop[dist] << "\t" << CVsqueezed_fromtop[dist] << "\t" << PADsum_incan_fromtop[dist] << "\t" << transmfull_arithm_incan_fromtop[dist] << "\t" << transmfilled_arithm_incan_fromtop[dist] << "\t" << transmfull_geom_incan_fromtop[dist] << "\t" << transmfilled_geom_incan_fromtop[dist] << endl;
    }
    
    output_vertical_ALS << "z" << "\t" << "voxtotal" << "\t" << "voxnona" << "\t" << "voxna"<< "\t" << "voxfilled" << "\t" << "voxempty" << "\t" << "voxblocked" << "\t" << "voxnotblocked"  << "\t" << "voxtotal_incan" << "\t" << "voxnona_incan" << "\t" << "voxna_incan" << "\t" << "voxfilled_incan" << "\t" << "voxempty_incan" << "\t" << "voxblocked_incan" << "\t" << "voxnotblocked_incan" << "\t" << "pulses" << "\t" << "pulsesfilled" << "\t" << "transmissions" << "\t" << "transmissionsfilled" << "\t" << "echos_incan" << "\t" << "pulses_incan" << "\t" << "pulsesfilled_incan" << "\t" << "transmissions_incan" << "\t" << "transmissionsfilled_incan" << "\t" << "transmfull_arithm" << "\t" << "transmfilled_arithm"  << "\t" << "PADsum_incan" << "\t" << "transmfull_arithm_incan"  << "\t" << "transmfilled_arithm_incan" <<  endl;
    

    for(int h = height_lowerlimit; h < height_max+1; h++){
        output_vertical_ALS << h << "\t" << voxtotal_ALS[h] << "\t" << voxnona_ALS[h] << "\t" << voxna_ALS[h] << "\t" << voxfilled_ALS[h] << "\t" << voxempty_ALS[h] << "\t" << voxblocked_ALS[h] << "\t" << voxnotblocked_ALS[h]  << "\t" << voxtotal_incan_ALS[h] << "\t" << voxnona_incan_ALS[h] << "\t" << voxna_incan_ALS[h] << "\t" << voxfilled_incan_ALS[h] << "\t" << voxempty_incan_ALS[h] << "\t" << voxblocked_incan_ALS[h] << "\t" << voxnotblocked_incan_ALS[h] << "\t" << pulses[h] << "\t" << pulsesfilled[h] << "\t" << transmissions[h] << "\t" << transmissionsfilled[h] << "\t" << echos_incan[h] << "\t" << pulses_incan[h] << "\t" << pulsesfilled_incan[h] << "\t" << transmissions_incan[h] << "\t" << transmissionsfilled_incan[h]  << "\t" << transmfull_arithm_ALS[h] << "\t" << transmfilled_arithm_ALS[h] << "\t" << PADsum_incan_ALS[h] << "\t" << transmfull_arithm_incan_ALS[h] << "\t" << transmfilled_arithm_incan_ALS[h] << endl;
    }

    output_verticalfromtop_ALS << "d_cnpy" << "\t" << "voxtotal_incan" << "\t" << "voxnona_incan" << "\t" << "voxna_incan" << "\t" << "voxfilled_incan" << "\t" << "voxempty_incan" << "\t" << "voxblocked_incan" << "\t" << "voxnotblocked_incan" << "\t" << "echos_incan" << "\t" << "pulses_incan" << "\t" << "pulsesfilled_incan" << "\t" << "transmissions_incan" << "\t" << "transmissionsfilled_incan" << "\t" << "PADsum_incan" << "\t" << "transmfull_arithm_incan" << "\t" << "transmfilled_arithm_incan" << endl;
    
    for(int dist = 0; dist < height_max+1; dist++){
        output_verticalfromtop_ALS << dist << "\t" << voxtotal_incan_ALS_fromtop[dist] << "\t" << voxnona_incan_ALS_fromtop[dist] << "\t" << voxna_incan_ALS_fromtop[dist] << "\t" << voxfilled_incan_ALS_fromtop[dist] << "\t" << voxempty_incan_ALS_fromtop[dist] << "\t" << voxblocked_incan_ALS_fromtop[dist] << "\t" << voxnotblocked_incan_ALS_fromtop[dist] << "\t" << echos_incan_fromtop[dist] << "\t" << pulses_incan_fromtop[dist] << "\t" << pulsesfilled_incan_fromtop[dist] << "\t" << transmissions_incan_fromtop[dist] << "\t" << transmissionsfilled_incan_fromtop[dist] << "\t" << PADsum_incan_ALS_fromtop[dist] << "\t" << transmfull_arithm_incan_ALS_fromtop[dist] << "\t" << transmfilled_arithm_incan_ALS_fromtop[dist] << endl;
    }
    
    output_field << "site" << "\t" << "x" << "\t" << "y" << "\t" << "chm_emp" << "\t" << "chm_random" << "\t" << "chm" << "\t" << "chm_ALS" << "\t" << "CA" << "\t" << "CV" << "\t" << "CVsqueezed" << "\t" << "empty_incan" << "\t" << "LAI" << "\t" << "transm" << "\t" << "agb_avg" << "\t" << "heightsweighted" << endl;
    for(int s = 0; s < sites; s++){
        output_field << s << "\t" << s%cols << "\t" << s/cols << "\t" << chm_empirical[s] << "\t" << chm_field_random[s] << "\t" << chm_field[s] << "\t" << fieldCHM_ALS[s] << "\t" << fieldCA[s] << "\t" << fieldCV[s] << "\t" << fieldCV_squeezed[s] << "\t" << fieldempty_incan[s] << "\t" << fieldLAI[s] << "\t" << fieldtransm[s] << "\t" << fieldAGB[s] << "\t" << fieldHeightsweighted[s] << endl;
    }

    /* output in two dimensions */
    /* we ignore height fields above canopy height (NA) */
    
    vector<float> CP_twodim(height_canopy_maxtwodim, 0.0);
    vector<float> CP_twodimcutoff(height_canopy_maxtwodim, 0.0);

    /* write in long format */
    for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
        for(int h = 0; h < height_canopy + 1; h++){
            
            float CVperc;
            float CVperccutoff;
            int voxtotal_incanopy = voxtotal_incan_twodim[height_canopy + h * height_canopy_max];
            if(voxtotal_incanopy > 0){
                CVperc = float(CV_twodim[height_canopy + h * height_canopy_max])/float(voxtotal_incanopy);
                CVperccutoff = float(CV_twodimcutoff[height_canopy + h * height_canopy_max])/float(voxtotal_incanopy);
            } else {
                CVperc = -1.0;
                CVperccutoff = -1.0;
            }
            
            CP_twodim[height_canopy + h * height_canopy_max] = CVperc;
            CP_twodimcutoff[height_canopy + h * height_canopy_max] = CVperccutoff;
        }
    }

    /* now create a density plot */
    /* we interpolate and extrapolate based on a minimum number of voxels */
    /* the general process is as following */
    /*  - we go to a particular canopy height and determine the number of pixels with that particular height */
    /*  - if the number is smaller than the minimum required number ("nbsamples_minimum"), we also draw in the two canopy heights below and above, and so forth, till enough pixels are sampled, TODO: nbsamples_minimum should be part of general input sheet */
    /*  - we then create an average density distribution across the sampled pixels */
    /* this metric automatically adapts itself (heights that occur more often are smoothed less, whereas poorly sampled heights will draw from considerable information from outside) */
    /* gaps in the height distribution will be interpolated, lack of information for very high or low heights will be extrapolated from the boundaries */
    /* to ensure even sampling, we first project each additional height's information onto the voxels of the sampled height, interpolating where information is missing */
    /* we then simply average across the distributions */
    
    vector<float> CP_twodim_fullextent(height_canopy_maxtwodim, 0.0);
    vector<float> CP_twodimcutoff_fullextent(height_canopy_maxtwodim, 0.0);
    int nbsamples_minimum = 1000;
    
    CreateFullCPmatrix(CP_twodim_fullextent, CP_twodim, voxtotal_incan_twodim, nbsamples_minimum, height_canopy_max);
    CreateFullCPmatrix(CP_twodimcutoff_fullextent, CP_twodimcutoff, voxtotal_incan_twodim, nbsamples_minimum, height_canopy_max);
    
    /* write everything to output file */
    output_twodim << "height_canopy" << "\t" << "h" << "\t" << "voxtotal_incan" << "\t" << "voxfilled_incan" << "\t" << "voxempty_incan" << "\t" << "CV" << "\t" << "CVsqueezed" << "\t" << "PADsum" << "\t" << "CVperc" << "\t" << "CVperc_fullextent" << "\t" << "CVperccutoff" << "\t" << "CVperccutoff_fullextent" << endl;
    
    /* write in long format */
    for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
        for(int h = 0; h < height_canopy + 1; h++){

            output_twodim << height_canopy << "\t" << h << "\t" << voxtotal_incan_twodim[height_canopy + h * height_canopy_max] << "\t" << voxfilled_incan_twodim[height_canopy + h * height_canopy_max] << "\t" << voxempty_incan_twodim[height_canopy + h * height_canopy_max] << "\t" << CV_twodim[height_canopy + h * height_canopy_max] << "\t" << CVsqueezed_twodim[height_canopy + h * height_canopy_max] << "\t" << PADsum_twodim[height_canopy + h * height_canopy_max] << "\t" << CP_twodim[height_canopy + h * height_canopy_max] << "\t" << CP_twodim_fullextent[height_canopy + h * height_canopy_max] << "\t" << CP_twodimcutoff[height_canopy + h * height_canopy_max] << "\t" << CP_twodimcutoff_fullextent[height_canopy + h * height_canopy_max] << endl;
        }
    }
    
    /* write packing coefficients to output file */
    output_input_cp << "paramID" << "\t" << "height_canopy" << "\t" << "height_within_canopy" << "\t" << "CP_height" << endl;
    for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
        for(int height_within_canopy = 0; height_within_canopy < height_canopy + 1; height_within_canopy++){
            
            output_input_cp << paramIDcurrent << "\t" << height_canopy << "\t" << height_within_canopy << "\t" << CP_twodim_fullextent[height_canopy + height_within_canopy * height_canopy_max] << endl;
        }
    }
    
    output_input_cpcutoff << "paramID" << "\t" << "height_canopy" << "\t" << "height_within_canopy" << "\t" << "CP_height" << endl;
    for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
        for(int height_within_canopy = 0; height_within_canopy < height_canopy + 1; height_within_canopy++){
            
            output_input_cpcutoff << paramIDcurrent << "\t" << height_canopy << "\t" << height_within_canopy << "\t" << CP_twodimcutoff_fullextent[height_canopy + height_within_canopy * height_canopy_max] << endl;
        }
    }
}

void CreateFullCPmatrix(vector<float> &CP_twodim_fullextent, vector<float> &CP_twodim, vector<int> &voxtotal_incan_twodim, int nbsamples_minimum, int height_canopy_max){
    for(int height_canopy = 0; height_canopy < height_canopy_max; height_canopy++){
        
        int nbsamples = 0;
        
        /* first get the number of pixels and write to array */
        vector<float> packing_height(height_canopy+1,0.0);
        int nbsamples_height = voxtotal_incan_twodim[height_canopy + height_canopy * height_canopy_max];
        
        /* if there is a sample, add to the array */
        int nbheights = 0;
        
        if(nbsamples_height > 0){
            nbheights++;
            nbsamples += nbsamples_height;
            for(int h = 0; h < height_canopy + 1; h++){
                packing_height[h] += CP_twodim[height_canopy + h * height_canopy_max];
            }
        }

        /* now, if we do not have enough samples, we add more height layers */
        int height_canopy_indexupper = height_canopy;
        int height_canopy_indexlower = height_canopy;
        
        int fullextent = 0;
        
        while(nbsamples < nbsamples_minimum && fullextent == 0){
            //if(height_canopy == 65) cout << height_canopy << " nbheights: " << nbheights << " nbsamples: " << nbsamples << endl;
            /* udpate the indices */
            height_canopy_indexupper += 1;
            height_canopy_indexlower -= 1;
            
            int extent_upper;
            if(height_canopy_indexupper >= height_canopy_max) extent_upper = 1;
            else extent_upper = 0;
            
            /* cannot be 0, because at zero, there is no internal structure */
            int extent_lower;
            if(height_canopy_indexlower <= 0) extent_lower = 1;
            else extent_lower = 0;
            
            /* if there is no more height layer to add, stop */
            if(extent_upper == 1 && extent_lower == 1) fullextent = 1;
            else {
                /* determine whether height layer can be used (i.e. is within the range ) */
                if(extent_upper == 0){
                    int nbsamples_upper = voxtotal_incan_twodim[height_canopy_indexupper + height_canopy_indexupper * height_canopy_max];
                    
                    /* now determine whether the height layer actually is sampled at all */
                    if(nbsamples_upper > 0){
                        /* now we convert the additional height layer to the scale of the original height layer */
                        vector<float> packing_upper(height_canopy+1,0.0);
                        vector<int> sampling_upper(height_canopy+1,0);
                        
                        for(int h = 0; h < height_canopy_indexupper + 1; h++){
                            float packing = CP_twodim[height_canopy_indexupper + h * height_canopy_max];
                            float percentage_height = float(h)/float(height_canopy_indexupper);
                            int height_backtransformed = int(percentage_height * float(height_canopy));
                            packing_upper[height_backtransformed] += packing;
                            sampling_upper[height_backtransformed]++;
                        }
                        
                        /* first we update layers that have been sampled */
                        for(int h = 0; h < height_canopy + 1; h++){
                            int sampling = sampling_upper[h];
                            if(sampling > 0){
                                /* if the particular voxel is sampled, take the mean */
                                float inv_sampling = 1.0/float(sampling);
                                packing_upper[h] *= inv_sampling;
                            }
                        }
                        
                        /* then we update the layers that need to be interpolated */
                        for(int h = 0; h < height_canopy + 1; h++){
                            int sampling = sampling_upper[h];
                            if(sampling == 0){
                                /* we loop in both directions, this works because we generally exclude the 0 layer */
                                /* both the uppermost layer (h_incanopy_upper == height_canopy) and the lowest layer (h_incanopy_lower == 0) are always filled */
                                int h_incanopy_upper = h;
                                int h_incanopy_lower = h;
                                int sampling_incanopy_upper = 0;
                                int sampling_incanopy_lower = 0;
                                
                                while(h_incanopy_upper < height_canopy + 1 && sampling_incanopy_upper == 0){
                                    h_incanopy_upper++;
                                    if(sampling_upper[h_incanopy_upper] > 0) sampling_incanopy_upper = 1;
                                }
                                
                                while(h_incanopy_lower > 0 && sampling_incanopy_lower == 0){
                                    h_incanopy_lower--;
                                    if(sampling_upper[h_incanopy_lower] > 0) sampling_incanopy_lower = 1;
                                }
                                
                                float packing = 0.5 * (packing_upper[h_incanopy_upper] + packing_upper[h_incanopy_lower]);
                                packing_upper[h] = packing;
                                sampling_upper[h] = 1;
                            }
                        }
                        
                        /* now add to the overall packing density */
                        for(int h = 0; h < height_canopy + 1; h++){
                            packing_height[h] += packing_upper[h];
                        }
                        
                        nbheights++;
                        nbsamples += nbsamples_upper;
                    }
                }
                
                /* determine whether height layer can be used (i.e. is within the range) */
                if(extent_lower == 0){
                    int nbsamples_lower = voxtotal_incan_twodim[height_canopy_indexlower + height_canopy_indexlower * height_canopy_max];
                    
                    /* now determine whether the height layer actually is sampled at all */
                    if(nbsamples_lower > 0){
                        /* now we convert the additional height layer to the scale of the original height layer */
                        vector<float> packing_lower(height_canopy+1,0.0);
                        vector<int> sampling_lower(height_canopy+1,0);
                        
                        for(int h = 0; h < height_canopy_indexlower + 1; h++){
                            float packing = CP_twodim[height_canopy_indexlower + h * height_canopy_max];
                            float percentage_height = float(h)/float(height_canopy_indexlower);
                            int height_backtransformed = int(percentage_height * float(height_canopy));
                            packing_lower[height_backtransformed] += packing;
                            sampling_lower[height_backtransformed]++;
                        }
                        
                        /* first we update layers that have been sampled */
                        for(int h = 0; h < height_canopy + 1; h++){
                            int sampling = sampling_lower[h];
                            if(sampling > 0){
                                /* if the particular voxel is sampled, take the mean */
                                float inv_sampling = 1.0/float(sampling);
                                packing_lower[h] *= inv_sampling;
                            }
                        }
                        
                        /* then we update the layers that need to be interpolated */
                        for(int h = 0; h < height_canopy + 1; h++){
                            int sampling = sampling_lower[h];
                            if(sampling == 0){
                                /* we loop in both directions, this works because we generally exclude the 0 layer */
                                /* both the uppermost layer (h_incanopy_upper == height_canopy) and the lowest layer (h_incanopy_lower == 0) are always filled */
                                int h_incanopy_upper = h;
                                int h_incanopy_lower = h;
                                int sampling_incanopy_upper = 0;
                                int sampling_incanopy_lower = 0;
                                
                                while(h_incanopy_upper < height_canopy + 1 && sampling_incanopy_upper == 0){
                                    h_incanopy_upper++;
                                    if(sampling_lower[h_incanopy_upper] > 0) sampling_incanopy_upper = 1;
                                }
                                
                                while(h_incanopy_lower > 0 && sampling_incanopy_lower == 0){
                                    h_incanopy_lower--;
                                    if(sampling_lower[h_incanopy_lower] > 0) sampling_incanopy_lower = 1;
                                }
                                
                                float packing = 0.5 * (packing_lower[h_incanopy_upper] + packing_lower[h_incanopy_lower]);
                                packing_lower[h] = packing;
                                sampling_lower[h] = 1;
                            }
                        }
                        
                        /* now add to the overall packing density */
                        for(int h = 0; h < height_canopy + 1; h++){
                            packing_height[h] += packing_lower[h];
                        }
                        
                        nbheights++;
                        nbsamples += nbsamples_lower;
                        
                    }
                }
            }
        }
        
        /* now calculate the overall average */
        float inv_nbheights;
        if(nbheights > 0) inv_nbheights = 1.0/float(nbheights);
        else inv_nbheights = 0.0;
        
        for(int h = 0; h < height_canopy + 1; h++){
            packing_height[h] *= inv_nbheights;
            float packing = packing_height[h];
            CP_twodim_fullextent[height_canopy + h * height_canopy_max] = packing;
        }
    }
}


void OutputInputCHMinventory(){
    /* output the final CHM and inventory*/
    
    output_input_CHM << "x\ty\tz";
    
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int site = col + row * cols;
            int height = chm_field[site];
            output_input_CHM << endl;
            output_input_CHM << col << "\t" << row << "\t" << height;
        }
    }
    
    output_input_inventory << "x\ty\tdbh\twsg";
    cout << "Outputting full inventory." << endl;
    
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int site = col + row * cols;
            if(T[site].t_age > 0.0){
                output_input_inventory << endl;
                output_input_inventory << col << "\t" << row << "\t" << T[site].t_dbh << "\t" << T[site].t_wsg;
            }
        }
    }
    
    output_input_inventorycutoff << "x\ty\tdbh\twsg";
    
    float dbh_cutoff_fromground_output;
    if(ConstructorStep == 0){
        dbh_cutoff_fromground_output = 0.1;
        cout << "Outputting inventory with lower limit for dbh. Assuming default cutoff of 10cm." << endl;
    } else {
        dbh_cutoff_fromground_output = dbh_cutoff_fromground;
        cout << "Outputting inventory with lower limit for dbh. Cutoff is: " << lround(dbh_cutoff_fromground_output * 100) << "cm." << endl;
    }
    
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int site = col + row * cols;
            float dbh = T[site].t_dbh;
            if(T[site].t_age > 0.0 && dbh >= dbh_cutoff_fromground_output){
                output_input_inventorycutoff << endl;
                output_input_inventorycutoff << col << "\t" << row << "\t" << dbh << "\t" << T[site].t_wsg;
            }
        }
    }
}

void UpdateLeaves(){
    
    nbtrees_carbonstarv = 0;
    
    /* create a leaves-on canopy */
    /* sort the trees according to height (descending) */
    /* create tree array */
    
    int nbtrees_update = 0;
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            nbtrees_update++;
        }
    }
    
    vector<int> trees_alive;
    trees_alive.reserve(nbtrees_update);
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            trees_alive.push_back(site);
        }
    }
    
    /* now order the array */
    for(int i = 0 ;i < nbtrees_update; i++){
        for(int tree = 0; tree < nbtrees_update-i-1; tree++){
            int site_current = trees_alive[tree];
            int site_next = trees_alive[tree+1];
            float height_current = T[site_current].t_Tree_Height;
            float height_next = T[site_next].t_Tree_Height;
            
            if(height_current < height_next){
                trees_alive[tree] = site_next;
                trees_alive[tree+1] = site_current;
            }
        }
    }
    
    /* clear the LAI3D field */
    /* this is crucial and has to be done each time to create comparable conditions that the Canopy Constructor improves upon */
    for(int h=height_max;h>=0;h--){
        for(int site=0;site<sites;site++){
            LAI3D[h][site] = 0.0;
        }
    }
    
    /* cycle through all trees, now ordered by height, and update the LAI field */
    for(int tree = 0 ;tree < nbtrees_update; tree++){
        int site_census = trees_alive[tree];
        
        /* get the tree dimensions */
        float dbh = T[site_census].t_dbh;
        float height = T[site_census].t_Tree_Height;
        float CR = T[site_census].t_Crown_Radius;
        float CD = T[site_census].t_Crown_Depth;
        int crown_displacement = T[site_census].t_CrownDisplacement;
        
        int site_crown = site_census + crown_displacement;
        
        float LAImax = T[site_census].t_LAImax;

        float LAIabove = UpdateLAIabove_eff(site_crown, height, CR, CD, dbh);

        
        //cout << site_census << " dbh: " << T[site_census].t_dbh << " Nmass: " << T[site_census].t_Nmass << " LAI: " << LAIabove << endl;
        /* now we update the LAI */
        /* if the tree experiences too much shade to even allocate a single leaf (LAIabove > LAImax), we set LAI = 0 without further checks */
        /* if the tree experiences some sun, we first calculate the maximum LAItree that it can have, which is limited by two factors: a) the effective LAI above it (LAIabove), and b) a physiological limitation, which here is modelled to lie somewhere between being able to optimally exploit light (LAI = LAImax) and just about surviving (LAI = LAImin). We summarize this through LAItree_max, which is the minimum of LAImax - LAIabove and LAImin + (LAImax - LAImin); a tree property that accounts for gaps and trees not being able to allocate all the leaves they potentially could allocate */
        /* We then allocate the LAItree_max, and calculate NPP for the whole tree. If NPP is < 0, the tree cannot survive within the current limits, so we set its LAI to zero */
        /* All trees with LAI zero (or NPP < 0) are essentially living from stored carbon and will be counted as "under carbon starvation" */
        
        
        if(LAIabove > LAImax){
            T[site_census].t_LAIabove = LAIabove;
            T[site_census].t_LAI = 0.0;
            T[site_census].UpdateDensity();
            
        } else {
            float LAItree_max_light = LAImax - LAIabove;
            float LAItree_max_physiology = LAImax;
            float LAI = minf(LAItree_max_light, LAItree_max_physiology);
            
            T[site_census].t_LAIabove = LAIabove;
            T[site_census].t_LAI = LAI;
            T[site_census].UpdateDensity();
            CalcLAItrial(0, site_crown, height, CR, CD, LAI,T[site_census].t_dbh);
        }
    }
    
    /* cycle through all trees again, now calculate the GPP/NPP */
    
    for(int tree = 0 ;tree < nbtrees_update; tree++){
        int site_census = trees_alive[tree];
        
        /* get the tree dimensions */
        float height = T[site_census].t_Tree_Height;
        float CR = T[site_census].t_Crown_Radius;
        float CD = T[site_census].t_Crown_Depth;
        int crown_displacement = T[site_census].t_CrownDisplacement;
        
        float dbh = T[site_census].t_dbh;
      
        float LAI = T[site_census].t_LAI;

//        float treeGPP_previous = T[site_census].t_GPP;
//        float treeNPP_previous = T[site_census].t_NPP;


        float treeGPP = 0.0;
        float treeNPP = 0.0;
        
        if(LAI > 0.0){
            float Nmass = T[site_census].t_Nmass;
            float Pmass = T[site_census].t_Pmass;
            float LMA = T[site_census].t_LMA;
            
            float wsg = T[site_census].t_wsg;
            
            float SLA=10000.0/LMA;
            float Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
            float Jmaxm=pow(10.0, minf((-1.50+0.41*log10(Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(Pmass*1000.0)+0.32*log10(SLA)))); //  this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
            float Vcmax=Vcmaxm*LMA;
            float Jmax=Jmaxm*LMA;
            float Narea = LMA * Nmass;
            float Parea = LMA * Pmass;
            
            float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2));
            
            /* now calculate GPP and NPP */
            GettreeNPPGPP(treeGPP, treeNPP, Vcmax, Jmax, Rdark, site_census, crown_displacement, wsg, height, CR, CD, dbh, LAI);

            if(treeNPP < 0.0) treeNPP = 0.0;
        }
        
        T[site_census].t_GPP = treeGPP;
        T[site_census].t_NPP = treeNPP;
        
        if(treeNPP <= 0.0 && dbh >= dbh_limitcarbon) nbtrees_carbonstarv++;

    }

    if(nbtrees_total_threshold > 0){
        carbonstarv = float(nbtrees_carbonstarv)/float(nbtrees_total_threshold);
    } else {
        carbonstarv = 0.0;
    }
    
}
                        
int GetNbTreesThreshold(vector<int> &list_trees_potential, vector<int> &list_trees_examined, vector<float> &heights_list_trees_potential, int site_tree, float dbh, float height, float CR, int col_crowncenter, int row_crowncenter, bool addtree, int &functioncalls){
    int nbtrees_threshold = 0;
    functioncalls++;
    
    /* we only enter the loop if it is possible that there are trees > threshold below this tree */
    float maxCR = maxCR_perheight[height];
    float maxdist = CR + maxCR;
    
    float maxarea = PI * maxdist * maxdist;                  // floor of crown_area to bound area accumulation
    int maxarea_int = int(maxarea);                                   // floor of crown_area to bound area accumulation
    maxarea_int = max(maxarea_int,1);                                 // minimum area of crown (1), important for radii below ~0.5
    maxarea_int = min(maxarea_int,7853);                              // maximum area of crown (radius 50)
    //maxarea_int = min(maxarea_int,1963);                              // maximum area of crown (radius 25)
    
    for(int i = 0; i < maxarea_int; i++){
        int site_relative = LookUp_neighbour_site[i];

        int row = row_crowncenter + site_relative/101 - 50;
        int col = col_crowncenter + site_relative%101 - 50;
        
//        int site_relative = LookUp_Crown_site[i];
//        
//        int row = row_crowncenter + site_relative/51 - 25;
//        int col = col_crowncenter + site_relative%51 - 25;
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
            int site_test = col + row * cols;
            
            //if(site_test == 75361) cout << "Recursion site: " << site_test << endl;
            
            //cout << site_tree << " within loop" << endl;
            
            /* the tree needs to climb several hurdles before being evaluated, i.e. is it alive, above a certain height threshold, and within the radius */
            
            /* we check whether the tree references itself, since crown site and tree site do not correspond */
            /* self-reference would be prohibited later as well through the height_test (the tree cannot be smaller than itself), but this is cleaner and slightly better performance-wise */
            if(site_test != site_tree && T[site_test].t_age > 0){
                
                float height_test = T[site_test].t_Tree_Height;
                if(height_test < height && height_test >= heightminimum_threshold){
                    float dbh_test = T[site_test].t_dbh ;
                    if(dbh > dbh_limitcarbon || dbh_test >= dbh_limitcarbon){
                    
                        /* only check if the tree has not been examined before */
                        int examined = 0;
                        for(int j=0; j < list_trees_examined.size(); j++){
                            if(list_trees_examined[j] == site_test) examined = 1;
                        }
                        if(examined == 0){
                            /* now we check whether tree is within radius */
                            /* this is the only check that can be performed many times, because we do not know whether a tree might be within another tree's radius */
                            
                            float CR_test = T[site_test].t_Crown_Radius;
                            int crown_displacement_test = T[site_test].t_CrownDisplacement;
                            int site_crown_test = site_test + crown_displacement_test;
                            int col_crowncenter_test = site_crown_test%cols;
                            int row_crowncenter_test = site_crown_test/cols;
                            
                            int distance_crowns_squared = (col_crowncenter_test - col_crowncenter) * (col_crowncenter_test - col_crowncenter) + (row_crowncenter_test - row_crowncenter) * (row_crowncenter_test - row_crowncenter);
                            int distance_crownoverlap_squared = int((CR + CR_test) * (CR + CR_test)) + 1; /* we add a bit of safety margin */
                            
                            //                        if(site_test == 104153) cout << "Recursion site: " << site_test << " dbh_test: " << dbh_test << " height_test: " << height_test << " CR_test: " << CR_test << " current site: " << site_tree << " height: " << height << " CR: " << CR  << " minimum threshold: " << heightminimum_threshold << " distance_trees_squared: " << distance_crowns_squared << " CR distance: " << distance_crownoverlap_squared << endl;
                            //
                            /* the trees need to be smaller than the current tree to be affected, taller than the minimum height threshold and have overlapping radii */
                            if(distance_crowns_squared <= distance_crownoverlap_squared){
                                //cout << "Recursion site: " << site_test << " col: " << col_crowncenter << " row: " << row_crowncenter << " col2: " << col_crowncenter_test << " row2: " << row_crowncenter_test << " dbh_test: " << dbh_test << " height_test: " << height_test << " CR_test: " << CR_test << " current site: " << site_tree << " height: " << height << " CR: " << CR  << " minimum threshold: " << heightminimum_threshold << " distance_trees_squared: " << distance_crowns_squared << " CR distance: " << distance_crownoverlap_squared << endl;
                                
                                bool addtree_test = 1;
                                //float dbh_test = T[site_test].t_dbh;
                                int nbtrees_threshold_test = GetNbTreesThreshold(list_trees_potential, list_trees_examined, heights_list_trees_potential, site_test, dbh_test, height_test, CR_test, col_crowncenter_test, row_crowncenter_test, addtree_test, functioncalls);
                                nbtrees_threshold += nbtrees_threshold_test;
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    list_trees_examined.push_back(site_tree);
    
    if(addtree == 1){
        if(dbh >= dbh_limitcarbon){
            nbtrees_threshold++;
        }
        
        if(nbtrees_threshold > 0){
            list_trees_potential.push_back(site_tree);
            heights_list_trees_potential.push_back(height);
        }
    }
    return(nbtrees_threshold);
}
    
    
void GetTreeList(int site_tree1, int site_tree2, float CR1_previous, float CR1_new, float CR2_previous, float CR2_new, float CD1_previous, float CD1_new, float CD2_previous, float CD2_new, float height1_previous, float height1_new, float height2_previous, float height2_new, int crown_displacement1_previous, int crown_displacement1_new, int crown_displacement2_previous, int crown_displacement2_new, vector<int> &list_trees){
    /* this is the version of the function where tree crowns are exchanged */
    /* we get a list of all trees that are potentially affected, including the changed trees with their new dimensions */
    /* the latter is important, since we will base our calculations on the new dimensions */
    
    /* in a first round, we exclude trees:
     a) that are the trees that have been exchanged (1-2 trees)
     b) that are smaller than the dbh threshold (typically 0.1m)
     c) whose crown base is above the maximum height of the changed voxels (> height_max), i.e. the height until which the tree or trees influence the voxel field (i.e. highest height where LAI was allocated)
     d) another potential criteria would be overlap with crown radius, but this is tricky, since changes in LAI will likely provoke trickle-down effects (i.e. LAI changes in crowns below that spread the area of influence further). At the moment we operate with the maximum height as area of influence, since this will incorporate the increase in trickle-down effects for larger trees, but that this works needs to be verified.
     */

    /* we create two separate lists that we merge later, after separately sorting them */
    vector<int> list_trees_potential1;
    vector<int> list_trees_examined1;
    vector<float> heights_list_trees_potential1;
    
    vector<int> list_trees_potential2;
    vector<int> list_trees_examined2;
    vector<float> heights_list_trees_potential2;
    
    int height_max1 = int(maxf(height1_previous, height1_new));
    int height_max2 = int(maxf(height2_previous, height2_new));
    
    list_trees_potential1.reserve(height_max1 * 10);
    list_trees_examined1.reserve(height_max1 * 25);
    heights_list_trees_potential1.reserve(height_max1 * 10);

    list_trees_potential2.reserve(height_max2 * 10);
    list_trees_examined2.reserve(height_max2 * 25);
    heights_list_trees_potential2.reserve(height_max2 * 10);
    
    float dbh1 = T[site_tree1].t_dbh;
    float dbh2 = T[site_tree2].t_dbh;
    
    int site_crown1_previous = site_tree1 + crown_displacement1_previous;
    int site_crown1_new = site_tree1 + crown_displacement1_new;
    int site_crown2_previous = site_tree2 + crown_displacement2_previous;
    int site_crown2_new = site_tree2 + crown_displacement2_new;

    int col_crowncenter1_previous = site_crown1_previous%cols;
    int col_crowncenter1_new = site_crown1_new%cols;
    int col_crowncenter2_previous = site_crown2_previous%cols;
    int col_crowncenter2_new = site_crown2_new%cols;

    int row_crowncenter1_previous = site_crown1_previous/cols;
    int row_crowncenter1_new = site_crown1_new/cols;
    int row_crowncenter2_previous = site_crown2_previous/cols;
    int row_crowncenter2_new = site_crown2_new/cols;

    bool addtree_previous = 0;
    bool addtree_new = 1;
    
    /* now get all the potentially affected trees for both sites */
    /* first for the old configuration */
    int functioncalls1_previous = 0;
    int functioncalls2_previous = 0;

    int trees1_previous = GetNbTreesThreshold(list_trees_potential1,list_trees_examined1, heights_list_trees_potential1, site_tree1, dbh1, height1_previous, CR1_previous, col_crowncenter1_previous, row_crowncenter1_previous, addtree_previous, functioncalls1_previous);
    int trees2_previous = GetNbTreesThreshold(list_trees_potential2, list_trees_examined2, heights_list_trees_potential2, site_tree2, dbh2, height2_previous, CR2_previous, col_crowncenter2_previous, row_crowncenter2_previous, addtree_previous, functioncalls2_previous);

    /* now for the new configuration */
    /* we avoid overlap with the old configuration by temporarily "turning off" the old values */
    T[site_tree1].t_age = 0;
    T[site_tree2].t_age = 0;
    
    int functioncalls1_new = 0;
    int functioncalls2_new = 0;

    int trees1_new = GetNbTreesThreshold(list_trees_potential1, list_trees_examined1, heights_list_trees_potential1, site_tree1, dbh1, height1_new, CR1_new, col_crowncenter1_new, row_crowncenter1_new, addtree_new, functioncalls1_new );
    int trees2_new = GetNbTreesThreshold(list_trees_potential2, list_trees_examined2, heights_list_trees_potential2, site_tree2, dbh2, height2_new, CR2_new, col_crowncenter2_new, row_crowncenter2_new, addtree_new, functioncalls2_new);
    /* reactivate */
    T[site_tree1].t_age = 1;
    T[site_tree2].t_age = 1;

    /* first we test for overlap */
    int merge = 0;
    for(int i = 0; i < list_trees_potential1.size(); i++){
        for(int j = 0; j < list_trees_potential2.size(); j++){
            if(list_trees_potential1[i] == list_trees_potential2[j]) merge = 1;
        }
    }
    
    vector<int> list_trees_potential;
    vector<float> heights_list_trees_potential;
    
    int trees_combined = int(list_trees_potential1.size() + list_trees_potential2.size());
    
    list_trees_potential.reserve(trees_combined);
    heights_list_trees_potential.reserve(trees_combined);
    
    if(merge == 1){
        /* if they are to be merged, we simply add them up and sort them together */
        for(int i = 0; i < list_trees_potential1.size();i++){
            list_trees_potential.push_back(list_trees_potential1[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential1[i]);
        }
        for(int i = 0; i < list_trees_potential2.size();i++){
            list_trees_potential.push_back(list_trees_potential2[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential2[i]);
        }
        
        for(int i = 0 ;i < list_trees_potential.size(); i++){
            for(int tree = 0; tree < list_trees_potential.size()-i-1; tree++){
                int site_current = list_trees_potential[tree];
                int site_next = list_trees_potential[tree+1];
                float height_current = heights_list_trees_potential[tree];
                float height_next = heights_list_trees_potential[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential[tree] = site_next;
                    list_trees_potential[tree+1] = site_current;
                    heights_list_trees_potential[tree] = height_next;
                    heights_list_trees_potential[tree+1] = height_current;
                }
            }
        }
    } else {
        /* otherwise we sort them separately and then add them up */
        for(int i = 0 ;i < list_trees_potential1.size(); i++){
            for(int tree = 0; tree < list_trees_potential1.size()-i-1; tree++){
                int site_current = list_trees_potential1[tree];
                int site_next = list_trees_potential1[tree+1];
                float height_current = heights_list_trees_potential1[tree];
                float height_next = heights_list_trees_potential1[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential1[tree] = site_next;
                    list_trees_potential1[tree+1] = site_current;
                    heights_list_trees_potential1[tree] = height_next;
                    heights_list_trees_potential1[tree+1] = height_current;
                }
            }
        }
        
        for(int i = 0 ;i < list_trees_potential2.size(); i++){
            for(int tree = 0; tree < list_trees_potential2.size()-i-1; tree++){
                int site_current = list_trees_potential2[tree];
                int site_next = list_trees_potential2[tree+1];
                float height_current = heights_list_trees_potential2[tree];
                float height_next = heights_list_trees_potential2[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential2[tree] = site_next;
                    list_trees_potential2[tree+1] = site_current;
                    heights_list_trees_potential2[tree] = height_next;
                    heights_list_trees_potential2[tree+1] = height_current;
                }
            }
        }
        
        /* now add them up */
        for(int i = 0; i < list_trees_potential1.size();i++){
            list_trees_potential.push_back(list_trees_potential1[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential1[i]);
        }
        for(int i = 0; i < list_trees_potential2.size();i++){
            list_trees_potential.push_back(list_trees_potential2[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential2[i]);
        }
    }

    
    /* now we have a list of trees that could be potentially influenced */
    /* the problem is that there are potentially a lot of duplicates in it, and we have to remove those */
    list_trees.reserve(list_trees_potential.size());
    
    int site_tree_previous = list_trees_potential[0];
    list_trees.push_back(site_tree_previous);
    
    for(int i = 1 ;i < list_trees_potential.size(); i++){
        int site_tree = list_trees_potential[i];
        if(site_tree != site_tree_previous){
            list_trees.push_back(site_tree);
            site_tree_previous = site_tree;
        }
    }
}


void GetTreeList(int site_original, int site_new, vector<int> &list_trees){
    /* this is the version of the function where a tree's position is changed */
    /* we get a list of all trees that are potentially affected, including the changed tree */
    /* the latter is important, since we will base our calculations on the new dimensions */
    /* in the case of shifting trees, since all dimensions remain the same, there is on need for extra information, however */
    
    /* in a first round, we exclude trees:
     a) that are the trees that have been exchanged (1-2 trees)
     b) that are smaller than the dbh threshold (typically 0.1m)
     c) whose crown base is above the maximum height of the changed voxels (> height_max), i.e. the height until which the tree or trees influence the voxel field (i.e. highest height where LAI was allocated)
     d) another potential criteria would be overlap with crown radius, but this is tricky, since changes in LAI will likely provoke trickle-down effects (i.e. LAI changes in crowns below that spread the area of influence further). At the moment we operate with the maximum height as area of influence, since this will incorporate the increase in trickle-down effects for larger trees, but that this works needs to be verified.
     */

    /* first get the tree dimensions */
    float dbh = T[site_original].t_dbh;
    float height = T[site_original].t_Tree_Height;
    float CR = T[site_original].t_Crown_Radius;
    int height_max = int(height);
    
    /* create tree array */
    /* we create two separate lists that we merge later, after separately sorting them */
    vector<int> list_trees_potential1;
    vector<int> list_trees_examined1;
    vector<float> heights_list_trees_potential1;
    
    vector<int> list_trees_potential2;
    vector<int> list_trees_examined2;
    vector<float> heights_list_trees_potential2;
    
    list_trees_potential1.reserve(height_max * 5);
    list_trees_examined1.reserve(height_max * 25);
    heights_list_trees_potential1.reserve(height_max * 5);
    
    list_trees_potential2.reserve(height_max * 5);
    list_trees_examined2.reserve(height_max * 25);
    heights_list_trees_potential2.reserve(height_max * 5);
    
    int col_crowncenter_original = site_original%cols;
    int row_crowncenter_original = site_original/cols;
    
    int col_crowncenter_new = site_new%cols;
    int row_crowncenter_new = site_new/cols;
    
    bool addtree_original = 0;
    bool addtree_new = 1;
    
    /* temporarily deactivate the tree */
    T[site_original].t_age = 0;
    
    int functioncalls_original = 0;
    int functioncalls_new = 0;

    GetNbTreesThreshold(list_trees_potential1, list_trees_examined1, heights_list_trees_potential1, site_original, dbh, height, CR, col_crowncenter_original, row_crowncenter_original, addtree_original, functioncalls_original);
    GetNbTreesThreshold(list_trees_potential2, list_trees_examined2, heights_list_trees_potential2, site_new, dbh, height, CR, col_crowncenter_new, row_crowncenter_new, addtree_new, functioncalls_new);
    
    T[site_original].t_age = 1;

    /* first we test for overlap */
    int merge = 0;
    for(int i = 0; i < list_trees_potential1.size(); i++){
        for(int j = 0; j < list_trees_potential2.size(); j++){
            if(list_trees_potential1[i] == list_trees_potential2[j]) merge = 1;
        }
    }
    
    vector<int> list_trees_potential;
    vector<float> heights_list_trees_potential;
    
    int trees_combined = int(list_trees_potential1.size() + list_trees_potential2.size());
    
    list_trees_potential.reserve(trees_combined);
    heights_list_trees_potential.reserve(trees_combined);
    
    if(merge == 1){
        /* if they are to be merged, we simply add them up and sort them together */
        for(int i = 0; i < list_trees_potential1.size();i++){
            list_trees_potential.push_back(list_trees_potential1[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential1[i]);
        }
        for(int i = 0; i < list_trees_potential2.size();i++){
            list_trees_potential.push_back(list_trees_potential2[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential2[i]);
        }
        
        for(int i = 0 ;i < list_trees_potential.size(); i++){
            for(int tree = 0; tree < list_trees_potential.size()-i-1; tree++){
                int site_current = list_trees_potential[tree];
                int site_next = list_trees_potential[tree+1];
                float height_current = heights_list_trees_potential[tree];
                float height_next = heights_list_trees_potential[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential[tree] = site_next;
                    list_trees_potential[tree+1] = site_current;
                    heights_list_trees_potential[tree] = height_next;
                    heights_list_trees_potential[tree+1] = height_current;
                }
            }
        }
    } else {
        /* otherwise we sort them separately and then add them up */
        for(int i = 0 ;i < list_trees_potential1.size(); i++){
            for(int tree = 0; tree < list_trees_potential1.size()-i-1; tree++){
                int site_current = list_trees_potential1[tree];
                int site_next = list_trees_potential1[tree+1];
                float height_current = heights_list_trees_potential1[tree];
                float height_next = heights_list_trees_potential1[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential1[tree] = site_next;
                    list_trees_potential1[tree+1] = site_current;
                    heights_list_trees_potential1[tree] = height_next;
                    heights_list_trees_potential1[tree+1] = height_current;
                }
            }
        }
        
        for(int i = 0 ;i < list_trees_potential2.size(); i++){
            for(int tree = 0; tree < list_trees_potential2.size()-i-1; tree++){
                int site_current = list_trees_potential2[tree];
                int site_next = list_trees_potential2[tree+1];
                float height_current = heights_list_trees_potential2[tree];
                float height_next = heights_list_trees_potential2[tree+1];
                
                if(height_current < height_next){
                    list_trees_potential2[tree] = site_next;
                    list_trees_potential2[tree+1] = site_current;
                    heights_list_trees_potential2[tree] = height_next;
                    heights_list_trees_potential2[tree+1] = height_current;
                }
            }
        }
        
        /* now add them up */
        for(int i = 0; i < list_trees_potential1.size();i++){
            list_trees_potential.push_back(list_trees_potential1[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential1[i]);
        }
        for(int i = 0; i < list_trees_potential2.size();i++){
            list_trees_potential.push_back(list_trees_potential2[i]);
            heights_list_trees_potential.push_back(heights_list_trees_potential2[i]);
        }
    }
    
    /* now we have a list of trees that could be potentially influenced */
    /* the problem is that there are still duplicates in it, and we have to remove those */
    list_trees.reserve(list_trees_potential.size());
    
    if(list_trees_potential.size() > 0){
        int site_tree_previous = list_trees_potential[0];
        list_trees.push_back(site_tree_previous);
        
        for(int i = 1 ;i < list_trees_potential.size(); i++){
            int site_tree = list_trees_potential[i];
            if(site_tree != site_tree_previous){
                list_trees.push_back(site_tree);
                site_tree_previous = site_tree;
            }
        }
    }
}

void GetTreeList(int site_tree, float CR_previous, float CR_new, float CD_previous, float CD_new, float height_previous, float height_new, float crown_displacement_previous, float crown_displacement_new, vector<int> &list_trees){
    /* this is the version of the function where a tree crown is only modified */
    /* we get a list of all trees that are potentially affected, including the changed tree with its new dimensions */
    /* the latter is important, since we will base our calculations on the new dimensions */
    
    /* in a first round, we exclude trees:
     a) that are the trees that have been exchanged (1-2 trees)
     b) that are smaller than the dbh threshold (typically 0.1m)
     c) whose crown base is above the maximum height of the changed voxels (> height_max), i.e. the height until which the tree or trees influence the voxel field (i.e. highest height where LAI was allocated)
     d) another potential criteria would be overlap with crown radius, but this is tricky, since changes in LAI will likely provoke trickle-down effects (i.e. LAI changes in crowns below that spread the area of influence further). At the moment we operate with the maximum height as area of influence, since this will incorporate the increase in trickle-down effects for larger trees, but that this works needs to be verified.
     */
    
    /* create tree array */
    vector<int> list_trees_potential;
    vector<int> list_trees_examined;
    vector<float> heights_list_trees_potential;
    
    int height_max = int(maxf(height_previous, height_new));
    
    list_trees_potential.reserve(height_max * 10);
    list_trees_examined.reserve(height_max * 50);
    heights_list_trees_potential.reserve(height_max * 10);
    
    float dbh = T[site_tree].t_dbh;
    
    int site_crown_previous = site_tree + crown_displacement_previous;
    int site_crown_new = site_tree + crown_displacement_new;

    int col_crowncenter_previous = site_crown_previous%cols;
    int col_crowncenter_new = site_crown_new%cols;

    int row_crowncenter_previous = site_crown_previous/cols;
    int row_crowncenter_new = site_crown_new/cols;

    bool addtree_previous = 0;
    bool addtree_new = 1;
    
    /* now get all the potentially affected trees for both sites */
    /* first for the old configuration */
    int functioncalls_previous = 0;
    int trees_previous = GetNbTreesThreshold(list_trees_potential, list_trees_examined, heights_list_trees_potential, site_tree, dbh, height_previous, CR_previous, col_crowncenter_previous, row_crowncenter_previous, addtree_previous, functioncalls_previous);

    /* now for the new configuration */
    /* we avoid overlap with the old configuration by temporarily "turning off" the old values */
    T[site_tree].t_age = 0;
    
    int functioncalls_new = 0;
    int trees_new = GetNbTreesThreshold(list_trees_potential, list_trees_examined, heights_list_trees_potential, site_tree, dbh, height_new, CR_new, col_crowncenter_new, row_crowncenter_new, addtree_new, functioncalls_new);

    /* reactivate */
    T[site_tree].t_age = 1;
    
    //cout << "Results from recursive algorithm at site: " << site_tree << ". Trees_prev: " << trees_previous << " Trees_new: " << trees_new << " total size of array: " << list_trees_potential.size() << endl;
    
    /* now we sort the list of trees  */
    for(int i = 0 ;i < list_trees_potential.size(); i++){
        for(int tree = 0; tree < list_trees_potential.size()-i-1; tree++){
            int site_current = list_trees_potential[tree];
            int site_next = list_trees_potential[tree+1];
            float height_current = heights_list_trees_potential[tree];
            float height_next = heights_list_trees_potential[tree+1];
            
            if(height_current < height_next){
                list_trees_potential[tree] = site_next;
                list_trees_potential[tree+1] = site_current;
                heights_list_trees_potential[tree] = height_next;
                heights_list_trees_potential[tree+1] = height_current;
            }
        }
    }
    
    /* now we have a list of trees that could be potentially influenced */
    /* the problem is that there are potentially a lot of duplicates in it, and we have to remove those */
    list_trees.reserve(list_trees_potential.size());
    
    int site_tree_previous = list_trees_potential[0];
    list_trees.push_back(site_tree_previous);
    
    for(int i = 1 ;i < list_trees_potential.size(); i++){
        int site_tree = list_trees_potential[i];
        if(site_tree != site_tree_previous){
            list_trees.push_back(site_tree);
            site_tree_previous = site_tree;
        }
    }
}

void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site1, float height1_new, float CR1_new, float CD1_new, int crown_displacement1_new, float Nmass1_new, float Pmass1_new, float LMA1_new, int site2, float height2_new, float CR2_new, float CD2_new, int crown_displacement2_new, float Nmass2_new, float Pmass2_new, float LMA2_new){
    /* this function updates the tree's leaves, i.e. modifying LAI3D as well as the respective tree variables */
    /* the two tree sites where change has happened are included in the normal iteration, just with a slight difference: we remove leaves at a different place than where we allocate (typically) */
    
    list_trees_properties.clear();
    
    /* remove all trees first, this is in analogy with the full routine */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        float dbh = T[site_tree_update].t_dbh;
        float height = T[site_tree_update].t_Tree_Height;
        float CR = T[site_tree_update].t_Crown_Radius;
        float CD = T[site_tree_update].t_Crown_Depth;
        float LAI = T[site_tree_update].t_LAI;
        int crown_displacement = T[site_tree_update].t_CrownDisplacement;
        int site_crown = site_tree_update + crown_displacement;
        CalcLAItrial(1, site_crown, height, CR, CD, LAI,dbh);
    }
    
    /* now run through all the selected trees, and update the LAI first */
    /* only update the LAI, to have a consistent approach with UpdateLeaves() */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        
        float dbh = T[site_tree_update].t_dbh;
        float height, CR, CD, crown_displacement;
        if(site_tree_update == site1){
            height = height1_new;
            CR = CR1_new;
            CD = CD1_new;
            crown_displacement = crown_displacement1_new;
        } else if(site_tree_update == site2){
            height = height2_new;
            CR = CR2_new;
            CD = CD2_new;
            crown_displacement = crown_displacement2_new;
        } else{
            height = T[site_tree_update].t_Tree_Height;
            CR = T[site_tree_update].t_Crown_Radius;
            CD = T[site_tree_update].t_Crown_Depth;
            crown_displacement = T[site_tree_update].t_CrownDisplacement;
        }
        float LAI_previous = T[site_tree_update].t_LAI;
        float LAImax = T[site_tree_update].t_LAImax;
        float LAIabove_previous = T[site_tree_update].t_LAIabove;
        int site_crown = site_tree_update + crown_displacement;

        float LAIabove = UpdateLAIabove_eff(site_crown, height, CR, CD, dbh);
 
        /* calculate the LAI */
        /* we only do this, if the LAIabove has changed, and by default for the updated trees (probably not necessary, as this will introduce changes by default) */
        float LAI;
        if((LAIabove != LAIabove_previous) || site_tree_update == site1 || site_tree_update == site2){
        //if((int(LAIabove * 1000.0) != int(LAIabove_previous * 1000.0)) || site_tree_update == site1 || site_tree_update == site2){
            if(LAIabove <= LAImax){
                float LAItree_max_light = LAImax - LAIabove;
                float LAItree_max_physiology = LAImax;
                LAI = minf(LAItree_max_light, LAItree_max_physiology);
            } else {
                LAI = 0.0;
            }
        } else {
            LAI = LAI_previous;
        }
        
        /* update the LAI field */
        CalcLAItrial(0, site_crown, height, CR, CD, LAI,dbh);

        list_trees_properties.push_back(LAIabove);
        list_trees_properties.push_back(LAI);
        list_trees_properties.push_back(0.0); // adding GPP pro forma
        list_trees_properties.push_back(0.0); // adding NPP pro forma
    }
    
    /* now update the GPP and NPP */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
    
        float dbh = T[site_tree_update].t_dbh;
        
        float height, CR, CD, crown_displacement;
        if(site_tree_update == site1){
            height = height1_new;
            CR = CR1_new;
            CD = CD1_new;
            crown_displacement = crown_displacement1_new;
        } else if(site_tree_update == site2){
            height = height2_new;
            CR = CR2_new;
            CD = CD2_new;
            crown_displacement = crown_displacement2_new;
        } else{
            height = T[site_tree_update].t_Tree_Height;
            CR = T[site_tree_update].t_Crown_Radius;
            CD = T[site_tree_update].t_Crown_Depth;
            crown_displacement = T[site_tree_update].t_CrownDisplacement;
        }

        float treeNPP_previous = T[site_tree_update].t_NPP;

        float treeGPP = 0.0;
        float treeNPP = 0.0;
        
        float LAI = list_trees_properties[1 + i * 4];

        /* we only need to update the GPP/NPP if LAI > 0.0, otherwise trivial: 0.0 */
        if(LAI > 0.0){
            float Nmass, Pmass, LMA;
            if(site_tree_update == site1){
                Nmass = Nmass1_new;
                Pmass = Pmass1_new;
                LMA = LMA1_new;
            } else if(site_tree_update == site2){
                Nmass = Nmass2_new;
                Pmass = Pmass2_new;
                LMA = LMA2_new;
            } else {
                Nmass = T[site_tree_update].t_Nmass;
                Pmass = T[site_tree_update].t_Pmass;
                LMA = T[site_tree_update].t_LMA;
            }
            
            float wsg = T[site_tree_update].t_wsg;
            float SLA=10000.0/LMA;
            float Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
            float Jmaxm=pow(10.0, minf((-1.50+0.41*log10(Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(Pmass*1000.0)+0.32*log10(SLA)))); // this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
            float Vcmax=Vcmaxm*LMA;
            float Jmax=Jmaxm*LMA;
            float Narea = LMA * Nmass;
            float Parea = LMA * Pmass;
            
            float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2));
            
            /* now calculate GPP and NPP */
            GettreeNPPGPP(treeGPP, treeNPP, Vcmax, Jmax, Rdark, site_tree_update, crown_displacement, wsg, height, CR, CD, dbh, LAI);

            if(treeNPP < 0.0) treeNPP = 0.0;
        }
        
        list_trees_properties[2 + i * 4] = treeGPP;
        list_trees_properties[3 + i * 4] = treeNPP;
        
        if(dbh >= dbh_limitcarbon){
            if(treeNPP_previous <= 0.0){
                if(treeNPP > 0.0){
                    nbtrees_carbonstarv--;

                } else {

                }
            } else {
                if(treeNPP <= 0.0){
                    nbtrees_carbonstarv++;
                } else {
                }
            }
        } else {
        }
    }

    if(nbtrees_total_threshold > 0){
        carbonstarv = float(nbtrees_carbonstarv)/float(nbtrees_total_threshold);
    } else {
        carbonstarv = 0.0;
    }
    
}


void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site_original, int site_new){
    /* this function updates the tree's leaves, i.e. modifying LAI3D as well as the respective tree variables */
    /* the two tree sites where change has happened are included in the normal iteration, just with a slight difference: we remove leaves at a different place than where we allocate (typically) */
    
    list_trees_properties.clear();
    
    /* remove all trees first, this is in analogy with the full routine */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        float dbh = T[site_tree_update].t_dbh;
        float height = T[site_tree_update].t_Tree_Height;
        float CR = T[site_tree_update].t_Crown_Radius;
        float CD = T[site_tree_update].t_Crown_Depth;
        int crown_displacement = T[site_tree_update].t_CrownDisplacement;
        int site_crown = site_tree_update + crown_displacement;
        float LAI = T[site_tree_update].t_LAI;
        CalcLAItrial(1, site_crown, height, CR, CD, LAI, dbh);
    }
    
    /* now run through all the selected trees, and update the LAI */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        
        float dbh = T[site_tree_update].t_dbh;
        float height = T[site_tree_update].t_Tree_Height;
        float CR = T[site_tree_update].t_Crown_Radius;
        float CD = T[site_tree_update].t_Crown_Depth;
        float crown_displacement = T[site_tree_update].t_CrownDisplacement;
        
        /* the reference to the tree is the old site (site_original), but the crown where we now allocate is the new site (site_new) */
        int site_crown;
        if(site_tree_update == site_original) site_crown = site_new + crown_displacement;
        else  site_crown = site_tree_update + crown_displacement;
        
        float LAImax = T[site_tree_update].t_LAImax;
        float LAIabove_previous = T[site_tree_update].t_LAIabove;
        float LAI_previous = T[site_tree_update].t_LAI;
        

        float LAIabove = UpdateLAIabove_eff(site_crown, height, CR, CD, dbh);

        /* calculate the LAI */
        /* two cases where we update: the tree has been moved (probably no update necessary, as this will induce changes by default) or there is change in the light environment above the crown */
        float LAI;
        if(LAIabove != LAIabove_previous || site_tree_update == site_original){
        //if((int(LAIabove * 1000.0) != int(LAIabove_previous * 1000.0)) || site_tree_update == site_original){
            if(LAIabove <= LAImax){
                float LAItree_max_light = LAImax - LAIabove;
                float LAItree_max_physiology = LAImax;
                LAI = minf(LAItree_max_light, LAItree_max_physiology);
                
                /* now update the LAI field */
                CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
            } else {
                LAI = 0.0;
            }
        } else {
            LAI = LAI_previous;
            CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
        }
        
        list_trees_properties.push_back(LAIabove);
        list_trees_properties.push_back(LAI);
        list_trees_properties.push_back(0.0); // adding GPP pro forma
        list_trees_properties.push_back(0.0); // adding NPP pro forma
    }
    
    
    /* now update the GPP and NPP */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        
        float height = T[site_tree_update].t_Tree_Height;
        float CR = T[site_tree_update].t_Crown_Radius;
        float CD = T[site_tree_update].t_Crown_Depth;
        float crown_displacement = T[site_tree_update].t_CrownDisplacement;
        float dbh = T[site_tree_update].t_dbh;
        
        float LAI = list_trees_properties[1 + i * 4];
        
        float treeNPP_previous = T[site_tree_update].t_NPP;
        
        float treeGPP = 0.0;
        float treeNPP = 0.0;
        
        /* we only need to update the GPP/NPP if LAI > 0.0, otherwise trivial: 0.0 */
        if(LAI > 0.0){
                
            float Nmass = T[site_tree_update].t_Nmass;
            float Pmass = T[site_tree_update].t_Pmass;
            float LMA = T[site_tree_update].t_LMA;
            
            float wsg = T[site_tree_update].t_wsg;
            
            float SLA=10000.0/LMA;
            float Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
            float Jmaxm=pow(10.0, minf((-1.50+0.41*log10(Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(Pmass*1000.0)+0.32*log10(SLA)))); // this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
            float Vcmax=Vcmaxm*LMA;
            float Jmax=Jmaxm*LMA;
            float Narea = LMA * Nmass;
            float Parea = LMA * Pmass;
            
            float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2));
            
            /* again, we need to differentiate: if we look at the tree that has been moved, then we need to use the new site */
            if(site_tree_update == site_original) GettreeNPPGPP(treeGPP, treeNPP, Vcmax, Jmax, Rdark, site_new, crown_displacement, wsg, height, CR, CD, dbh, LAI);

            else GettreeNPPGPP(treeGPP, treeNPP, Vcmax, Jmax, Rdark, site_tree_update, crown_displacement, wsg, height, CR, CD, dbh, LAI);
            
            if(treeNPP < 0.0) treeNPP = 0.0;
        }
        
        list_trees_properties[2 + i * 4] = treeGPP;
        list_trees_properties[3 + i * 4] = treeNPP;
        
        if(dbh >= dbh_limitcarbon){
            if(treeNPP_previous <= 0.0){
                if(treeNPP > 0.0) nbtrees_carbonstarv--;
            } else {
                if(treeNPP <= 0.0) nbtrees_carbonstarv++;
            }
        }
    }


    if(nbtrees_total_threshold > 0){
        carbonstarv = float(nbtrees_carbonstarv)/float(nbtrees_total_threshold);
    } else {
        carbonstarv = 0.0;
    }
}
void UpdateLeaves_fromList(vector<int> &list_trees, vector<float> &list_trees_properties, int site_tree, float height_new, float CR_new, float CD_new, int crown_displacement_new){
    /* this function updates the tree's leaves, i.e. modifying LAI3D as well as the respective tree variables */
    /* the two tree sites where change has happened are included in the normal iteration, just with a slight difference: we remove leaves at a different place than where we allocate (typically) */
    
    list_trees_properties.clear();
    
    /* remove all trees first, this is in analogy with the full routine */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        float dbh = T[site_tree_update].t_dbh;
        float height = T[site_tree_update].t_Tree_Height;
        float CR = T[site_tree_update].t_Crown_Radius;
        float CD = T[site_tree_update].t_Crown_Depth;
        int crown_displacement = T[site_tree_update].t_CrownDisplacement;
        int site_crown = site_tree_update + crown_displacement;
        float LAI = T[site_tree_update].t_LAI;
        CalcLAItrial(1, site_crown, height, CR, CD, LAI, dbh);
    }

    /* now run through all the selected trees and update LAI */
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        
        float dbh = T[site_tree_update].t_dbh;
        float height, CR, CD, crown_displacement;
        if(site_tree_update == site_tree){
            height = height_new;
            CR = CR_new;
            CD = CD_new;
            crown_displacement = crown_displacement_new;
        } else {
            height = T[site_tree_update].t_Tree_Height;
            CR = T[site_tree_update].t_Crown_Radius;
            CD = T[site_tree_update].t_Crown_Depth;
            crown_displacement = T[site_tree_update].t_CrownDisplacement;
        }
        
        int site_crown = site_tree_update + crown_displacement;
        
        float LAImax = T[site_tree_update].t_LAImax;
        float LAIabove_previous = T[site_tree_update].t_LAIabove;
        float LAI_previous = T[site_tree_update].t_LAI;
        float LAIabove = UpdateLAIabove_eff(site_crown, height, CR, CD, dbh);
        
        /* calculate the LAI */
        /* we only do this, if the LAIabove has changed, and by default for the updated trees (probably not necessary, as this will introduce changes by default) */
        float LAI;
        if((LAIabove != LAIabove_previous) || site_tree_update == site_tree){
        //if((int(LAIabove * 1000.0) != int(LAIabove_previous * 1000.0)) || site_tree_update == site_tree){
            if(LAIabove <= LAImax){
                float LAItree_max_light = LAImax - LAIabove;
                float LAItree_max_physiology = LAImax;
                LAI = minf(LAItree_max_light, LAItree_max_physiology);
                
                /* Update LAI field */
                CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
            } else {
                LAI = 0.0;
            }
        } else {
            LAI = LAI_previous;
            CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
        }

        list_trees_properties.push_back(LAIabove);
        list_trees_properties.push_back(LAI);
        list_trees_properties.push_back(0.0); // adding GPP pro forma
        list_trees_properties.push_back(0.0); // adding NPP pro forma
    }
    
    /* now update the GPP and NPP */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree_update = list_trees[i];
        
        float height, CR, CD, crown_displacement;
        if(site_tree_update == site_tree){
            height = height_new;
            CR = CR_new;
            CD = CD_new;
            crown_displacement = crown_displacement_new;
        } else {
            height = T[site_tree_update].t_Tree_Height;
            CR = T[site_tree_update].t_Crown_Radius;
            CD = T[site_tree_update].t_Crown_Depth;
            crown_displacement = T[site_tree_update].t_CrownDisplacement;
        }
        
        float dbh = T[site_tree_update].t_dbh;
        
        float LAI = list_trees_properties[1 + i * 4];
        
        float treeNPP_previous = T[site_tree_update].t_NPP;
        
        float treeGPP = 0.0;
        float treeNPP = 0.0;
        
        if(LAI > 0.0){
            float Nmass = T[site_tree_update].t_Nmass;
            float Pmass = T[site_tree_update].t_Pmass;
            float LMA = T[site_tree_update].t_LMA;
            
            float wsg = T[site_tree_update].t_wsg;
            
            float SLA=10000.0/LMA;
            float Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
            float Jmaxm=pow(10.0, minf((-1.50+0.41*log10(Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(Pmass*1000.0)+0.32*log10(SLA)))); // this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
            float Vcmax=Vcmaxm*LMA;
            float Jmax=Jmaxm*LMA;
            float Narea = LMA * Nmass;
            float Parea = LMA * Pmass;
            
            float Rdark = (1.3893 + (0.0728 * Narea) + (0.0015 * Parea) + (0.0095 * Vcmax) - (0.0358 * 26.2));
            
            GettreeNPPGPP(treeGPP, treeNPP, Vcmax, Jmax, Rdark, site_tree_update, crown_displacement, wsg, height, CR, CD, dbh, LAI);

            //cout << site_tree_update << " site_crown: " << site_crown << " LAImax: " << LAImax << " LAImaxabove: " << LAImaxabove << " LAImin: " << LAImin << " LAIabove: " << LAIabove << " LAItree_max_light: " << LAItree_max_light << " LAItree_max_physiology: " << LAItree_max_physiology << " LAItree_max: " << LAItree_max << " NPP: " << treeNPP << " LAI: " << LAI << endl;
            
            if(treeNPP < 0.0) treeNPP = 0.0;
        }
        
        list_trees_properties[2 + i * 4] = treeGPP;
        list_trees_properties[3 + i * 4] = treeNPP;
 
        if(dbh >= dbh_limitcarbon){
            if(treeNPP_previous <= 0.0){
                if(treeNPP > 0.0) nbtrees_carbonstarv--;
            } else {
                if(treeNPP <= 0.0) nbtrees_carbonstarv++;
            }
        }
    }
    
    if(nbtrees_total_threshold > 0){
        carbonstarv = float(nbtrees_carbonstarv)/float(nbtrees_total_threshold);
    } else {
        carbonstarv = 0.0;
    }
}
    
        
// TODO: these could likely be merged into a common function
void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site_tree1, float height1_new, float CR1_new, float CD1_new, int crown_displacement1_new, int site_tree2, float height2_new, float CR2_new, float CD2_new, int crown_displacement2_new){
    /* this is a function to reverse to the previous state */
    /* first we get the tree properties back, and then we recalculate the voxel field */
    /* tree order does not matter, since we are not recalculating, but simply replacing the original values */
    /* this includes also the two trees that are swapped as a trial, since all their properties are conserved */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree = list_trees[i];
        
        int crown_displacement = T[site_tree].t_CrownDisplacement;
        int site_crown = site_tree + crown_displacement;
        float dbh = T[site_tree].t_dbh;
        float height = T[site_tree].t_Tree_Height;
        float CR = T[site_tree].t_Crown_Radius;
        float CD = T[site_tree].t_Crown_Depth;
        float LAI_previous = T[site_tree].t_LAI;
        float LAI_new = list_trees_properties[1 + i * 4];
        
        if(site_tree == site_tree1){
            CalcLAItrial(1, site_tree1 + crown_displacement1_new, height1_new, CR1_new, CD1_new, LAI_new, dbh);
            CalcLAItrial(0, site_crown, height, CR, CD, LAI_previous, dbh);
        }
        else if(site_tree == site_tree2){
            CalcLAItrial(1, site_tree2 + crown_displacement2_new, height2_new, CR2_new, CD2_new, LAI_new,dbh);
            CalcLAItrial(0, site_crown, height, CR, CD, LAI_previous,dbh);
        }
        else{
            float LAIdiff = LAI_previous - LAI_new;
            
            CalcLAItrial(0, site_crown, height, CR, CD, LAIdiff, dbh);
            //CalcLAItrial(1, site_crown, height, CR, CD, LAI_new);
            //CalcLAItrial(0, site_crown, height, CR, CD, LAI_previous);
        }
    }
}

void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site_original, int site_new){
    /* this is a function to reverse to the previous state */
    /* first we get the tree properties back, and then we recalculate the voxel field */
    /* tree order does not matter, since we are not recalculating, but simply replacing the original values */
    /* this includes also the two trees that are swapped as a trial, since all their properties are conserved */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree = list_trees[i];
    
        int crown_displacement = T[site_tree].t_CrownDisplacement;
        int site_crown = site_tree + crown_displacement;
        float dbh = T[site_tree].t_dbh;
        float height = T[site_tree].t_Tree_Height;
        float CR = T[site_tree].t_Crown_Radius;
        float CD = T[site_tree].t_Crown_Depth;
        float LAI_new = list_trees_properties[1 + i * 4];
        
        /* we remove leaves */
        /* the moved tree is always saved at the original site */
        if(site_tree == site_original){
            CalcLAItrial(1, site_new + crown_displacement, height, CR, CD, LAI_new,dbh);
        }
        else{
            CalcLAItrial(1, site_crown, height, CR, CD, LAI_new, dbh);
        }
        
        /* now add the old ones again */
        float LAI = T[site_tree].t_LAI;
        CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
    }
}
  
void ReverseLeaves(vector<int> &list_trees, vector<float> &list_trees_properties,int site, float height_new, float CR_new, float CD_new, int crown_displacement_new){
    /* this is a function to reverse to the previous state */
    /* first we get the tree properties back, and then we recalculate the voxel field */
    /* tree order does not matter, since we are not recalculating, but simply replacing the original values */
    /* this includes also the two trees that are swapped as a trial, since all their properties are conserved */
    
    for(int i = 0; i < list_trees.size(); i++){
        int site_tree = list_trees[i];
        
        int crown_displacement = T[site_tree].t_CrownDisplacement;
        int site_crown = site_tree + crown_displacement;
        float dbh = T[site_tree].t_dbh;
        float height = T[site_tree].t_Tree_Height;
        float CR = T[site_tree].t_Crown_Radius;
        float CD = T[site_tree].t_Crown_Depth;
        float LAI_new = list_trees_properties[1 + i * 4];
        
        if(site_tree == site){
            CalcLAItrial(1, site_tree + crown_displacement_new, height_new, CR_new, CD_new, LAI_new, dbh);
        } else {
            CalcLAItrial(1, site_crown, height, CR, CD, LAI_new, dbh);
        }
    
        /* now add the old ones again */
        float LAI = T[site_tree].t_LAI;
        CalcLAItrial(0, site_crown, height, CR, CD, LAI, dbh);
    }
}

void FitTree(int site_tree, float height, float CR, float CD, int &success, vector<int> &CV_cumulatedemp, vector<int> &CV_cumulated){
    success = 1;
    trees_sampled++;

    /* now assess whether crown fits into canopy, and make sure we do not add in a biased fashion */
    vector<int> CVtree_layers(height_max+1,0);
    CalcCV_layers(CVtree_layers, site_tree, height, CR, CD);

    int CV_net = 0;

    for(int h = 0; h < height_max+1; h++){
        int CV_cumulatedemp_h = CV_cumulatedemp[h];
        int CV_cumulated_h = CV_cumulated[h];
        int CVtree_layers_h = CVtree_layers[h];
        
        /* only consider the areas where the tree actually contributes */
        if(CVtree_layers_h > 0){
            int CV_available_h =  CV_cumulatedemp_h - CV_cumulated_h;
            int CV_net_h;
            
            if(CV_available_h < 0){
                /* if the respective height layer is already overfilled, then all added CV voxels are counted as negative, plus the ones that were there before */
                CV_net_h = CV_available_h - CVtree_layers_h;
            } else {
                /* if the height layer still has space, we subtract the new CV voxels from it */
                /* if what remains is positive, this means that we can add all voxels, if it is negative, the rest is counted as negative too */
                int CV_available_after_h = CV_available_h - CVtree_layers_h;
                if(CV_available_after_h > 0) CV_net_h = CVtree_layers_h;
                else CV_net_h = CV_available_after_h + (CVtree_layers_h + CV_available_after_h); // this is the negative contribution (CV_available_after_h, negative according to condition), plus what remains from the original added volume
            }
            CV_net += CV_net_h;
        }
    }

    if(CV_net < 0){
        success = 0;
        //cout << "NO SUCCESS CV_net " << endl;
    } else {
        for(int h = 0; h < height_max+1; h++){
            CV_cumulated[h] += CVtree_layers[h];
        }
    }
}
    
int CalcCV_generalized(float height, float CR, float CD){
    int volume = 0;
    
    int crown_base = int(height - CD),
    crown_top = int(height);
    if(CD > 3.0){
        float crownshell_base = height - CD + 1.5;
        float crown_extent = height - crownshell_base;
        float crown_slope = CR * (1.0 - shape_crown) / (crown_extent - 0.5);
        /* area for the innermost (and innermost) cylinder: this is the part of the crown that keeps pushing out and creates new layers */
        int crown_extent_toplayer = int(crown_extent - 0.5);                 /* the extent gives the number of full layers both upwards and downwards not including the base layer (i.e. half of it, ~ 0.5) */
        float radius_innermost = CR - crown_slope * float(crown_extent_toplayer);
        
        int nblayers_fromtop = 4;
        
        for(int layer_fromtop = 0; layer_fromtop < nblayers_fromtop; layer_fromtop++){
            
            /* calculate the innermost/uppermost cylinder */
            int crown_area_previous = 0;
            
            /* crown area */
            int crown_area_innermost_int = GetCrownIntarea(radius_innermost);
            
            /* height at which we allocate */
            int height_innermost = crown_top - layer_fromtop;
            
            for(int i = 0; i < crown_area_innermost_int; i++){
                volume++;
            }
            
            crown_area_previous = crown_area_innermost_int;
            
            /* now loop through the outer crown cylinders */
            int height_toplayer = int(crownshell_base + 0.5 + float(crown_extent_toplayer)) - layer_fromtop;
            int height_baselayer = int(crownshell_base + 1.5) - layer_fromtop;
        
            if(layer_fromtop == 3){
                for(int h = height_innermost-1; h >= height_baselayer; h--){
                    for(int i = 0; i < crown_area_innermost_int; i++){
                        volume++;
                    }
                }
            }

            for(int h = height_toplayer; h >= height_baselayer; h--){
                /* calculating the radius of the current layer depending on the respective slopes */
                float radius_height = CR - crown_slope * (h - height_baselayer);    /* for the lowest layer, i.e. h == height_baselayer, radius = CR */
                /* crown area */
                int crown_area_int = GetCrownIntarea(radius_height);
                
                for(int i = crown_area_previous; i < crown_area_int; i++){
                    volume++;
                }
                if(layer_fromtop < 3) crown_area_previous = crown_area_int;
            }
        }
    } else {
        int crown_area_int = GetCrownIntarea(CR);
        for(int h=crown_top;h>=crown_base;h--){
            for(int i = 0; i < crown_area_int; i++){
                volume++;
            }
        }
    }
    return(volume);
}
    
void CalcCV_layers(vector<int> &CVtree_layers, int site_crown, float height, float CR, float CD){
    /* we will calculate this numerically, as if we looped through the voxel field */
    /* this could, however, be significantly faster, by simply adding up canopy areas and deducting the gap area */
    
    int crown_base = int(height - CD),
    crown_top = int(height),
    row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
      
    float fraction_filled_target = 1.0;
    int shell_fromtop = 0;
    int voxel_volume = 1;
    int layers_filled = min(crown_top - crown_base,3);

    LoopLayerUpdateCrownStatistic_template(row_crowncenter, col_crowncenter, height, CR, CD, fraction_filled_target, shell_fromtop, layers_filled, GetRadiusSlope, voxel_volume, CVtree_layers, ModifyNoModification_int, AddupVolumeLayers);
}

int PlaceTree(float height_test, float CR_test, float CD_test, vector<int> &sites_free){

    /* get a random default site */
    int nbsites_free = int(sites_free.size());
    int index_site_test = int(gsl_rng_uniform_int(gslrand, nbsites_free));
    int site_test = sites_free[index_site_test];
    
    /* check other sites, if other conditions are enforced */
    int nbsites_testedmax = sites;
    int nbsites_tested = 0;
    int success = 0;
    
    float avg_negdiff_test;
    if(flag_Prefitting == 1) avg_negdiff_test = GetAvgNegdiff(site_test, height_test, CR_test, CD_test);
    else avg_negdiff_test = 0.0;
        
    /* we now run through sites that potentially could fit the tree and pick the one that creates the best fit */
    while(nbsites_tested < nbsites_testedmax && success == 0){
        nbsites_tested++;
        
        if(flag_PreventCrownpiercing == 1){
            int ingrowth = TreeIngrowth(site_test, height_test, CR_test);
            if(ingrowth == 0){
                success = 1;
            }
        } else {
            success = 1;
        }
        
        if(flag_Prefitting == 1 && success == 1){
            /* we check whether tree fits into canopy (on average) */
            /* get the sum of negative deviations from the empirical canopy (with a boundary layer, determined by the factor_fit) */
            if(avg_negdiff_test < 0.0){
                success = 0;
            }
        }
            
        if(success == 0){
            int index_site_testnew = int(gsl_rng_uniform_int(gslrand, nbsites_free));
            int site_testnew = sites_free[index_site_testnew];
                
            if(flag_Prefitting == 1){
                float avg_negdiff_testnew = GetAvgNegdiff(site_testnew, height_test, CR_test, CD_test);
                /* if the difference is negative, but closer to zero than before, then we update the site */
                if(avg_negdiff_testnew > avg_negdiff_test){
                    site_test = site_testnew;
                    avg_negdiff_test = avg_negdiff_testnew;
                }
            }
        }
    }
    return(site_test);
}

    
void InitialiseClimate(int flag_climate, int &error_climate){
    /*** Initialization of environmental variables ***/
    /*************************************************/
    
    if(flag_climate == 1){
        if(flag_climate) sprintf(inputfile_climate,"%s",bufi_climate);
        fstream InClim(inputfile_climate, ios::in);
        
        if(InClim){
            error_climate = 0;
            cout << "Reading in file: " << inputfile_climate << endl;
            
            for(int line=0;line<4;line++) InClim.getline(buffer,128,'\n');
            
            /* normalized daily variation in PPFD, VPD, T */
            for (int i=0; i<=23; i++) InClim >> daily_light[i];
            InClim.getline(buffer,128,'\n');
            for (int i=0; i<=23; i++) InClim >> daily_vpd[i];
            InClim.getline(buffer,128,'\n');
            for (int i=0; i<=23; i++) InClim >> daily_T[i];
            
            InClim.getline(buffer,128,'\n');
            InClim.getline(buffer,128,'\n');
            InClim.getline(buffer,128,'\n');
            
            /* monthly averages */

            for (int i=0; i<iterperyear; i++) InClim >> Temperature[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> DailyMaxTemperature[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> NightTemperature[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> Rainfall[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> WindSpeed[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> MaxIrradiance[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> MeanIrradiance[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> SaturatedVapourPressure[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> VapourPressure[i];
            InClim.getline(buffer,128,'\n');
     
            for (int i=0; i<iterperyear; i++) InClim >> VapourPressureDeficit[i];
            InClim.getline(buffer,128,'\n');
            
            for (int i=0; i<iterperyear; i++) InClim >> DailyVapourPressureDeficit[i];
            InClim.getline(buffer,128,'\n');
        
            for (int i=0; i<iterperyear; i++) InClim >> DailyMaxVapourPressureDeficit[i];
            InClim.getline(buffer,128,'\n');
            
            /* choose average conditions */
            Wmax_avg = Tmax_avg = VPDmax_avg = Tnight_avg = temp_avg = 0.0;
            for (int i=0; i<iterperyear; i++) {
                Wmax_avg += MaxIrradiance[i]*2.25;
                Tmax_avg += DailyMaxTemperature[i];
                VPDmax_avg += DailyMaxVapourPressureDeficit[i];
                Tnight_avg += NightTemperature[i];
                temp_avg += Temperature[i];
            }
            Wmax_avg *= 0.08333333;
            Tmax_avg *= 0.08333333;
            VPDmax_avg *= 0.08333333;
            Tnight_avg *= 0.08333333;
            temp_avg *= 0.08333333;
            
            cout << "Climate variables from climate input sheet: Wmax_avg: " << Wmax_avg << " Tmax_avg: " << Tmax_avg << " VPDmax_avg: " << VPDmax_avg << " Tnight_avg: " << Tnight_avg << " temp_avg: " << temp_avg << endl;
            
        } else {
            error_climate = 2;
            cout << "ERROR. Climate file could not be read. Canopy Constructor will exit." << endl;
        }

        InClim.close();
    } else {
        error_climate = 1;
        nbsteps_carbonbalance = 0;
        cout << "WARNING. No climate file was provided. Canopy Constructor defaults to mode without carbon balance optimization." << endl;
    }
}

void InitialiseCP(int flag_cp, int &error_cp){

    if(flag_cp == 1){
        if(flag_cp) sprintf(inputfile_cp,"%s",bufi_cp);
        fstream Incp(inputfile_cp, ios::in);
        int height_canopy_max = height_max + 1;

        if(Incp){
   
            error_cp = 0;
            cout << "Reading from file: " << inputfile_cp << endl;
        
            string line;
            getline(Incp,line);
            istringstream firstlinestream(line);

            //vector<string> variable_names{"paramID","height_canopy","height_within_canopy","CP_height"};
            //int nb_variables = int(variable_names.size());
            
            string variable_names[4] = {"paramID","height_canopy","height_within_canopy","CP_height"};
            int nb_variables = 4;
            
            vector<int> variable_index;
        
            string variable_name;
                    
            while(firstlinestream >> variable_name){
                int index = -1;
                for(int i = 0; i < nb_variables; i++){
                    if(variable_name == variable_names[i]) index = i;
                }
                if(index == -1){
                    cout << "Ignoring unknown variable: " << variable_name << ". Variables should be one of:";
                    for(int i = 0; i < nb_variables; i++){
                        cout << "\t" << variable_names[i];
                    }
                    cout << endl;
                }
                variable_index.push_back(index);
            }


            int linenb = 0;
            int linenb_error = 0;
            
            /* first iteration just checks whether the file is good */
            while(getline(Incp, line) && linenb_error == 0){
                istringstream linestream(line);
                linenb++;
                
                string variable_value;
                
                int paramID = -1, height_canopy = -1, height_within_canopy = -1;
                float CP_height = -1.0;
                
                int v = 0;
                while(linestream >> variable_value){
                    bool quiet;
                    if(linenb < 5) quiet = 0;
                    else quiet = 1;
                    
                    int index = variable_index[v];
                    
                    if(index >= 0){
                        string variable_name = variable_names[index];
     
                        if(variable_name == "paramID") SetParameter(variable_name, variable_value, paramID, 0, INT_MAX, -1, quiet);
                        else if(variable_name == "height_canopy") SetParameter(variable_name, variable_value, height_canopy, 0, height_max, -1, quiet);
                        else if(variable_name == "height_within_canopy") SetParameter(variable_name, variable_value, height_within_canopy, 0, height_max, -1, quiet);
                        else if(variable_name == "CP_height") SetParameter(variable_name, variable_value, CP_height, 0.0f, 5.0f, -1.0f, quiet);
                    }
                    v++;
                }
                
                int parameterset = -1;
                for(int i = 0; i < nb_parametersets; i++){
                    int paramIDtest = int(lround(parameters_detailed[0 + i * nb_parametersdetailed]));
                    if(paramID == paramIDtest) parameterset = i;
                }
                
                if(parameterset >= 0 && height_canopy >= 0 && height_within_canopy >= 0 && CP_height >= 0.0){
                    CP_twodimemp_params[parameterset + nb_parametersets * (height_canopy + height_within_canopy * height_canopy_max)] = CP_height;
                } else {
                    linenb_error = linenb;
                }
            }
            
            if(linenb == 0){
                error_cp = 1;
                cout << "WARNING. Crown packing input was empty. Options reduced to Step 0 or Step 1." << endl;
            } else if(linenb_error > 0){
                error_cp = 2;
                cout << "ERROR. Line: " << linenb_error << " of crown packing file could not be read in. Canopy Constructor will exit. Please make sure that parameterset IDs correspond to IDs of the detailed parameter file, and that all heights and packing densities are greater 0. Canopy Constructor will exit." << endl;
            }
        } else {
            error_cp = 2;
              cout << "ERROR. Crown packing input could not be read. Canopy Constructor will exit." << endl;
        }
        Incp.close();
    } else {
        error_cp = 1;
        cout << "WARNING. No crown packing input was provided. Options reduced to Step 0 or Step 1." << endl;
    }
}
        

void SupplementTraits(float dbh, int species_label, float &height, float &CR, float &CD, float &Nmass, float &Pmass, float &LMA, float &wsg){
    
    /* if no values are provided for key traits (values <= 0), this function complements traits from the species level traits, with random variation */
    float log_dbh = log10(dbh*100.0);
    int bin_tree = int(log_dbh/dbhlog_binsize);
    float sigma_heightbin_tree = sigma_heightbin[bin_tree];
    float sigma_CRbin_tree = sigma_CRbin[bin_tree];
    
    if(!(height > 0.0)){
        float deviation_height = float(gsl_ran_gaussian(gslrand, sigma_heightbin_tree));
        height= b_height * dbh / (dbh + a_height) * exp(deviation_height);
        height= minf(height, height_max-1);
        height= maxf(height, 1.0);   // minimum 1m height!
    }
    
    if(!(CD > 0.0)){
        float deviation_CD = float(gsl_ran_gaussian(gslrand, sigma_CD));
        CD = (a_CD + b_CD * height) * exp(deviation_CD);
        CD = minf(CD, height * 0.5);
    }

    if(!(CR > 0.0)){
        float deviation_CR = float(gsl_ran_gaussian(gslrand, sigma_CRbin_tree));
        CR = exp(a_CR + b_CR * log(dbh) + deviation_CR);
        CR = minf(CR, 25.0);
    }

    /* the other traits and photosynthetic quantities */
    if(!(Nmass > 0.0) || !(Pmass > 0.0) || !(LMA > 0.0)){
        double deviation_N, deviation_P, deviation_LMA;
        if(covariance_status == 0){
            deviation_N = gsl_ran_gaussian(gslrand, sigma_N);
            deviation_P = gsl_ran_gaussian(gslrand, sigma_P);
            deviation_LMA = gsl_ran_gaussian(gslrand, sigma_LMA);
        } else {
            gsl_ran_multivariate_gaussian(gslrand, mu_N_P_LMA, mcov_N_P_LMA, variation_N_P_LMA);
            deviation_N = gsl_vector_get(variation_N_P_LMA, 0);
            deviation_P = gsl_vector_get(variation_N_P_LMA, 1);
            deviation_LMA = gsl_vector_get(variation_N_P_LMA, 2);
        }
        float sp_mean_LMA = mean_LMA;
        float sp_mean_Nmass = mean_N;
        float sp_mean_Pmass = mean_P;
        
        Nmass = sp_mean_Nmass * float(exp(deviation_N));
        Pmass = sp_mean_Pmass * float(exp(deviation_P));
        LMA = sp_mean_LMA * float(exp(deviation_LMA));
    }
    float deviation_wsg = float(gsl_ran_gaussian(gslrand, sigma_wsg));
    float sp_mean_wsg = mean_wsg;
    wsg = sp_mean_wsg + deviation_wsg;
}

void InitialiseInventory(int flag_inventory, int &error_inventory){
    /* Error handling: error == 0 means no error, error == 1 means nonfatal error, i.e. no file is provided and default options are used, error == 2 is fatal error (file was inconsistent with other parameter/data sets provided or nonexistent) */
    
    if(flag_inventory == 1){
        sprintf(inputfile_inventory,"%s",bufi_inventory);
        fstream InInventory(inputfile_inventory, ios::in);
        
        if(InInventory){
            error_inventory = 0;
            cout << "Reading in file " << inputfile_inventory << "\n";
            
            string line;
            getline(InInventory,line);
            istringstream firstlinestream(line);
            
            /* define potential variables, the first three are mandatory */
            //vector<string> variable_names{"x","y","species_label","dbh","height","CR","CD","LMA","Nmass","Pmass","wsg"};
            //int nb_variables = int(variable_names.size());
            
            string variable_names[11] = {"x","y","species_label","dbh","height","CR","CD","LMA","Nmass","Pmass","wsg"};
            int nb_variables = 11;
            int nb_variables_fromfile = 0;
            
            vector<int> variable_index;
            string variable_name;
            
            while(firstlinestream >> variable_name){
                int index = -1;
                for(int i = 0; i < nb_variables; i++){
                     if(variable_name == variable_names[i]){
                         index = i;
                         nb_variables_fromfile++;
                     }
                }
                if(index == -1){
                    cout << "Ignoring unknown variable: " << variable_name << ". Variables must be one of:";
                    for(int i = 0; i < nb_variables; i++){
                        cout << "\t" << variable_names[i];
                    }
                    cout << endl;
                }
                variable_index.push_back(index);
            }
            
            trees_sitesnspecies_empirical.reserve(2*sites);
            trees_traits_empirical.reserve(8*sites);
            
            dbh_cutoff_fromground = dbh_limittop;

            int linenb = 0;
            int linenb_error = 0;
            int  sites_outside = 0, sites_duplicate = 0, trees_accepted = 0;
            
            // restricting to data sets with a maximum number of values corresponding to no. of sites
            while(getline(InInventory, line) && linenb_error == 0){
                istringstream linestream(line);
                linenb++;
                 
                /* we update the values from the file */
                int col = -1, row = -1, species_label = -1;
                float dbh = -1.0, height = -1.0, CR = -1.0, CD = -1.0, LMA = -1.0, Nmass = -1.0, Pmass = -1.0, wsg = -1.0;
                
                int v = 0;
                string variable_value;
                while(linestream >> variable_value){
                    bool quiet;
                    if(linenb < 5) quiet = 0;
                    else quiet = 1;
                    int index = variable_index[v];
                    
                    if(index >= 0){
                        string variable_name = variable_names[index];
                        
                        if(variable_name == "x") SetParameter(variable_name, variable_value, col, 0, INT_MAX, -1, quiet);
                        else if(variable_name == "y") SetParameter(variable_name, variable_value, row, 0, INT_MAX, -1, quiet);
                        else if(variable_name == "species_label") SetParameter(variable_name, variable_value, species_label, 0, INT_MAX, 1, quiet);
                        else if(variable_name == "dbh") SetParameter(variable_name, variable_value, dbh, 0.0f, dbh_limittop, -1.0f, quiet);
                        else if(variable_name == "height") SetParameter(variable_name, variable_value, height, 0.0f, float(height_max), 0.0f, quiet);
                        else if(variable_name == "LMA") SetParameter(variable_name, variable_value, LMA, 0.0f, 10000.0f, 0.0f, quiet);
                        else if(variable_name == "Nmass") SetParameter(variable_name, variable_value, Nmass, 0.0f, 1.0f, 0.0f, quiet);
                        else if(variable_name == "Pmass") SetParameter(variable_name, variable_value, Pmass, 0.0f, 1.0f, 0.0f, quiet);
                        else if(variable_name == "wsg") SetParameter(variable_name, variable_value, wsg, 0.0f, 1.5f, 0.0f, quiet);
                    }
                    v++;
                }
                
                if(dbh > 0.0){
                    /* now standardize the location */
                    col -= float(mincol_absolute);
                    row -= float(minrow_absolute);
                   
                    if((col >= 0) && (col < cols) && (row >= 0) && (row < rows)){
                       
                        int site = col + row * cols;
                        
                        int site_occupied = 0;
                        for(int i = 0; i < trees_sitesnspecies_empirical.size()/2; i++){
                            if(trees_sitesnspecies_empirical[0 + i * 2] == site) site_occupied = 1;
                        }
                        
                        if(!site_occupied){
                            trees_accepted++;
                            dbh_cutoff_fromground = minf(dbh_cutoff_fromground,dbh);
                            
                            trees_sitesnspecies_empirical.push_back(site);
                            trees_sitesnspecies_empirical.push_back(species_label);
                            trees_traits_empirical.push_back(dbh);
                            trees_traits_empirical.push_back(height);
                            trees_traits_empirical.push_back(CR);
                            trees_traits_empirical.push_back(CD);
                            trees_traits_empirical.push_back(Nmass);
                            trees_traits_empirical.push_back(Pmass);
                            trees_traits_empirical.push_back(LMA);
                            trees_traits_empirical.push_back(wsg);
                        } else {
                            sites_duplicate++;
                        }
                    } else {
                        sites_outside++;
                    }
                } else {
                    linenb_error = linenb;
                }
            }
            
            if(linenb == 0){
                error_inventory = 1;
                dbh_cutoff_fromground = dbh_limitbottom;
                cout << "WARNING. Inventory input was empty. Canopy Constructor defaults to drawing trees from distribution." << endl;
            } else if (linenb_error > 0){
                error_inventory = 2;
                cout << "ERROR. Line: " << linenb_error << " of inventory file could not be read in. Canopy Constructor will exit." << endl;
            } else {
                cout << "Successfully read inventory file from input. Read " << linenb << " lines, out of which " << trees_accepted << " were accepted. " << sites_duplicate << " were removed due to duplicate sites, and " << sites_outside << " due to location outside of the coordinate frame." << endl;
                if(float(sites_outside)/float(trees_accepted) >= 0.1) cout << "WARNING. This corresponds to " << 100.0 * float(sites_outside)/float(trees_accepted) << " % of trees outside the coordinate frame." << endl;
            }
        } else {
            error_inventory = 2;
            cout << "ERROR. Inventory input file could not be read. Canopy Constructor will exit." << endl;
        }
        InInventory.close();
    } else {
        error_inventory = 1;
        dbh_cutoff_fromground = dbh_limitbottom;
        cout << "WARNING. No inventory file was provided. Canopy Constructor defaults to drawing trees from distribution." << endl;
    }
}

void CanopyDistanceFull(float &dissimilarity, int &sum_abserror){
    vector<int> hist_simulated(height_max+1, 0);
    vector<int> hist_empirical(height_max+1, 0);
    
    int diff_cumulated = 0;
    int height_count = 0;
    
    for(int site = 0; site < sites; site++){
        int height_emp = chm_empirical[site];
        int height_sim = chm_field[site];
        /* we exclude all voxels where both the empirical height is lower than the cutoff and the simulated height is also lower than the empirical one */
        /* this means that we do not assume that the simulation can fill up low canopy sections, and do not penalise it for such deviations */
        /* is is, however, penalised for deviations in the other direction, e.g. where it overestimates the canopy height */
        //if(height_emp > height_cutoff_fromground){
        
        hist_empirical[height_emp]++;
        hist_simulated[height_sim]++;
        int diff = height_sim - height_emp;
        diff_cumulated += diff;
        if(diff < 0) diff = -diff;
        sum_abserror += diff;
        height_count++;
    }
    
    /* normalizing the histograms: height count depends a lot on the chosen tree crown size, and a bias towards smaller sizes will increase the chance of overlap */
    /* one way to deal with this is by normalizing with the whole measured field and/or using nonoverlapping area, instead of overlapping area */
    if(height_count > 0){
        int hbin_size = 1;
        
        for(int hbin = 0; hbin <= height_max + 1; hbin += hbin_size){
            float percentagebin_emp = 0.0, percentagebin_sim = 0.0;
            
            for(int h = hbin; h < min(hbin + hbin_size,height_max+1); h++){
                /* we normalize the histogram */
                percentagebin_emp += float(hist_empirical[h])/float(height_count);
                percentagebin_sim += float(hist_simulated[h])/float(height_count);
            }
            float percentagebin_error = percentagebin_sim - percentagebin_emp;
            if(percentagebin_error < 0.0) percentagebin_error *= (-1.0);
            dissimilarity += percentagebin_error*0.5;
        }
    } else {
        sum_abserror = height_max;
        dissimilarity = 1.0;
    }
};
                   
                   
void CanopyDevMetrics(int site_crown, int extent, int remove){
    //cout << "CalcDevMetrics" << endl;
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    int min_row = max(0, row_crowncenter - extent);
    int min_col = max(0, col_crowncenter - extent);
    int max_row = min(rows - 1, row_crowncenter + extent);
    int max_col = min(cols - 1, col_crowncenter + extent);
    
    int sum_abserror_local = 0;

    for(int row = min_row; row <= max_row; row++){
        for(int col = min_col; col <= max_col; col++){
            
            //cout << site_crown << " col: " << col << " row: " << row << endl;
            int site = col + row * cols;
            int height_emp = chm_empirical[site];
            int height_sim = chm_field[site];
            
            int diff = height_sim - height_emp;
            int abserror;
            if(diff < 0) abserror = -diff;
            else abserror = diff;
            
            if(remove == 0){
                hist_CHM_sim[height_sim]++;
                sum_abserror_local += abserror;
            } else {
                hist_CHM_sim[height_sim]--;
                sum_abserror_local -= abserror;
            }
        }
    }
    
    sum_abserror += sum_abserror_local;
}
   
float CalcDissimilarity(int hbin_size){
    float dissimilarity = 0.0;
    float inv_sites = 1.0/float(sites);
        
    for(int hbin = 0; hbin <= height_max + 1; hbin += hbin_size){
        float percentagebin_emp = 0.0, percentagebin_sim = 0.0;
        
        for(int h = hbin; h < min(hbin + hbin_size,height_max+1); h++){
            /* we normalize the histogram */
            percentagebin_emp += float(hist_CHM_emp[h]) * inv_sites;
            percentagebin_sim += float(hist_CHM_sim[h]) * inv_sites;
        }
        float percentagebin_error = percentagebin_sim - percentagebin_emp;
        if(percentagebin_error < 0.0) percentagebin_error *= (-1.0);
        dissimilarity += percentagebin_error*0.5;
    }

    return(dissimilarity);
}
   
float CanopyMeanheight(){
   int sumheight = 0;
   
   for(int h = 0; h < height_max + 1; h++) sumheight += hist_CHM_sim[h] * h;
   
   float meanheight = float(sumheight)/float(sites);
   return(meanheight);
}
                
void UpdateCandistFull(){
    sum_abserror = 0;
    
    for(int h = 0; h < height_max+1; h++){
        hist_CHM_sim[h] = 0;
    }
    
    for(int site = 0; site < sites; site++){
        int height_sim = chm_field[site];
        hist_CHM_sim[height_sim]++;
    }
    
    /* we also initialise the values that we keep track of */
    for(int site = 0; site < sites; site++){
        int height_sim = chm_field[site];
        int height_emp = chm_empirical[site];
        
        /* provision, in case trees below a certain cutoff should not be fitted (e.g. when understory is not filled up) */
        if(height_emp < height_cutoff_fromground && height_sim < height_emp){
            height_sim = height_emp;
        }
        
        int diff = height_sim - height_emp;
        int abserror;
        if(diff < 0) abserror = -diff;
        else abserror = diff;
        
        sum_abserror += abserror;
    }
}

                   
int CalcCrownDisplacement(int site_tree, int height){
    float crowndisplacement_maximum = float(height) * crowndisplacement_factor;
    
    int col_displacement = int(lround(gsl_ran_flat(gslrand, -crowndisplacement_maximum,crowndisplacement_maximum)));
    int row_displacement = int(lround(gsl_ran_flat(gslrand, -crowndisplacement_maximum,crowndisplacement_maximum)));
    
    /* make sure, we do not allocate the crown center outside the plot */
    int col_tree, row_tree;
    col_tree = site_tree%cols;
    row_tree = site_tree/cols;
    
    /* we calculate the respective distances to the edges */
    int coldist_upper = cols - 1 - col_tree;
    int coldist_lower = -col_tree;
    int rowdist_upper = rows - 1 - row_tree;
    int rowdist_lower = -row_tree;
    
    /* now limit the displacement distances */
    col_displacement = min(coldist_upper, col_displacement);
    col_displacement = max(coldist_lower, col_displacement);
    row_displacement = min(rowdist_upper, row_displacement);
    row_displacement = max(rowdist_lower, row_displacement);
    
    /* now put it together */
    int crowndisplacement = col_displacement + row_displacement * cols;
    
    return(crowndisplacement);
}

float DrawPowerlaw_truncated(float slope, float dbh_lowerlimit, float dbh_upperlimit){
    /* we use a method, suggested by wolfram-alpha: http://mathworld.wolfram.com/RandomNumber.html */
    double unif = gsl_rng_uniform(gslrand);
    float dbh_test = pow(((pow(dbh_upperlimit,(slope+1)) - pow(dbh_lowerlimit,(slope+1)))*unif + pow(dbh_lowerlimit,(slope+1))),(1/(slope+1)));
    return(dbh_test);
}

float DrawExponential_truncated(float rate_exp, float dbh_lowerlimit, float dbh_upperlimit){
    /* we draw from exponential distribution, limited to the given range */
    float dbh_test = 0.0;
    int success_dbh = 0;
    while(success_dbh == 0){
        dbh_test = gsl_ran_exponential(gslrand, rate_exp);
        if(dbh_test >= dbh_lowerlimit && dbh_test < dbh_upperlimit) success_dbh = 1;
    }
    return(dbh_test);
}

float DrawWeibull_truncated(float shape, float scale, float dbh_lowerlimit, float dbh_upperlimit){
    /* brute-force drawing from distribution and comparing to */
    float dbh_test = 0.0;
    int success_dbh = 0;
    while(success_dbh == 0){
        dbh_test = float(gsl_ran_flat(gslrand, dbh_lowerlimit, dbh_upperlimit));
        float ratio_ab = shape/scale;
        float factor_wb = pow((dbh_test*100.0)/scale, shape - 1.0);
        float exponent_wb = -pow((dbh_test*100.0)/scale, shape);
        float probability =  ratio_ab * factor_wb * exp(exponent_wb);
        double unif = gsl_rng_uniform(gslrand);
        if(unif < probability) success_dbh = 1;
    }
    return(dbh_test);
}
    
int FillupTrees(vector<int> &sites_free){
    /* we assume a distribution for trees underneath the cutoff and fill up until both trees above and below the cutoff are equivalent */
    int nb_belowcutoff = 0;
    int nb_abovecutoff = 0;
    
    /* this is a vector to describe the tree distribution below the cutoff */
    /* we ensure that there are not more trees than sites */
    int nbtrees_max_belowcutoff = sites-trees_fromground;
    
    /* first we need to draw trees underneath the cutoff value */
    
    /* determine trees within range where the distribution is supposed to hold */
    int nbtrees_max_abovecutoff = 0;
    float dbh_limittop_fillup;
    
    if(flag_powerlaw == 1){
        dbh_limittop_fillup = dbh_transition;
        cout << "Assuming power law with slope: " << a_sdd << " and upper cutoff: " << dbh_limittop_fillup << endl;
    } else {
        dbh_limittop_fillup = dbh_limittop;
        cout << "Assuming Weibull distribution with shape: " << a_sdd << " and scale: " << b_sdd << endl;
    }
    
    for(int s = 0; s < sites; s++){
        if(T[s].t_age > 0.0){
            float dbh = T[s].t_dbh;
            if(dbh > dbh_limitbottom && dbh < dbh_limittop_fillup) nbtrees_max_abovecutoff++;
        }
    }
    
    /* we then need to multiply this number with the tree dimensions: dbh, height, CR, CD, var_height, var_CR, var_CD */
    /* access like in a two dimensionsal array with [dim][tree], e.g. dim + tree_id * dimensions */
   
    while(nb_abovecutoff < nbtrees_max_abovecutoff && nb_belowcutoff < nbtrees_max_belowcutoff){

        float dbh_test;
        if(flag_powerlaw == 1) dbh_test = DrawPowerlaw_truncated(a_sdd, dbh_limitbottom, dbh_limittop_fillup); // if a power law is assumed, we limit the matching to b_sdd
        else dbh_test = DrawWeibull_truncated(a_sdd, b_sdd, dbh_limitbottom, dbh_limittop_fillup);
        
        if(dbh_test >= dbh_limitbottom && dbh_test < dbh_cutoff_fromground){
            nb_belowcutoff++;
            int species_label = 1;

            /* before filling in the traits, update them */
            /* update the tratis */
            float height = 0.0, CR = 0.0, CD = 0.0, Nmass = 0.0, Pmass = 0.0,LMA = 0.0,wsg = 0.0;
            
            SupplementTraits(dbh_test, species_label, height, CR, CD, Nmass, Pmass, LMA, wsg);
            
            /* determine a free spot where to place the tree*/
            int site = PlaceTree(height, CR, CD, sites_free);
            
            T[site].InitialiseTree(site, dbh_test, height, CR, CD, Nmass, Pmass, LMA, wsg, species_label);
            FillVoxel(1, site + T[site].t_CrownDisplacement, height, CR, CD);
            
            /* update the free sites vector */
            int index_site = 0;
            for(int i = 1; i < sites_free.size(); i++){
                int site_test = sites_free[i];
                if(site_test == site) index_site = i;
            }
            sites_free.erase(sites_free.begin()+index_site);
            
        } else {
            //cout << "DBH_test: " << dbh_test << endl;
            nb_abovecutoff++;
        }
    }
    return(nb_belowcutoff);
}
    
void UpdateCorrelation(){
    cout << "\nResetting covariance matrix." << endl;
    /* first set array to zero */
    for(int bin = 0; bin < dbhlog_binnb; bin++){
        correlation_structure[bin][0] = 0.0;
        correlation_structure[bin][1] = 0.0;
        correlation_structure[bin][2] = 0.0;
        correlation_structure[bin][3] = 0.0;
        correlation_structure[bin][4] = 0.0;
        correlation_structure[bin][5] = 0.0;
    }
    
    int nbtrees = 0;
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            nbtrees++;
            
            float dbh = T[site].t_dbh;
            float log_dbh = log10(dbh * 100.0);
            int bin_tree = int(log_dbh/dbhlog_binsize);
            
            float dev_height = T[site].t_dev_height;
            float dev_CR = T[site].t_dev_CR;
            
            /* we save: number of trees in bin, sum of deviations, sum of squared deviations and product of deviations */
            correlation_structure[bin_tree][0] += 1.0;
            correlation_structure[bin_tree][1] += dev_height;
            correlation_structure[bin_tree][2] += dev_height * dev_height;
            correlation_structure[bin_tree][3] += dev_CR;
            correlation_structure[bin_tree][4] += dev_CR * dev_CR;
            correlation_structure[bin_tree][5] += dev_CR * dev_height;
        }
    }
    cout << "Covariance matrix created for " << nbtrees << " trees." << endl;
}

void UpdateDraw(){
 
    int nbtrees = int(trees_draw.size()/5);
    
    cout << "Updating tree dimensions (allometries + random variation) for " << nbtrees << " trees in total." << endl;
    
    for(int tree = 0; tree < nbtrees; tree++){
        float dbh = trees_draw[0 + tree * 5];
        float log_dbh = log10(dbh * 100.0);
        int bin_tree = int(log_dbh/dbhlog_binsize);

        /* compute corresponding dimensions */
        float sigma_heightbin_tree = sigma_heightbin[bin_tree];
        float sigma_CRbin_tree = sigma_CRbin[bin_tree];
        
        float dev_height_new = float(gsl_ran_gaussian(gslrand, sigma_heightbin_tree));
        float dev_CR_new = float(gsl_ran_gaussian(gslrand, sigma_CRbin_tree));
        float dev_CD_new = float(gsl_ran_gaussian(gslrand, sigma_CD));
        
        trees_draw[1 + tree * 5] = dev_height_new;
        trees_draw[2 + tree * 5] = dev_CR_new;
        trees_draw[3 + tree * 5] = dev_CD_new;
    }
    cout << "Tree dimensions are updated" << endl;
}

/* This is the main function to construct the virtual canopy and improve its fit */
/* Trees are obtained via field inventories or by drawing from distributions */
/* An initial canopy is constructed */
/* Then the initial canopy is optimized according to several criteria */
    
/* This function works differently if there is a known ground dbh distribution */
void ConstructCanopy(){
    
    /* we make sure that all trees are considered alive, before testing whether they fit into the canopy */
    /* otherwise progressive loss of trees */
    
    /* random initialisation of ground trees */
    cout << "\n####################################################################";
    cout << "\n##### Creating virtual inventory and filling up initial canopy #####";
    cout << "\n####################################################################" << endl;

    ConstructInitialCanopy();
    
    /* empty the tree vectors, since everything is now saved in the Tree class */
    trees_traits_empirical.clear();
    trees_sitesnspecies_empirical.clear();
    
    /* calculate summary statistics and save them for final output */
    CalcSumstatCanopy(sum_abserror_random, dissimilarity_random, carbonstarv_random, nbingrowth_random, nbmaxovl_random, nbmaxovl10_random, CA_exterior_random, CA_full_random);
    SaveCHM(chm_field_random);
    
    if(ConstructorStep > 0){
        cout << "\n###############################################################";
        cout << "\n##### Updating tree arrays and starting fitting procedure #####";
        cout << "\n###############################################################" << endl;
        cout << "Create randomized and non-randomized tree and site arrays and begin the comparison between empirical and simulated canopies" << endl;
        UpdateArrays();
        
        /* set precision for output streams */
        cout << fixed;
        cout << setprecision(3);
        
        /* create vector for random reshuffling of sites */
        vector<int> sites_random;
        sites_random.reserve(sites);
        for(int site = 0; site < sites; site++){
            sites_random.push_back(site);
        }
        
        cout << "Initial statistics: Dissimilarity (%): " << 100.0 *  CalcDissimilarity(1) << " MAE: " << sum_abserror/float(sites) << " Carbonstarv (%): " << 100.0 * carbonstarv << " Meanheight: " << CanopyMeanheight() << " CA_exterior (%): " << 100.0 * CA_exterior/float(CA_full) << endl;
        
        /* first for the CHM distribution */
        cout << "\n########################################################";
        cout << "\n##### Creating best-fit canopy height distribution #####";
        cout << "\n########################################################" << endl;
        cout << "Optimizing the overlap of empirical and simulated canopy height models (non-spatial)" << endl;
        int iteration_minimum = 0;
        int dissimilarity_flag = 1, mae_flag = 0, carbonstarv_flag = 0;

        FitCanopy(iteration_minimum, nbsteps_dissimilarity, dissimilarity_flag, mae_flag, carbonstarv_flag, sites_random);
       
        int nbingrowth_dist, nbmaxovl_dist, nbmaxovl10_dist, CA_exterior_dist, CA_full_dist;
        CalcSumstatCanopy(sum_abserror_dist, dissimilarity_dist, carbonstarv_dist, nbingrowth_dist, nbmaxovl_dist, nbmaxovl10_dist, CA_exterior_dist , CA_full_dist);
        
        /* then the carbon balance */
        cout << "\n############################################";
        cout << "\n##### Creating best-fit carbon balance #####";
        cout << "\n############################################" << endl;
        cout << "Creating the most viable virtual canopy, i.e. with the lowest percentage of trees in carbon starvation" << endl;
        iteration_minimum = iteration_minimum + nbsteps_dissimilarity;
        dissimilarity_flag = 0; mae_flag = 0; carbonstarv_flag = 1;
        
        FitCanopy(iteration_minimum, nbsteps_carbonbalance, dissimilarity_flag, mae_flag, carbonstarv_flag, sites_random);
        
        int nbingrowth_physiology, nbmaxovl_physiology, nbmaxovl10_physiology, CA_exterior_physiology, CA_full_physiology;
        CalcSumstatCanopy(sum_abserror_physiology, dissimilarity_physiology, carbonstarv_physiology, nbingrowth_physiology, nbmaxovl_physiology, nbmaxovl10_physiology, CA_exterior_physiology, CA_full_physiology);
            
        /* then the spatial optimization of the CHM */
        cout << "\n########################################################";
        cout << "\n##### Creating best-fit spatial crown distribution #####";
        cout << "\n########################################################" << endl;
        cout << "Optimizing the overlap of empirical and simulated canopy height models in space" << endl;
        iteration_minimum = iteration_minimum + nbsteps_carbonbalance;
        dissimilarity_flag = 0; mae_flag = 1; carbonstarv_flag = 0;
        
        FitCanopy(iteration_minimum, nbsteps_mae, dissimilarity_flag, mae_flag, carbonstarv_flag, sites_random);
        int nbingrowth_spatial, nbmaxovl_spatial, nbmaxovl10_spatial, CA_exterior_spatial, CA_full_spatial;
        CalcSumstatCanopy(sum_abserror_spatial, dissimilarity_spatial, carbonstarv_spatial, nbingrowth_spatial, nbmaxovl_spatial, nbmaxovl10_spatial, CA_exterior_spatial, CA_full_spatial);
        
        /* finally, the overall combined fit */
        cout << "\n#############################################";
        cout << "\n##### Creating combined best-fit canopy #####";
        cout << "\n#############################################" << endl;
        cout << "Combining all metrics" << endl;
        iteration_minimum = iteration_minimum + nbsteps_mae;

        if(nbsteps_dissimilarity > 0) dissimilarity_flag = 1;
        else mae_flag = 0;
        if(nbsteps_carbonbalance > 0) carbonstarv_flag = 1;
        else mae_flag = 0;
        if(nbsteps_mae > 0)mae_flag = 1;
        else mae_flag = 0;

        /* update the min and max values */
        /* for the maximum values, we choose the maximum of the two other fitting procedures */
        dissimilarity_max = dissimilarity_random;
        carbonstarv_max = carbonstarv_random;
        sum_abserrormax = sum_abserror_random;
        
        dissimilarity_min = dissimilarity_dist;
        carbonstarv_min = carbonstarv_physiology;
        sum_abserrormin = sum_abserror_spatial;

        FitCanopy(iteration_minimum, nbsteps_combined, dissimilarity_flag, mae_flag, carbonstarv_flag, sites_random);
        
        cout << "\n#####################################";
        cout << "\n##### Finalizing the simulation #####";
        cout << "\n#####################################" << endl;
        
        iteration_final = iteration_minimum + nbsteps_combined;
        rate_final = float(nbsuccess)/float(nbattempt);
        CalcSumstatCanopy(sum_abserror_final, dissimilarity_final, carbonstarv_final, nbingrowth_final, nbmaxovl_final, nbmaxovl10_final, CA_exterior_final, CA_full_final);
    }
}
 
void ConstructInitialCanopy(){
    
    /* define the free space on the plot */
    vector<int> sites_free;
    sites_free.reserve(sites);
    for(int site = 0; site < sites; site++){
        sites_free.push_back(site);
    }
    
    if(ConstructorStep == 2){
        /* This is a function that first calculates the theoretical volume available for trees, across height layers */
        /* And then, from this we deduce the underlying tree distribution */
        CA_exterior = 0;
        CA_full = 0;

        /* Vectors for virtual crown volume per height layer */
        int CV_totalemp = 0;
        vector<int> CV_cumulatedemp(height_max + 1, 0);
        
        /* Transform the CHM into a CV distribution */
        GetCVcumulatedemp(CV_cumulatedemp, CV_totalemp);
        
        int nbtrees = int(trees_sitesnspecies_empirical.size()/2);
        vector<int> trees_random;
        for(int tree = 0; tree < nbtrees; tree++){
            trees_random.push_back(tree);
        }
        random_shuffle(trees_random.begin(), trees_random.end());

        /* a few metrics to keep track of the fitting */
        trees_sampled = 0;
        trees_fitted = 0;
        trees_fitted10 = 0;
        
        int trees_sampled_max = 2 * sites;
        int canopy_filled = 0;
        
        vector<int> CV_cumulated(height_max + 1, 0);
        
        int trees_powerlaw = 0, trees_exponential = 0, trees_weibull = 0;
        
        while(trees_sampled < trees_sampled_max && canopy_filled == 0){
            
            int trees_fitted_beforedraw = trees_fitted;

            for(int index_tree = 0; index_tree < nbtrees; index_tree++){
                
                int flag_distribution = 0;  /* 0 if from the field inventory, 1 if power law, 2 if exponential, 3 if weibull */
                float dbh, height, CR, CD, Nmass, Pmass, LMA, wsg;
                int species_label;
                if(trees_random.size() > 0){
                    int tree = trees_random[index_tree];
                    
                    species_label = trees_sitesnspecies_empirical[1 + tree * 2];
                    dbh = trees_traits_empirical[0 + tree * 8]; height = trees_traits_empirical[1 + tree * 8]; CR = trees_traits_empirical[2 + tree * 8]; CD = trees_traits_empirical[3 + tree * 8]; Nmass = trees_traits_empirical[4 + tree * 8]; Pmass = trees_traits_empirical[5 + tree * 8]; LMA = trees_traits_empirical[6 + tree * 8]; wsg = trees_traits_empirical[7 + tree * 8];
                } else {
                    species_label = 1;
                    if(flag_powerlaw == 1){
                        double unif = gsl_rng_uniform(gslrand);
                        if(unif < fraction_drawspowerlaw){
                            flag_distribution = 1;
                            dbh = DrawPowerlaw_truncated(a_sdd, dbh_limitbottom, dbh_transition); // if a power law is assumed, we limit the matching to b_sdd
                        } else {
                            flag_distribution = 2;
                            float rate_exp = -1.0/b_sdd;
                            dbh = DrawExponential_truncated(rate_exp, dbh_transition, dbh_limittop);
                        }
                    } else {
                        flag_distribution = 3;
                        dbh = DrawWeibull_truncated(a_sdd, b_sdd, dbh_limitbottom, dbh_limittop);
                    }
                    height = 0.0; CR = 0.0; CD = 0.0; Nmass = 0.0; Pmass = 0.0; LMA = 0.0; wsg = 0.0;
                }
                
                SupplementTraits(dbh, species_label, height, CR, CD, Nmass, Pmass, LMA, wsg);

                /* place tree */
                int site = PlaceTree(height, CR, CD, sites_free);
                
                int success;
                
                FitTree(site, height, CR, CD, success, CV_cumulatedemp, CV_cumulated);
                
                if(success == 1){
                    /* update statistics */
                    trees_fitted++;
                    if(dbh >= 0.1) trees_fitted10++;
                    if(flag_distribution == 1) trees_powerlaw++;
                    if(flag_distribution == 2) trees_exponential++;
                    if(flag_distribution == 3) trees_weibull++;
                    
                    T[site].InitialiseTree(site, dbh, height, CR, CD, Nmass, Pmass, LMA, wsg, species_label);
                    FillVoxel(1, site + T[site].t_CrownDisplacement, height, CR, CD);
                    
                    /* remove site from vector of free sites */
                    int index_site = 0;
                    for(int i = 1; i < sites_free.size(); i++){
                        int site_tested = sites_free[i];
                        if(site_tested == site) index_site = i;
                    }
                    sites_free.erase(sites_free.begin()+index_site);
                    
                }
            }
            /* if the whole draw was rejected (for trees above 10cm), then stop fitting (potentially losing ~ 1-5 trees) */
            if(trees_fitted == trees_fitted_beforedraw) canopy_filled = 1;
        }
        
        cout << "Trees fitted into canopy: " << trees_fitted << endl;
        if(trees_powerlaw > 0) cout << "Drawn from power law: " << trees_powerlaw << endl;
        if(trees_exponential > 0) cout << "Drawn from exponential tail of power law: " << trees_exponential << endl;
        if(trees_weibull > 0) cout << "Drawn from 2-parameter Weibull distribution: " << trees_weibull << endl;
        int CV_total = 0;
        for(int h = 0; h < (height_max + 1); h++){
            CV_total += CV_cumulated[h];
        }
        cout << "Total crown volume (simulated): " << CV_totalemp << endl;
    } else {
        /* in the case of simple ground fitting, we only need to copy the vectors that have been read in empirically and update missing traits */
        
        int trees_powerlaw = 0, trees_exponential = 0, trees_weibull = 0;
        
        if(trees_sitesnspecies_empirical.size() > 0){
            cout << "Supplementing inventory with traits." << endl;

            for(int tree = 0; tree < trees_sitesnspecies_empirical.size()/2; tree++){
                int site = trees_sitesnspecies_empirical[0 + tree * 2];
                
                /* update traits from species information, if necessary */
                int species_label = trees_sitesnspecies_empirical[1 + tree * 2];
                
                /* use local variables for clarity's sake */
                float dbh = trees_traits_empirical[0 + tree * 8], height = trees_traits_empirical[1 + tree * 8], CR = trees_traits_empirical[2 + tree * 8], CD = trees_traits_empirical[3 + tree * 8], Nmass = trees_traits_empirical[4 + tree * 8], Pmass = trees_traits_empirical[5 + tree * 8], LMA = trees_traits_empirical[6 + tree * 8], wsg = trees_traits_empirical[7 + tree * 8];
                SupplementTraits(dbh, species_label, height, CR, CD, Nmass, Pmass, LMA, wsg);
                
                /* update array from local variables */
                T[site].InitialiseTree(site, dbh, height, CR, CD, Nmass, Pmass, LMA, wsg, species_label);
                FillVoxel(1, site + T[site].t_CrownDisplacement, height, CR, CD);
                
                /* update the free sites vector */
                int index_site = 0;
                for(int i = 1; i < sites_free.size(); i++){
                    int site_tested = sites_free[i];
                    if(site_tested == site) index_site = i;
                }
                sites_free.erase(sites_free.begin()+index_site);
            }
        } else {
            cout << "Drawing trees from distribution and supplementing them with traits." << endl;
            for(int tree = 0; tree < nb_sdd; tree++){
                float dbh;
                if(flag_powerlaw == 1){
                    double unif = gsl_rng_uniform(gslrand);
                    if(unif < fraction_drawspowerlaw){
                        trees_powerlaw++;
                        dbh = DrawPowerlaw_truncated(a_sdd, dbh_limitbottom, dbh_transition); // if a power law is assumed, we limit the matching to b_sdd
                    } else {
                        trees_exponential++;
                        float rate_exp = -1.0/b_sdd;
                        dbh = DrawExponential_truncated(rate_exp, dbh_transition, dbh_limittop);
                    }
                } else {
                    trees_weibull++;
                    dbh = DrawWeibull_truncated(a_sdd, b_sdd, dbh_limitbottom, dbh_limittop);
                }
                float height = 0.0, CR = 0.0, CD = 0.0, Nmass = 0.0, Pmass = 0.0, LMA = 0.0, wsg = 0.0;
                int species_label = 1;
                
                SupplementTraits(dbh, species_label, height, CR, CD, Nmass, Pmass, LMA, wsg);
                
                int site = PlaceTree(height, CR, CD, sites_free);

                T[site].InitialiseTree(site, dbh, height, CR, CD, Nmass, Pmass, LMA, wsg, species_label);
                FillVoxel(1, site + T[site].t_CrownDisplacement, height, CR, CD);
                
                /* update the free sites vector */
                int index_site = 0;
                for(int i = 1; i < sites_free.size(); i++){
                    int site_tested = sites_free[i];
                    if(site_tested == site) index_site = i;
                }
                sites_free.erase(sites_free.begin()+index_site);
            }
        }
        cout << "Successfully created " << sites - sites_free.size()  << " virtual trees." << endl;
        if(trees_powerlaw > 0) cout << "Drawn from power law: " << trees_powerlaw << endl;
        if(trees_exponential > 0) cout << "Drawn from exponential tail of power law: " << trees_exponential << endl;
        if(trees_weibull > 0) cout << "Drawn from 2-parameter Weibull distribution: " << trees_weibull << endl;
    }

    /* now we fill up underneath threshold */
    cout << "Supplementing trees below trunk diameter threshold of: " << round(100.0 * dbh_cutoff_fromground) << "cm " << endl;
    int nb_belowcutoff = 0;
    if(dbh_cutoff_fromground > 0.01){
        nb_belowcutoff = FillupTrees(sites_free);
    }
    cout << "Number of supplemented trees: " << nb_belowcutoff << endl;
}
 
void GetCVcumulatedemp(vector<int> &CV_cumulatedemp, int &CV_totalemp){
    cout << "Inferring crown volume from crown packing densities and CHM" << endl;
    vector<int> CV_cumulated_twodimemp((height_max + 1) * (height_max + 1), 0);
    
    /* The chm frequency distribution */
    vector<int> chm_frequency(height_max + 1, 0);

    for(int s=0;s<sites;s++){
       int height_canopy = min(chm_empirical[s],height_max);
       chm_frequency[height_canopy]++;
    }
    /* now loop over the crown volume array and multiply the CV with the number of times the corresponding canopy height occurs */
    for(int height_canopy = 0; height_canopy < (height_max) + 1; height_canopy++){
       int frequency = chm_frequency[height_canopy];
       for(int h = 0; h <= height_canopy; h++){
           CV_cumulated_twodimemp[height_canopy + h * (height_max + 1)] = int(lround(CP_twodimemp[height_canopy + h * (height_max + 1)] * float(frequency)));
       }
    }

    cout << "Calculating crown volume per height layer" << endl;
    /* now sum across the height layers */
    for(int h = 0; h < (height_max + 1); h++){
        for(int height_canopy = 0; height_canopy < (height_max + 1); height_canopy++){
           CV_cumulatedemp[h] += CV_cumulated_twodimemp[height_canopy + h * (height_max + 1)];
           CV_totalemp += CV_cumulated_twodimemp[height_canopy + h * (height_max + 1)];
        }
    }

    cout << "Total crown volume (empirical): " << CV_totalemp << endl;
}

void SaveCHM(vector<int> &chm_field_safe){
    chm_field_safe.clear();
    chm_field_safe.reserve(sites);
    
    for(int site = 0; site < sites; site++){
        chm_field_safe[site] = chm_field[site];
    }
}

void CalcSumstatCanopy(int &sum_abserror_safe, float &dissimilarity_safe, float &carbonstarv_safe, int &nbingrowth_safe, int &nbmaxovl_safe, int &nbmaxovl10_safe, int &CA_exterior_safe, int &CA_full_safe){
    /* update the covariation within bins */
    UpdateCorrelation();

    /* update rudimentary tree statistics */
    nbtrees_total = 0;
    nbtrees_total_threshold = 0;
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            nbtrees_total++;
            if(T[site].t_dbh >= dbh_limitcarbon) nbtrees_total_threshold++;
        }
    }
    
    /* update the canopy height statistics and save them */
    UpdateCHM();
    UpdateCandistFull();
    sum_abserror_safe = sum_abserror;
    dissimilarity_safe = CalcDissimilarity(1);

    /* update the carbon starvation statistics and save them */
    //UpdateLeaves();
    carbonstarv_safe = carbonstarv;
    
    /* update the ingrowth of trees */
    nbingrowth_safe = 0;
    nbmaxovl_safe = 0;
    nbmaxovl10_safe = 0;
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            int treeingrowth = TreeIngrowthFromabove(site + T[site].t_CrownDisplacement, T[site].t_Tree_Height, T[site].t_Crown_Radius);
            nbingrowth_safe += treeingrowth;
            
            float overlap_tree = GetVoxeloverlap(site + T[site].t_CrownDisplacement, T[site].t_Tree_Height, T[site].t_Crown_Radius, T[site].t_Crown_Depth, 0);
            if(overlap_tree > maxoverlap_percentage){
                nbmaxovl_safe++;
                if(T[site].t_dbh >= 0.1) nbmaxovl10_safe++;
            }
        }
    }
    
    /* update the exterior area (i.e. tree crowns reaching outside the plot) */
    CA_exterior = 0;
    CA_full = 0;
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            int CA_exterior_tree = CalcCrownAreaExterior(site + T[site].t_CrownDisplacement, T[site].t_Crown_Radius);
            CA_exterior += CA_exterior_tree;
            int CA_tree = GetCrownIntarea(T[site].t_Crown_Radius);
            CA_full += CA_tree;
            //if(T[site].t_dbh < 0.01) cout << "WARNING" << endl;
        }
    }
    
    CA_exterior_safe = CA_exterior;
    CA_full_safe = CA_full;
    
    cout << "Nb of trees piercing other trees: " << nbingrowth_safe << "\nNb of trees exceeding max overlap: " << nbmaxovl_safe << "\nNb of trees exceeding max overlap (>=10cm): " << nbmaxovl10_safe << "\nCA_full " << CA_full_safe << "\nCA_exterior: " << CA_exterior_safe << endl;
}
    
void UpdateArrays(){
    trees_dbhsorted.clear();
    trees_dbhsorted.reserve(nbtrees_total);

    trees_dbhlogbins.clear();
    trees_dbhlogbins.reserve(dbhlog_binnb);

    for(int bin = 0; bin < dbhlog_binnb; bin++){
        trees_dbhlogbins.push_back(0);
    }

    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0 && T[site].t_dbh >= dbh_cutoff_fromground){
            trees_dbhsorted.push_back(site);
            float log_dbh = log10(T[site].t_dbh * 100.0);
            int bin_tree = int(log_dbh/dbhlog_binsize);
            trees_dbhlogbins[bin_tree]++;
        }
    }
    
    /* and sort the other according to dbh */
    for(int i = 0; i < trees_dbhsorted.size(); i++){
        for(int rank = 0; rank < trees_dbhsorted.size()-i-1; rank++){
            int site_current = trees_dbhsorted[rank];
            int site_next = trees_dbhsorted[rank+1];
            float dbh_current = T[site_current].t_dbh;
            float dbh_next = T[site_next].t_dbh;

            if(dbh_current > dbh_next){
                trees_dbhsorted[rank] = site_next;
                trees_dbhsorted[rank+1] = site_current;
            }
        }
    }
    
    /* for physiology submodule: get maximum crown radii for each group*/
    maxCR_perheight.reserve(height_max+1);
    for(int h = 0; h < height_max+1; h++){
        maxCR_perheight.push_back(0.0);
    }
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            int height = int(T[site].t_Tree_Height);
            float CR = T[site].t_Crown_Radius;

            if(CR > maxCR_perheight[height]) maxCR_perheight[height] = CR;
        }
    }
    
    for(int h = 0; h < height_max+1; h++){
        for(int h_lower = 0; h_lower < h; h_lower++){
            float maxCR_lower = maxCR_perheight[h_lower];
            float maxCR = maxCR_perheight[h];
            
            if(maxCR_lower > maxCR) maxCR_perheight[h] = maxCR_lower;
        }
    }
    
    heightminimum_threshold = float(height_max);
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0 && T[site].t_dbh >= dbh_limitcarbon){
            float height = T[site].t_Tree_Height;
            if(height < heightminimum_threshold) heightminimum_threshold = height;
        }
    }
}

void FitCanopy(int iteration_minimum, int nbsteps, int dissimilarity_flag, int mae_flag, int carbonstarv_flag, vector<int> &sites_random){

    int iteration_maximum = iteration_minimum + nbsteps;
    
    cout << "Nb of steps: " << nbsteps << endl;
    /* This is the core fitting loop */
    for(int iteration = iteration_minimum; iteration < iteration_maximum; iteration++){
        nbattempt = 0;
        nbsuccess = 0;

        /* shuffle the order in which tree crowns are modified */
        random_shuffle(sites_random.begin(), sites_random.end());
        
        if(ConstructorStep == 2){
            for(int i = 0; i < sites_random.size(); i++){
                int site_census = sites_random[i];
                if(T[site_census].t_age > 0){
#ifdef exchangepositions_step2
                    float dbh = T[site_census].t_dbh;
                       if(dbh >= 0.01){
                           FindBetterPositions(site_census, mae_flag, dissimilarity_flag, carbonstarv_flag);
                       }
                       else {
                           /* do not shift position of the smallest trees when optimizing the carbon starvation percentage, probably needs to be changed */
                           //FindBetterPosition(site_census, mae_flag, dissimilarity_flag, carbonstarv_flag);
                       }
#else
                    FindBetterPosition(site_census, mae_flag, dissimilarity_flag, carbonstarv_flag);
#endif
                }
            }
        } else {
            for(int i = 0; i < sites_random.size(); i++){
                int site_census = sites_random[i];
                if(T[site_census].t_age > 0){
                    float dbh = T[site_census].t_dbh;
                    if(dbh >= dbh_cutoff_fromground){
                        FindBetterCrown(site_census, mae_flag, dissimilarity_flag, carbonstarv_flag);
                    }
                    else {
                        /* do not shift position of the smallest trees when optimizing the carbon starvation percentage, probably needs to be changed */
                        FindBetterPosition(site_census, mae_flag, dissimilarity_flag, carbonstarv_flag);
                    }
                }
            }
        }
        
        float rate_success = float(nbsuccess)/float(nbattempt);
        float meanheight = CanopyMeanheight();
        if(carbonstarv_flag == 1) UpdateLeaves();
       
        cout << iteration << " Tries: " << nbattempt << " Successes: " << nbsuccess << " Rate (%): " << 100.0 * rate_success << " Dissimilarity (%): " << 100.0 *  CalcDissimilarity(1) << " MAE: " << float(sum_abserror)/float(sites) << " Carbonstarv (%): " << 100.0 * carbonstarv << " Meanheight: " << meanheight << " CA_exterior (%): " << 100.0 * CA_exterior/float(CA_full);
        cout << endl;
 
//        float dissimilarity_test = 0.0;
//        int sum_abserror_test = 0;
//        //UpdateCHM();
//        CanopyDistanceFull(dissimilarity_test, sum_abserror_test);
//
//        if(carbonstarv_flag == 1){
//            UpdateLeaves();
//            float carbonstarv_test = carbonstarv;
//            cout << " Recalc: " << dissimilarity_test << " | " << sum_abserror_test << " | " << carbonstarv_test;
//        } else {
//            cout << " Recalc: " << dissimilarity_test << " | " << sum_abserror_test;
//        }
    }
}

    
void UpdateCHM(){
    /* compute the chm_field */
    for(int site = 0; site < sites; site++){
        int height_canopy = 0;
        
        for(int h=0;h<(height_max+1);h++){
            if(Voxel3D[h][site]>0) height_canopy = h;
        }
        
        chm_field[site] = height_canopy;
    }
}
    
void UpdateCHMradius(int site_crown, float crown_radius){
    /* compute the chm_field */
    
    //cout << "Update CHM" << endl;
    
    int row_crowncenter = site_crown/cols,
    col_crowncenter = site_crown%cols;
    
    float crown_area = PI * crown_radius * crown_radius;     // floor of crown_area to bound area accumulation
    int crown_area_int = int(crown_area);   // floor of crown_area to bound area accumulation
    crown_area_int = max(crown_area_int,1);                                 // minimum area of crown (1), important for radii below ~0.5
    crown_area_int = min(crown_area_int,1963);
    
    /* alternative: if(crown_area_int > 1963) UpdateCHM(voxel); */
    
    for(int i = 0; i < crown_area_int; i++){
        int site_relative = LookUp_Crown_site[i];
    
        int row = row_crowncenter + site_relative/51 - 25;
        int col = col_crowncenter + site_relative%51 - 25;
        
        //cout << site_crown << " col: " << col << " row: " << row << endl;
        
        if(row >= 0 && row < rows && col >= 0 && col < cols){
            int site=col+cols*row;
            int height_canopy = height_max+1;
            int voxfull = 0;
            while(height_canopy > 0 && voxfull == 0){
                height_canopy--;
                if(Voxel3D[height_canopy][site] > 0) voxfull = 1;
            }

            chm_field[site] = height_canopy;
        }
    }
}

void OutputSumstatFinal(){

    /*###########################*/
    /*#### Canopy statistics ####*/
    /*###########################*/
    
    /* output_chm */
    
    output_chm << parameterset << "\t" << paramIDcurrent << "\t" << iteration_final;
    
    for(int h = 0; h < height_max + 1; h++){
        output_chm << "\t" << hist_CHM_sim[h];
    }
    
    output_chm << endl;
    
    float mean_height_sim = 0.0;
    float sd_height_sim = 0.0;

    float mean_heightsquared_sim = 0.0;

    for(int h=height_cutoff_fromground;h<height_max+1;h++){
        float density_sim = float(hist_CHM_sim[h])/float(sites);

        /* other summary statistics, no need for smoothed versions */
        mean_height_sim += density_sim * float(h);
        mean_heightsquared_sim += density_sim * float(h*h);
    }
    
    sd_height_sim = sqrt(mean_heightsquared_sim - mean_height_sim * mean_height_sim);
    
    /* column names have been added at initialisation of empirical field */
    
    output_sumstat << parameterset << "\t" << paramIDcurrent << "\t" << iteration_final << "\t" << height_cutoff_fromground << "\t" << mean_height_emp << "\t" << sd_height_emp << "\t" << mean_height_sim << "\t" << sd_height_sim << "\t" << dissimilarity_random << "\t" << dissimilarity_spatial << "\t" << dissimilarity_dist << "\t" << dissimilarity_physiology << "\t" << dissimilarity_final << "\t" << sum_abserror_random << "\t" << sum_abserror_spatial << "\t" << sum_abserror_dist << "\t" << sum_abserror_physiology << "\t" << sum_abserror_final << "\t" << carbonstarv_random << "\t" << carbonstarv_spatial << "\t" << carbonstarv_dist << "\t" << carbonstarv_physiology << "\t" << carbonstarv_final << "\t" << rate_spatial << "\t" << rate_dist << "\t" << rate_physiology << "\t" << rate_final << "\t" << nbingrowth_random << "\t" << nbingrowth_final << "\t" << nbmaxovl_random << "\t" << nbmaxovl_final << "\t" << nbmaxovl10_random << "\t" << nbmaxovl10_final << "\t" << CA_exterior_random << "\t" << CA_full_random << "\t" << CA_exterior_final << "\t" << CA_full_final << "\t";
    
    
    /*##############################*/
    /*###### BIOMASS AND GPP #######*/
    /*##############################*/
    
    /* we also quantify the AGB */
    float AGB = 0.0;
    float AGBcrown_goodman = 0.0;
    
    /* the number of trees that do not receive enough light */
    int nbtrees_shaded = 0;

    /* get the AGB of trees above cutoff */
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0 && T[site].t_dbh > dbh_cutoff_fromground){
            float mass_cylinder = T[site].t_wsg*T[site].t_Tree_Height*T[site].t_dbh*T[site].t_dbh; // the mass of the stem if assumed to be cylindric. This is converted into real biomass using tapering factors, etc.
            AGB += 0.0673*pow(mass_cylinder*10000, 0.976);
            AGBcrown_goodman += 0.001 * exp(12.45 + 0.86 * log(mass_cylinder) + 0.5 * log(T[site].t_Crown_Radius));  // from Goodman et al. 2014
            if(T[site].t_LAIabove > T[site].t_LAImax) nbtrees_shaded++;
        }
    }
    
    /* now calculate per ha and convert from kg to tons (factor 0.001) */
    float invha = 10000.0/float(sites);
    AGB *= invha * 0.001;
    AGBcrown_goodman *= invha * 0.001;
    
    /* add sumstats for GPP and NPP */
    float GPPcum = 0.0;
    float NPPcum = 0.0;
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            GPPcum += T[site].t_GPP;
            NPPcum += T[site].t_NPP;
        }
    }
    
    /* calculate per ha and convert to MgC (factor 10e-6)*/
    GPPcum *= invha * 10e-6;
    NPPcum *= invha * 10e-6;
    
    output_sumstat << AGB << "\t" << AGBcrown_goodman << "\t" << GPPcum << "\t" << NPPcum << endl;
}
    
void OutputAGBtrees_finescale(){
    /* tiling of the plot area */
    int cols_quarterha = cols/50;
    int rows_quarterha = rows/50;
    int sites_quarterha = cols_quarterha * rows_quarterha;
    
    vector<float> AGBquarterha_10(sites_quarterha,0.0);
    vector<float> AGBquarterha_full(sites_quarterha, 0.0);
    
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int site = col + row * cols;
            if(T[site].t_age > 0){
                int colquarterha = col/50;
                int rowquarterha = row/50;
                int sitequarterha = colquarterha + rowquarterha * cols_quarterha;
                
                float mass_cylinder = T[site].t_wsg*T[site].t_Tree_Height*T[site].t_dbh*T[site].t_dbh; // the mass of the stem if assumed to be cylindric. This is converted into real biomass using tapering factors, etc.
                float agb = 0.0673*pow(mass_cylinder*10000, 0.976);
                AGBquarterha_full[sitequarterha] += agb;
                
                float dbh = T[site].t_dbh;
                if(dbh >= 0.1){
                    AGBquarterha_10[sitequarterha] += agb;
                }
            }
        }
    }
    
    for(int site_quarterha = 0; site_quarterha < sites_quarterha; site_quarterha++){
        int col_quarterha = site_quarterha%cols_quarterha;
        int row_quarterha = site_quarterha/cols_quarterha;
        
        /* transform to absolute coordinates, add 50m to take the centre of the quarterha */
        int col_quarterha_abs = col_quarterha * 50 + mincol_absolute + 25;
        int row_quarterha_abs = row_quarterha * 50 + minrow_absolute + 25;
        
        float AGBquarterha_site = 0.001 * AGBquarterha_full[site_quarterha];
        float AGBquarterha_10_site = 0.001 * AGBquarterha_10[site_quarterha];
        
        output_agbquarterha << parameterset << "\t" << paramIDcurrent << "\t" << site_quarterha << "\t" << col_quarterha << "\t" << row_quarterha << "\t" << col_quarterha_abs << "\t" << row_quarterha_abs << "\t" << AGBquarterha_site << "\t" << AGBquarterha_10_site << endl;
    }
}
    
void OutputAGBtrees_distributions(){

    int sdd_binnumber_max = int(dbh_limittop * 10.0);
    vector<int> SDD_count(sdd_binnumber_max,0);
    output_sdd <<  parameterset << "\t" << paramIDcurrent;

    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int site = col + row * cols;
            if(T[site].t_age > 0){
                float dbh = T[site].t_dbh;
                int SDD_bin = int(dbh * 10.0);
                SDD_count[SDD_bin]++;
            }
        }
    }
    
    int trees_total = 0;
    int trees10_total = 0;
    for(int SDD_bin = 0; SDD_bin < sdd_binnumber_max; SDD_bin++){
        output_sdd << "\t" << SDD_count[SDD_bin];
        trees_total += SDD_count[SDD_bin];
        if(SDD_bin > 0) trees10_total += SDD_count[SDD_bin];
    }
    output_sdd << "\t" << trees_total << "\t" << trees10_total << endl;
}

void OutputSumstatEmp(){
    /* This is just a companion function to OutputSumstatFinal() */
    /* it writes the header line for the output and one line for the empirical canopy */
    
    mean_height_emp = 0.0;
    sd_height_emp = 0.0;
    
    float mean_heightsquared_emp = 0.0;
 
    for(int h=height_cutoff_fromground;h<height_max+1;h++){
        float density_emp = float(hist_CHM_emp[h])/float(sites);
        
        /* other summary statistics, no need for smoothed versions */
        mean_height_emp += density_emp * float(h);
        mean_heightsquared_emp += density_emp * float(h*h);
        
    }
    
    sd_height_emp = sqrt(mean_heightsquared_emp - mean_height_emp * mean_height_emp);
    
    /*################*/
    /*#### OUTPUT ####*/
    /*################*/
    
    /* at the top we add column names and a first row summarizing the empirical data */
    /* the header line */
    output_sumstat << "parameterset" << "\t" << "paramIDcurrent" << "\t" << "iteration_final" << "\t" << "height_cutoff_fromground" << "\t" << "mean_height_emp" << "\t" << "sd_height_emp" << "\t" << "mean_height_sim" << "\t" << "sd_height_sim" << "\t" << "dissimilarity_random" << "\t" << "dissimilarity_spatial" << "\t" << "dissimilarity_dist" << "\t" << "dissimilarity_physiology" << "\t" << "dissimilarity_final" << "\t" << "sum_abserror_random" << "\t" << "sum_abserror_spatial" << "\t" << "sum_abserror_dist" << "\t" << "sum_abserror_physiology" << "\t" << "sum_abserror_final" << "\t" << "carbonstarv_random" << "\t" << "carbonstarv_spatial" << "\t" << "carbonstarv_dist" << "\t" << "carbonstarv_physiology" << "\t" << "carbonstarv_final" << "\t" << "rate_spatial" << "\t" << "rate_dist" << "\t" << "rate_physiology" << "\t" << "rate_final" << "\t" << "nbingrowth_random" << "\t" << "nbingrowth_final" << "\t" << "nbmaxovl_random" << "\t" << "nbmaxovl_final" << "\t" << "nbmaxovl10_random" << "\t" << "nbmaxovl10_final" << "\t" << "CA_exterior_random" << "\t" << "CA_full_random" << "\t" << "CA_exterior_final" << "\t" << "CA_full_final" << "\t" << "AGB" << "\t" << "AGBcrown_goodman" << "\t" << "GPPcum" << "\t" << "NPPcum" << endl;
    
    output_chm << "nbparam" << "\t" << "ID" << "\t" << "iter";
    
    for(int h = 0; h < height_max + 1; h++){
        output_chm << "\t" << h;
    }
    
    output_chm << endl;
    
    output_chm << -1 << "\t" << -1 << "\t" << -1;
    
    for(int h = 0; h < height_max + 1; h++){
         output_chm << "\t" << hist_CHM_emp[h];
    }
    
    output_chm << endl;
    
    /* the header for the other output files */
    output_agbquarterha << "nbparam" << "\t" << "ID" << "\t" << "site_unit" << "\t" << "col_unit" << "\t" << "row_unit" << "\t" << "col_unit_abs" << "\t" << "row_unit_abs" << "\t" << "AGB_unit" << "\t" << "AGB_unit_10" << endl;
   
    output_sdd << "nbparam" << "\t" << "ID";
    int sdd_binnumber_max = int(dbh_limittop * 10.0);
    for(int i = 0; i < sdd_binnumber_max; i++){
        output_sdd << "\t" << i;
    }
    output_sdd << "\t" << "trees_total" << "\t" << "trees10_total" << endl;
}

void OutputTroll(){
    cout << "\n##################################################################";
    cout << "\n##### Creating input files for TROLL forest growth simulator #####";
    cout << "\n##################################################################" << endl;
    cout << "Writing input files for transfer of Canopy Constructor canopies into TROLL model" << endl;
    
    /* this recreates an input sheet for troll straight from the Canopy Constructor */
    /* this ensure highest possible compatibility */
    //output_troll.setf(ios::fixed,ios::floatfield);
    output_troll.precision(6);
    
    output_troll << "############################################" << endl;
    output_troll << "##### Parameter file for TROLL program #####" << endl;
    output_troll << "############################################" << endl;
    output_troll << "### GENERAL PARAMETERS ###" << endl;
    output_troll << cols << "\t/* cols # nb of columns */" << endl;
    output_troll << rows << "\t/* rows # nb of rows  */" << endl;
    output_troll << "6000\t/* nbiter # total nb of timesteps */" << endl;
    output_troll << "12\t/* number of iteration per year */" << endl;
    output_troll << "1\t/* NV # vertical nb of cells (nb per m) */" << endl;
    output_troll << "1\t/* NH # horizontal nb of cells (nb per m) */" << endl;
    output_troll << "4\t/* nbout # Number of outputs */" << endl;
    output_troll << 1 << "\t/* numesp # Number of species */" << endl;
    output_troll << "0.05\t/* p # light incidence param (diff through turbid medium) */" << endl;
    output_troll << endl;
    output_troll << "###	Characters shared by species ###" << endl;
    output_troll << klight << "\t/* klight # light attenuation in the canopy Beer-Lambert */" << endl;
    output_troll << phi << "\t/* phi # quantum yield (in micromol C/micromol photon) */" << endl;
    output_troll << g1 << "\t/* parameter g1 of Medlyn et al ’s stomatal conductance model */" << endl;
    output_troll << "0.10\t/* vC  # variance of the flexion moment */" << endl;
    output_troll << "0.005\t/* DBH0 # initial dbh (m) */" << endl;
    output_troll << "0.95\t/* H0 # initial height (m) */" << endl;
    output_troll << "0.3\t/* CR_min # minimum crown radius (in m) */" << endl;
    output_troll << a_CR << "\t/* CR_a # CR log intercept */" << endl;
    output_troll << b_CR << "\t/* CR_b # CR log slope */" << endl;
    output_troll << "0.3\t/* de0 # initial crown depth(in m) */" << endl;
    output_troll << shape_crown << "\t/* shape_crown # fraction of crown radius at crown top */" << endl;
    output_troll << "1\t/* dens # leaf density (m^2/m^2) */" << endl;
    output_troll << "0.35\t/* fraction of biomass allocated to above ground wood (branch turnover+stem) */" << endl;
    output_troll << "0.25\t/* fraction of biomass allocated to canopy (leaves + reproductive organs + twigs) */" << endl;
    output_troll << "50000\t/* constant used to scale total seed rain per hectare across species (in next computation) */" << endl;
    output_troll << "10\t/* nb of seeds produced and dispersed by each mature tree when SEEDTRADEOFF is not defined */" << endl;
    output_troll << sigma_height << "\t/* sigma_height # intraspecific variation in tree height (lognormal) */" << endl;
    output_troll << sigma_CR << "\t/* sigma_CR # intraspecific variation in crown radius (lognormal) */" << endl;
    output_troll << sigma_CD << "\t/* sigma_CD # intraspecific variation in crown depth (lognormal) */" << endl;
    output_troll << sigma_P << "\t/* sigma_P # intraspecific variation in leaf phosphorus (lognormal) */" << endl;
    output_troll << sigma_N << "\t/* sigma_N # intraspecific variation in leaf nitrogen (lognormal) */" << endl;
    output_troll << sigma_LMA << "\t/* sigma_LMA # intraspecific variation in LMA (lognormal) */" << endl;
    output_troll << sigma_wsg << "\t/* sigma_wsg # intraspecific variation in wood specific gravity */" << endl;
    output_troll << "0.05\t/* sigma_dmax # intraspecific variation in maximum diameter */" << endl;
    output_troll << corr_CR_height << "\t/* corr_CR_height # correlation coefficient between crown radius and tree height */" << endl;
    output_troll << corr_N_P << "\t/* corr_N_P # correlation coefficient between leaf nitrogen and leaf phosphorus */" << endl;
    output_troll << corr_N_LMA << "\t/* corr_N_LMA # correlation coefficient between leaf nitrogen and LMA */" << endl;
    output_troll << corr_P_LMA << "\t/* corr_P_LMA # correlation coefficient between leaf phosphorus and LMA */" << endl;
    output_troll << "30\t/* leafdem_resolution # resolution of leaf demography model */" << endl;
    output_troll << "1\t/* p_tfsecondary # probability of secondary treefall */" << endl;
    output_troll << "0\t/* hurt_decay # parameter determining how tree damages are repaired */" << endl;
    output_troll << crown_gap_fraction << "\t/* crown_gap_fraction # fraction of gaps in the crown */" << endl;
    output_troll << "0.03\t/* minimal death rate */" << endl;
    output_troll << "0.03\t/* m1 (slope of death rate) */" << endl;
    output_troll << Cair << "\t/* atmospheric CO2 concentration in micromol/mol */" << endl;
    output_troll << "1\t/* LL parameterizations: Reich empirical, Kitajima model, and Kitajima model with leaf plasticity (0,1,2) */" << endl;
    output_troll << "2\t/* dynamic LA regulation: off, 1.0, 0.75, or 0.5 (0,1,2,3) */" << endl;
    output_troll << "2\t/* sapwood parameterizations: constant thickness (0.04), Fyllas percentage, Fyllas lower limit (0,1,2) */" << endl;
    output_troll << "0\t/* agb parameterization without or including crown (0,1) */" << endl;
    output_troll << "0\t/* excess biomass into seeds after maturation (0,1) */" << endl;
    output_troll << endl;
    output_troll << "### Species description ###" << endl;
    output_troll << "species\tLMA\tNmass\tPmass\twsg\tdmax\thmax\tah\tseedvolume\tFreg\tsp_label" << endl;

    output_troll << "Community_mean"<< "\t" << mean_LMA << "\t" << mean_N << "\t" << mean_P << "\t" << mean_wsg << "\t" << dbh_limittop << "\t" << b_height << "\t" << a_height << "\t" << 1.0 << "\t" << 1.0 << "\t" << 1 << endl;
}
    
void OutputTrees(){
    /* and now the trees */
    output_trees << "col\trow\tdbh\theight\tCR\tCD\tvar_height\tvar_CR\tvar_CD\tCrownDisplacement\tNmass\tPmass\tLMA\twsg\tleafarea\tspecies_label\tspecies\tLAIabove\tLAImax\tLAI\tdens_top\tdens_belowtop\tdens_belowtop\tGPP\tNPP\tCrownVolume\tbin\tagb";
    
    for(int i = 1; i < 4; i++){
         float radiusfactor = powf(2.0,float(i));
         output_trees << "\tfactor" << radiusfactor;
    }
     
    output_trees << "\tcanopydistance" << endl;
    
    for(int site = 0; site < sites; site++){
        if(T[site].t_age > 0){
            
            /* determine leafarea of the tree */
            /* in theory this is just LAI * crown area, but in practice, since we have heterogeneities inside the crown, we need to account for them */
            float dbh = T[site].t_dbh;
            float log_dbh = log10(dbh * 100.0);
            int bin_tree = int(log_dbh/dbhlog_binsize);
            
            float CR = T[site].t_Crown_Radius;
            
            /* we calculate the leafarea */
            /* upwards deviation of crown radius converge asymptotically to a gap fraction of 100%, lower cutoff ensures that trees don't exceed their LAImax */
           
            float CR_mean = exp(a_CR + b_CR * log(dbh));
            float multiplier_CR = CR/CR_mean;
            float fraction_filled_general = 1.0 - crown_gap_fraction;
            float fraction_filled_target = minf(fraction_filled_general/(multiplier_CR * multiplier_CR),1.0);

            float crown_area_filled = GetCrownAreaFilled(CR, fraction_filled_target);
            
            float LAI = T[site].t_LAI;
            float leafarea = LAI * (crown_area_filled);
            
            T[site].CrownVolume();
            
            float mass_cylinder = T[site].t_wsg*T[site].t_Tree_Height*T[site].t_dbh*T[site].t_dbh; // the mass of the stem if assumed to be cylindric. This is converted into real biomass using tapering factors, etc.
            float agb = 0.0673*pow(mass_cylinder*10000, 0.976);
            
            output_trees << T[site].t_site%cols << "\t" << T[site].t_site/cols << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth <<"\t" << T[site].t_dev_height << "\t" << T[site].t_dev_CR << "\t" << T[site].t_dev_CD  << "\t" << T[site].t_CrownDisplacement << "\t" << T[site].t_Nmass << "\t" << T[site].t_Pmass << "\t" << T[site].t_LMA << "\t" << T[site].t_wsg << "\t" << leafarea << "\t" << T[site].t_species_label << "\t" << "Community_mean" << "\t" << T[site].t_LAIabove << "\t" << T[site].t_LAImax << "\t" << T[site].t_LAI << "\t" << T[site].t_dens_top << "\t" << T[site].t_dens_belowtop << "\t" << T[site].t_dens_belowtop << "\t" << T[site].t_GPP << "\t" << T[site].t_NPP << "\t" << T[site].t_CrownVolume << "\t" << bin_tree << "\t" << agb;
            
            /* add canopy exposure */
            int site_crown = site + T[site].t_CrownDisplacement;
            for(int i = 1; i < 4; i++){
                float radiusfactor = powf(2.0,float(i));
                float canopyexposure = GetCrownExposure(false,site_crown, T[site].t_Tree_Height, T[site].t_Crown_Radius, T[site].t_Crown_Depth, radiusfactor);
                output_trees << "\t" << canopyexposure;
            }
            
            /* and distance to the canopy above */
            float canopydistance = GetCanopyStatus(site_crown, T[site].t_Tree_Height, T[site].t_Crown_Radius, T[site].t_Crown_Depth);
            output_trees << "\t" << canopydistance;
        
            output_trees << endl;
        }
    }
}

void FreeMem(){
    delete [] T;

    for (int h=0; h<(height_max+1); h++) {
        delete [] LAI3D[h];
        delete [] Voxel3D[h];
        delete [] transmittance_simulated[h];
        delete [] ALS_sampling[h];
        delete [] ALS_echos[h];
    }
    delete [] LAI3D;
    delete [] Voxel3D;
    
    delete [] transmittance_simulated;
    delete [] ALS_sampling;
    delete [] ALS_echos;

    delete [] chm_empirical;
    delete [] hist_CHM_emp;
    delete [] hist_CHM_empcompcumulative;
    
#ifdef CAMPAIGN_REPEATED
    delete [] chm_empirical_rep;
    delete [] hist_CHM_emp_rep;
    delete [] hist_CHM_empcompcumulative_rep;
#endif
    
    for(int bin=0;bin < dbhlog_binnb; bin++) {
        delete [] correlation_structure[bin];
    }
    
    delete [] correlation_structure;

    delete [] Temperature;
    delete [] DailyMaxTemperature;
    delete [] NightTemperature;
    delete [] Rainfall;
    delete [] WindSpeed;
    delete [] MaxIrradiance;
    delete [] MeanIrradiance;
    delete [] SaturatedVapourPressure;
    delete [] VapourPressure;
    delete [] VapourPressureDeficit;
    delete [] DailyVapourPressureDeficit;
    delete [] DailyMaxVapourPressureDeficit;
    
    delete [] LookUp_T;
    delete [] LookUp_KmT;
    delete [] LookUp_VPD;
    delete [] LookUp_flux_absorption;
    delete [] LookUp_flux;
    delete [] LookUp_Rday;
    delete [] LookUp_JmaxT;
    delete [] LookUp_Rstem;
    delete [] LookUp_GammaT;
    delete [] LookUp_Rnight;
    delete [] LookUp_VcmaxT;
}
