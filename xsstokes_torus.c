/* stokes_torus v2.0 - polarized reflection from an axially symmetric surface
 *      of an optically thick circular or elliptical torus illuminated
 *      by an (un)polarised power law,
 *
 *        - an example of a polarisation subroutine for XSPEC using tables 
 *          computed with torus_integrator.py code
 *
 *    References: Podgorny et al. (2024), Podgorny (2025)
 *
 * par1 ... PhoIndex - power-law photon index of the primary flux
 * par2 ... cos_incl - cosine of the observer inclination (1.-pole, 0.-disc),
 *                    with XSET you can recalculate to inclination in degrees
 * par3 ... trTheta - transformed torus half-opening angle between 0 and 1
 *                    with XSET you can recalculate to Theta in degrees
 * par4 ... rhorhoin - the skew in the units of rho_in of the inner walls of the torus
 *                    between 1 and 2, if this is parameter is set to -1,
 *                    the geometry is fixed to a perfectly circular torus
 *                    (then keep rho frozen)
 * par5 ... xi0 - the ionization parameter at the inner-most ring of the torus
 *                in the equatorial plane
 *                     xi0 >= 5 - it uses the local reflection STOKES tables
 *                           for partial ionization
 *                     5 > xi0 > 0 - it assumes fully neutral reflection and uses
 *                           the fully neutral local reflection STOKES tables
 *                     xi0 < 0 - it assumes fully ionized reflection with
 *                           100% albedo and uses Chandrasekhar's diffuse
 *                           reflection prescription; in this case
 *                           the model is renormalized to X-ray luminosity between
 *                           10^(-1.1) and 10^(2.4) keV
 *                           via L_X = -xi0 * n_H0 * (rho_in)^2 where
 *                           n_H0 = 10^{15} cm^{-3} as for the table reprocessing,
 *                           rho_in is the distance between the torus and the center;
 *                           with XSET you can calculate L_{2-10 keV} luminosity
 *                           in erg/s and rho_in in pc, provided the norm value,
 *                           expected n_H0 value and D_Mpc distance
 *                           to the source (variables  "NORMVAL", "D_MPC"
 *                           in Mpc and "NH0"  in cm^{-3}), if additionally the central
 *                           object's mass is provided in XSET (variable "MASS"
 *                           in M_sun), then rho_in in r_g is calculated
 * par6 ... beta - the power-law index for the ionization parameter profile
 *                 across the surface with distance from the source; this
 *                 parameter is a proxy for density profile; keep it frozen
 *                 at 2 for a constant density across the surface; do not currently
 *                 interpolate in this parameter, use only -2, 2, 6 values and
 *                 do not switch between them within one XSPEC session to prevent
 *                 excessive memory allocation
 * par7 ... NpNr - a switch between
 *                    = 1 - reflection + primary
 *                    = 0 - reflection only
 * par8 ... pol_deg  - anisotropy and intrinsic polarisation degree of
 *                    primary radiation
 *                    >= 0 - polarisation degree of isotropic central emission
 *                    = -1 - a particular prescription for anisotropy and
 *                           polarisation of a slab corona (similar to
 *                           a typical lamp-post polarisation), keep par8
 *                           frozen, then par9 (chi) is not used
 * par9 ... chi - intrinsic polarisation angle (in degrees, -90 < chi < 90)
 *		       of primary radiation, the orientation is degenarate by 
 *		       180 degrees
 * par10 ... pos_ang - orientation of the system (in degrees, -90 < pos_ang < 90),
 *                    the orientation is degenarate by 180 degrees
 * par11 ... zshift  - overall redshift
 * par12 ... Stokes  - what should be stored in photar() array, i.e. as output
 *                    = -1 - the output is defined according to the XFLT0001 
 *                           keyword of the SPECTRUM extension of the data file,
 *                           where "Stokes:0" means photon number density flux,
 *                           "Stokes:1" means Stokes parameter Q devided by 
 *                           energy and "Stokes:2" means Stokes parameter U 
 *                           devided by energy
 *                    =  0 - array of photon number density flux per bin
 *                          (array of Stokes parameter I devided by energy)
 *                           with the polarisation computations switched off
 *                          (unpolarised primary is assumed)
 *                    =  1 - array of photon number density flux per bin
 *                          (array of Stokes parameter I devided by energy),
 *                           with the polarisation computations switched on
 *                    =  2 - array of Stokes parameter Q devided by energy
 *                    =  3 - array of Stokes parameter U devided by energy
 *                    =  4 - array of Stokes parameter V devided by energy
 *                    =  5 - array of degree of polarization
 *                    =  6 - array of polarization angle psi=0.5*atan(U/Q)
 *                    =  7 - array of "Stokes" angle
 *                           beta=0.5*asin(V/sqrt(Q*Q+U*U+V*V))
 *                    =  8 - array of Stokes parameter Q devided by I
 *                    =  9 - array of Stokes parameter U devided by I
 *                    = 10 - array of Stokes parameter V devided by I
 * par13 ... norm - equal to (1e2*rho_in/D_Mpc)^2*(n_H0/1e15) where rho_in is
 *                  the physical distance between the center and the torus in pc,
 *                  D_Mpc is the distance to the source in Mpc, and n_H0 is the
 *                  expected density at the inner edge of the torus in cm^{-3};
 *                  with XSET you can obtain 2-10 keV luminosity in erg/s
 *                  for a given norm (re-enter the norm value as XSET variable
 *                  "NORMVAL") and D_Mpc in Mpc (variable "D_MPC"); and additionally
 *                  rho_in in pc for a given expected n_H0 in cm^{-3} (variable
 *                  "NH0") and rho_in in r_g, provided additionally the central
 *                  object's mass (variable "MASS" in M_sun)
 *
 ******************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>

/*******************************************************************************
*******************************************************************************/
#ifdef OUTSIDE_XSPEC

#define IFL    1
#define NPARAM 12
#define NE     300
#define E_MIN  0.1
#define E_MAX  100.


int main() {

int stokesxitorus(const double *ear, int ne, const double *param, int ifl,
               double *photar, double *photer, const char* init);

double ear[NE+1], photar[NE], photer[NE], param[NPARAM];
char   initstr[0] = "";
int    ie;

param[ 0] = 2.0;        // PhoIndex
param[ 1] = 0.775;      // cos_incl
param[ 2] = 0.5;        // trTheta
param[ 3] = 1.5;        // rho/rho_in
param[ 4] = 100.;        // xi0
param[ 5] = 2.;        // beta
param[ 6] = 0.;        // NpNr
param[ 7] = 0.;         // pol_deg
param[ 8] = 0.;         // chi
param[ 9] = 0.;         // pos_ang
param[ 10] = 0.;         // zshift
param[ 11] = 1.;         // STOKES

for(ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

stokesxitorus(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define REFSPECTRAneu1 "tables/stokes-n-iso-UNPOL-torus.fits" // isotropic fully neutral UNPOLARISED
#define REFSPECTRAneu2 "tables/stokes-n-iso-VRPOL-torus.fits" // isotropic fully neutral VERTICALLY POLARISED
#define REFSPECTRAneu3 "tables/stokes-n-iso-45DEG-torus.fits" // isotropic fully neutral DIAGONALLY POLARISED
#define REFSPECTRAneu4 "tables/stokes-n-slab-torus.fits" // slab fully neutral
#define REFSPECTRAion1 "tables/stokes-i-iso-UNPOL-torus.fits" // isotropic fully ionized UNPOLARISED
#define REFSPECTRAion2 "tables/stokes-i-iso-VRPOL-torus.fits" // isotropic fully ionized VERTICALLY POLARISED
#define REFSPECTRAion3 "tables/stokes-i-iso-45DEG-torus.fits" // isotropic fully ionized DIAGONALLY POLARISED
#define REFSPECTRAion4 "tables/stokes-i-slab-torus.fits" // slab fully ionized

#define cREFSPECTRAneu1 "tables/stokes-n-iso-UNPOL-ctorus.fits" // isotropic fully neutral UNPOLARISED circular torus
#define cREFSPECTRAneu2 "tables/stokes-n-iso-VRPOL-ctorus.fits" // isotropic fully neutral VERTICALLY POLARISED circular torus
#define cREFSPECTRAneu3 "tables/stokes-n-iso-45DEG-ctorus.fits" // isotropic fully neutral DIAGONALLY POLARISED circular torus
#define cREFSPECTRAneu4 "tables/stokes-n-slab-ctorus.fits" // slab fully neutral circular torus
#define cREFSPECTRAion1 "tables/stokes-i-iso-UNPOL-ctorus.fits" // isotropic fully ionized UNPOLARISED circular torus
#define cREFSPECTRAion2 "tables/stokes-i-iso-VRPOL-ctorus.fits" // isotropic fully ionized VERTICALLY POLARISED circular torus
#define cREFSPECTRAion3 "tables/stokes-i-iso-45DEG-ctorus.fits" // isotropic fully ionized DIAGONALLY POLARISED circular torus
#define cREFSPECTRAion4 "tables/stokes-i-slab-ctorus.fits" // slab fully ionized circular torus

#define VISIBILITY_FILE "visibility_line.txt" // Text file containing Theta_limit for each rho/rho_in
#define cVISIBILITY_FILE "visibility_line_c.txt" // Text file containing Theta_limit
#define MAX_SIZE_c 40000  // Size of the visibility_line_c.txt text file
#define MAX_SIZE 4000000  // Size of the visibility_line.txt text file
#define MAX_VIS_THETA  2000
#define MAX_VIS_RHO     256
#define MAX_LINE_LEN  32768
#define VIS_EPS       1.0e-12

#define PI   3.14159265358979
#define E0 0.079432823
#define EC 251.18864
#define ERG 6.241509e8
#define CLIGHT 29979245800.
#define MSUN 1.9891e33
#define GCONST 6.6743e-8
#define PCCM 3.08567758128e18
#define NPAR 6
#define NPARc 5
#define NPARni 7

extern int    xs_write(char* wrtstr, int idest);
extern float  DGFILT(int ifl, const char* key);
extern void   FPMSTR(const char* value1, const char* value2);
extern char*  FGMSTR(char* dname);
extern void   tabintxflt(float* ear, int ne, float* param, const int npar, 
                         const char* filenm, const char **xfltname, 
                         const float *xfltvalue, const int nxflt,
                         const char* tabtyp, float* photar, float* photer);

static int get_optional_xset_double(const char *name, double *value)
{
    char name_buf[128];
    char *pname;
    char *xset_string;
    char *ptr;
    char *endptr;

    if (name == NULL || value == NULL) {
        fprintf(stderr, "get_optional_xset_double: null input.\n");
        return -1;
    }

    snprintf(name_buf, sizeof(name_buf), "%s", name);
    pname = name_buf;

    xset_string = FGMSTR(pname);

    /*
     * Depending on XSPEC/version/wrapper behavior, an unset string may come
     * back as NULL or as an empty/blank string. Treat both as "not set".
     */
    if (xset_string == NULL) {
        return 0;
    }

    ptr = xset_string;
    while (*ptr != '\0' && isspace((unsigned char)*ptr)) {
        ptr++;
    }

    if (*ptr == '\0') {
        return 0;
    }

    errno = 0;
    *value = strtod(ptr, &endptr);

    while (*endptr != '\0' && isspace((unsigned char)*endptr)) {
        endptr++;
    }

    if (endptr == ptr || *endptr != '\0' || errno == ERANGE || !isfinite(*value)) {
        fprintf(stderr,
                "Invalid XSPEC xset variable %s = '%s'. Please set a numeric value or leave it unset.\n",
                name, xset_string);
        return -1;
    }

    return 1;
}


int stokesxitorus(const double *ear, int ne, const double *param, int ifl,
            double *photar, double *photer, const char* init) {


FILE   *fw;
//static int first = 1;
int status = 0;
static char   xsdir[255]="";
static char   pname[128]="XSDIR", ptrue_Theta_out[128] = "Theta_degrees", pinc_degrees[128] = "inc_degrees";
static char   pRefl[128] = "RF", prho_in_out[128] = "rho_in_pc", prho_in_out_r_g[128] = "rho_in_r_g", pL_out[128] = "L";
static char refspectra[16][255], visibility_path[255], visibility_path_c[255];

// - if set try XSDIR directory, otherwise look in the working directory
// Initialize refspectra elements and visibility file path
{
    char *xsdir_tmp = FGMSTR(pname);

    if (xsdir_tmp == NULL) {
        xsdir[0] = '\0';
    } else {
        snprintf(xsdir, sizeof(xsdir), "%s", xsdir_tmp);
    }
}
//sprintf(xsdir, "%s", FGMSTR(pname));
if (strlen(xsdir) == 0) {
    strcpy(refspectra[0], REFSPECTRAneu1);
    strcpy(refspectra[1], REFSPECTRAneu2);
    strcpy(refspectra[2], REFSPECTRAneu3);
    strcpy(refspectra[3], REFSPECTRAneu4);
    strcpy(refspectra[4], REFSPECTRAion1);
    strcpy(refspectra[5], REFSPECTRAion2);
    strcpy(refspectra[6], REFSPECTRAion3);
    strcpy(refspectra[7], REFSPECTRAion4);
    strcpy(refspectra[8], cREFSPECTRAneu1);
    strcpy(refspectra[9], cREFSPECTRAneu2);
    strcpy(refspectra[10], cREFSPECTRAneu3);
    strcpy(refspectra[11], cREFSPECTRAneu4);
    strcpy(refspectra[12], cREFSPECTRAion1);
    strcpy(refspectra[13], cREFSPECTRAion2);
    strcpy(refspectra[14], cREFSPECTRAion3);
    strcpy(refspectra[15], cREFSPECTRAion4);

    strcpy(visibility_path, VISIBILITY_FILE);
    strcpy(visibility_path_c, cVISIBILITY_FILE);
} else {
    if (xsdir[strlen(xsdir) - 1] == '/') {
        sprintf(refspectra[0], "%s%s", xsdir, REFSPECTRAneu1);
        sprintf(refspectra[1], "%s%s", xsdir, REFSPECTRAneu2);
        sprintf(refspectra[2], "%s%s", xsdir, REFSPECTRAneu3);
        sprintf(refspectra[3], "%s%s", xsdir, REFSPECTRAneu4);
        sprintf(refspectra[4], "%s%s", xsdir, REFSPECTRAion1);
        sprintf(refspectra[5], "%s%s", xsdir, REFSPECTRAion2);
        sprintf(refspectra[6], "%s%s", xsdir, REFSPECTRAion3);
        sprintf(refspectra[7], "%s%s", xsdir, REFSPECTRAion4);
        sprintf(refspectra[8], "%s%s", xsdir, cREFSPECTRAneu1);
        sprintf(refspectra[9], "%s%s", xsdir, cREFSPECTRAneu2);
        sprintf(refspectra[10], "%s%s", xsdir, cREFSPECTRAneu3);
        sprintf(refspectra[11], "%s%s", xsdir, cREFSPECTRAneu4);
        sprintf(refspectra[12], "%s%s", xsdir, cREFSPECTRAion1);
        sprintf(refspectra[13], "%s%s", xsdir, cREFSPECTRAion2);
        sprintf(refspectra[14], "%s%s", xsdir, cREFSPECTRAion3);
        sprintf(refspectra[15], "%s%s", xsdir, cREFSPECTRAion4);
        sprintf(visibility_path, "%s%s", xsdir, VISIBILITY_FILE);
        sprintf(visibility_path_c, "%s%s", xsdir, cVISIBILITY_FILE);
    } else {
        sprintf(refspectra[0], "%s/%s", xsdir, REFSPECTRAneu1);
        sprintf(refspectra[1], "%s/%s", xsdir, REFSPECTRAneu2);
        sprintf(refspectra[2], "%s/%s", xsdir, REFSPECTRAneu3);
        sprintf(refspectra[3], "%s/%s", xsdir, REFSPECTRAneu4);
        sprintf(refspectra[4], "%s/%s", xsdir, REFSPECTRAion1);
        sprintf(refspectra[5], "%s/%s", xsdir, REFSPECTRAion2);
        sprintf(refspectra[6], "%s/%s", xsdir, REFSPECTRAion3);
        sprintf(refspectra[7], "%s/%s", xsdir, REFSPECTRAion4);
        sprintf(refspectra[8], "%s/%s", xsdir, cREFSPECTRAneu1);
        sprintf(refspectra[9], "%s/%s", xsdir, cREFSPECTRAneu2);
        sprintf(refspectra[10], "%s/%s", xsdir, cREFSPECTRAneu3);
        sprintf(refspectra[11], "%s/%s", xsdir, cREFSPECTRAneu4);
        sprintf(refspectra[12], "%s/%s", xsdir, cREFSPECTRAion1);
        sprintf(refspectra[13], "%s/%s", xsdir, cREFSPECTRAion2);
        sprintf(refspectra[14], "%s/%s", xsdir, cREFSPECTRAion3);
        sprintf(refspectra[15], "%s/%s", xsdir, cREFSPECTRAion4);
        sprintf(visibility_path, "%s/%s", xsdir, VISIBILITY_FILE);
        sprintf(visibility_path_c, "%s/%s", xsdir, cVISIBILITY_FILE);
    }
}

int    i, j, ie, stokes;
double pol_deg, pos_ang, chi, p0, chi0;
const char*   xfltname = "Stokes";
float  xfltvalue;
float  Svector[3][ne];
float  Smatrix[9][ne];
float  fl_paramc[NPARc]={(float) param[0], (float) param[1], (float) param[2], (float) param[4], (float) param[10]};
float  fl_param[NPAR]={(float) param[0], (float) param[1], (float) param[2], (float) param[4], (float) param[3], (float) param[10]};
float  fl_paramni[NPARni]={(float) param[0], (float) param[1], (float) param[2], (float) param[5], (float) param[3], 1., (float) param[10]};
static char xi_files[8][3][255];
static const char *dummy_xfltname = "none";
static float dummy_xfltvalue = 0.0f;
const char*  tabtyp="add";
float  fl_ear[ne+1], fl_photer[ne];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne], 
       qar_final[ne], uar_final[ne], Dnorm, Dnormprim, refl_ratio, flux_refl, flux_prim;
double pamin, pamax, pa2min, pa2max, D, normval, L, L_X, mass, Emid;
double zzshift;
double norm_energy;
double photar1[ne];

FILE *file;
double vis_Theta_c[MAX_SIZE_c], vis_inc_c[MAX_SIZE_c];
static double vis_theta[MAX_VIS_THETA];
static double vis_rho[MAX_VIS_RHO];
static double vis_inc[MAX_VIS_THETA][MAX_VIS_RHO];
static double vis_inc_r0[MAX_VIS_THETA];
char line[MAX_LINE_LEN];
char *tok;
char *endptr;
int ntheta;
int nrho;
int it;
int ir;
int irho;
int found_theta;
int D_status;
int nH0_status;
int norm_status;
int M_status;
double r0;        /* assign this before using the interpolation block */
double r1;
double r2;
double wr;
double th1;
double th2;
double z1;
double z2;
double x0, x1, y1, x2, y2;
int count = 0;
double trTheta, mue_tot, Theta_limit, min_Theta, max_Theta, true_Theta;
double xival, betaval, rhoval, redshift, NpNr, Eh, El, gam, a0, rho_in, rho_in_cm, rho_in_pc, rho_in_r_g, n_H0;
char true_Theta_out[32], inc_degrees[32], Refl[32], rho_in_out[32], rho_in_out_r_g[32], L_out[32];

gam = (float) param[0];
mue_tot = (float) param[1];
trTheta = (float) param[2];
rhoval = (float) param[3];
xival = (float) param[4];
betaval = (float) param[5];
NpNr = (float) param[6];
pol_deg = param[7];
chi = param[8]/180.*PI;
pos_ang = param[9]/180.*PI;
redshift = (float) param[10];
stokes = (int) param[11];

if (strlen(xsdir) == 0) {
    snprintf(xi_files[0][0], sizeof(xi_files[0][0]),
             "tables/stokes-xi-iso-UNPOL-torus-beta2.fits");
    snprintf(xi_files[0][1], sizeof(xi_files[0][1]),
             "tables/stokes-xi-iso-UNPOL-torus-beta-2.fits");
    snprintf(xi_files[0][2], sizeof(xi_files[0][2]),
             "tables/stokes-xi-iso-UNPOL-torus-beta6.fits");
    snprintf(xi_files[1][0], sizeof(xi_files[1][0]),
             "tables/stokes-xi-iso-VRPOL-torus-beta2.fits");
    snprintf(xi_files[1][1], sizeof(xi_files[1][1]),
             "tables/stokes-xi-iso-VRPOL-torus-beta-2.fits");
    snprintf(xi_files[1][2], sizeof(xi_files[1][2]),
             "tables/stokes-xi-iso-VRPOL-torus-beta6.fits");
    snprintf(xi_files[2][0], sizeof(xi_files[2][0]),
             "tables/stokes-xi-iso-45DEG-torus-beta2.fits");
    snprintf(xi_files[2][1], sizeof(xi_files[2][1]),
             "tables/stokes-xi-iso-45DEG-torus-beta-2.fits");
    snprintf(xi_files[2][2], sizeof(xi_files[2][2]),
             "tables/stokes-xi-iso-45DEG-torus-beta6.fits");
    snprintf(xi_files[3][0], sizeof(xi_files[3][0]),
             "tables/stokes-xi-slab-torus-beta2.fits");
    snprintf(xi_files[3][1], sizeof(xi_files[3][1]),
             "tables/stokes-xi-slab-torus-beta-2.fits");
    snprintf(xi_files[3][2], sizeof(xi_files[3][2]),
             "tables/stokes-xi-slab-torus-beta6.fits");
    snprintf(xi_files[4][0], sizeof(xi_files[4][0]),
             "tables/stokes-xi-iso-UNPOL-ctorus-beta2.fits");
    snprintf(xi_files[4][1], sizeof(xi_files[4][1]),
             "tables/stokes-xi-iso-UNPOL-ctorus-beta-2.fits");
    snprintf(xi_files[4][2], sizeof(xi_files[4][2]),
             "tables/stokes-xi-iso-UNPOL-ctorus-beta6.fits");
    snprintf(xi_files[5][0], sizeof(xi_files[5][0]),
             "tables/stokes-xi-iso-VRPOL-ctorus-beta2.fits");
    snprintf(xi_files[5][1], sizeof(xi_files[5][1]),
             "tables/stokes-xi-iso-VRPOL-ctorus-beta-2.fits");
    snprintf(xi_files[5][2], sizeof(xi_files[5][2]),
             "tables/stokes-xi-iso-VRPOL-ctorus-beta6.fits");
    snprintf(xi_files[6][0], sizeof(xi_files[6][0]),
             "tables/stokes-xi-iso-45DEG-ctorus-beta2.fits");
    snprintf(xi_files[6][1], sizeof(xi_files[6][1]),
             "tables/stokes-xi-iso-45DEG-ctorus-beta-2.fits");
    snprintf(xi_files[6][2], sizeof(xi_files[6][2]),
             "tables/stokes-xi-iso-45DEG-ctorus-beta6.fits");
    snprintf(xi_files[7][0], sizeof(xi_files[7][0]),
             "tables/stokes-xi-slab-ctorus-beta2.fits");
    snprintf(xi_files[7][1], sizeof(xi_files[7][1]),
             "tables/stokes-xi-slab-ctorus-beta-2.fits");
    snprintf(xi_files[7][2], sizeof(xi_files[7][2]),
             "tables/stokes-xi-slab-ctorus-beta6.fits");
} else if (xsdir[strlen(xsdir) - 1] == '/') {
    snprintf(xi_files[0][0], sizeof(xi_files[0][0]),
             "%stables/stokes-xi-iso-UNPOL-torus-beta2.fits", xsdir);
    snprintf(xi_files[0][1], sizeof(xi_files[0][1]),
             "%stables/stokes-xi-iso-UNPOL-torus-beta-2.fits", xsdir);
    snprintf(xi_files[0][2], sizeof(xi_files[0][2]),
             "%stables/stokes-xi-iso-UNPOL-torus-beta6.fits", xsdir);
    snprintf(xi_files[1][0], sizeof(xi_files[1][0]),
             "%stables/stokes-xi-iso-VRPOL-torus-beta2.fits", xsdir);
    snprintf(xi_files[1][1], sizeof(xi_files[1][1]),
             "%stables/stokes-xi-iso-VRPOL-torus-beta-2.fits", xsdir);
    snprintf(xi_files[1][2], sizeof(xi_files[1][2]),
             "%stables/stokes-xi-iso-VRPOL-torus-beta6.fits", xsdir);
    snprintf(xi_files[2][0], sizeof(xi_files[2][0]),
             "%stables/stokes-xi-iso-45DEG-torus-beta2.fits", xsdir);
    snprintf(xi_files[2][1], sizeof(xi_files[2][1]),
             "%stables/stokes-xi-iso-45DEG-torus-beta-2.fits", xsdir);
    snprintf(xi_files[2][2], sizeof(xi_files[2][2]),
             "%stables/stokes-xi-iso-45DEG-torus-beta6.fits", xsdir);
    snprintf(xi_files[3][0], sizeof(xi_files[3][0]),
             "%stables/stokes-xi-slab-torus-beta2.fits", xsdir);
    snprintf(xi_files[3][1], sizeof(xi_files[3][1]),
             "%stables/stokes-xi-slab-torus-beta-2.fits", xsdir);
    snprintf(xi_files[3][2], sizeof(xi_files[3][2]),
             "%stables/stokes-xi-slab-torus-beta6.fits", xsdir);
    snprintf(xi_files[4][0], sizeof(xi_files[4][0]),
             "%stables/stokes-xi-iso-UNPOL-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[4][1], sizeof(xi_files[4][1]),
             "%stables/stokes-xi-iso-UNPOL-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[4][2], sizeof(xi_files[4][2]),
             "%stables/stokes-xi-iso-UNPOL-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[5][0], sizeof(xi_files[5][0]),
             "%stables/stokes-xi-iso-VRPOL-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[5][1], sizeof(xi_files[5][1]),
             "%stables/stokes-xi-iso-VRPOL-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[5][2], sizeof(xi_files[5][2]),
             "%stables/stokes-xi-iso-VRPOL-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[6][0], sizeof(xi_files[6][0]),
             "%stables/stokes-xi-iso-45DEG-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[6][1], sizeof(xi_files[6][1]),
             "%stables/stokes-xi-iso-45DEG-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[6][2], sizeof(xi_files[6][2]),
             "%stables/stokes-xi-iso-45DEG-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[7][0], sizeof(xi_files[7][0]),
             "%stables/stokes-xi-slab-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[7][1], sizeof(xi_files[7][1]),
             "%stables/stokes-xi-slab-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[7][2], sizeof(xi_files[7][2]),
             "%stables/stokes-xi-slab-ctorus-beta6.fits", xsdir);
} else {
    snprintf(xi_files[0][0], sizeof(xi_files[0][0]),
             "%s/tables/stokes-xi-iso-UNPOL-torus-beta2.fits", xsdir);
    snprintf(xi_files[0][1], sizeof(xi_files[0][1]),
             "%s/tables/stokes-xi-iso-UNPOL-torus-beta-2.fits", xsdir);
    snprintf(xi_files[0][2], sizeof(xi_files[0][2]),
             "%s/tables/stokes-xi-iso-UNPOL-torus-beta6.fits", xsdir);
    snprintf(xi_files[1][0], sizeof(xi_files[1][0]),
             "%s/tables/stokes-xi-iso-VRPOL-torus-beta2.fits", xsdir);
    snprintf(xi_files[1][1], sizeof(xi_files[1][1]),
             "%s/tables/stokes-xi-iso-VRPOL-torus-beta-2.fits", xsdir);
    snprintf(xi_files[1][2], sizeof(xi_files[1][2]),
             "%s/tables/stokes-xi-iso-VRPOL-torus-beta6.fits", xsdir);
    snprintf(xi_files[2][0], sizeof(xi_files[2][0]),
             "%s/tables/stokes-xi-iso-45DEG-torus-beta2.fits", xsdir);
    snprintf(xi_files[2][1], sizeof(xi_files[2][1]),
             "%s/tables/stokes-xi-iso-45DEG-torus-beta-2.fits", xsdir);
    snprintf(xi_files[2][2], sizeof(xi_files[2][2]),
             "%s/tables/stokes-xi-iso-45DEG-torus-beta6.fits", xsdir);
    snprintf(xi_files[3][0], sizeof(xi_files[3][0]),
             "%s/tables/stokes-xi-slab-torus-beta2.fits", xsdir);
    snprintf(xi_files[3][1], sizeof(xi_files[3][1]),
             "%s/tables/stokes-xi-slab-torus-beta-2.fits", xsdir);
    snprintf(xi_files[3][2], sizeof(xi_files[3][2]),
             "%s/tables/stokes-xi-slab-torus-beta6.fits", xsdir);
    snprintf(xi_files[4][0], sizeof(xi_files[4][0]),
             "%s/tables/stokes-xi-iso-UNPOL-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[4][1], sizeof(xi_files[4][1]),
             "%s/tables/stokes-xi-iso-UNPOL-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[4][2], sizeof(xi_files[4][2]),
             "%s/tables/stokes-xi-iso-UNPOL-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[5][0], sizeof(xi_files[5][0]),
             "%s/tables/stokes-xi-iso-VRPOL-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[5][1], sizeof(xi_files[5][1]),
             "%s/tables/stokes-xi-iso-VRPOL-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[5][2], sizeof(xi_files[5][2]),
             "%s/tables/stokes-xi-iso-VRPOL-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[6][0], sizeof(xi_files[6][0]),
             "%s/tables/stokes-xi-iso-45DEG-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[6][1], sizeof(xi_files[6][1]),
             "%s/tables/stokes-xi-iso-45DEG-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[6][2], sizeof(xi_files[6][2]),
             "%s/tables/stokes-xi-iso-45DEG-ctorus-beta6.fits", xsdir);
    snprintf(xi_files[7][0], sizeof(xi_files[7][0]),
             "%s/tables/stokes-xi-slab-ctorus-beta2.fits", xsdir);
    snprintf(xi_files[7][1], sizeof(xi_files[7][1]),
             "%s/tables/stokes-xi-slab-ctorus-beta-2.fits", xsdir);
    snprintf(xi_files[7][2], sizeof(xi_files[7][2]),
             "%s/tables/stokes-xi-slab-ctorus-beta6.fits", xsdir);
}

if(stokes == -1){
  xfltvalue = DGFILT(ifl, xfltname);
  if (xfltvalue == 0. || xfltvalue == 1. || xfltvalue == 2.){
    stokes = 1 + (int) xfltvalue;
  }
  else {
    xs_write("sttorus: no or wrong information on data type (counts, q, u)", 5);
    xs_write("sttorus: stokes = par8 = 0 (i.e. counts) will be used", 5);
    stokes=0;
  }
}

if (betaval != -2 && betaval != 2 && betaval != 6 && xival >= 5){
    xs_write("sttorus: for partially ionized tables, other beta values than -2, 2, 6 currently cannot be used", 5);
    xs_write("sttorus: beta = par6 = 2 (standard value) will be used", 5);
    betaval = 2;
}

if (xival == 0){
    xs_write("sttorus: xi0 = 0 cannot be used", 5);
    xs_write("sttorus: xi0 = par5 = 1 (standard value for cold reflection) will be used", 5);
    xival = 1.;
}

//Let's read and interpolate the FITS tables that we will need using internal
//XSPEC routine tabintxflt
//Note that we do not use errors here
for(ie = 0; ie <= ne; ie++) fl_ear[ie] = (float) ear[ie];
if(stokes){//we use polarised tables
//  if (first) {
// The status parameter must always be initialized.
  status = 0;
  if (rhoval == -1.){
      if (xival >= 5){
          if (pol_deg == -1){
                if (betaval == -2){
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][1], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }else if (betaval == 2){
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][0], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }else{
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][2], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
              if (betaval == -2){
                for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4+i][1], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }else if (betaval == 2){
                  for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4+i][0], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }else{
                  for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4+i][2], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }
              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }
      }else if (xival >= 0. && xival < 5.){
          if (pol_deg == -1){
                for (i = 0; i <= 2; i++){
                      xfltvalue = (float) i;
                      tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[8+3], &xfltname,
                                &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
              for (i = 0; i <= 2; i++)
                for (j = 0; j <= 2; j++){
                  xfltvalue = (float) j;
                  tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[8+i], &xfltname,
                            &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                  }

              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }

      }else{
        if (pol_deg == -1){
                for (i = 0; i <= 2; i++){
                      xfltvalue = (float) i;
                      tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[12+3], &xfltname,
                                &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
              for (i = 0; i <= 2; i++)
                for (j = 0; j <= 2; j++){
                  xfltvalue = (float) j;
                  tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[12+i], &xfltname,
                            &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                  }

            //  first = 0;
            //  }


              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }
      }

  }else{
      if (xival >= 5){
          if (pol_deg == -1){
               if (betaval == -2){
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][1], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }else if (betaval == 2){
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][0], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }else{
                    for (i = 0; i <= 2; i++){
                        xfltvalue = (float) i;
                        tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][2], &xfltname,
                                    &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                    }
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
                if (betaval == -2){
                for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[i][1], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }else if (betaval == 2){
                  for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[i][0], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }else{
                  for (i = 0; i <= 2; i++)
                    for (j = 0; j <= 2; j++){
                    xfltvalue = (float) j;
                    tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[i][2], &xfltname,
                                    &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                    }
              }

              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }
      }else if (xival >= 0. && xival < 5.){
          if (pol_deg == -1){
                for (i = 0; i <= 2; i++){
                      xfltvalue = (float) i;
                      tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[3], &xfltname,
                                &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
              for (i = 0; i <= 2; i++)
                for (j = 0; j <= 2; j++){
                  xfltvalue = (float) j;
                  tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[i], &xfltname,
                            &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                  }

            //  first = 0;
            //  }


              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }

      }else{
        if (pol_deg == -1){
                for (i = 0; i <= 2; i++){
                      xfltvalue = (float) i;
                      tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[4+3], &xfltname,
                                &xfltvalue, 1, tabtyp, Svector[i], fl_photer);
                }
                for(ie = 0; ie < ne; ie++) {
                    far[ie] = Svector[0][ie];
                    qar[ie] = Svector[1][ie];
                    uar[ie] = Svector[2][ie];
                    var[ie] = 0.;
                }
            }else{
              for (i = 0; i <= 2; i++)
                for (j = 0; j <= 2; j++){
                  xfltvalue = (float) j;
                  tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[4+i], &xfltname,
                            &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);
                  }

            //  first = 0;
            //  }


              //Let's perform the transformation to initial primary polarisation degree and angle
              for(ie = 0; ie < ne; ie++) {
                for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
                for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

                far[ie] = Smatrix[0][ie] +
                                    pol_deg * ( Smatrix[3][ie] * cos(2.*(chi)) +
                                                Smatrix[6][ie] * sin(2.*(chi)) );
                qar[ie] = Smatrix[1][ie] +
                                    pol_deg * ( Smatrix[4][ie] * cos(2.*(chi))+
                                                Smatrix[7][ie] * sin(2.*(chi)) );
                uar[ie] = Smatrix[2][ie] +
                                    pol_deg * ( Smatrix[5][ie] * cos(2.*(chi))+
                                                Smatrix[8][ie] * sin(2.*(chi)) );
                var[ie] = 0.;

                // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
                // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
                // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
                // var[ie] = 0.;
              }
          }
      }
  }
}else{
  // The status parameter must always be initialized.
  status = 0;
  if (rhoval == -1.){
      if (xival >= 5.){
          if (pol_deg == -1){
              if (betaval == -2){
                  xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][1], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }else if (betaval == 2){
                  xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][0], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }else{
                   xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[7][2], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }
              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              if (betaval == -2){
                xfltvalue = 0.;
                tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4][1], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              }else if (betaval == 2){
                xfltvalue = 0.;
                tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4][0], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
              }else{
                 xfltvalue = 0.;
                 tabintxflt(fl_ear, ne, fl_paramc, NPARc, xi_files[4][2], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
              }

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }else if (xival >= 0. && xival < 5.){
          if (pol_deg == -1){

              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[11], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[8], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }else{
          if (pol_deg == -1){

              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[15], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[12], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }

  }else{
      if (xival >= 5.){
          if (pol_deg == -1){
              if (betaval == -2){
                  xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][1], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }else if (betaval == 2){
                  xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][0], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }else{
                   xfltvalue = 0.;
                  tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[3][2], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
                }
              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              if (betaval == -2){
                xfltvalue = 0.;
                tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[0][1], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              }else if (betaval == 2){
                xfltvalue = 0.;
                tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[0][0], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
              }else{
                 xfltvalue = 0.;
                 tabintxflt(fl_ear, ne, fl_param, NPAR, xi_files[0][2], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);
              }
              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }else if (xival >= 0. && xival < 5.){
          if (pol_deg == -1){

              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[3], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[0], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }else{
          if (pol_deg == -1){

              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[7], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }else{
              xfltvalue = 0.;
              tabintxflt(fl_ear, ne, fl_paramni, NPARni, refspectra[4], &xfltname, &xfltvalue,
                        1, tabtyp, Smatrix[0], fl_photer);

              for(ie = 0; ie < ne; ie++){
                far[ie] = Smatrix[0][ie];
              }
          }
      }

  }
}

// let's transform trTheta back to true_theta

// The status parameter must always be initialized.
status = 0;
if (rhoval == -1){
  file = fopen(visibility_path_c, "rt");
  if (file == NULL) {
  perror("Error opening file");
  return -1;
  }

  while (fscanf(file, "%lf %lf", &vis_Theta_c[count], &vis_inc_c[count]) == 2) {
  count++;
  }
  fclose(file);

  if (mue_tot > 1.) mue_tot = 1.;
  if (mue_tot < -1.) mue_tot = -1.;
  x0 = acos(mue_tot) / PI * 180.0;  // Interpolation point

  // Inline linear interpolation
  for (i = 0; i < count - 1; i++) {
  x1 = vis_inc_c[i];
  y1 = vis_Theta_c[i];
  x2 = vis_inc_c[i + 1];
  y2 = vis_Theta_c[i + 1];
  if (x0 >= x1 && x0 <= x2) {
      Theta_limit = y1 + (x0 - x1) * (y2 - y1) / (x2 - x1);
      break;
  }
  }

  // Ensure Theta_limit is set
  if (i == count - 1) {
  fprintf(stderr, "Interpolation point is out of range.\n");
  return -1;
  }

  /////min_Theta = fmax(25.0, Theta_limit);
  min_Theta = fmax(0., Theta_limit);
  max_Theta = 90.0;
  true_Theta = trTheta * (max_Theta - min_Theta) + min_Theta;

  sprintf(true_Theta_out, "%12.6f", true_Theta);
  FPMSTR(ptrue_Theta_out, true_Theta_out);

  sprintf(inc_degrees, "%12.6f", x0);
  FPMSTR(pinc_degrees, inc_degrees);

}else{
  file = fopen(visibility_path, "rt");  /* visibility_path_c should point to visibility_line.txt */
  if (file == NULL) {
      perror("Error opening file");
      return -1;
  }

  /* ---------- Read rho header row ---------- */

  if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "Error reading header from visibility table.\n");
      fclose(file);
      return -1;
  }

  nrho = 0;

  /* First token is something like Theta\rho; skip it. */
  tok = strtok(line, " \t\r\n");
  if (tok == NULL) {
      fprintf(stderr, "Visibility table header is empty.\n");
      fclose(file);
      return -1;
  }

  while ((tok = strtok(NULL, " \t\r\n")) != NULL) {
      if (nrho >= MAX_VIS_RHO) {
          fprintf(stderr, "Too many rho columns in visibility table; increase MAX_VIS_RHO.\n");
          fclose(file);
          return -1;
      }

      vis_rho[nrho] = strtod(tok, &endptr);
      if (endptr == tok) {
          fprintf(stderr, "Could not parse rho value in visibility table header: %s\n", tok);
          fclose(file);
          return -1;
      }

      nrho++;
  }

  if (nrho < 2) {
      fprintf(stderr, "Visibility table needs at least two rho columns for interpolation.\n");
      fclose(file);
      return -1;
  }

  /* ---------- Read Theta rows and inclination grid ---------- */

  ntheta = 0;

  while (fgets(line, sizeof(line), file) != NULL) {
      tok = strtok(line, " \t\r\n");

      /* Skip blank lines. */
      if (tok == NULL) {
          continue;
      }

      if (ntheta >= MAX_VIS_THETA) {
          fprintf(stderr, "Too many Theta rows in visibility table; increase MAX_VIS_THETA.\n");
          fclose(file);
          return -1;
      }

      vis_theta[ntheta] = strtod(tok, &endptr);
      if (endptr == tok) {
          fprintf(stderr, "Could not parse Theta value in visibility table row %d: %s\n",
                  ntheta, tok);
          fclose(file);
          return -1;
      }

      for (ir = 0; ir < nrho; ir++) {
          tok = strtok(NULL, " \t\r\n");
          if (tok == NULL) {
              fprintf(stderr,
                      "Visibility table row %d has too few columns: expected %d rho values.\n",
                      ntheta, nrho);
              fclose(file);
              return -1;
          }

          vis_inc[ntheta][ir] = strtod(tok, &endptr);
          if (endptr == tok) {
              fprintf(stderr,
                      "Could not parse table value at Theta row %d, rho column %d: %s\n",
                      ntheta, ir, tok);
              fclose(file);
              return -1;
          }
      }

      ntheta++;
  }

  fclose(file);

  if (ntheta < 2) {
      fprintf(stderr, "Visibility table needs at least two Theta rows for interpolation.\n");
      return -1;
  }

  /* ---------- Set interpolation points ---------- */

  if (mue_tot > 1.) mue_tot = 1.;
  if (mue_tot < -1.) mue_tot = -1.;
  x0 = acos(mue_tot) / PI * 180.0;  /* inclination in degrees */
  r0 = rhoval;

  /* ---------- Find bracketing rho columns ---------- */

  irho = -1;

  for (ir = 0; ir < nrho - 1; ir++) {
      r1 = vis_rho[ir];
      r2 = vis_rho[ir + 1];

      if ((r0 >= r1 && r0 <= r2) || (r0 >= r2 && r0 <= r1)) {
          irho = ir;
          break;
      }
  }

  if (irho < 0) {
      fprintf(stderr,
              "rho interpolation point is out of range: r0 = %.12g, table range = %.12g to %.12g\n",
              r0, vis_rho[0], vis_rho[nrho - 1]);
      return -1;
  }

  r1 = vis_rho[irho];
  r2 = vis_rho[irho + 1];

  if (fabs(r2 - r1) <= VIS_EPS) {
      wr = 0.0;
  } else {
      wr = (r0 - r1) / (r2 - r1);
  }

  /* ---------- Interpolate table values linearly in rho at each Theta row ---------- */

  for (it = 0; it < ntheta; it++) {
      z1 = vis_inc[it][irho];
      z2 = vis_inc[it][irho + 1];

      vis_inc_r0[it] = z1 + wr * (z2 - z1);
  }

  /* ---------- Invert inclination-vs-Theta curve to get Theta_limit ---------- */

  found_theta = 0;
  Theta_limit = 0.0;

  /* Exact lower endpoint. */
  if (fabs(x0 - vis_inc_r0[0]) <= VIS_EPS) {
      Theta_limit = vis_theta[0];
      found_theta = 1;
  } else {
      for (it = 0; it < ntheta - 1; it++) {
          z1 = vis_inc_r0[it];
          z2 = vis_inc_r0[it + 1];

          th1 = vis_theta[it];
          th2 = vis_theta[it + 1];

          /*
          * Handle flat table segments. If x0 lies exactly on a flat segment,
          * choose the lower-Theta edge of that segment.
          */
          if (fabs(z2 - z1) <= VIS_EPS) {
              if (fabs(x0 - z1) <= VIS_EPS) {
                  Theta_limit = th1;
                  found_theta = 1;
                  break;
              } else {
                  continue;
              }
          }

          if ((x0 >= z1 && x0 <= z2) || (x0 >= z2 && x0 <= z1)) {
              Theta_limit = th1 + (x0 - z1) * (th2 - th1) / (z2 - z1);
              found_theta = 1;
              break;
          }
      }
  }

  if (!found_theta) {
      fprintf(stderr,
              "Inclination interpolation point is out of range: x0 = %.12g deg at r0 = %.12g\n",
              x0, r0);
      return -1;
  }

  /* ---------- Continue as before ---------- */

  min_Theta = fmax(0.0, Theta_limit);
  max_Theta = 90.0;
  true_Theta = trTheta * (max_Theta - min_Theta) + min_Theta;

  sprintf(true_Theta_out, "%12.6f", true_Theta);
  FPMSTR(ptrue_Theta_out, true_Theta_out);

  sprintf(inc_degrees, "%12.6f", x0);
  FPMSTR(pinc_degrees, inc_degrees);

}

// Let's calculate L_{2-10keV} and rho_in
rho_in = rho_in_cm = rho_in_pc = 0.;
D_status = get_optional_xset_double("D_MPC", &D);
norm_status = get_optional_xset_double("NORMVAL", &normval);

if (D_status < 0 || norm_status < 0) {
    return -1; /* D_MPC or NORMVAL was set but invalid */
}
if (D_status == 1 && norm_status == 1) {
    if (xival > 0.) {
        L_X = 1e15 * normval * D * D * xival * PCCM * PCCM / 1e4;
    } else {
        L_X = - 1e15 * normval * D * D * xival * PCCM * PCCM / 1e4;
    }
    if (gam != 2.) {
        L = L_X / (pow(EC, 2. - gam) - pow(E0, 2. - gam)) *
                    (pow(10., 2. - gam) - pow(2., 2. - gam));
    } else {
        L = L_X / log(EC / E0) * log(10. / 2.);
    }
    snprintf(L_out, sizeof(L_out), "%12.6e", L);
    FPMSTR(pL_out, L_out);

    nH0_status = get_optional_xset_double("NH0", &n_H0);
    if (nH0_status < 0) {
        return -1;  /* n_H0 was set but invalid */
    }
    if (nH0_status == 1) {
        if (xival > 0.) {
            rho_in_pc = sqrt(L_X / xival / n_H0 / PCCM / PCCM);
        } else {
            rho_in_pc = sqrt(- L_X / xival / n_H0 / PCCM / PCCM);
        }
        rho_in_cm = rho_in_pc * PCCM;
        snprintf(rho_in_out, sizeof(rho_in_out), "%12.6e", rho_in_pc);
        FPMSTR(prho_in_out, rho_in_out);

        M_status = get_optional_xset_double("MASS", &mass);
        if (M_status < 0) {
            return -1;  /* MASS was set but invalid */
        }
        if (M_status == 1) {
            rho_in_r_g = rho_in_cm * CLIGHT * CLIGHT / GCONST / mass / MSUN;
            snprintf(rho_in_out_r_g, sizeof(rho_in_out_r_g), "%12.6e", rho_in_r_g);
            FPMSTR(prho_in_out_r_g, rho_in_out_r_g);
        }else{
            FPMSTR(prho_in_out_r_g, "not set");
        }
    }else{
        FPMSTR(prho_in_out, "not set");
        FPMSTR(prho_in_out_r_g, "not set");
    }
}else{
    FPMSTR(prho_in_out, "not set");
    FPMSTR(prho_in_out_r_g, "not set");
    FPMSTR(pL_out, "not set");
}

// The reflection tables have to be multiplied by 10^29 for storage purposes,
// but the distance ratio (to be squared) is divided by 10^8: 29-2*8=13
Dnorm = 1e13;// * pow(distratio, 2.);
// For the primary we assume nH = 10^15 cm-3 in reflection tables
// and no storage factors, then 15-2*8 = -1
Dnormprim = 1e-1;// * pow(distratio, 2.);

// Let's normalize reflected spectra for a distant observer
// and convert to flux
if (xival >= 5.){
    for(ie = 0; ie < ne; ie++){
      far[ie] *= Dnorm;
      if (stokes) {
        qar[ie] *= Dnorm;
        uar[ie] *= Dnorm;
        var[ie] *= Dnorm;
      }
    }
}else if (xival >= 0. && xival < 5.){
    for(ie = 0; ie < ne; ie++){
      // neutral reflection tables have a forgotten norm factor 76.394372684,
      // which is 1 / (\Delta\mu_e * \Delta\Phi_e)
      far[ie] *= Dnorm * xival * 76.394372684;
      if (stokes) {
        qar[ie] *= Dnorm * xival * 76.394372684;
        uar[ie] *= Dnorm * xival * 76.394372684;
        var[ie] *= Dnorm * xival * 76.394372684;
      }
    }
}else{
      for(ie = 0; ie < ne; ie++){
      far[ie] *= -Dnorm * xival / 5e4;
      if (stokes) {
        qar[ie] *= -Dnorm * xival / 5e4;
        uar[ie] *= -Dnorm * xival / 5e4;
        var[ie] *= -Dnorm * xival / 5e4;
      }
    }
}


// Let's add primary flux to the solution and redshift it
refl_ratio=-1.;
if ( true_Theta > x0 && NpNr != 0 ) {

  if (pol_deg == -1.){
      a0 = 1.73*mue_tot*(1.+5.3*mue_tot-0.2*mue_tot*mue_tot)/(1.+1.3*mue_tot+4.4*mue_tot*mue_tot);
  }else{
      a0 = 1.;
  }

  zzshift = 1.0 / (1.0 + redshift);

  if (zzshift <= 0.0) {
      fprintf(stderr, "Invalid redshift: z = %.12g gives 1/(1+z) <= 0.\n", redshift);
      return -1;
  }

  if (gam != 2.0) {
      norm_energy = (pow(EC, 2.0 - gam) - pow(E0, 2.0 - gam)) / (2.0 - gam);
  } else {
      norm_energy = log(EC / E0);
  }

  for (ie = 0; ie < ne; ie++) {
      El = fmax(ear[ie],     zzshift * E0);
      Eh = fmin(ear[ie + 1], zzshift * EC);

      if (Eh <= El) {
          photar1[ie] = 0.0;
      } else if (gam != 1.0) {
          photar1[ie] =
              pow(zzshift, gam) *
              (pow(Eh, 1.0 - gam) - pow(El, 1.0 - gam)) /
              (1.0 - gam) / norm_energy;
      } else {
          photar1[ie] =
              pow(zzshift, gam) *
              log(Eh / El) / norm_energy;
      }
  }

  for ( ie = 0; ie < ne; ie++ ){
    if (xival >= 0){
      photar1[ie] *= NpNr * Dnormprim * xival * ERG * a0 / 4. / PI;
    }else{
      photar1[ie] *= - NpNr * Dnormprim * xival * ERG * a0 / 4. / PI;
    }
  }

  flux_refl = flux_prim = 0.0;

  if (stokes) {
      if (pol_deg == -1.0) {
          p0 = 0.064 * (1.0 - mue_tot) *
              (1.0 + 16.3 * mue_tot + 6.2 * mue_tot * mue_tot) /
              (1.0 + 8.2 * mue_tot - 2.1 * mue_tot * mue_tot);
          chi0 = 0.0;
      } else {
          p0 = pol_deg;
          chi0 = chi;
      }
  }

  for (ie = 0; ie < ne; ie++) {
      Emid = 0.5 * (ear[ie] + ear[ie + 1]);

        if (photar1[ie] != 0.0) {
            flux_refl += far[ie] * Emid;
            flux_prim += photar1[ie] * Emid;
        }

      far[ie] += photar1[ie];

      if (stokes && photar1[ie] != 0.0) {
          qar[ie] += photar1[ie] * p0 * cos(2.0 * chi0);
          uar[ie] += photar1[ie] * p0 * sin(2.0 * chi0);
      }
  }

  if (flux_prim > 0.0) {
      refl_ratio = flux_refl / flux_prim;
  } else {
      refl_ratio = -1.0;
  }
}
snprintf(Refl, sizeof(Refl), "%e", refl_ratio);
FPMSTR(pRefl, Refl);




/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
fw = fopen("parameters.txt", "w");
fprintf(fw, "PhoIndex        %12.6f\n", param[0]);
fprintf(fw, "cos_incl     %12.6f\n", param[1]);
fprintf(fw, "trTheta         %12.6f\n", param[2]);
fprintf(fw, "rhorhoin   %12.6f\n", param[3]);
fprintf(fw, "xi0          %12.6f\n", param[4]);
fprintf(fw, "beta         %12.6f\n", param[5]);
fprintf(fw, "NpNr          %12.6f\n", param[6]);
fprintf(fw, "pol_deg        %12.6f\n", param[7]);
fprintf(fw, "chi         %12.6f\n", param[8]);
fprintf(fw, "pos_ang        %12.6f\n", param[9]);
fprintf(fw, "zshift      %12.6f\n", param[10]);
fprintf(fw, "Stokes      %12d\n", (int) param[11]);
fprintf(fw, "true_Theta      %12.6f\n", true_Theta);
fprintf(fw, "inc_degrees      %12.6f\n", x0);
fprintf(fw, "rho_in      %12.6f\n", rho_in_pc);
fprintf(fw, "rho_in_r_g      %12.6f\n", rho_in_r_g);
fprintf(fw, "n_H0      %12.6f\n", n_H0);
fprintf(fw, "refl_ratio      %12.6f\n", refl_ratio);
fclose(fw);
#endif
/******************************************************************************/

// interface with XSPEC
if (!stokes) for (ie = 0; ie < ne; ie++) photar[ie] = far[ie];
else {
//   let's change the orientation of the system 
  if(pos_ang != 0.)
    for( ie=0; ie<ne; ie++ ){
      qar_final[ie] = qar[ie] * cos(2*pos_ang) - uar[ie] * sin(2*pos_ang);
     uar_final[ie] = uar[ie] * cos(2*pos_ang) + qar[ie] * sin(2*pos_ang);
    }
   else
  for( ie=0; ie<ne; ie++ ){
      qar_final[ie] = qar[ie];
      uar_final[ie] = uar[ie];
    }
  pamin = 1e30;
  pamax = -1e30;
  pa2min = 1e30;
  pa2max = -1e30;
  for (ie = ne - 1; ie >= 0; ie--) {
    pd[ie] = sqrt(qar_final[ie] * qar_final[ie] + uar_final[ie] * uar_final[ie] 
                  + var[ie] * var[ie]) / (far[ie] + 1e-99);
    pa[ie] = 0.5 * atan2(uar_final[ie], qar_final[ie]) / PI * 180.;
    if (ie < (ne - 1)) {
      while ((pa[ie] - pa[ie + 1]) > 90.) pa[ie] -= 180.;
      while ((pa[ie + 1] - pa[ie]) > 90.) pa[ie] += 180.;
    }
    if (pa[ie] < pamin) pamin = pa[ie];
    if (pa[ie] > pamax) pamax = pa[ie];
    pa2[ie] = 0.5 * asin(var[ie] / sqrt(qar_final[ie] * qar_final[ie] 
                         + uar_final[ie] * uar_final[ie] + var[ie] * var[ie] 
                         + 1e-99)) / PI * 180.;
    if (ie < (ne - 1)) {
      while ((pa2[ie] - pa2[ie + 1]) > 90.) pa2[ie] -= 180.;
      while ((pa2[ie + 1] - pa2[ie]) > 90.) pa2[ie] += 180.;
    }
    if (pa2[ie] < pa2min) pa2min = pa2[ie];
    if (pa2[ie] > pa2max) pa2max = pa2[ie];
  }
  #ifdef OUTSIDE_XSPEC
  fw = fopen("stokes.dat", "w");
  #endif
  for (ie = 0; ie < ne; ie++) {
    if ((pamax + pamin) > 180.) pa[ie] -= 180.;
    if ((pamax + pamin) < -180.) pa[ie] += 180.;
    if ((pa2max + pa2min) > 180.) pa2[ie] -= 180.;
    if ((pa2max + pa2min) < -180.) pa2[ie] += 180.;
    #ifdef OUTSIDE_XSPEC
    fprintf(fw,
      "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
      0.5 * (ear[ie] + ear[ie+1]), far[ie] / (ear[ie+1] - ear[ie]), 
      qar_final[ie] / (ear[ie+1] - ear[ie]), 
      uar_final[ie] / (ear[ie+1] - ear[ie]), 
      var[ie] / (ear[ie+1] - ear[ie]), pd[ie], pa[ie], pa2[ie]);
  #endif
//interface with XSPEC..........................................................
    if (stokes ==  1) photar[ie] = far[ie];
    if (stokes ==  2) photar[ie] = qar_final[ie];
    if (stokes ==  3) photar[ie] = uar_final[ie];
    if (stokes ==  4) photar[ie] = var[ie];
    if (stokes ==  5) photar[ie] = pd[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  6) photar[ie] = pa[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  7) photar[ie] = pa2[ie] * (ear[ie + 1] - ear[ie]);
    if (stokes ==  8) photar[ie] = qar_final[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
    if (stokes ==  9) photar[ie] = uar_final[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
    if (stokes == 10) photar[ie] = var[ie] / (far[ie]+1e-99) * (ear[ie + 1] - ear[ie]);
  }
  #ifdef OUTSIDE_XSPEC
  fclose(fw);
  #endif
}


return 0;
}
