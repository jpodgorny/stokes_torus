/* xsstokes - polarized reflection from an axially symmetric surface illuminated 
 *	    isotropically by an (un)polarised power law,
 *        - an example of a polarisation subroutine for XSPEC using tables 
 *          computed with torus_integrator.py code
 * 
 *
 * par1 ... PhoIndex - power-law photon index of the primary flux
 * par2 ... cos_incl - cosine of the observer inclination (1.-pole, 0.-disc)
 * par3 ... trTheta - transformed torus half-opening angle between 0 and 1
 * par4 ... B - taking into account region below equator in integration
 *			= 1  - yes
 *			= 0  - no
 * par5 ... poldeg  - intrinsic polarisation degree of primary radiation
 * par6 ... chi - intrinsic polarisation angle (in degrees, -90 < chi < 90)
 *		       of primary radiation, the orientation is degenarate by 
 *		       180 degrees
 * par7 ... pos_ang - orientation of the system (-90 < pos_ang < 90), 
 *                    the position angle (in degrees) of the system 
 *                    rotation axis with direction up,
 *                    the orientation is degenarate by 180 degrees
 * par8 ... zshift  - overall Doppler shift
 * par9 ... Stokes  - what should be stored in photar() array, i.e. as output
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
 * par10 ... norm - if the distance to the source D [pc] is known, this
 *                  parameter allows to check consistency of the best-fit value
 *                  with the expected illuminating total flux from the central
 *                  source in the equatorial plane at the distance of 1 pc from
 *                  the center as F_in [erg / cm^2 / s] = 3.98*10^(-15)*D^2*norm,
 *                  which can be also estimated e.g. with a simultaneously used
 *                  comptonization model
 *
 ******************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
*******************************************************************************/
#ifdef OUTSIDE_XSPEC

#define IFL    1
#define NPARAM 5
#define NE     200
#define E_MIN  1.
#define E_MAX  100.

int main() {

int stokesnitorus(const double *ear, int ne, const double *param, int ifl,
               double *photar, double *photer, const char* init);

double ear[NE+1], photar[NE], photer[NE], param[NPARAM];
char   initstr[0] = "";
int    ie;

param[ 0] = 2.0;        // PhoIndex
param[ 1] = 0.775;      // cos_incl
param[ 2] = 0.5;        // trTheta
param[ 3] = 1.;         // B
param[ 4] = 0.;         // poldeg
param[ 5] = 0.;         // chi
param[ 6] = 0.;         // pos_ang
param[ 7] = 0.;         // zshift
param[ 8] = 1.;         // Stokes

for(ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

stokesnitorus(ear, NE, param, IFL, photar, photer, initstr);
return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/

#define REFSPECTRA1 "stokes-neutral-iso-UNPOL-torus.fits" // UNPOLARISED
#define REFSPECTRA2 "stokes-neutral-iso-HRPOL-torus.fits" // HORIZONTALLY POLARISED
#define REFSPECTRA3 "stokes-neutral-iso-45DEG-torus.fits" // DIAGONALLY POLARISED
#define VISIBILITY_FILE "visibility_line.txt" // Text file containing Theta_limit
#define MAX_SIZE 40000  // Size of the visibility line text file

#define PI   3.14159265358979
#define NPAR 5

extern int    xs_write(char* wrtstr, int idest);
extern float  DGFILT(int ifl, const char* key);
extern void   FPMSTR(const char* value1, const char* value2);
extern char*  FGMSTR(char* dname);
extern void   tabintxflt(float* ear, int ne, float* param, const int npar, 
                         const char* filenm, const char **xfltname, 
                         const float *xfltvalue, const int nxflt,
                         const char* tabtyp, float* photar, float* photer);

int stokesnitorus(const double *ear, int ne, const double *param, int ifl,
            double *photar, double *photer, const char* init) {


FILE   *fw;
//static int first = 1;
int status = 0;
static char   xsdir[255]="";
static char   pname[128]="XSDIR", ptrue_Theta_out[128] = "true_Theta", pinc_degrees[128] = "inc_degrees";
static char refspectra[3][255], visibility_path[255];

// - if set try XSDIR directory, otherwise look in the working directory
// Initialize refspectra elements and visibility file path
sprintf(xsdir, "%s", FGMSTR(pname));
if (strlen(xsdir) == 0) {
    strcpy(refspectra[0], REFSPECTRA1);
    strcpy(refspectra[1], REFSPECTRA2);
    strcpy(refspectra[2], REFSPECTRA3);
    strcpy(visibility_path, VISIBILITY_FILE);
} else {
    if (xsdir[strlen(xsdir) - 1] == '/') {
        sprintf(refspectra[0], "%s%s", xsdir, REFSPECTRA1);
        sprintf(refspectra[1], "%s%s", xsdir, REFSPECTRA2);
        sprintf(refspectra[2], "%s%s", xsdir, REFSPECTRA3);
        sprintf(visibility_path, "%s%s", xsdir, VISIBILITY_FILE);
    } else {
        sprintf(refspectra[0], "%s/%s", xsdir, REFSPECTRA1);
        sprintf(refspectra[1], "%s/%s", xsdir, REFSPECTRA2);
        sprintf(refspectra[2], "%s/%s", xsdir, REFSPECTRA3);
        sprintf(visibility_path, "%s/%s", xsdir, VISIBILITY_FILE);
    }
}
    
int    i, j, ie, stokes;
double pol_deg, pos_ang, chi;
const char*   xfltname = "Stokes";
float  xfltvalue;
float  Smatrix[9][ne];
float  fl_param[NPAR]={(float) param[0], (float) param[1], (float) param[2], (float) param[3], (float) param[7]};
const char*  tabtyp="add";
float  fl_ear[ne+1], fl_photer[ne];
double far[ne], qar[ne], uar[ne], var[ne], pd[ne], pa[ne], pa2[ne], 
       qar_final[ne], uar_final[ne];
double pamin, pamax, pa2min, pa2max;

FILE *file;
double vis_Theta[MAX_SIZE], vis_inc[MAX_SIZE];
double x0, y, x1, y1, x2, y2;
int count = 0, k;
double trTheta, mue_tot, Theta_limit, min_Theta, max_Theta, true_Theta;
char true_Theta_out[32], inc_degrees[32];

mue_tot = (float) param[1];
trTheta = (float) param[2];
pol_deg = param[4];
chi = param[5]/180.*PI;
pos_ang = param[6]/180.*PI;
stokes = (int) param[8];
if(stokes == -1){
  xfltvalue = DGFILT(ifl, xfltname);
  if (xfltvalue == 0. || xfltvalue == 1. || xfltvalue == 2.){
    stokes = 1 + (int) xfltvalue;
  }
  else {
    xs_write("stokes: no or wrong information on data type (counts, q, u)", 5);
    xs_write("stokes: stokes = par8 = 0 (i.e. counts) will be used", 5);
    stokes=0;
  }
}

//Let's read and interpolate the FITS tables that we will need using internal
//XSPEC routine tabintxflt
//Note that we do not use errors here
for(ie = 0; ie <= ne; ie++) fl_ear[ie] = (float) ear[ie];
if(stokes){//we use polarised tables
//  if (first) {
// The status parameter must always be initialized.
  status = 0;

  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 2; j++){
      xfltvalue = (float) j;
      tabintxflt(fl_ear, ne, fl_param, NPAR, refspectra[i], &xfltname, 
                 &xfltvalue, 1, tabtyp, Smatrix[i*3+j], fl_photer);  
      }

//  first = 0;
//  }


  //Let's perform the transformation to initial primary polarisation degree and angle
  for(ie = 0; ie < ne; ie++) {
    for(j=0; j<=2; j++) Smatrix[j+3][ie] -= Smatrix[j][ie];
    for(j=0; j<=2; j++) Smatrix[j+6][ie] -= Smatrix[j][ie];

    far[ie] = Smatrix[0][ie] +
                        pol_deg * ( -Smatrix[3][ie] * cos(2.*(chi)) +
                                    Smatrix[6][ie] * sin(2.*(chi)) );
    qar[ie] = Smatrix[1][ie] +
                        pol_deg * ( -Smatrix[4][ie] * cos(2.*(chi))+
                                    Smatrix[7][ie] * sin(2.*(chi)) );
    uar[ie] = Smatrix[2][ie] +
                        pol_deg * ( -Smatrix[5][ie] * cos(2.*(chi))+
                                    Smatrix[8][ie] * sin(2.*(chi)) );
    var[ie] = 0.;   
                              
    // far[ie] = ( 1. + pol_deg ) * Smatrix[0][ie] - pol_deg * Smatrix[3][ie];
    // qar[ie] = ( 1. + pol_deg ) * Smatrix[1][ie] - pol_deg * Smatrix[4][ie];
    // uar[ie] = ( 1. + pol_deg ) * Smatrix[2][ie] - pol_deg * Smatrix[5][ie];
    // var[ie] = 0.;   
  }
}else{//we just use unpolarised counts
//  if (first) {

  // The status parameter must always be initialized.
  status = 0;

  xfltvalue = 0.;
  tabintxflt(fl_ear, ne, fl_param, NPAR, refspectra[0], &xfltname, &xfltvalue, 
             1, tabtyp, Smatrix[0], fl_photer);  

//  first = 0;
//  }

  for(ie = 0; ie < ne; ie++){
    far[ie] = Smatrix[0][ie];
  }
}


// let's transform trTheta back to true_theta



// The status parameter must always be initialized.
status = 0;


file = fopen(visibility_path, "rt");
if (file == NULL) {
perror("Error opening file");
return -1;
}

while (fscanf(file, "%lf %lf", &vis_Theta[count], &vis_inc[count]) == 2) {
count++;
}
fclose(file);



x0 = acos(mue_tot) / PI * 180.0;  // Interpolation point

// Inline linear interpolation
for (i = 0; i < count - 1; i++) {
x1 = vis_inc[i];
y1 = vis_Theta[i];
x2 = vis_inc[i + 1];
y2 = vis_Theta[i + 1];
if (x0 >= x1 && x0 <= x2) {
    Theta_limit = y1 + (x0 - x1) * (y2 - y1) / (x2 - x1);
    break;
}
}

// Ensure Theta_limit is set; you may need additional logic for extrapolation.
if (i == count - 1) {
fprintf(stderr, "Interpolation point is out of range.\n");
return -1;
}

min_Theta = fmax(25.0, Theta_limit);
max_Theta = 90.0;
true_Theta = trTheta * (max_Theta - min_Theta) + min_Theta;

sprintf(true_Theta_out, "%12.6f", true_Theta);
FPMSTR(ptrue_Theta_out, true_Theta_out);

sprintf(inc_degrees, "%12.6f", x0);
FPMSTR(pinc_degrees, inc_degrees);

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
fw = fopen("parameters.txt", "w");
fprintf(fw, "PhoIndex        %12.6f\n", param[0]);
fprintf(fw, "cos_incl     %12.6f\n", param[1]);
fprintf(fw, "trTheta         %12.6f\n", param[2]);
fprintf(fw, "B          %12.6f\n", param[3]);
fprintf(fw, "poldeg        %12.6f\n", param[4]);
fprintf(fw, "chi         %12.6f\n", param[5]);
fprintf(fw, "pos_ang        %12.6f\n", param[6]);
fprintf(fw, "zshift      %12.6f\n", param[7]);
fprintf(fw, "Stokes      %12d\n", (int) param[8]);
fprintf(fw, "true_Theta      %12.6f\n", true_Theta);
fprintf(fw, "inc_degrees      %12.6f\n", x0);
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
  fw = fopen("stokes.dat", "w");
  for (ie = 0; ie < ne; ie++) {
    if ((pamax + pamin) > 180.) pa[ie] -= 180.;
    if ((pamax + pamin) < -180.) pa[ie] += 180.;
    if ((pa2max + pa2min) > 180.) pa2[ie] -= 180.;
    if ((pa2max + pa2min) < -180.) pa2[ie] += 180.;
    fprintf(fw,
      "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
      0.5 * (ear[ie] + ear[ie+1]), far[ie] / (ear[ie+1] - ear[ie]), 
      qar_final[ie] / (ear[ie+1] - ear[ie]), 
      uar_final[ie] / (ear[ie+1] - ear[ie]), 
      var[ie] / (ear[ie+1] - ear[ie]), pd[ie], pa[ie], pa2[ie]);
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
  fclose(fw);
}


return 0;
}
