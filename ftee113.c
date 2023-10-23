/* Program ftee v. 1.13
				written 11 - 16 May, 2010
				UPDATES:
				corrected 09 May     2012 - RYAN - This corrected the volume (and therefore eU) calculation for the ellipsoid
				updated 8 January 2013 - added code to allow ftee to be run from any location in the user PATH
				v. 1.03(1)  9-30-2013 bug in reading input file: was reading sample[i] not &sample[i] address!!! Somehow it worked more or less...
				       Aargh! _And_ several other variables not referenced as an address!! 
				v. 103(2) 9-30-2013 
				       BUT.... the problem was an extra input field in the first read from file (leftover from old version??).
				       Who knows how the code worked with that error, plus the variable/address errors.
				v. 103(3) 10-1-2013. Doh. Error was in incorrect# of read parameters. Ignore all the b.s. I just wrote about
					   reading the variables wrong: since the strings are arrays of char, they are pointers and so we read them without &
					   Cleaned up warnings displayed using -Wall compiler flag.
				v. 1100(b) 27 Jan 2014. Tried to implement openMP to speed things up a bit
				     seems to be working! about 2X improvement with two cores. only test was no openMP sample data, and test against ftee1033
				     using that data.
				v. 1100 14 April, 2015. Accepted openMP tweak. Added calculation of ppm values for U, Th, and Sm, plus Ft-equivalent radius for bulk Ft
				v. 1100 18 November, 2015.			     
				     Finished version 1.10 - got U, Th, Sm ppm working as well as Ft-equivalent radius. Needs testing, folks.
				     Modified reporting to also give condensed table of just sample info
				     Note that (obviously) ppm values obtained using calculated grain volume and mineral density
				v. 1110 3 December 2019:
					 Updated eU calculation after Cooperdock et al. (2019) to be U + .238Th + 0.0012Sm .   Changes are coefficient from 
					 0.235 to 0.238 and addition of Sm term (older Ftee code ignored Sm).
					 NOT TESTED!!!!  
				v. 112 4 February 2020:
					 Now reporting volume and mass of grain as well as density used. Also prophylactic bug fix in that
					 these new reported calculations should work even if user has mixed minerals in their input file.
					 (needed to put sample density into an array).
					 Check for goodness just with a few old files from older versions. NOT extensively tested.
				v. 113 September 2023:
					 In prep for moving this to GitHub, checked Apple Silicon compilation and changed
					 output file suffix to "txt" from "xls", which doesn't really work well. Also fixed
					 small reporting bug in printing out number of cores found.

	IMPORTANT: This code It gives the same value as concrete examples from Hourigan et al., and checks out against 
	           the analytical solution for a sphere. However, it has NOT been _exhaustively_ checked against all combinations of inputs.
	           
	           
	AUTHOR:	Dr. Peter K. Zeitler
			Department of Earth and Environmental Sciences
			1 W. Packer Avenue
			Lehigh University
			Bethlehem, Pennsylvania  18015-3188	U.S.A.
			phone: (610) 758-3671  fax: (610) 758-3677
			internet: peter.zeitler@lehigh.edu

------------------------------------------------------------------------------------------------------
    PURPOSE:  Calculate alpha-loss correction for U-Th-Sm/He dating, using Monte-Carlo model for several
    		  symmetrical geometries for apatite, zircon, titanite/sphene, monazite and rutile. The
    		  geometries are: ellipsoid (includes sphere); cylinder and flattened cylinder (pinacoidal only);
    		  and tetragonal prism (pinacoidal or with optional pyramidal terminations; includes cube).
    		  
			  See the manual for details about the code, and using it.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include  <stdbool.h>
#include  <string.h>
#include  <unistd.h>
#include <omp.h>
#include <locale.h>

#define	pi        3.14159265

// following block is for random-number-generator routines
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

// function prototypes
void init_genrand(unsigned long s);  // declare these for the ones we're using
double genrand_res53(void);

int in_ellipsoid(double z, double x, double y);
int in_cylinder(double z, double x, double y);
int in_tetragonal(double z, double x, double y);
double rectx (double z, double diagonal_length);
double recty (double z, double diagonal_length);

// global variables (declare these outside of main() to be able to share them with subroutines -- crude, I know).
double length[501], width1[501], width2[501], tip1[501], tip2[501]; // storage for input data
double z_end, z_start, z_end1, z_end2;
double rlength, rwidth1, rwidth2;
double zposn, xposn, yposn, xposn_n, yposn_n, zposn_n;
double xspan, yspan, zspan;
double diag_angle, diag;
long isample;
long size;

int main(int argc, char *argv[])
{
// i-o related things, including files, handling argv[] stuff from command line, and more input data
	FILE *ofp, *ifp;
	size_t ibytes;
	char *infile = NULL;
	char *outfile = NULL;
	char *suffix = ".txt";
	char sample[501][50],mineral[501][20],calibration[501][20],geometry[501][20];
	long ngrains[501];
	long report_flag[501];  // for use in reporting whole-sample output for multi-grain samples
	long grain_count; // for use in cleaning up the reporting just sample data, not grain data

// timing things for assessing openMP
	clock_t started, ended;
	double par_started, par_ended;
	double elapsed;
	long nthreads;

// working variables
	double ft_group[501],ft_sample[501],eu_sample[501], uval[501], thval[501],smval[501];
	double U_ppm[501],Th_ppm[501],Sm_ppm[501],Req[501],Req_group[501];
	double volume[501], moles238[501],moles235[501],moles232[501],moles147[501];	
	double stopping238,stopping235,stopping232,stopping147;
	double heprod238, heprod235, heprod232, heprod147, heprod_total;
	double density,sample_density[501];

	double sampling_density; // use to report mean density of sampling that monte-carlo routine does per micron

	long i, ejections;
	long nsamples;
	long check;

	double inclination, azimuth;
	
	double ftee238, ftee235, ftee232, ftee147, ftee;
	double tot_volume;
	double eff, effprime, Ress, Ress238, Ress235, Ress232, Ress147;  // for calculating Ft-equivalent radius

// declare following variables as long long to be sure that we don't hit max_int if someone chooses a really large number
// of monte-carlo iterations

	long long decays_total, decays238_retained, decays235_retained, decays232_retained, decays147_retained, ipoint, imontes;
	long long expected_decays;

// to help us find ourselves...
	char path[255];

	getcwd(path,(size_t)size);
	chdir(path);
	
// ***** Start up 

printf("\n");
printf("+++++++++++++++++++++++++++  PROGRAM ftee 1.13 +++++++++++++++++++++++++++\n");
printf("\n");

// ***** Handle command-line argument (expect one, the input file name) 	
if(argc == 1) // there is only one value (which is name of this program), and we need two
{
		puts("ftee 1.13:");
		puts("USAGE: [./]ftee inputfilename\n");
		exit(1);
}
else  // assume that the second argument (argv[1]) is the name. Strip off characters after delimiter to build output filename of type .txt
{
	ibytes=strcspn(argv[1],"."); // get length of string before dot-delimiter
	infile = malloc(strlen(argv[1]) + 1); // allocate memory for inermediate buffer (will be variable hence malloc)	 
	if (!infile)   //make sure memory allocation worked
	{
		printf("failed allocating memory for input file... exiting\n");
		exit(0);	
	}

	strncpy(infile,argv[1],ibytes);   // write characters before delimiter to buffer

	outfile = malloc(strlen(infile) + strlen(suffix) + 1); // allocate memory for output file name (will be variable hence malloc)	 
	if (!outfile)   //make sure memory allocation worked
	{
		printf("failed allocating memory for input file... exiting\n");
		exit(0);	
	}

strcpy(outfile, infile);  // write the filename from character buffer to output string

strcat(outfile,suffix);  // append the suffix to build the full output filename

// ***** Get inputs from file

	ifp = fopen(argv[1], "r"); //open input file

	if (!ifp)  // give up if ifp is NULL - this means file not found
	{
		printf("unable to open input file... exiting\n");
		exit(0);
	}
	fscanf(ifp,"%ld%lld%ld",&nsamples,&imontes,&ejections);

	assert((nsamples >= 1)&&(nsamples <= 500)); // gotta have one sample at least, cannot have more than 500 because of array declaration (could be increased)
	assert((imontes >= 1)&&(imontes <= 1000000000)); // gotta have one carlo iteration (assume), get suspicious if user wants more than a billion
	assert((ejections >= 1)&&(ejections <= 10000)); // gotta have one ejection but get nervous about a request for too many
    setlocale(LC_ALL, "");
	printf("<<<< INPUTS >>>>\n\n");
	printf("  nsamples: %ld   Monte-Carlo iterations: %'lld   Ejections at each site: %ld\n\n",nsamples,imontes,ejections);

	for (i = 1; i <= nsamples; i++)
	{
		fscanf(ifp, "%s%ld%lf%lf%lf%s%s%s%lf%lf%lf%lf%lf", sample[i],&ngrains[i],&uval[i],&thval[i],&smval[i],mineral[i],calibration[i],geometry[i],&length[i],&width1[i],&width2[i],&tip1[i],&tip2[i]);
		assert((uval[i] > 0.0) || (thval[i] > 0.0));
		assert ((length[i] > 0.0) && (width1[i] > 0.0) && (width2[i] > 0.0));
	}
	fclose(ifp);
}

// ***** Echo inputs from file back to screen
for (i = 1; i <= nsamples; i++)
{
	printf("  %ld  %s  %ld %7.3f %7.3f %7.3f %s %s %s %6.2f %6.2f %6.2f %5.2f %5.2f\n",i,sample[i],ngrains[i],uval[i],thval[i],smval[i],mineral[i],calibration[i],geometry[i],length[i],width1[i],width2[i],tip1[i],tip2[i]);
}
printf("\n");
	if (imontes < 100000)
	{
		printf("WARNING: rather meager number of monte-carlo iterations requested - are you sure??\n");
		printf("         please check your inputs and type 1 to continue, 0 to stop: ");
		scanf("%ld",&check);
		if (check != 1) exit(0);
		printf("\n");
	}
	if (ejections > 500)
	{
		printf("WARNING: more than 500 ejections per Monte-Carlo node. This code may run very slowly or \n");
		printf("         you may be undersampling, using too few monte-carlo trials taken as a tradeoff - are you sure??\n");
		printf("         please check your inputs and type 1 to continue, 0 to stop: ");
		scanf("%ld",&check);
		if (check != 1) exit(0);
		printf("\n");
	}

// ***** OK, we now have what we hope is good input in the relevant arrays. Start the main loop through individual samples

// first set up timing for overall run through multiple samples	
		par_started = omp_get_wtime();
		started = clock();

// initialize Nishimura-Matsumoto random-number routines
		init_genrand(-abs(time(NULL)));	
			
printf("<<<< WORKING >>>>\n\n");
for (isample=1;isample <= nsamples;isample++)  // MAIN SAMPLE PROCESSING LOOP!!

// ***** First, assign the stopping values. Note that in theory we would want to parse the input and complain if a user switches calibration
// or mineral type _within_ the processing of grains from one sample, but for now, we just let the person bungle if they screw up.
// Actually first thing we do is make sure we find a mineral type we understand and can handle
// and we also check to make sure we have a calibration type and a geometry we can handle

{
	if ((strncasecmp(mineral[isample], "apatite",20) != 0) && (strncasecmp(mineral[isample], "zircon",20) != 0) && (strncasecmp(mineral[isample], "titanite",20) != 0) && (strncasecmp(mineral[isample], "sphene",20) != 0) && (strncasecmp(mineral[isample], "monazite",20) != 0) && (strncasecmp(mineral[isample], "rutile",20) != 0))
		{
			printf("mineral either not recognized or not yet implemented -- ending...\n");
			exit(1);	
		}
	if ((strncasecmp(calibration[isample], "farley",20) != 0) && (strncasecmp(calibration[isample], "ketcham",20) != 0) && (strncasecmp(calibration[isample], "oldfarley",20) != 0))
		{
			printf("calibration not recognized -- ending...\n");
			exit(1);	
		}
	if ((strncasecmp(geometry[isample], "ellipsoid",20) != 0) && (strncasecmp(geometry[isample], "cylinder",20) != 0) && (strncasecmp(geometry[isample], "tetragonal",20) != 0))
		{
			printf("geometry either not recognized or not yet implemented -- ending...\n");
			exit(1);	
		}
	if (strncasecmp(mineral[isample], "apatite",20) == 0)
	{
		density = 3.19;
		if (strncasecmp(calibration[isample],"farley",20) == 0)
		{	

			stopping238 = 19.3;
			stopping235 = 22.4;
			stopping232 = 22.1;
			stopping147 = 5.93;   // usually not done back in the day; so either set Sm to zero or use this new value (difference will be miniscule)
		}
		else
		{
			stopping238 = 18.81;
			stopping235 = 21.80;
			stopping232 = 22.25;
			stopping147 = 5.93;
		}
		if (strncasecmp(calibration[isample],"oldfarley",20) == 0) // bury old values in here for debugging, comparison, and historical purposes
		{	
			stopping238 = 19.68;
			stopping235 = 22.83;
			stopping232 = 22.46;
		}
	}
	if (strncasecmp(mineral[isample], "zircon",20) == 0)
	{
		density = 4.65;
		if (strncasecmp(calibration[isample],"farley",20) == 0)
		{
			stopping238 = 16.97;
			stopping235 = 19.64;
			stopping232 = 19.32;
			stopping147 = 0.0;
		}
		else
		{
			stopping238 = 15.55;
			stopping235 = 18.05;
			stopping232 = 18.43;
			stopping147 = 4.76;
		}
	}
	if (strncasecmp(mineral[isample], "titanite",20) == 0)
	{
		density = 3.53;
		if (strncasecmp(calibration[isample],"farley",20) == 0)
		{
			stopping238 = 18.12;
			stopping235 = 21.01;
			stopping232 = 20.68;
			stopping147 = 0.0;
		}
		else
		{
			stopping238 = 17.46;
			stopping235 = 20.25;
			stopping232 = 20.68;
			stopping147 = 5.47;
		}
	}
	if (strncasecmp(mineral[isample], "sphene",20) == 0)
	{
		density = 3.53;
		if (strncasecmp(calibration[isample],"farley",20) == 0)
		{
			stopping238 = 18.12;
			stopping235 = 21.01;
			stopping232 = 20.68;
			stopping147 = 0.0;
		}
		else
		{
			stopping238 = 17.46;
			stopping235 = 20.25;
			stopping232 = 20.68;
			stopping147 = 5.47;
		}
	}
	if (strncasecmp(mineral[isample], "monazite",20) == 0)
	{
		density = 5.26;
		if (strncasecmp(calibration[isample],"farley",20) == 0)
		{
			stopping238 = 16.6;
			stopping235 = 19.3;
			stopping232 = 19.0;
			stopping147 = 0.0;
		}
		else
		{
			stopping238 = 16.18;
			stopping235 = 18.74;
			stopping232 = 19.11;
			stopping147 = 4.98;
		}
	}
	if (strncasecmp(mineral[isample], "rutile",20) == 0)
	{
		density = 4.25;
		if (strncasecmp(calibration[isample],"farley",20) == 0) //just use Ketcham values for now
		{
			stopping238 = 15.30;
			stopping235 = 17.76;
			stopping232 = 18.14;
			stopping147 = 4.77;
		}
		else
		{
			stopping238 = 15.30;
			stopping235 = 17.76;
			stopping232 = 18.14;
			stopping147 = 4.77;
		}
	}
	
// Keep record of density for sample (in case user does run with different minerals!)
	sample_density[isample] = density;

// ***** Now set up the coordinates and problem space with parameters neede do describe grain

// get grain radii for use in geometry-based routines
	rlength = length[isample]/2.0;
	rwidth1 = width1[isample]/2.0;
	rwidth2 = width2[isample]/2.0;
		
// determine the full dimensions of the problem space
// in this non-gridded version, we use the grain's maximum limits in x, y, and z as the size of the space (we dont need to leave room for grid elements to handle alpha vector)
	zspan = length[isample];
	xspan = width1[isample];
	yspan = width2[isample];

// get U,Th, and Sm amounts in moles for later use in determining production values for Ft weighting

	moles238[isample] = uval[isample]/238.029/1e9*137.88/138.88;
	moles235[isample] = uval[isample]/238.029/1e9*1/138.88; 
	moles232[isample] = thval[isample]/232.039/1e9;	
	moles147[isample] = smval[isample]/150.36/1e9*0.1499;
	
// get He production to use in weighting component Ft's for grain
// Note that for multigrain samples all we can use as an estimate are the bulk-sample measurements

	heprod238 = 8*1.551*moles238[isample];  // for these weights we just need relative production so use lamda*1e10
	heprod235 = 7*9.849*moles235[isample];
	heprod232 = 6*0.4948*moles232[isample];
	heprod147 = 1*0.06539*moles147[isample];
	heprod_total = heprod238 + heprod235 + heprod232 + heprod147;

// ***** Finally, do the actual calculations. Each geometry has a separate set-up block and then its own monte-carlo routine (because I am too lazy to compact this)

// ********** START   TETRAGONAL PRISM WITH OPTIONAL PYRAMIDAL TIPS *********************************

	if (strncasecmp(geometry[isample], "tetragonal",20) == 0)
	{	
		// IMPORTANT: Because of possibly unequal pyramidal tips, this geometry does not use symmetry for the z-axis.
		
		// We break the problem into three z-zones: two possibly different pyramidal regions and the rectangular prism
		// We don't have to fret with complex trigonometry for the pyramidal zones because as we move along the z-axis, we simply
		// have a rectangular area whose bounds can be determined using two straight lines drawn from the crystal edges.
		
		z_start = 0.0;
		z_end1 = tip1[isample];
		z_end2 = length[isample] - tip2[isample];
		
		expected_decays = ((tip1[isample]+tip2[isample])/3 + length[isample]-tip1[isample] - tip2[isample])/length[isample] * imontes * ejections;  // ratio, volume of (pyramidal) tetragonal prism inscribed in confining retangular space (net formula!)

		diag = sqrt(rwidth1*rwidth1 + rwidth2*rwidth2);  // half-length of prism diagonal -- we use these to get x and y bounds on pyramid as a function of z
		diag_angle  = atan(width2[isample]/width1[isample]);  // angle diagonal makes with prism face 
		
// get volume for use in determining eU
		volume[isample] = 4*rwidth1*rwidth2*((length[isample] - tip1[isample] - tip2[isample]) + (tip1[isample] + tip2[isample])/3);

// ***** MAIN MONTE CARLO LOOP *****
	
		decays238_retained = decays235_retained = decays232_retained = decays147_retained = 0;
		decays_total = 0;
		
#pragma omp parallel private(zposn,xposn,yposn,zposn_n,xposn_n,yposn_n,i,azimuth,inclination)  num_threads(omp_get_num_procs())

#pragma omp for reduction(+:decays238_retained,decays235_retained,decays232_retained,decays147_retained,decays_total) nowait

		for (ipoint = 1; ipoint<=imontes; ipoint++)
		{
		nthreads = omp_get_num_threads();
//			zposn = zspan*((double)rand())/RAND_MAX;  // this statement// left in to show how to use rand(), which generates integer values and needs rescaling
			zposn = zspan * genrand_res53();
			xposn = xspan/2.0 - xspan * genrand_res53();
			yposn = yspan/2.0 - yspan * genrand_res53();

// now consider doing the alpha-ejection simulation; cheat a bit by using each selected position and direction of alpha-vector to calculate ejection(s) for all four nuclides
			if (in_tetragonal(zposn,xposn,yposn) == 1)  // call function that determines location; if we are in crystal or on surface, we do alpha-ejection(s)
			{
				for (i = 1; i <= ejections; i++)  // for each monte-carlo position, do ejections number of random alpha ejections (see notes at start of code)
				{
					decays_total = decays_total + 1;
	//				inclination = pi * genrand_res53(); // old (bad) randomization of just angles
					inclination = acos(1.0 - 2.0 * genrand_res53());
					azimuth = 2*pi * genrand_res53();
	
					xposn_n = xposn + stopping238 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping238 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping238 * cos(inclination);
					
					if (in_tetragonal(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays238_retained = decays238_retained + 1;
					}
					xposn_n = xposn + stopping235 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping235 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping235 * cos(inclination);
					if (in_tetragonal(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays235_retained = decays235_retained + 1;
					}
					xposn_n = xposn + stopping232 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping232 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping232 * cos(inclination);
					if (in_tetragonal(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays232_retained = decays232_retained + 1;
					}
					xposn_n = xposn + stopping147 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping147 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping147 * cos(inclination);
					if (in_tetragonal(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays147_retained = decays147_retained + 1;
					}
				} // end of ejection loop
			}
		} // end of monte-carlo loop		
	} // end of tetragonal-prism case

// ********** START   ELLIPSOID (INCLUDES SPHERE) *********************************
		
	if (strncasecmp(geometry[isample], "ellipsoid",20) == 0)
	{	
		// ellipsoid also handles case of sphere (all dimensions equal)
		volume[isample] = 0.66666667*pi*rwidth1*rwidth2*length[isample];
		expected_decays = 0.523598775598 * imontes * ejections;  // ratio, volume of ellipsoid inscribed in confining space
		
		z_start = 0.0;
		z_end1 = zspan;
		
	// *****MAIN MONTE CARLO LOOP *****
	
		decays238_retained = decays235_retained = decays232_retained = decays147_retained = 0;
		decays_total = 0;
		
#pragma omp parallel private(ipoint,zposn,xposn,yposn,zposn_n,xposn_n,yposn_n,i,azimuth,inclination)  num_threads(omp_get_num_procs())
#pragma omp for reduction(+:decays238_retained,decays235_retained,decays232_retained,decays147_retained,decays_total) schedule(dynamic,1) nowait

		for (ipoint = 1; ipoint<=imontes; ipoint++)
		{	
		nthreads = omp_get_num_threads();
			zposn = zspan/2.0 - zspan * genrand_res53(); // for ellipsoid we need z to range from -zmax to + zmax
			xposn = xspan/2.0 - xspan * genrand_res53(); // or we end up modeling only a hemisphere!
			yposn = yspan/2.0 - yspan * genrand_res53();	

			if (in_ellipsoid(zposn,xposn,yposn) == 1)  // if we are in crystal or on surface, do an alpha-ejection
			{
				for (i = 1; i <= ejections; i++)  // for each monte-carlo position, do ejections number of random alpha ejections (see notes at start of code)
				{
					decays_total = decays_total + 1;
	//				inclination = pi * genrand_res53(); // old (bad) randomization of just angles
					inclination = acos(1.0 - 2.0 * genrand_res53());
					azimuth = 2*pi * genrand_res53();
	
					xposn_n = xposn + stopping238 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping238 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping238 * cos(inclination);
					
					if (in_ellipsoid(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays238_retained = decays238_retained + 1;
					}
					xposn_n = xposn + stopping235 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping235 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping235 * cos(inclination);
					if (in_ellipsoid(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays235_retained = decays235_retained + 1;
					}
					xposn_n = xposn + stopping232 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping232 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping232 * cos(inclination);
					if (in_ellipsoid(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays232_retained = decays232_retained + 1;
					}
					xposn_n = xposn + stopping147 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping147 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping147 * cos(inclination);
					if (in_ellipsoid(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays147_retained = decays147_retained + 1;
					}
				} // end ejection loop
			}
		} // end of monte-carlo loop
	} // end of ellipsoid case

// ********** START   CYLINDER  *********************************
	
	if (strncasecmp(geometry[isample], "cylinder",20) == 0)
	{	
		// this simple cylinder case could be extended to handle tips following approach used with tetragonal geometry, or modified with some fussing to do hex prism with tips (in which case inputs will need to diverge from Farley measuring convention)
		volume[isample] = pi * rwidth1 * rwidth2 * length[isample];
		
		z_start = 0.0;
		z_end1 = zspan;
		
		expected_decays = 0.785398163397 * imontes * ejections;  // ratio, volume of (elliptical) cylinder inscribed in confining space

	// *****MAIN MONTE CARLO LOOP *****
	
		decays238_retained = decays235_retained = decays232_retained = decays147_retained = 0;
		decays_total = 0;
		
#pragma omp parallel private(ipoint,zposn,xposn,yposn,zposn_n,xposn_n,yposn_n,i,azimuth,inclination)  num_threads(omp_get_num_procs())
#pragma omp for reduction(+:decays238_retained,decays235_retained,decays232_retained,decays147_retained,decays_total) schedule(dynamic,1) nowait				
		for (ipoint = 1; ipoint<=imontes; ipoint++)
		{	
		nthreads = omp_get_num_threads();
			zposn = zspan * genrand_res53();
			assert ((zposn >= 0.0)&&(zposn <=zspan));  //remove once code working (avoid extra statements in MC loop!)
			xposn = xspan/2.0 - xspan * genrand_res53();
			yposn = yspan/2.0 - yspan * genrand_res53();	

// now do the alpha-ejection simulation; cheat a bit by using each selected position and alpha-vector to calculate ejection for all four nuclides
			if (in_cylinder(zposn,xposn,yposn) == 1)  // if we are in crystal or on surface, do an alpha-ejection
			{
				for (i = 1; i <= ejections; i++)  // for each monte-carlo position, do ejections number of random alpha ejections (see notes at start of code)
				{
					decays_total = decays_total + 1;
	//				inclination = pi * genrand_res53(); // old (bad) randomization of just angles
					inclination = acos(1.0 - 2.0 * genrand_res53());
					azimuth = 2*pi * genrand_res53();
	
					xposn_n = xposn + stopping238 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping238 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping238 * cos(inclination);
					
					if (in_cylinder(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays238_retained = decays238_retained + 1;
					}
					xposn_n = xposn + stopping235 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping235 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping235 * cos(inclination);
					if (in_cylinder(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays235_retained = decays235_retained + 1;
					}
					xposn_n = xposn + stopping232 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping232 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping232 * cos(inclination);
					if (in_cylinder(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays232_retained = decays232_retained + 1;
					}
					xposn_n = xposn + stopping147 * sin(inclination) * cos(azimuth);
					yposn_n = yposn + stopping147 * sin(inclination) * sin(azimuth);
					zposn_n = zposn + stopping147 * cos(inclination);
					if (in_cylinder(zposn_n,xposn_n,yposn_n) == 1)
					{
						decays147_retained = decays147_retained + 1;
					}
				} // end ejection loop
			}
		} // end of monte-carlo loop
	} // end of cylinder case	

// ********** FINISHED WITH GRAIN-GEOMETRY BLOCKS  *********************************

// We're done with the monte-carlo simulation for the component Ft's. Now do wrap-up calculations

// As a check report decays simulated and compare to expected decays if the monte carlo simulation uniformly distributed points across the
// problem space (expected decays should scale with fraction of problem space that is on or within the grain).

	printf("%3ld: Points in grain: %'lld out of %'lld requested (within %3.2f%%).\n",isample,decays_total,expected_decays,((float)(expected_decays - decays_total)/decays_total)*100.0);
	printf("     Retained decays: %'lld  %'lld  %'lld  %'lld   (8-5-2-7)\n",decays238_retained,decays235_retained,decays232_retained,decays147_retained);

	ftee238 = ((double) decays238_retained)/decays_total;
	ftee235 = ((double) decays235_retained)/decays_total;
	ftee232 = ((double) decays232_retained)/decays_total;
	ftee147 = ((double) decays147_retained)/decays_total;

//  here do bulk grain Ft, weighting by the previously calculated relative helium production rates
	ft_sample[isample] = (heprod238 * ftee238 + heprod235 * ftee235 + heprod232 * ftee232 + heprod147*ftee147)/heprod_total;

sampling_density = cbrt(((double)decays_total/ejections)/volume[isample]);

// Now find radius equivalent sphere for each grain. We are using weighted mean values for each decay chain. Let's try this:
// Determine Rs for each component ft, and then calculate a weighted mean Rs based on production.
// Need to use Newton-Raphson to invert the Ft polynomial for spheres.

	Ress238 = stopping238; // eff is quite sharply non linear with multiple roots. Plot shows
	                        // that minimum at negative values is about half a stopping distance
	                        //  higher values can also fail so need to keep starting seed lower.
	for (i=0; i<=1000;i++)
	{
		eff = 1 - ftee238 - 0.75*stopping238/Ress238 + pow(stopping238,3)/16./pow(Ress238,3);
		effprime = 0.75*stopping238/pow(Ress238,2) - 3./16.*pow(stopping238,3)/pow(Ress238,4);
		Ress238 = Ress238 - eff/effprime;
	}
	Ress235 = stopping235;
	for (i=0; i<=1000;i++)
	{
		eff = 1 - ftee235 - 0.75*stopping235/Ress235 + pow(stopping235,3)/16./pow(Ress235,3);
		effprime = 0.75*stopping235/pow(Ress235,2) - 3./16.*pow(stopping235,3)/pow(Ress235,4);
		Ress235 = Ress235 - eff/effprime;
	}
	Ress232 = stopping232;
	for (i=0; i<=1000;i++)
	{
		eff = 1 - ftee232 - 0.75*stopping232/Ress232 + pow(stopping232,3)/16./pow(Ress232,3);
		effprime = 0.75*stopping232/pow(Ress232,2) - 3./16.*pow(stopping232,3)/pow(Ress232,4);
		Ress232 = Ress232 - eff/effprime;
	}
	Ress147 = stopping147;
	for (i=0; i<=1000;i++)
	{
		eff = 1 - ftee147 - 0.75*stopping147/Ress147 + pow(stopping147,3)/16./pow(Ress147,3);
		effprime = 0.75*stopping147/pow(Ress147,2) - 3./16.*pow(stopping147,3)/pow(Ress147,4);
		Ress147 = Ress147 - eff/effprime;
	}
	
	Req[isample] = (heprod238 * Ress238 + heprod235 * Ress235 + heprod232 * Ress232 + heprod147 * Ress147)/heprod_total;

	printf("     Component Ft's: %6.4f %6.4f %6.4f %6.4f - grain Ft: %6.4f (%3.1f samples per micron)\n",ftee238,ftee235,ftee232,ftee147,ft_sample[isample],sampling_density);
	printf("\n");

} //  *************** END MAIN SAMPLE-PROCESSING ROUTINE: now do clean-up, bulk Ft, eU, other things

// Almost done! Now, where relevant, determine the bulk-sample Ft value for multi-grain samples, and calculate samples eU
// we're also going to do something fugly and for multigrain samples calculate a volume-weighted Ft-equivalent sphere (on the grounds that bigger grains
// would contribute more (we've already weighted Rs by the production for each decay chain, but in multigrain samples we have no way to know the parent
// values for each grain). See papers by Rich Ketcham for a more full discussion and rigorous approach.

	for (isample=1;isample<=nsamples;isample++)
	{
		eu_sample[isample] = U_ppm[isample] = Th_ppm[isample] = Sm_ppm[isample] = 0.0;  //initialize eU and ppm arrays
	}

	i = 1;
	do
	{	
		ftee = Ress = 0.0;
		tot_volume = 0.0;

		for (isample=i;isample<=i+ngrains[i]-1;isample++)
		{
			ftee = ftee + volume[isample] * ft_sample[isample];  // we'd like to weight by alpha production but for bulk multigrain analysis we can't just do that
			tot_volume = tot_volume + volume[isample];
			Ress = Ress + volume[isample] * Req[isample];
		}
		ftee = ftee/tot_volume;
		Ress = Ress/tot_volume;

		U_ppm[i] = uval[i]*1e9/tot_volume/sample_density[i];	
		Th_ppm[i] = thval[i]*1e9/tot_volume/sample_density[i];	
		Sm_ppm[i] = smval[i]*1e9/tot_volume/sample_density[i];
		
		eu_sample[i] = U_ppm[i] + 0.238 * Th_ppm[i] + 0.0012 * Sm_ppm[i];
		
		for (isample=i;isample<=i+ngrains[i]-1;isample++)
		{
			ft_group[isample] = ftee;
			Req_group[isample] = Ress;
		}
		i = i + ngrains[i];
	
	} while (i <= nsamples);

	// Parse inputs and set flags for writing data

		report_flag[1] = 1;
		grain_count = 1;
	for (i=2;i<=nsamples;i++)
	{
		if ((ngrains[i] == ngrains[i-1])&&(grain_count<ngrains[i]))
		{
			report_flag[i] = 0;	
			grain_count++;
		}
		else  //moved to new sample
		{
			report_flag[i] = 1;
			grain_count = 1;
		}

	}

	ofp = fopen(outfile, "w");

	fprintf(ofp,"Analysis\tFt-sample\tU (ppm)\tTh (ppm)\tSm (ppm)\tRs(Ft) (um)\teU-sample\tFt-grain\tGrains\tU (ng)\tTh (ng)\tSm (ng)\tMineral\tCalibration\tGeometry\tLength (um)\tWidth1 (um)\tWidth2 (um)\tTip1 (um)\tTip2 (um)\tvolume (um3)\tmass (ng)\tdensity\n");
	printf("<<<< SUMMARY RESULTS >>>>\n\n");
	printf("\t     Sample              Ft       Rs(Ft)    U (ppm)     Th (ppm)    Sm (ppm)    eU (ppm)    vol (µm3)   mass (ng)\n");
	printf("\t --------------------  -------   --------  ---------   ----------  ----------  ----------  ----------  ----------\n");


	for (i=1;i<=nsamples;i++)
	{
		if (report_flag[i] == 1)
		{
			fprintf(ofp,"%s\t%5.3f\t%6.1f\t%6.1f\t%6.1f\t%6.2f\t%6.1f\t%6.4f\t%2ld\t%7.4f\t%7.4f\t%7.4f\t%s\t%s\t%s\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sample[i],ft_group[i],U_ppm[i],Th_ppm[i],Sm_ppm[i],Req_group[i],eu_sample[i],ft_sample[i],ngrains[i],uval[i],thval[i],smval[i],mineral[i],calibration[i],geometry[i],length[i],width1[i],width2[i],tip1[i],tip2[i]);
			printf("  %ld\t  %-15s\t%5.3f     %6.2f    %7.2f      %7.2f     %7.2f     %7.2f   %8.0f      %6.1f \n",i,sample[i],ft_group[i],Req_group[i],U_ppm[i],Th_ppm[i],Sm_ppm[i],eu_sample[i],volume[i],volume[i]*sample_density[i]*1.0e-3);
		}
	}

	fprintf(ofp,"\t \nComplete Grain Information\n");

	for (i=1;i<=nsamples;i++)
	{
		if (report_flag[i] == 1)
		{
			fprintf(ofp,"%s\t%6.4f\t%6.1f\t%6.1f\t%6.1f\t%6.2f\t%6.1f\t%6.4f\t%2ld\t%7.4f\t%7.4f\t%7.4f\t%s\t%s\t%s\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%8.0f\t%6.1f\t%4.2f\n",sample[i],ft_group[i],U_ppm[i],Th_ppm[i],Sm_ppm[i],Req[i],eu_sample[i],ft_sample[i],ngrains[i],uval[i],thval[i],smval[i],mineral[i],calibration[i],geometry[i],length[i],width1[i],width2[i],tip1[i],tip2[i],volume[i],volume[i]*sample_density[i]*1.0e-3,sample_density[i]);
		}
		else
		{
			fprintf(ofp,"%s\t\t\t\t\t%6.2f\t\t%6.4f\t%ld\t\t\t\t\t%s\t%s\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n",sample[i],Req[i],ft_sample[i],ngrains[i],calibration[i],geometry[i],length[i],width1[i],width2[i],tip1[i],tip2[i]);	
		}
	}
	fclose(ofp);

// report rough performance (system clock and walltime for all cores)
	par_ended = omp_get_wtime();	
	ended = clock();
	elapsed = (double) (ended - started);

	printf("\n<<<< FINISHED >>>>\n");

	printf("\n  Done! Total run time was %6.2fs; total CPU time was %6.2fs\n",par_ended - par_started,elapsed/CLOCKS_PER_SEC);
	printf("        .... tried to use all %d cores found.\n",omp_get_num_procs());	

	printf("\n");
	printf("+++++++++++++++++++++++++++  ftee 1.13 finished  +++++++++++++++++++++++++++\n");
	printf("\n");
		
	return 0;
}              // END of MAIN()

// ********** Functions follow (first those related to grain geometry, and then those that belong to Mersenne Twister random-number generation)

int in_ellipsoid(double z, double x, double y)
{
	extern double rlength, rwidth1, rwidth2;
	double ellipsoid_arg;

	ellipsoid_arg = x*x/rwidth1/rwidth1 + y*y/rwidth2/rwidth2 + z*z/rlength/rlength;
	if (ellipsoid_arg <= 1.0)  // we must be on surface or in grain
	{
		return 1;		
	}
	else
	{
		return 0;
	}	
}

int in_cylinder(double z, double x, double y)
{
	extern double rlength, rwidth1, rwidth2;
	extern double zspan;

	double cylinder_arg;

	if ((z < 0.0)||(z > zspan))  // we can't be in grain -- we're beyond the z-limits
	{
		return 0;
	}
	else  // we could be in grain as far as zposn is concerned
	{
		cylinder_arg = x*x/rwidth1/rwidth1 + y*y/rwidth2/rwidth2;
		if (cylinder_arg <= 1.0)  //  we are inside cylinder or on its surface
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

int in_tetragonal(double z, double x, double y)
{
	extern double tip1[501], tip2[501];
	extern double z_end1, z_end2;
	extern long isample;
	extern double rwidth1, rwidth2;
	extern double zspan;

	if ((z < 0.0)||(z > zspan))  // can't be in grain -- beyond tips
	{
		return 0;
	}

	if (z < z_end1)  // we are in region of pyramidal tip1
	{
		if (fabs(x) <= rectx(z,tip1[isample]))  //  we are on or inside x-range of grain
		{
			if (fabs(y) <= recty(z,tip1[isample])) // we are on or inside grain
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			return 0;
		}
	}
	else
	{
		if (z > z_end2)  // we are in region of pyramidal tip2
		{
			if (fabs(x) <= rectx(zspan - z,tip2[isample]))  //  we are on or inside x-range of grain
			{
				if (fabs(y) <= recty(zspan - z,tip2[isample])) // we are on or inside grain
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return 0;
			}
		}
		else  // we are in tetragonal-prism section of grain
		{
			if (fabs(x) <= rwidth1)  //  we are on or inside x-range of grain
			{
				if (fabs(y) <= rwidth2) // we are on or inside grain
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return 0;
			}
		}
	}
}

double rectx (double z, double tiplength)
{
	extern double diag_angle, diag;
	double coord;
	
	if (tiplength <= 0.0) // should have already been handled by z-code, but to be sure trap no-pyramid case
	{
		coord = 0.0;
	}
	else
	{
		coord = fabs(z/tiplength * diag * cos(diag_angle));
	}
	return coord;
}

double recty (double z, double tiplength)
{
	extern double diag_angle, diag;
	double coord;
	
	if (tiplength <= 0.0) // should have already been handled by z-code, but to be sure trap no-pyramid case
	{
		coord = 0.0;
	}
	else
	{
		coord = fabs(z/tiplength * diag * sin(diag_angle));	
	}
	return coord;
}

// all of the following code is from Nishimura and Matsumoto as follows:
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


