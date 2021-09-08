/******************* uGMRT filterbank converter *******************/
/*
What?
-------------
All-in-one filterbank code to convert full polar or total power
GMRT data to sigproc filterbank format file. Currently written 
to convert only LSB data (highest frequency channel first).

Why?
-------------
GMRT full polarisation data is written in T1_C1[PQRS], T1_C2[PQRS], 
T1_C3[PQRS] ... format. While sigproc accepts data in a format which
is T1_P[C1...CN], T1_Q[C1...CN], T1_R[C1...CN], T1_S[C1...CN]. This
code essentially written to perform this conversion. 

After the conversion to filterbank, if asked, the code can convert
the filterbank file to psrfits using  DSPSR and also perform  RFI
excision using iterative_cleaner.py extracted from coast_guard 
pipeline (https://github.com/larskuenkel/iterative_cleaner).

------
Compilation:

gcc -o process_gmrt_psr process_gmrt_psr.c 
------
Usage:

process_gmrt_psr <input_file> {options}

For detailed instruction of usage, see the help page (-h option).
/******************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// Defining all the functions up front.
const char* mjd_calc(FILE *tfile);
struct get_levels {
float RRavg, RLavg, LLavg, LRavg;
};
typedef struct get_levels Struct;
Struct calc_levels(char *inputfile, int nch, int npol);
void send_string(FILE *outfile, char *string);
void send_int(FILE *outfile, char *name, int number);
void send_double(FILE *outfile, char *name, double number);
const char* get_srcdetails(FILE *parfile);
void run_dspsr (char *filename, char *parfile, int nbins, double tsub, int threads);
void run_clean(char *filename);
float seek_filesize(char *fpin);
void usage();



int main(int argc, char *argv[]) 
{
	if(argc < 7) 
	{
		usage();
		return 0;
	}

	FILE *infile, *outfile, *parfile;
	int i=1, j, nch=1024, npol=4, outpol=4, nbins=128, threads=1;
	int equal_poln=0, dspsr=0, rclean=0, m=5, rfi=0, sideband=-1;
	char input[500], output[500], psr[20], parafile[100], mjdstr[100];
	double fch1=500.0, bw=200.0, mjd, tsmpl, tsub=10.0;
	// Looping over the input arguments and reads them in.
	for (i=1;i<argc;i++)
	{
		if (fopen(argv[i],"rb") != NULL)
		{
			strcpy(input,argv[i]);
		} else if ((strcmp(argv[i], "-n")) == 0)
		{
			nch = atoi(argv[++i]);
		} else if ( (strcmp(argv[i], "-p")) == 0)
		{
			npol = atoi(argv[++i]);
		} else if ( (strcmp(argv[i], "-op")) == 0)
		{
			outpol = atoi(argv[++i]);
		} else if ( (strcmp(argv[i], "-f")) ==0 )
		{
			fch1 = atof(argv[++i]);
		} else if ( (strcmp(argv[i], "-bw")) ==0 )
		{
			bw = atof(argv[++i]);
		} else if ( (strcmp(argv[i], "-sb")) ==0 )
		{
			sideband = atoi(argv[++i]);
		} else if ( (strcmp(argv[i], "-mjd")) ==0 )
		{
			mjd = atof(argv[++i]);
		} else if ( (strcmp(argv[i], "-ts")) ==0 )
		{
			tsmpl = atof(argv[++i]);
		} else if ( (strcmp(argv[i], "-E")) ==0 )
		{
			strcpy(parafile,argv[++i]);
 		} else if ( (strcmp(argv[i], "-eq")) ==0 )
		{
			equal_poln = 1;			
		} else if ( (strcmp(argv[i], "-dspsr")) ==0 )
		{
			dspsr = 1;
		} 
		else if ( (strcmp(argv[i], "-b")) ==0 )
		{
			nbins = atoi(argv[++i]);
		} else if ( (strcmp(argv[i], "-tsub")) ==0 )
		{
			tsub = atof(argv[++i]);
		} else if ( (strcmp(argv[i], "-t")) ==0 )
		{
			threads = atoi(argv[++i]);
		} 
		else if ( (strcmp(argv[i], "-clean")) ==0 )
		{
			rclean = 1;
		} 
		else if ( (strcmp(argv[i], "-m")) ==0 )
		{
			m = atoi(argv[++i]);
		} 
		else if ( (strcmp(argv[i], "-rfi")) ==0 )
		{
			rfi = atoi(argv[++i]);
		} 
		else if ( (strcmp(argv[i], "-h")) ==0 )
		{
			usage();
		}
	}
	if (npol == outpol)
	{
		printf("\nWriting out the data in the input polarisation format\n");
	} if (outpol == 2)
	{
		printf("\nWriting out only RR,LL channels\n");
	} if (npol == 1)
	{
		outpol = 1;
		printf("\nWriting out the data in the input polarisation format\n");
	}
	// Defines the buffer sizes
	short buf[npol*nch], data4p[npol][nch], data1p[nch];
	//int intdata4p[npol][nch];
	// Opens the input file
	infile = fopen(input,"rb");
	if ((infile == NULL))
	{
		printf("\nThe input file is empty!!\n");
	}
	// Estimates the file size
	float infile_size = seek_filesize(input);
	float Lfactor[4] = {0.0, 0.0, 0.0, 0.0};
	if (outpol == 4)
	{
		if (equal_poln == 1)
		{
			Struct pl;
			pl = calc_levels(input, nch, npol);
			//Lfactor[0] = (pl.RRavg);			
			//Lfactor[1] = (pl.RLavg);
			//Lfactor[2] = (pl.LLavg);
			//Lfactor[3] = (pl.LRavg);
			Lfactor[1] = (pl.RRavg - pl.RLavg);
			Lfactor[2] = (pl.RRavg - pl.LLavg);
			Lfactor[3] = (pl.RRavg - pl.LRavg);
		}
	}
	printf("%f, %f, %f, %f\n",Lfactor[0], Lfactor[1], Lfactor[3], Lfactor[3]);
	//printf("%f %f %f %f\n",pl.RRavg, pl.RLavg, pl.LLavg, pl.LRavg);
	// If MJD is not given, it looks for the GMRT hdr
	// file and calculates the MJD from it using the
	// mjd_calc() function.
	if (mjd <= 1.0)
	{
		char tstampfile[500];
		FILE *tfile;
		strcpy(tstampfile,input);
		strcat(tstampfile,".hdr");
		tfile = fopen(tstampfile,"r");
		strcpy(mjdstr,mjd_calc(tfile));
		mjd = atof(mjdstr);
		//printf("%f\n",mjd);
	}	

	// Tries to open the parameter file and read the
	// pulsar name, RA and DEC for writing to the 
	// sigproc header.
	char *srcdetails, rfistat[20];
	char rajstr[100], decjstr[100];
	srcdetails = malloc(sizeof(char)*200);
	parfile = fopen(parafile,"r");
	if ((parfile != NULL))
	{
		strcpy(srcdetails,get_srcdetails(parfile));
		const char* aa = strtok(srcdetails, " ");
		strcpy(psr,aa);
		const char* bb = strtok(NULL ," ");
		strcpy(rajstr,bb);
		const char* cc = strtok(NULL ," ");
		strcpy(decjstr,cc);
	} else if (parfile == NULL)
	{
	strcpy(psr, "J0123+4567");
	strcpy(rajstr, "012345.67");
	strcpy(decjstr, "012345.67");	
	}
	fclose(parfile);

	// Checks the RFI mitigation status of the file
	// for giving the proper output file name.	
	if (rfi == 0)
	{
		strcpy(rfistat,"norfix");
	} else if (rfi == 1)
	{
		strcpy(rfistat,"gptool");
	} else if (rfi == 2)
	{
		strcpy(rfistat,"rfiClean");
	}

	// Creates the output filename in the regular format
	sprintf(output, "%s_%.6f_%d_%s.fil", psr, mjd, (int)fch1,rfistat);

	// Opens the output file for writing
	outfile = fopen(output, "wb");

	// Generates the sigproc header and writes it to
	// the output file.
	int machine_id = 14;
	int telescope_id = 7;
	int data_type = 1;
	int nbeams = 1;
	int ibeam = 1;
	int nbits = 16;
	
	send_string(outfile, "HEADER_START");
	send_int(outfile, "machine_id", machine_id);
	send_int(outfile, "telescope_id", telescope_id);
	send_int(outfile, "data_type", data_type);
	send_double(outfile, "fch1", fch1);
	send_double(outfile, "foff", (bw*sideband/nch));
	send_int(outfile, "nchans", nch);
	send_int(outfile, "nbeams", nbeams);
	send_int(outfile, "ibeam", ibeam);
	send_int(outfile, "nbits", nbits);
	send_double(outfile, "tstart", mjd);
	send_double(outfile, "tsamp", tsmpl);
	send_int(outfile, "nifs", outpol);
	send_string(outfile,"source_name");
	send_string(outfile,psr);
	send_double(outfile, "src_raj", atof(rajstr));
	send_double(outfile, "src_dej", atof(decjstr));
	send_string(outfile, "HEADER_END");
	// Beginning of finterbank conversion
	int progress=0;
	//float RRa, RLa, LLa, LRa;
	while (!feof(infile))
	{
		// Loop handling 4 polar data
		if (npol == 4) 
		{
			// Reads in data to buffer
			fread(buf, 2, nch*4, infile);
			for (i=0;i<npol;i++)
			{
				for (j=0;j<nch;j++)
				{
					data4p[i][j] = ((float)buf[(j*4)+i] + Lfactor[i]);
				}
			}
			// Writing out each polarisation channel
			// separately. Order is RR, LL, RL, LR.
			//int sum, k;
			//for (k=0;k<nch;k++)
			//{
			//	sum+=intdata4p[0][k];
			//}
			//printf("%d\n",sum/nch);
			if (outpol == 4)
			{
				//data4p[1] = ((float)*data4p[1]) / 988.;
				fwrite(&data4p[0],2,nch,outfile);
				fwrite(&data4p[2],2,nch,outfile);
				fwrite(&data4p[1],2,nch,outfile);
				fwrite(&data4p[3],2,nch,outfile);
			}
			if (outpol == 2)
			{
				
				//printf("%d %f %f %f %f\n",progress, (float)*data4p[0], (float)*data4p[1], (float)*data4p[2], (float)*data4p[3]);
				fwrite(&data4p[0],2,nch,outfile);
				fwrite(&data4p[2],2,nch,outfile);
				//fwrite(&data4p[1],2,nch,outfile);
				//fwrite(&data4p[3],2,nch,outfile);
			}

		}
		// Loop handling total power data
		else if (npol == 1)
		{
			int16_t buf[nch];
			// Reads in data to buffer
			fread(buf, 2, nch, infile);
			for (i=0;i<nch;i++)
			{
				data1p[i] = buf[i];
			}
			// Writes out the data in the same order
			// as the input.
  			fwrite(&data1p,2,nch,outfile);
		}
		progress++;
		printf("\r%.2f%% converted to filterbank",
				100.0*((float)(sizeof(buf)*progress)/infile_size));
	} // end of filterbank conversion
	printf("\n");
//	printf("\nhere: %f %f %f %f\n",RRa/progress, RLa/progress, LLa/progress, LRa/progress);
//	printf("\nhere: %f %f %f %f\n",RRa, RLa, LLa, LRa);
	
	fclose(infile);
	fclose(outfile);
	
	char fileroot[200];
	strcpy(fileroot,output);
	fileroot[strlen(fileroot)-4] = 0;

	// Running dspsr, if asked.
	if (dspsr == 1)	run_dspsr(fileroot, parafile, nbins, tsub, threads);

	// Running iterative_cleaner.py, if asked.
	if (dspsr == 0 && rclean == 1) 
	{
		printf("No fits files to run clean\n");
		exit(0);
	} else if (dspsr == 1 && rclean == 1) run_clean(fileroot);

} // End of main loop

///////////////////////////////////////////////////////////////////////////////


// Function for checking the input file size.
float seek_filesize(char *inpfile)
{
	FILE *fpin;
	fpin = fopen(inpfile,"rb");
	fseek(fpin, 0L, SEEK_END);
	float fpin_size = ftell(fpin);
	fclose(fpin);
	return(fpin_size);
}

Struct calc_levels(char *inputfile, int nch, int npol)
{
	Struct s;
	FILE *fpin1;
	int i, j, nsample;
	short buf1[npol*nch], data14p[npol][nch];
	float RRa, RLa, LLa, LRa;
	fpin1 = fopen(inputfile,"rb");
	while (!feof(fpin1))
	{	
		fread(buf1, 2, nch*4, fpin1);
		for (i=0;i<npol;i++)
		{
			for (j=0;j<nch;j++)
			{
				data14p[i][j] = buf1[(j*4)+i];
			}
		}
		RRa += (float)*data14p[0];
		RLa += (float)*data14p[1];
		LLa += (float)*data14p[2];
		LRa += (float)*data14p[3];
		nsample++;
	}
	fclose(fpin1);
	//printf("%d\n", nsample);
	s.RRavg = RRa / (float) nsample;
	s.RLavg = RLa / (float) nsample;
	s.LLavg = LLa / (float) nsample;
	s.LRavg = LRa / (float) nsample;
	return s;
}
	

// Function that runs DSPSR on the filterbank file.
void run_dspsr (char *filename, char *parfile, int nbins, double tsub, int threads)
{
	char dsprun[500];
	sprintf(dsprun, "dspsr -b %d -E %s -L %.1f -A -t %d -e fits -O %s %s.fil", 
			nbins, parfile, tsub, threads, filename, filename);
	printf("\n\nNow running DSPSR\n");
	printf("%s\n", dsprun);
	system(dsprun);
}

// Function to clean the fits file created by DSPSR 
// using iterative_cleaner.py
void run_clean(char *filename)
{
	char clrun[500];
	sprintf(clrun, "iterative_cleaner.py %s.fits -o %s.ar",filename, filename);
	printf("\n\nNow cleaning the fits file\n");
	printf("%s\n",clrun);
	system(clrun);
	
}

// Copied the below three functions from sigproc's
// "send_stuff.c" to create the header.
void send_string(FILE *outfile, char *string)
{
	int len;
	len = strlen(string);
	fwrite(&len, sizeof(int), 1, outfile);
	fwrite(string, sizeof(char), len, outfile);
}
void send_int(FILE *outfile, char *name, int number)
{
	send_string(outfile, name);
	fwrite(&number, sizeof(int), 1, outfile);
}
void send_double(FILE *outfile, char *name, double number)
{
	send_string(outfile, name);
	fwrite(&number, sizeof(double), 1, outfile);
}

// Function to extract pulsar name, RAJ and DECJ
// for writing it to the sigproc header.
const char* get_srcdetails(FILE *parfile)
{
	char line[300], tmp[30];
	char  psr[20], raj[30], decj[30];
	while (fgets(line, sizeof(line), parfile))
	{
		const char* c1 = strtok(line, " ");
		const char* c2 = strtok(NULL, " ");
		if (strcmp(c1, "PSRJ")==0){
			strcpy(psr,c2);
			strtok(psr, "\n");
		} else if (strcmp(c1, "RAJ")==0){
			strcpy(raj,c2);
		}else if (strcmp(c1, "DECJ")==0){
			strcpy(decj,c2);
		}

	}
	const char* r1 = strtok(raj,":");
	strcpy(tmp,r1);
	const char* r2 = strtok(NULL,":");
	strcat(tmp,r2);
	const char* r3 = strtok(NULL,":");
	strcat(tmp,r3);
	strcpy(raj,tmp);
	const char* d1 = strtok(decj,":");
	strcpy(tmp,d1);
	const char* d2 = strtok(NULL,":");
	strcat(tmp,d2);
	const char* d3 = strtok(NULL,":");
	strcat(tmp,d3);
	strcpy(decj,tmp);

	char * readout;
	readout = malloc(sizeof(char)*200);
	sprintf(readout,"%s %s %s",psr, raj, decj);

	return readout;
}


// Calculates the MJD of the observation from the GMRT
// standard header file. Part of the code is taken from
// https://www.atnf.csiro.au/computing/software/gipsy/sub/julianday.c
const char* mjd_calc(FILE *tfile)
{
char line[300], datestr[100], ddstr[10], mnstr[10], yystr[10];
char timestr[100], hhstr[10], mmstr[10], ssstr[20];

	while (fgets(line, sizeof(line), tfile))
	{
		const char* c1 = strtok(line, " ");
		const char* c2 = strtok(NULL, " ");
		if (strcmp(c1, "Date:")==0){
			strcpy(datestr,c2);
			strtok(datestr, "\n");
		} else if (strcmp(c1, "IST")==0){
			const char* c3 = strtok(NULL, " ");
			strcpy(timestr,c3);
		}

	}
	const char* r1 = strtok(datestr,":");
	strcpy(ddstr,r1);
	const char* r2 = strtok(NULL,":");
	strcpy(mnstr,r2);
	const char* r3 = strtok(NULL,":");
	strcpy(yystr,r3);
	const char* c1 = strtok(timestr,":");
	strcpy(hhstr,c1);
	const char* c2 = strtok(NULL,":");
	strcpy(mmstr,c2);
	const char* c3 = strtok(NULL,":");
	strcpy(ssstr,c3);
	
	// Start of the copied function from the above link.
	double A, B;
	double JD;
	double Gregorian = 1582.0 + 10.0/12.0 + 15.0/365.25;
	double day, month, year;
	double hour, mins, secs, dayfrac, mjd;
   
	hour  = atof(hhstr);
	mins  = atof(mmstr);
	secs  = atof(ssstr);
	day   = atof(ddstr);
	month = atof(mnstr);
	year  = atof(yystr);

	dayfrac = ((hour + (mins/60.) + (secs/3600))-5.5)/24.0;


	if (month < 3.0) 
	{
		year  -= 1.0; 
		month += 12.0; 
	}
	if ( (year + month/12.0 + day/365.25) >= Gregorian )
	{
		A = (double) ( (int)(year/100.0) );
		B = 2 - A + (double) ( (int)(A/4.0) );
	} else 
	{
		B = 0;
	}
	if (year >= 0.0) {
		JD = (double) ( (int)( 365.25 * year ) ) + 
			(double) ( (int)( 30.6001 * (month + 1.0) ) ) +
		day + 1720994.5 + B;
	} else {
		JD = (double) ( (int)( 365.25 * year - 0.75 ) ) + 
				(double) ( (int)( 30.6001 * (month + 1.0) ) ) +
		day + 1720994.5 + B;                       
	}
	// End of the copied function from the above link.
	
	mjd = JD+dayfrac-2400000.5;
	char * readout;
	readout = malloc(sizeof(char)*200);
	sprintf(readout,"%.16lf",mjd);
   return(readout);
}


// Prints help
void usage() 
{

	printf( "process_gmrt_psr - Processes  the  GMRT filterbank  data and\n"
			"                   and  converts it  to  different  formats.\n"
			"                   Currently supports conversion  to sigproc\n"
			"                   and psrfits format. There is an option to\n"
			"                   clean  the data  in psrfits  format using\n"
			"                   iterative_cleaner.py\n"
			"\nUsage: process_gmrt_psr {filename} -dspsr -clean {options}\n"
			"\n-h\t: Prints this help page\n"
			"\nfilterbank options\n\n"
			"-n\t: number of channels\n"
			"-p\t: number of polarisations in the input data\n"
			"-op\t: number of polarisations in the output data\n"
			"-f\t: Frequency of the highest channel in MHz (def: 500.0)\n"
			"-bw\t: Recording bandwidth in MHz (def: 200.0)\n"
			"-sb\t: Sideband sense. -1 for LSB; 1 for USB (def: -1)\n"
			"-mjd\t: Timestamp of the first sample in MJD\n"
			"-ts\t: Sampling time (in secs)\n"
			"-E\t: Parameter file\n"
			"-rfi\t: Specify if the data is cleaned of RFI\n"
			"\t (0: norfix[default]; 1: gptool; 2: rfiClean)\n"
			"-eq\t: Equalise the different polarisation bands.\n"
			"\t (Experimental and only works with -op=4 option)\n"
			"\nDSPSR options\n\n"
			"-dspsr\t: Process filterbank file with DSPSR (def: False)\n"
			"-b\t: Number of bins to fold (def: 128)\n"
			"-tsub\t: Sub-integration length in seconds (def: 10.0)\n"
			"-t\t: Number of thread to use (def: 1)\n"
			"\niterative_clean.py options\n\n"
			"-clean\t: Clean fits file produced by DSPSR (def: False)\n"
			"-m\t: Number of iterations (def: 5)\n"
		 );
}

