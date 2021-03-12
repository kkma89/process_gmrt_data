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


Compilation:
------------
gcc -o process_gmrt_psr process_gmrt_psr.c 

Usage:
------------
process_gmrt_psr <input_file> {options}

For detailed instruction of usage, see the help page (-h option).

Detailed help:
------------
filterbank options

-n      : number of channels

-p      : number of polarisations

-f      : Frequency of the highest channel in MHz (def: 500.0)

-bw     : Recording bandwidth in MHz (def: 200.0)

-sb     : Sideband sense. -1 for LSB; 1 for USB (def: -1)

-mjd    : Timestamp of the first sample in MJD

-ts	: Sampling time (in secs)

-E	: Parameter file

-rfi    : Specify if the data is cleaned of RFI
            (0: norfix[default]; 1: gptool; 2: rfiClean)

DSPSR options

-dspsr : Process filterbank file with DSPSR (def: False)

-b     : Number of bins to fold (def: 128)

-tsub  : Sub-integration length in seconds (def: 10.0)

-t     : Number of thread to use (def: 1)

iterative_clean.py options

-clean : Clean fits file produced by DSPSR (def: False)

-m     : Number of iterations (def: 5)
