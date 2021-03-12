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
