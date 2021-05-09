/*
qLD - High performance computation of Linkage disequilibrium
Copyright (C) 2020  C. Theodoris, N. Alachiotis

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "zlib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <pthread.h>
#include <math.h>

#define INFILENAMESIZE 1000
#define STRINGLENGTH 1000
#define VCF_HLENGTH 9 // number of fields in the VCF header line
#define MAX_CHROM_NAME_VCF 100
#define MAX_STATES_VCF 5
#define MAXAFLENGTH 9

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3

#define STATESALL 8
#define ZERO '0'
#define ONE  '1'
#define GAP '-'
#define AD 'A'
#define CY 'C'
#define GU 'G'
#define TH 'T'
#define UN 'N'
#define ad 'a'
#define cy 'c'
#define gu 'g'
#define th 't'

#define min(a,b)     ( (a) > (b) ? (b) : (a) )
#define max(a,b)     ( (a) <= (b) ? (b) : (a) )

double gettime(void);
int getNextWord(gzFile fp, 
                char * word, 
                int *readEOL, 
                int *readEOF, 
                int *wordLength);
int getNextLine(gzFile fpIn, 
                char ** line, 
                int *readEOL, 
                int *readEOF, 
                int *lineLength);
void writeLine(gzFile fpOut,char **line);
int countLines(gzFile fpIn, char **line, int* lineLength);
int getWordFromString(char* line, 
                      char ** word, 
                      int *readEOL, 
                      int *wordLength, 
                      int *index);
int sortList(char* inputListName,int** list, char** line, int *lineLength);

gzFile createHeaderFile(gzFile fpIn, 
                        char* outputPathName, 
                        char ** headerLine1, 
                        char ** headerLine2, 
                        char **line, 
                        int* lineLength, 
                        char **word, 
                        int* wordLength);
int createPart(gzFile fpIn, 
               char* outputPathName, 
               char ** headerLine1, 
               char ** headerLine2, 
               char **line, 
               int* lineLength, 
               char **word, 
               int* wordLength, 
               int lines,
               int size,
               int* snipSize, 
               char* alignmentId, 
               int * maxcount,
               int *counter, 
               int *posMin, 
               int* posMax);
int createPartIndexW(gzFile fpIn, 
                     char* outputPathName, 
                     char ** headerLine1, 
                     char ** headerLine2, 
                     char **line, 
                     int* lineLength, 
                     char **word, 
                     int* wordLength, 
                     int lines, 
                     int size,
                     int* snipSize, 
                     char* alignmentId, 
                     int * maxcount,
                     int *counter, 
                     int *posMin, 
                     int* posMax, 
                     int Wmin, 
                     int Wmax);
int createPartPosW(gzFile fpIn, 
                   char* outputPathName, 
                   char ** headerLine1, 
                   char ** headerLine2, 
                   char **line, 
                   int* lineLength, 
                   char **word, 
                   int* wordLength, 
                   int lines, 
                   int size,
                   int* snipSize, 
                   char* alignmentId, 
                   int * maxcount,
                   int *counter, 
                   int *posMin, 
                   int* posMax, 
                   int posWmin, 
                   int posWmax);
int createPartPosL(gzFile fpIn, 
                   char* outputPathName, 
                   char ** headerLine1, 
                   char ** headerLine2, 
                   char **line, 
                   int* lineLength, 
                   char **word, 
                   int* wordLength, 
                   int lines, 
                   int size,
                   int* snipSize, 
                   char* alignmentId, 
                   int * maxcount,
                   int *counter, 
                   int *posMin, 
                   int* posMax, 
                   int ** inList,
                   int listSize);
int createSingleFileIndexW(gzFile fpIn, 
                           char* outputPathName, 
                           int *counter, 
                           char **line, 
                           int* lineLength, 
                           char **word, 
                           int* wordLength, 
                           int lines, 
                           int Wmin, 
                           int Wmax);
int createSingleFilePosW(gzFile fpIn, 
                         char* outputPathName, 
                         int *counter, 
                         char **line, 
                         int* lineLength, 
                         char **word, 
                         int* wordLength, 
                         int lines,
                         int posWmin, 
                         int posWmax);
int createSingleFilePosL(gzFile fpIn, 
                         char* outputPathName, 
                         int *counter, 
                         char **line, 
                         int* lineLength, 
                         char **word, 
                         int* wordLength, 
                         int lines, 
                         int ** inList,
                         int listSize);
