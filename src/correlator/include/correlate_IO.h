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
#ifndef CIO_H
#define CIO_H

#include "header.h"
#include "correlate_IOx64.h"
#include "correlate_MEM.h"

char *getDir(char *output, char *input);

void getFilter(char *sampleListFile, sample_t *sampleList);

int sample_isValid(sample_t *sampleList, char *sample, int *valid_count, int ploidy);

void readHeaderFile(char* inputPathName,
                    char ** headerLine1,
                    char ** headerLine2,
                    char* alignmentId,
                    int* snipsPerFile,
                    int* snipSize,
                    int* totalSnips,
                    int *posMin,
                    int *posMax);

void makeValidList(sample_t *sampleList,
                   char *headerLine2,
                   valid_t *validData);

void findFiles(char* inputPathName,
        char* alignmentId,
        int snipsPerFile,
        int totalSnips,
        int posWmin,
        int posWmax,
        char*** filesList,
        int* filesListNum,
        int mdf);

void writeResultsHeader(outFileType fpOut);

void writeResults(outFileType fpOut,
                  unsigned int *table_A_posIndex,
                  char **tableA_IDindex,
                  unsigned int *tableA_bitcount,
                  int tableASize,
                  unsigned int *table_B_posIndex,
                  char **tableB_IDindex,
                  unsigned int *tableB_bitcount,
                  int tableBSize,
                  int snp_size,
                  ResultDataType *results,
                  float r2limit,
                  int posWset2);

void writeStateCount(outFileType fpOut,
        unsigned int *table_A_posIndex,
        char **tableA_IDindex,
        int tableASize,
        unsigned int *table_B_posIndex,
        char **tableB_IDindex,
        int tableBSize,
        uint8_t *STtable,
        int posWset2);

void countStates_mode(FILE *fpOut,
                      table_x64 *tableData,
                      table_x64 *tableData2,
                      int posWset2);
#endif
