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
#ifndef STRC_H
#define STRC_H

#include "header.h"

typedef struct preprocess_t
{
    char *input;
    char *output;
    char **sampleList;
    int sampleListSize;
    int ploidy;
    int impute;
}pre_t;

typedef struct header_t
{
    char *alignmentID;
    char *headerLine1;
    char *headerLine2;
    int *valid_mask;
    int valid_count;
    int valid_mask_size;
    int snipsPerFile;
    int snipSize;
    int totalSnips;
    char **filesList;
    int filesListNum;
}header_t;

typedef struct table_x64{
    inputDataType_x64 *SNPtable;
    unsigned int *POStable;
    unsigned int *BCtable;
    unsigned int *BCtable_x32;
    char **IDtable;
    unsigned int tableSize;
    unsigned int tableIndex;
    unsigned int compSize;
    unsigned int compIndex;
}table_x64;

typedef struct helper_s{
    int ploidy;
    char *snipChar;
    int snipCharLength;
    int snipCharIndex;
    char *snipPart;
    char *line;
    int lineLength;
    char *word;
    int wordLength;
}helper_t;

#endif
