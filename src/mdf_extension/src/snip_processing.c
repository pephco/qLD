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
#include "../include/header.h"
#include "../include/io.h"
#include "../include/snip_processing.h"

int isValidVCFBase(char input)
{
    if(input=='A' || input=='C' ||
       input=='G' || input=='T' ||
       input=='N' || input=='a' ||
       input=='c' || input=='g' ||
       input=='t' || input=='n' ||
       input == '.')
        return 1;
    return 0;
}

int scanStateVector(char * stateVector, char X)
{
    int i;
    for(i=0;i < MAX_STATES_VCF;i++)
        if(stateVector[i] == X)
            return 1;
    return 0;
}

int getStatesREF(char * string, char * stateVector, int line)
{
    int i, j, index=0, elen=0, slen=strlen(string);
    char CharToStore = 'X';

    for(j=0;j < MAX_STATES_VCF;j++)
        stateVector[j] = CharToStore;

    for(i=0;i < slen;i++)
    {
        if(string[i] != ',')
        {
            if(string[i] == '<')
            {
                for(j=0;j < MAX_STATES_VCF;j++)
                    stateVector[j]='X';
                return 0;
            }
            assert(isValidVCFBase(string[i]) == 1);
            elen++;
            CharToStore=string[i];
        }
        if(string[i] == ',' || i == slen-1)
        {
            if(elen == 1)
            {
                if(scanStateVector(stateVector, CharToStore) == 0)
                {
                    stateVector[index++]=CharToStore;
                    elen=0;
                }
                else
                {
                    fprintf(stderr, "\n\n ERROR: Nucleotide %c (field REF) in position \
                            %d appears twice.\n\n", CharToStore, line);
                    exit(1);
                }
            }
            else
            {
                for(j=0;j < MAX_STATES_VCF;j++)
                    stateVector[j]='X';
                return 0;
            }
        }
    }
    return 1;
}

int getStatesNum(char * stateVector)
{
    int i, states=0;

    for(i=0;i < MAX_STATES_VCF;i++)
    {
        if(stateVector[i] != 'X')
            states++;
        else
            return states;
    }
    return states;
}


int getStatesALT(char * string, char * stateVector, int line)
{
    int i, j, index=0, elen=0, slen=strlen(string);

    char CharToStore='X';

    for(j=0;j < MAX_STATES_VCF;j++)
        if(stateVector[j] == 'X')
        {
            index=j;
            break;
        }

    assert(index != 0);

    for(i=0;i < slen;i++)
    {
        if(string[i] != ',')
        {
            if(string[i] == '<')
            {
                for(j=0;j < MAX_STATES_VCF;j++)
                    stateVector[j]='X';
                return 0;
            }
            assert(isValidVCFBase(string[i]) == 1);
            elen++;
            CharToStore=string[i];
        }
        if(string[i] == ',' || i == slen-1)
        {
            if(elen == 1)
            {
                if(scanStateVector(stateVector, CharToStore) == 0)
                {
                    stateVector[index++]=CharToStore;
                    elen=0;
                }
                else
                {
                    fprintf(stderr, "\n\n ERROR: Nucleotide %c (field ALT) in position %d\
                            appears twice (in REF or ALT).\n\n", CharToStore, line);
                    exit(1);
                }
            }
            else
            {
                for(j=0;j < MAX_STATES_VCF;j++)
                    stateVector[j]='X';
                return 0;
            }
        }
    }
    return 1;
}

float getValueAF(char * string, int line)
{
    int i, j, len=strlen(string), k=0;
    char AF_s[MAXAFLENGTH+1];
    AF_s[0]=0;
    float AF=-1.0, freqs[1000];

    if(strcmp(".",string) == 0 || len < 4)
        return -1;
    for(i=0;i < len-3;i++)
    {
        if(string[i] == 'A' && string[i+1] == 'F' && string[i+2] == '=')
        {
            if(i == 0)
                break;
            else
                if(string[i-1] == ';')
                    break;
        }
    }
    if(i == len-3)
        return -1.0;
    i=i+3;
    for(j=0;j < MAXAFLENGTH;j++)
    {
        if(string[i]==',')
        {
            assert(AF_s[0] != 0);
            freqs[k]=atof(AF_s);
            j=0;
            ++i;
            ++k;
        }
        if(string[i] == ';' || i == len)
            break;
        if(string[i] == '.' || (string[i] >= 48 && string[i] <= 57) ||
                string[i] == 'e' || string[i] == '-' || string[i] == '+')
        {
            AF_s[j]=string[i];
            i++;
            AF_s[j+1]=0;
        }
        else
        {
            fprintf(stderr, "\n\n ERROR: Invalid character (%c) in AF of field INFO in \
                    position %d.\n\n", string[i], line );
            exit(1);
        }
    }
    if(strlen(AF_s) == 0)
    {
        fprintf(stderr, "\n\n ERROR: AF in field INFO in position %d has no value.\n\n",
                line);
        assert(strlen(AF_s) != 0);
    }
    if(k == 0)
        AF = atof(AF_s);
    else
    {
        freqs[k] = atof(AF_s);
        ++k;
        AF=0.;
        for(i=0; i < k; ++i)
        {
            if(freqs[i] > 0. && freqs[i] < 1.)
            {
                AF=freqs[i];
                break;
            }
        }
    }
    if(AF < 0.0 || AF > 1.0)
    {
        fprintf(stderr, "\n\n ERROR: AF (%f) in line %d should be >= 0.0 and <= 1.0.\n\n",
                AF, line);
        assert(AF >= 0.0 && AF <= 1.0);
    }
    return AF;
}

int checkVTisSNP(char * string)
{
    int i, len=strlen(string);

    if(strcmp(".",string)==0 || len < 4)
        return -1;
    for(i=0;i < len-3;i++)
    {
        if(string[i] == 'V' && string[i+1] == 'T' && string[i+2] == '=')
        {
            if(i == 0)
                break;
            else
                if(string[i-1] == ';')
                    break;
        }
    }
    if(i == len-3)
        return -1;
    i=i+3;
    assert(i <= len-3);
    if(string[i] == 'S' && string[i+1] == 'N' && string[i+2] == 'P')
        return 1;
    else
        return 0;
}

int getGTpos(char * string)
{
    int i, len=strlen(string), GTposition=0;

    for(i=0;i < len-1;i++)
    {
        if(string[i] == ':')
            GTposition ++;

        if(string[i] == 'G' && string[i+1] == 'T')
        {
            if(i == 0)
            {
                if(i+2 == len)
                    return GTposition;
                else
                    if(string[i+2] == ':')
                        return GTposition;
            }
            else
            {
                if(string[i-1] == ':')
                {
                    if(i+2 == len)
                        return GTposition;
                    else
                        if(string[i+2] == ':')
                            return GTposition;
                }
            }
        }
    }
    GTposition=-1;
    return GTposition;
}

int getGTfield(char * string, int GTpos)
{
    int i=0, pos=GTpos, len=strlen(string), j=0, counter=0;
    assert(pos != -1);

    for(i=0;i < len;i++)
    {
        if(string[i] != ':')
        {
            if(pos == 0)
            {
                string[j++]=string[i];

                if(string[i] == '|' || string[i] == '/')
                    counter++;
            }
        }
        else
        {
            pos--;
            if(pos == -1)
            {
                string[j]=0;
                break;
            }
        }
    }
    assert(strlen(string) > 0);
    return counter+1;
}

void dataShuffleKnuth(char * data, int startIndex, int endIndex)
{
    if(startIndex == endIndex)
        return;

    int i, index;
    char tmp;

    for (i=endIndex; i > startIndex; i--)
    {
        index=startIndex + (rand() % (i - startIndex + 1));

        tmp=data[index];
        data[index]=data[i];
        data[i]=tmp;
    }
}

void getGTdata(char * string, 
               char * stateVector, 
               int statesTotal, 
               char * sampleData)
{
    int i, j=0, index=0, start=0, end=0, len=strlen(string);
    for(i=0;i < len;i++)
    {
        if(string[i] >= 48 && string[i] <= 57)
        {
            index=string[i]-48;

            assert(index < statesTotal);
            assert(j+1 < 4);
            sampleData[j++]=stateVector[index];
        }
        else
        {
            if(string[i] == '.')
            {
                sampleData[j++]='N';
            }

            if(string[i] == '/')
            {
                end++;
            }

            if(string[i] == '|')
            {

                assert(end < 4);
                dataShuffleKnuth(sampleData, start, end);
                start=j;
                end=j;
            }
        }
    }

    assert(end < 4);
    dataShuffleKnuth(sampleData, start, end);
}

void processSampleVCF(helper_t *helperData,
                      int GTpos,
                      char *stateVector,
                      int statesALT)
{
    int dataSize=getGTfield(helperData->word,GTpos);
    
    if(dataSize > helperData->ploidy)
    {
        fprintf(stderr,"ERROR: Detected data bigger than the ploidy of input\n");
        exit(1);
    }
    getGTdata(helperData->word, 
              stateVector, 
              statesALT, 
              helperData->snipPart);
    helperData->snipPart[dataSize]='\0';

    if(helperData->snipCharIndex + dataSize >= helperData->snipCharLength)
    {
        helperData->snipCharLength+=BUFFER_INCR;
        helperData->snipChar=(char *)realloc_buff((void *)helperData->snipChar,
                                                    helperData->snipCharLength, "char");
    }

    helperData->snipCharIndex+=dataSize;

    strcat(helperData->snipChar,helperData->snipPart);
}

int mapCharToInt(char a)
{
    if(a == ZERO)   return 0;
    if(a == ONE)    return 1;
    if(a == GAP || a == UN) return 2;
    if(a == AD)     return 3;
    if(a == CY)     return 4;
    if(a == GU)     return 5;
    if(a == TH)     return 6;
    assert(0);
    return 7;
}

int getDataType(int * states)
{
    int i,accum=0;

    if(states[STATESALL-1] != 0) // Invalid character was found
        return -1;
    for(i=3;i < STATESALL;i++)
        accum+=states[i];
    if((states[0] != 0 || states[1] != 0) && accum == 0) // BIN
    {
        if(states[2] != 0) // with gaps, Ns or both
            return 3;
        return 2;
    }
    accum=0;
    for(i=0;i < 2;i++)
        accum+=states[i];
    if((states[3] != 0 || states[4] != 0 || states[5] != 0 || states[6] != 0) &&
            accum == 0) // DNA
    {
        if(states[2] != 0) // with gaps, Ns or both
            return 5;
        return 4;
    }

    //den eixe case gia to ola missing....
    if(states[2] != 0)
        return 6;

    return 0;
}

int determineStates(int snipSize, char * line, int *states)
{
    int i;
    for(i=0;i < snipSize;i++)
        states[mapCharToInt(line[i])]++;
    return getDataType(states);
}

int removeNonPolymorphicSiteBIN(char* line, int snipSize, int filterOut)
{
    int i, dif, rm;
    int setVec[STATESALL];

    if(filterOut == 0)
    {
        for(i=0;i < STATESALL;i++)
            setVec[i]=0;
        rm=1;
        for(i=0;i < snipSize;i++)
        {
            setVec[mapCharToInt(line[i])]+=1;
            dif=setVec[0] != 0 ? 1 : 0;
            dif+=setVec[1] != 0 ? 1 : 0;
            if(dif == 2)
            {
                rm=0;
                break;
            }
        }
    }
    else
    {
        for(i=0;i < STATESALL;i++)
            setVec[i]=0;
        rm=1;
        for(i=0;i < snipSize;i++)
            setVec[mapCharToInt(line[i])]+=1;
        dif=setVec[0]!=0?1:0;
        dif+=setVec[1]!=0?1:0;
        if(dif == 2)
        {
            rm=0;
            if(setVec[0] == 1 || setVec[1] == 1)
                rm=1;
        }
    }
    if(setVec[2] != 0)
        rm=1;
    if(rm == 1)
    {
#ifdef VERBOSE
        printf("Discarded site\n");
#else
        ;
#endif
    }
    return rm;
}

int removeNonPolymorphicSiteDNA(char* line, int snipSize, int filterOut)
{
    int i, dif, rm;
    int setVec[STATESALL];

    if(filterOut == 0)
    {
        for(i=0;i < STATESALL;i++)
        {
            setVec[i]=0;
        }
        rm=1;
        for(i=0;i < snipSize;i++)
        {
            setVec[mapCharToInt(line[i])]+=1;
            dif=setVec[3] != 0 ? 1 : 0;
            dif+=setVec[4] != 0 ? 1 : 0;
            dif+=setVec[5] != 0 ? 1 : 0;
            dif+=setVec[6] != 0 ? 1 : 0;

            if(dif >= 2)
            {
                rm=0;
                break;
            }
        }
    }
    else
    {
        for(i=0;i < STATESALL;i++)
        {
            setVec[i]=0;
        }

        rm=1;

        for(i=0;i < snipSize;i++)
        {
            setVec[mapCharToInt(line[i])]+=1;
        }

        dif=setVec[3] != 0 ? 1 : 0;
        dif+=setVec[4] != 0 ? 1 : 0;
        dif+=setVec[5] != 0 ? 1 : 0;
        dif+=setVec[6] != 0 ? 1 : 0;

        if(dif >= 2)
        {
            rm=0;

            if(setVec[3] == 1 || setVec[4] == 1 || setVec[5] == 1 || setVec[6] == 1)
                rm=1;
        }
    }
    if(setVec[2] != 0)
    {
        rm=1;
    }
    if(rm == 1)
    {
#ifdef VERBOSE
        printf("Discarded site: %d %d %d %d %d\n",setVec[2], setVec[3],
                                                  setVec[4], setVec[5], setVec[6]);
        fflush(stdout);
#else
        ;
#endif
    }
    return rm;
}

int removeNonPolymorphicSite(char* line, int snipSize,int states, int filterOut)
{
    if(states == 2 || states == 3)
        return removeNonPolymorphicSiteBIN(line, snipSize, filterOut);

    if(states == 4 || states == 5)
        return removeNonPolymorphicSiteDNA(line, snipSize, filterOut);

    //den eixe case gia ola missing...
    if(states == 6)
        return 1;
    return 0;
}

void switchValues(int * input, char * sortedStates, int index0, int index1)
{
    int tmp;
    char tmpC;

    if(input[index1] > input[index0])
    {
        tmp=input[index0];
        tmpC=sortedStates[index0];
        input[index0]=input[index1];
        sortedStates[index0]=sortedStates[index1];
        input[index1]=tmp;
        sortedStates[index1]=tmpC;
    }
}

void sort(int * freqs, char * sortedStates)
{
    switchValues(freqs, sortedStates, 0, 1);
    switchValues(freqs, sortedStates, 2, 3);
    switchValues(freqs, sortedStates, 0, 2);
    switchValues(freqs, sortedStates, 1, 3);
    switchValues(freqs, sortedStates, 1, 2);
}

int minorToMajor(int* inp, char* outp)
{
    char major, alpha[4]={'A', 'C', 'G', 'T'};

    int total=0, rv, i;

    sort(inp,alpha);

    assert(inp[0] >= inp[1] && inp[1] >= inp[2] && inp[2] >= inp[3]);

    if(inp[2]+inp[3] == 0)
        return -1;

    total=inp[0]+inp[1];

    for(i=0; i < 2; ++i)
    {
        rv=rand() % total + 1;
        major=(rv <= inp[0]) ? alpha[0] : alpha[1];
        outp[i]=alpha[i+2];
        outp[i+2]=major;
    }
    return 1;
}

void countStates(int *numStates, int *gap, int *states)
{
    int i;
    if (*gap == 0)
        if (states[2] != 0)
            *gap=1;
    *numStates=0;
    states[2]=0;
    for(i=0;i < STATESALL;i++)
    {
        if(states[i] != 0)
        {
            states[i]=0;
            (*numStates)++;
        }
    }
}

char getCharBIN(char input, char state0, char state1)
{
    if(input == GAP)
        return GAP;
    if(input == UN)
        return UN;
    if(input == state0)
        return ZERO;
    if(input == state1)
        return ONE;
    return 'X';
}

void convertAlignmentBIN(int snipSize, char* line)
{
    int i;
    char state0, state1, tmp;

    state0='x';
    state1='x';

    for(i=0;i < snipSize;i++)
    {
        tmp=line[i];

        if(state0 == 'x' && (tmp != GAP && tmp != UN))
            state0=tmp;
        if(state0 != 'x' && state1 == 'x' && tmp != state0 && tmp != GAP && tmp != UN)
            state1=tmp;
        line[i]=getCharBIN(tmp,state0,state1);
    }
}

void binaryDeduction(int statesIn, int snipSize, char* line)
{
    if(statesIn == 2 || statesIn == 3)
        return;

    int i, states[STATESALL], withGaps=0, siteStates=0;
    char replaceMinorStates[4];

    for(i=0;i < 8;++i)
        states[i]=0;

    for(i=0;i < snipSize; i++)
    {
        int map=mapCharToInt(line[i]);

        if(map >= 3 && map <= 6)
            states[mapCharToInt(line[i])-3]++;
    }

    if(minorToMajor(states,replaceMinorStates) == 1)
    {
        for(i=0;i < snipSize;i++)
        {
            if(line[i] == replaceMinorStates[0])
                line[i]=replaceMinorStates[2];

            if(line[i] == replaceMinorStates[1])
                line[i]=replaceMinorStates[3];
        }
    }

    for(i=0;i < STATESALL;++i)
        states[i]=0;

    for(i=0;i < snipSize;i++)
        states[mapCharToInt(line[i])]++;

    countStates(&siteStates, &withGaps, states);
    assert(siteStates > 0);

    if(siteStates == 2)
    {
        convertAlignmentBIN(snipSize,line);
    }
}

void compressSnip_x64(table_x64 *tableData,
                     char * snipChar,                   // snipChar[0],
                     int size_snp)                      // snipSize,
{
    unsigned int i;
    int bitcount=0, compIndex=tableData->tableIndex*tableData->compSize;
    unsigned int remainder=0;
    inputDataType_x64 temp=0u;
    uint8_t newbit=0;

    for(i=0; i < (unsigned int)size_snp; i++)
    {
        uint8_t sC=(uint8_t)snipChar[i];
        newbit=sC & 1;
        temp= temp << 1;
        temp|=newbit;
        bitcount++;

        if(bitcount == (sizeof(inputDataType_x64)*8))
        {
            tableData->SNPtable[compIndex++]=temp;
            bitcount=0;
            temp=0;
        }
        tableData->BCtable[tableData->tableIndex]+=newbit;
    }

    remainder=sizeof(inputDataType_x64)*8-bitcount;

    if(remainder < sizeof(inputDataType_x64)*8)
        temp<<=remainder;

    tableData->SNPtable[compIndex]=temp;
}
