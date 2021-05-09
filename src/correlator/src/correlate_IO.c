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
#include "../include/correlate_IO.h"
#include "../include/read_file.h"
#include "../include/fast_print.h"

// Used to ensure that no missing "/" occur
char *getDir(char *output, char *input)
{
    char tmp[INFILENAMESIZE];
    strcpy(tmp, input);
    strcpy(output, tmp);
    if(output[strlen(output)-1] != '/')
        strcat(output,"/");
    return output;
}

// Get samples from sample list
void getFilter(char *sampleListFile, sample_t *sampleList)
{
    if(!sampleListFile)
        return;
    FILE *fp=fopen(sampleListFile, "r");
    char line[IDLENGTH];
    while((fscanf(fp, "%s", line)) != EOF)
    {
        if(sampleList->sampleListIndex >= sampleList->sampleListSize-1)
        {
            sampleList->sampleListSize+=BUFFER_INCR;
            sampleList->sampleList=(char **)realloc_buff(
                                (void *)sampleList->sampleList,
                                sampleList->sampleListSize,
                                "char*");
            for(int i=sampleList->sampleListSize-BUFFER_INCR; i < sampleList->sampleListSize; i++)
                sampleList->sampleList[i]=(char *)malloc(IDLENGTH*sizeof(char));
        }
        strcpy(sampleList->sampleList[sampleList->sampleListIndex], line);
        sampleList->sampleListIndex++;
    }
    fclose(fp);
}

// Check if sample is included in the sample list
int sample_isValid(sample_t *sampleList, char *sample, int *valid_count, int ploidy)
{
    if(!sampleList->sampleListIndex)
    {
        (*valid_count)+=ploidy;
        return 1;
    }
    for(int i=0; i < sampleList->sampleListIndex; i++)
    {
        if((strcmp(sampleList->sampleList[i],sample)) == 0)
        {
            (*valid_count)+=ploidy;
            return 1;
        }
    }
    return 0;
}

// Read VCF header file and parse info
void readHeaderFile(char* inputPathName,
        char ** headerLine1,
        char ** headerLine2,
        char* alignmentId,
        int* snipsPerFile,
        int* snipSize,
        int* totalSnips,
        int *posMin,
        int *posMax)
{
    int eol=0, eof=0, status, headerFound=0, index=0;
    char alignmentId_temp[STRINGLENGTH];
    char *line=(char *)malloc(STRINGLENGTH*sizeof(char));
    strcpy(line,"\0");
    assert(line);
    int lineLength=STRINGLENGTH;
    char *word=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(word);
    strcpy(word,"\0");
    int wordLength=STRINGLENGTH;
    char headerFile[STRINGLENGTH];
    char tmp_buff[STRINGLENGTH];
    DIR *dp;
    struct dirent *entry;
    struct stat statbuf;

    if((dp=opendir(inputPathName)) == NULL)
    {
        fprintf(stderr,"cannot open directory: %s\n", inputPathName);
        assert(0);
    }
    // Iterate through files
    while((entry=readdir(dp)) != NULL)
    {
        strcpy(tmp_buff,inputPathName);
        strcat(tmp_buff,entry->d_name);
        stat(tmp_buff,&statbuf);
        // Ensure it is indeed a file
        if(S_ISREG(statbuf.st_mode))
        {
            // Check for header in name
            if(strstr(entry->d_name, "header") != NULL)
            {
                // Duplicate checking
                if(headerFound == 0)
                {
                    // Get alignmentID from header
                    strcpy(headerFile,tmp_buff);
                    sscanf(entry->d_name,"%[^'_']",alignmentId);
                    sprintf(headerFile,"%s%s",inputPathName,entry->d_name);
//#ifdef VERBOSE
//                    fprintf(threadData[0].threadLog, "Found Header file %s\n",headerFile);
//#endif
                    headerFound=1;
                }
                else
                {
#ifdef VERBOSE
                    printf("NOTICE - Multiple Header files \
found!! : %s\nNOTICE - Processing first header file\n",entry->d_name);
#else
                    ;
#endif
                }
            }
        }
    }
    closedir(dp);

    // No header found
    if(headerFound == 0)
    {
#ifdef VERBOSE
        fprintf(stderr,"ERROR: Could not find Header file in %s\n",inputPathName);
#endif
        exit(1);
    }

    char * headerFields[VCF_HLENGTH];
    int fieldInd=0;
    int VCFsamples=0;
    headerFields[0] = "#CHROM";
    headerFields[1] = "POS";
    headerFields[2] = "ID";
    headerFields[3] = "REF";
    headerFields[4] = "ALT";
    headerFields[5] = "QUAL";
    headerFields[6] = "FILTER";
    headerFields[7] = "INFO";
    headerFields[8] = "FORMAT";

#ifdef VERBOSE
    printf("alignment %s found\n",alignmentId);
#endif

    inFileType fpHeader=FOPEN(headerFile,"r");
    status=getNextLine(fpHeader, &line, &eol, &eof, &lineLength);
    assert(line);
    headerLine1[0]=(char *)realloc_buff((void *)headerLine1[0], lineLength, "char");

    // Parse header file for info
    if(status == 1 && eol == 1)
    {
        strcpy(headerLine1[0], line);
        while(eof == 0)
        {
            status=getNextLine(fpHeader, &line, &eol, &eof, &lineLength);
            if(status == 1 && eol == 1)
            {
                if(strlen(line) >= 6)
                {
                    if((line)[0] == '#' && \
                            (line)[1] == 'C' && \
                            (line)[2] == 'H' && \
                            (line)[3] == 'R' && \
                            (line)[4] == 'O' && \
                            (line)[5] == 'M')
                    {
                        headerLine2[0]=(char *)realloc_buff((void *)headerLine2[0],
                                                            lineLength, "char");
                        strcpy(headerLine2[0],line);
                        index=0;
                        while(getWordFromString(line,&word,&eol,&wordLength,&index)==1)
                        {
                            // Count entries before samples
                            if(fieldInd<VCF_HLENGTH)
                                if(strcmp(headerFields[fieldInd],word)!=0)
                                {
                                    fprintf(stderr,"\n\n ERROR: VCF header field %s is \
                                            missing.\n\n",headerFields[fieldInd]);
                                    assert(0);
                                }
                            fieldInd++;
                            if(fieldInd>=VCF_HLENGTH)
                                break;

                            if(eol == 1)
                                assert(0);
                        }
                        while(getWordFromString(line,&word,&eol,&wordLength,&index) == 1)
                        {
                            // Count sample names
                            VCFsamples++;
                            if(eol == 1)
                                break;
                        }
                        assert(VCFsamples!=0);
                        status=getNextLine(fpHeader, &line, &eol, &eof, &lineLength);
                        if(status == 1 && eol==1)
                        {
                            // Get info stored in header
                            sscanf(line,"%s %d %d %d %d %d", alignmentId_temp,
                                    snipsPerFile, snipSize, totalSnips, posMin, posMax);
                            if((strcmp(alignmentId_temp,alignmentId))!=0)
                            {
#ifdef VERBOSE
                                fprintf(stderr, "Different internal and \
external ID: %s != %s\n", alignmentId_temp, alignmentId);
#endif
                                assert(0);
                            }
                        }
                        else
                        {
                            fprintf(stderr,"ERROR: Wrong header format\n");
                            assert(0);
                        }
                    }
                }
            }
        }
    }
    else if(eof==1)
    {
        fprintf(stderr,"ERROR: Empty header file\n");
        assert(0);
    }
    FCLOZE(fpHeader);
    free(line);
    free(word);
}

// Create valid samples map from sample list
void makeValidList(sample_t *sampleList,
                   char *headerLine1,
                   valid_t *validData)
{
    int eol=0, index=0, valid_index=0;
    char *word=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(word);
    strcpy(word,"\0");
    int wordLength=STRINGLENGTH;
    int fields=0;
    while(getWordFromString(headerLine1,&word,&eol,&wordLength,&index) == 1)
    {
        if(++fields >= VCF_HLENGTH)
            break;
    }
    validData->valid_count=0;
    while(getWordFromString(headerLine1,&word,&eol,&wordLength,&index) == 1)
    {
        validData->validList[valid_index++]=sample_isValid(sampleList,
                                                           word,
                                                           &(validData->valid_count),
                                                           sampleList->ploidy);
        if(eol == 1)
            break;
    }
    assert((validData->valid_count <= sampleList->ploidy*sampleList->sampleListIndex) || \
           (sampleList->sampleListIndex == 0));
    free(word);
}

// Check the file names and keep the ones that have snips
// included in the provided pos windows
void findFiles(char *inputPathName,
        char *alignmentId,
        int snipsPerFile,
        int totalSnips,
        int posWmin,
        int posWmax,
        char ***filesList,
        int *filesListNum,
        int mdf)
{
    DIR* dp;
    char file_found[INFILENAMESIZE];
    int tmpPosMinIndex=0, tmpMin, tmpMax, tmpPosMin;
    int tmpPosMax, prevTmpPosMax=0, prevTmpMax=1, eop=0;
    struct dirent *entry;
    struct stat statbuf;
    char tmpId[STRINGLENGTH];

    while(eop == 0)
    {
        if((dp=opendir(inputPathName)) == NULL)
        {
            fprintf(stderr,"cannot open directory: %s\n", inputPathName);
            assert(0);
        }
    
        // Iterate through files
        while((entry=readdir(dp)) != NULL)
        {
            strcpy(file_found,inputPathName);
            strcat(file_found,entry->d_name);
            stat(file_found,&statbuf);
            
            // Ensure it is indeed a file
            if(S_ISREG(statbuf.st_mode))
            {
                // Ensure it is the alignment and format we need
                if(strstr(entry->d_name, alignmentId) != NULL &&
                   strstr(entry->d_name, "header") == NULL &&
                   (mdf?strstr(entry->d_name, ".mdf"):strstr(entry->d_name, ".vcf.gz")))
                {
                    sscanf(entry->d_name,"%[^_]_%d_%d_%d_%d.",tmpId, &tmpMin, &tmpMax,
                            &tmpPosMin, &tmpPosMax);
                    if(prevTmpPosMax == 0)
                        prevTmpPosMax=tmpPosMax;
                    if(tmpMin == prevTmpMax)
                    {
                        // Ensure it is within the pos window and other validity checks
                        if((((tmpPosMin <= posWmin &&
                              posWmin <= tmpPosMax &&
                              tmpPosMinIndex == 0) ||
                             (prevTmpPosMax < posWmin &&
                              posWmin < tmpPosMin &&
                              tmpPosMinIndex == 0) ||
                             (tmpMin <= tmpPosMinIndex  &&
                              tmpPosMinIndex <= tmpMax)) &&
                            eop == 0) &&
                           (((tmpMax - tmpMin + 1)==snipsPerFile) ||
                            tmpMax == totalSnips))
                        {

                            // Add it to the files list
                            filesListNum[0]++;
                            filesList[0]=(char**)realloc_buff((void *)filesList[0],
                                                       filesListNum[0],"char*");
                            filesList[0][filesListNum[0]-1]=(char *)malloc(
                                                     sizeof(char)*INFILENAMESIZE);
                            sprintf(filesList[0][filesListNum[0]-1],"%s%s",inputPathName,
                                    entry->d_name);
                            if(tmpPosMax < posWmax)
                                tmpPosMinIndex=tmpMax+1;
                            else
                                eop=1;
                        }
                        prevTmpPosMax=tmpPosMax;
                        prevTmpMax=tmpMax+1;
                    }
                }
            }
       }
       closedir(dp);
    }
}

void writeResultsHeader(outFileType fpOut)
{
    FPRINT(fpOut,"POS2\tPOS1\tID2\tID1\tfreq[2]\tfreq[1]\tLD\n");
}

// Print the report
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
        int posWset2)
{
    int i,j;
    //printf("%d %d\n", tableASize,tableBSize);
    //FPRINT(fpOut," W2 pos\t- W1 pos\t| W2 ID\t- W1 ID\t| LD value\n");
    int cnt = 0;

    for(i=0;i<tableBSize;i++)
    {
        if(posWset2 == 0) cnt++;
        for(j=cnt;j<tableASize;j++)
        {
            if(results[i*tableASize+j] >= r2limit && results[i*tableASize+j] < 100)
            {
                //FPRINT(fpOut,"%u\t%u\t%s\t%s\t%.3f\t%.3f\t%1.5e\n",
                //       table_B_posIndex[i],table_A_posIndex[j],tableB_IDindex[i],
                //       tableA_IDindex[j],1.0f*tableB_bitcount[i]/snp_size,
                //       1.0f*tableA_bitcount[j]/snp_size, results[i*tableASize+j]);

                fprintf(fpOut, "%u\t%u\t",table_B_posIndex[i],table_A_posIndex[j]);
                fputs(tableB_IDindex[i], fpOut);
                fputc('\t', fpOut);
                fputs(tableA_IDindex[j], fpOut);
                fputc('\t', fpOut);
                fprint_f(fpOut,1.0f*tableB_bitcount[i]/snp_size,9); // Fast-print
                fputc('\t', fpOut);
                fprint_f(fpOut,1.0f*tableA_bitcount[j]/snp_size,9); //  >>
                fputc('\t', fpOut);
                fprint_f(fpOut,results[i*tableASize+j],9);          //  >>
                fputc('\n',fpOut);
            }
        }
    }
}

#ifdef COUNTSTATES
// Mode that writes number of states in a file. Different mode than Correlate.
void writeStateCount(outFileType fpOut,
        unsigned int *table_A_posIndex,
        char **tableA_IDindex,
        int tableASize,
        unsigned int *table_B_posIndex,
        char **tableB_IDindex,
        int tableBSize,
        uint8_t *STtable,
        int posWset2)
{
    int i,j;
    //printf("%d %d\n", tableASize,tableBSize);
    //FPRINT(fpOut," W2 pos\t- W1 pos\t| W2 ID\t- W1 ID\t| LD value\n");
    int cnt = 0;

    FPRINT(fpOut,"POS2\tPOS1\tID2\tID1\tNo_STATES\n");

    for(i=0;i<tableBSize;i++)
    {
        if(posWset2 == 0) cnt++;
        for(j=cnt;j<tableASize;j++)
        {
            //FPRINT(fpOut,"%u\t%u\t%s\t%s\t%.3f\t%.3f\t%1.5e\n",
            //       table_B_posIndex[i],table_A_posIndex[j],tableB_IDindex[i],
            //       tableA_IDindex[j],1.0f*tableB_bitcount[i]/snp_size,
            //       1.0f*tableA_bitcount[j]/snp_size, results[i*tableASize+j]);

            fprintf(fpOut,"%u\t%u\t",table_B_posIndex[i],table_A_posIndex[j]);
            fputs(tableB_IDindex[i], fpOut);
            fputc('\t', fpOut);
            fputs(tableA_IDindex[j], fpOut);
            fputc('\t', fpOut);
            fprintf(fpOut,"%u\t",STtable[j*tableBSize+i]);
            fputc('\n',fpOut);
        }
    }
}

// Counts number of states per SNP 
void countStates_mode(FILE *fpOut,
                      table_x64 *tableData,
                      table_x64 *tableData2,
                      int posWset2)
{
    int pos;
    uint8_t *STtable=(uint8_t*)malloc(tableData->tableSize*
                                      tableData2->tableSize*
                                      sizeof(uint8_t));
    assert(STtable);
    for(uint32_t indA=0; indA < (uint32_t)tableData->tableIndex; indA++)
    {
        for(uint32_t indB=0; indB < (uint32_t)tableData2->tableIndex; indB++)
        {
            pos=indA*tableData2->tableIndex+indB;
            STtable[pos]=0;
            for(uint32_t compD=0; compD < (uint32_t)tableData->compSize; compD++)
            {
                for(int bitP=0; bitP < 64; bitP++)
                {
                    if(STtable[pos] == 15)
                        break;
                    // STtable OR 1 in position "bitA*2 + bitB" where
                    // states are: A=0,B=0 -> STtable[0]
                    //             A=0,B=1 -> STtable[1]
                    //             A=1,B=0 -> STtable[2]
                    //             A=1,B=1 -> STtable[3]
                    STtable[pos]|=(1 << (((tableData->SNPtable[indA*tableData->compSize+compD] >> bitP) & 1)*2
                                        +((tableData2->SNPtable[indB*tableData->compSize+compD] >> bitP) & 1)));
                }
                STtable[pos]=__builtin_popcount(STtable[pos]);
            }
        }
    }
    writeStateCount(fpOut,
                    tableData->POStable,
                    tableData->IDtable,
                    tableData->tableIndex,
                    tableData2->POStable,
                    tableData2->IDtable,
                    tableData2->tableIndex,
                    STtable,
                    posWset2);
   free(STtable);
}
#endif
