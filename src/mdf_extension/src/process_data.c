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
#include "../include/mem.h"
#include "../include/io.h"
#include "../include/process_data.h"
#include "../include/snip_processing.h"
#include "../include/read_file.h"

//this function reads lines from the files, creates the compressed snips
//and fills the tables
int writeToTable_x64(pre_t *preData,
                     header_t *headerData,
                     table_x64 *tableData,
                     helper_t *helperData,
                     int pos,
                     int *skippedLines,
                     int *counter)
{
    int eol=0, index=0;
    static int procVCF_count=0;
    unsigned int i;
    int elementIndex=-1, lineSkipped=0, s, rnps_flag=0;
    char stateVector[MAX_STATES_VCF];
    int position=-1, statesREF=0, statesALT=0, VTisSNP=0, GTpos=-1;
    int sampleIndex=-1;
    float AF=0.5;

    int statesNum=0, states[STATESALL];
    for(int i=0; i < STATESALL; ++i)
        states[i]=0;
    
    int wordLength=STRINGLENGTH;
    char ID[IDLENGTH];

    while(getWordFromString(helperData->line,
                            &(helperData->word),
                            &eol,
                            &wordLength,
                            &index) == 1)
    {
        elementIndex++;
        switch(elementIndex)
        {
            case 0:
                // CHROM
                //if(strcmp(VCF_alignment_name,*word)!=0)
                //{
                //strncpy(VCF_alignment_name, *word, MAX_CHROM_NAME_VCF);
                //assert(strlen(VCF_alignment_name)!=0);
                //}
                break;

            case 1: // POS
                position=atoi(helperData->word);
                break;

            case 2: // ID
                if(helperData->word)
                {
                    strcpy(ID,helperData->word);
                    if(ID[strlen(ID)] != '\0')
                        ID[strlen(ID)]='\0';
                }
                else
                {
                    fprintf(stderr,"ERROR: word not found at writeToTable_x64\n");
                    exit(1);
                }
                break;

            case 3: // REF
                if(strlen(helperData->word) > 1)
                {
                    lineSkipped=1;
                    break;
                }
                getStatesREF(helperData->word, stateVector, pos);
                statesREF=getStatesNum(stateVector);

                if(statesREF == 0)
                    lineSkipped=1;
                break;

            case 4: // ALT
                if(strlen(helperData->word) > 1)
                {
                    lineSkipped=1;
                    break;
                }
                getStatesALT(helperData->word, stateVector, pos);
                statesALT=getStatesNum(stateVector);

                if(statesALT != 0 && stateVector[statesALT-1] == '.')
                    for( s = statesALT - 1; s < statesALT-1 + statesREF; ++s)
                        stateVector[s]=stateVector[s-statesREF];
                if(statesALT == 0)
                    lineSkipped=1;
                break;

            case 5: // QUAL
                break;

            case 6: // FILTER
                if(strcmp("PASS", helperData->word) != 0 &&
                   strcmp(".", helperData->word) != 0)
                {
                    lineSkipped=1;
                }
                break;

            case 7: // INFO
                if(!strcmp(".", helperData->word))
                {
                    AF=.9;
                    break;
                }
                AF=getValueAF((helperData->word), pos); // returns -1.0 if AF not found
                if(AF == 0.0 || AF == 1.0)
                {
                    lineSkipped=1;
                    break;
                }
                VTisSNP=checkVTisSNP(helperData->word); // returns -1.0 if VT not found
                if(VTisSNP == 0)
                    lineSkipped=1;
                break;

            case 8: // FORMAT
                GTpos=getGTpos(helperData->word); // returns -1 if GT not found
                if(GTpos == -1)
                    lineSkipped=1;
                break;

            default:
                break;

        }
        if(elementIndex == VCF_HLENGTH-1 || eol == 1 || lineSkipped == 1)
            break;
    }
    if (position == -1)
    {
        return -1;
    }

    if(lineSkipped != 1)
    {
        sampleIndex=0;
        helperData->snipChar[0]='\0';
        // zero-out for new word
        while(getWordFromString(helperData->line,
                                &(helperData->word),
                                &eol,
                                &wordLength,
                                &index) == 1)
        {
            assert(sampleIndex < headerData->valid_mask_size);

            if(headerData->valid_mask[sampleIndex] == 1)
            {
                processSampleVCF(helperData,
                                 GTpos,
                                 stateVector,
                                 statesALT);
                procVCF_count++;
            }
            sampleIndex++;

            
            if(eol == 1)
                break;
        }
        helperData->snipCharIndex=0;
        if(tableData->tableIndex == 0)
            headerData->snipSize=strlen(helperData->snipChar);

        if(strlen(helperData->snipChar) != (size_t)headerData->snipSize)
        {
            fprintf(stderr, "\n\n ERROR: There are %lu nucleotides in position %d. "
                            "Expected %d.\n\n", strlen(helperData->snipChar),
                            position, headerData->snipSize);
            assert(strlen(helperData->snipChar) != (size_t)headerData->snipSize);
        }
        statesNum=determineStates(headerData->snipSize, helperData->snipChar, &(states[0]));
        
        //-----------------------------------------------------------
        if(statesNum == 3 && preData->impute)
        {
            //printf("Line is: %s\n",helperData->snipChar);
            //printf("zeros:%d ones:%d\n",states[0],states[1]); 
            double tHold=((double)states[0]) / 
                         ((double)(states[0] + states[1])); 
            //printf("Threshold: %f\n",tHold);
            //printf("Primary: %c Secondary: %c\n",stateVector[0], stateVector[1]);
                
            for(int smp=0; smp < headerData->snipSize; smp++)
            {
                if(helperData->snipChar[smp] == '.')
                {
                    double rVal=((double)rand()) / ((double)RAND_MAX);
                    //printf("Rand Value: %f\n",rVal);
                    helperData->snipChar[smp]=((rVal<=tHold)?stateVector[0]:stateVector[1]);
                }
            }

            //printf("Line is now: %s\n",helperData->snipChar);
        }

        if(statesNum == 5 && preData->impute)
        {
            //printf("Line is: %s\n",helperData->snipChar);
            //printf("A: %d, C: %d, G: %d, T: %d\n", states[3],
            //                                       states[4],
            //                                       states[5],
            //                                       states[6]); 
            double cumm=(double)(states[3] + states[4] + states[5] + states[6]);
            double tAC=(double)states[3] / cumm;
            double tCG=tAC + ((double)states[4] / cumm);
            double tGT=tCG + ((double)states[5] / cumm);
            //printf("T1: %f  T2: %f  T3: %f\n", tAC, tCG,  tGT);
            //printf("Primary: %c Secondary: %c\n",stateVector[0], stateVector[1]);
                
            for(int smp=0; smp < headerData->snipSize; smp++)
            {
                if(helperData->snipChar[smp] == 'N')
                {
                    double rVal=((double)rand()) / ((double)RAND_MAX);
                    //printf("Rand Value: %f\n",rVal);
                    if(rVal <= tAC)
                    {
                        helperData->snipChar[smp]='A';
                        //printf("Pos[%d]: A\n",smp);
                    }
                    else if(rVal <= tCG)
                    {
                        helperData->snipChar[smp]='C';
                        //printf("Pos[%d]: C\n",smp);
                    }
                    else if(rVal <= tGT)
                    {
                        helperData->snipChar[smp]='G';
                        //printf("Pos[%d]: G\n",smp);
                    }
                    else
                    {
                        helperData->snipChar[smp]='T';
                        //printf("Pos[%d]: T\n",smp);
                    }
                }
            }

            //printf("Line is now: %s\n",helperData->snipChar);
        }
        //-----------------------------------------------------------

        if(statesNum == -1)
        {
            fprintf(stderr, "\n\n ERROR: Empty alignment (no states found).\n\n");
            exit(1);
        }

        rnps_flag=removeNonPolymorphicSite(helperData->snipChar,
                                           headerData->snipSize,
                                           statesNum,
                                           0);
        if(!rnps_flag)
        {
            (*counter)++;
            //if(tableData->tableIndex >= 2)
            //{
            //    assert((unsigned int)position >=
            //            tableData->POStable[0][tableData->tableIndex-2]);
            //}
            if(tableData->tableIndex >= tableData->tableSize-1)
            {
                tableData->tableSize+=BUFFER_INCR;
                tableData->SNPtable=(inputDataType_x64 *)realloc_buff(
                                            (void *)tableData->SNPtable,
                                            tableData->tableSize*tableData->compSize,
                                            "inputDataType_x64");

                //print_alloc_size("SNPtable",tableData->tableSize*tableData->compSize);

                tableData->POStable=(unsigned int *)realloc_buff(
                                            (void *)tableData->POStable,
                                            tableData->tableSize,
                                            "unsigned int");
                tableData->BCtable=(unsigned int *)realloc_buff(
                                        (void *)tableData->BCtable,
                                        tableData->tableSize,
                                        "unsigned int");
                tableData->IDtable=(char **)realloc_buff(
                                        (void *)tableData->IDtable,
                                        tableData->tableSize,
                                        "char*");
                for(i=tableData->tableSize-BUFFER_INCR; i < tableData->tableSize; i++)
                    tableData->IDtable[i]=(char *)malloc(IDLENGTH*sizeof(char));
            }

            tableData->POStable[tableData->tableIndex]=position;
            tableData->BCtable[tableData->tableIndex]=0;
            strcpy(tableData->IDtable[tableData->tableIndex],ID);
            binaryDeduction(statesNum, headerData->snipSize, (helperData->snipChar));
            compressSnip_x64(tableData, helperData->snipChar, headerData->snipSize);
            assert(tableData->BCtable[tableData->tableIndex] <=
                   (unsigned int)headerData->snipSize);

            tableData->tableIndex++;
        }
        else
            (*skippedLines)++;
    }
    else
        (*skippedLines)++;

    return 0;
}

//function that reads lines from the files
void readTable_x64(pre_t *preData,
                   header_t *headerData,
                   table_x64 *tableData,
                   helper_t *helperData)
{
    int pos, eol, eof, status, counter=0, index, prev_val=0;
    int skippedLines=0;
    inFileType fpIn;
    FILE *fpOut;
    char temp_name[INFILENAMESIZE];

    char *writeBuffer=(char *)malloc(FILEBUFFERSIZE*sizeof(char));
    assert(writeBuffer); 

    for(int i=0; i < headerData->filesListNum; i++)
    {
        eol=0;
        eof=0;
        fpIn = FOPEN(headerData->filesList[i], "rb");

        if(fpIn == NULL)
        {
            fprintf(stderr,"ERROR: Failed to open file %s\n", headerData->filesList[i]);
            exit(1);
        }
        status=getNextLine(fpIn, &(helperData->line), &eol, &eof, &(helperData->lineLength));

        if(status == 1 && eol == 1)
        {
            if(strcmp(helperData->line, headerData->headerLine1))
            {
                fprintf(stderr,"ERROR: Header 1 different\n");
                exit(1);
            }
            else
            {
                status=getNextLine(fpIn,
                                   &(helperData->line),
                                   &eol,
                                   &eof,
                                   &(helperData->lineLength));
                if(status == 1 && eol==1)
                {
                    if(strcmp(helperData->line, headerData->headerLine2))
                    {
                        fprintf(stderr,"Header 2 different\n");
                        exit(1);
                    }
                    else
                    {
                        do
                        {
                            status=getNextLine(fpIn,
                                               &(helperData->line),
                                               &eol,
                                               &eof,
                                               &(helperData->lineLength));
                            if(status == 1 && eol == 1)
                            {
                                index=0;
                                status=getWordFromString(helperData->line,
                                                         &(helperData->word),
                                                         &eol,
                                                         &(helperData->wordLength),
                                                         &index); //chrom
                                assert(status);
                                status=getWordFromString(helperData->line,
                                                         &(helperData->word),
                                                         &eol,
                                                         &(helperData->wordLength),
                                                         &index); //pos
                                assert(status);
                                pos=atoi(helperData->word);

                                writeToTable_x64(preData,
                                                 headerData,
                                                 tableData,
                                                 helperData,
                                                 pos,
                                                 &skippedLines,
                                                 &counter);
                            }
                            else
                                pos=0;
                        }while(!eof);
                    }
                }
            }
        }
        else
        {
            fprintf(stderr,"ERROR: File %s empty\n", headerData->filesList[i]);
            exit(1);
        }
        FCLOZE(fpIn);
        strip_ext(headerData->filesList[i]);
        strip_ext(headerData->filesList[i]);
        strcpy(temp_name,basename(headerData->filesList[i]));
        strcpy(headerData->filesList[i],preData->output);
        strcat(headerData->filesList[i],temp_name);
        strcat(headerData->filesList[i],".mdf");
        printf("%s\n",headerData->filesList[i]);
        fpOut=fopen(headerData->filesList[i], "w");
        assert(fpOut);
        setvbuf(fpOut, writeBuffer, _IOLBF, FILEBUFFERSIZE);
        if(fpOut == NULL)
        {
            fprintf(stderr,"ERROR: Failed to open file %s\n", headerData->filesList[i]);
            exit(1);
        }

        fprintf(fpOut,"##fileformat=MDFv0.1\n##SNPs:%d\n#CHROM\tPOS\tID\tBC",
                tableData->tableIndex - prev_val);
        for(int k=0; k < tableData->compSize; k++)
            fprintf(fpOut,"\tcompSNP_%d",k);
        for(int j=prev_val; j < tableData->tableIndex; j++)
        {
            fprintf(fpOut,"\n%s\t%d\t%s\t%u",headerData->alignmentID,
                                             tableData->POStable[j],
                                             tableData->IDtable[j],
                                             tableData->BCtable[j]);
            for(int k=0; k < tableData->compSize; k++)
            {
                fprintf(fpOut,"\t%lu",tableData->SNPtable[j*tableData->compSize+k]);
            }
        }
        prev_val=counter;
        fflush(fpOut);
        fclose(fpOut);
    }
    free(writeBuffer);
}
