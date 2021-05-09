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
#include "../include/io.h"
#include "../include/mem.h"
#include "../include/read_file.h"

int getValueFromFilename(const void* filepath)
{
    size_t len;
    char *pdest, *file=NULL;
    int val;

    pdest=strrchr(*(const char**)filepath, '/');
    if(pdest == NULL)
        pdest=(char *)filepath;
    else
        pdest++;

    len=strlen(pdest);
    file=malloc(len+1);
    strcpy(file, pdest);
    sscanf(file, "%*[^_]_%d_%*[^$]", &val);
    free(file);
    return val;
}

int myCompare(const void* a, const void* b)
{
    int a_val=0,b_val=0;
    a_val=getValueFromFilename(a);
    b_val=getValueFromFilename(b);
    return a_val - b_val;
}

// Function to sort the array
void sort_files(char** arr, int n)
{
    // calling qsort function to sort the array
    // with the help of Comparator
    qsort(arr, n, sizeof(const char*), myCompare);
}

header_t *init_header_struct(void)
{
    header_t *headerData=(header_t *)malloc(sizeof(header_t));
    assert(headerData);

    headerData->alignmentID=(char *)malloc(IDLENGTH*sizeof(char));
    assert(headerData->alignmentID);
    headerData->headerLine1=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerData->headerLine1);
    headerData->headerLine2=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerData->headerLine2);
    headerData->valid_mask=(int *)malloc(sizeof(int));
    assert(headerData->valid_mask);
    headerData->valid_mask_size=0;
    headerData->snipsPerFile=0;
    headerData->snipSize=0;
    headerData->totalSnips=0;
    headerData->filesList=(char **)malloc(sizeof(char *));
    assert(headerData->filesList);
    headerData->filesList[0]=NULL;
    headerData->filesListNum=0;

    return headerData;
}

void print_header_struct(header_t *headerData)
{
    printf("ID: %s\n"
           "snipsPerFile: %d\n"
           "snipSize: %d\n"
           "totalSnips: %d\n",headerData->alignmentID,
                              headerData->snipsPerFile,
                              headerData->snipSize,
                              headerData->totalSnips);
    printf("%d Files:\n",headerData->filesListNum);
    for(int i=0; i < headerData->filesListNum; i++)
    {
        printf("\t%s\n", headerData->filesList[i]);
    }
}

void free_header_struct(header_t *headerData)
{
    free(headerData->alignmentID);
    free(headerData->headerLine1);
    free(headerData->headerLine2);
    free(headerData->valid_mask);
    for(int i=0; i < headerData->filesListNum; i++)
        free(headerData->filesList[i]);
    free(headerData->filesList);
    free(headerData);
}

char *getDir(char *output, char *input)
{
    char tmp[INFILENAMESIZE];
    strcpy(tmp, input);
    strcpy(output, tmp);
    if(output[strlen(output)-1] != '/')
        strcat(output,"/");
    return output;
}

void getFilter(char *filterFile, char ***wordList, int *wordListSize)
{
    if(!filterFile)
        return;
    FILE *fp=fopen(filterFile, "r");
    if(!fp)
    {
        fprintf(stderr,"ERROR: Sample List file \"%s\" doesn't exist\n",filterFile);
    }
    char line[STRINGLENGTH];
    char **temp;
    *wordListSize=0;
    while((fscanf(fp, "%s", line)) != EOF)
    {
        *wordListSize+=1;
        assert(*wordList);
        temp=(char **)realloc(*wordList,sizeof(char *)*(*wordListSize));
        assert(temp);
        *wordList=temp;
        (*wordList)[*wordListSize-1]=(char *)malloc(sizeof(char)*STRINGLENGTH);
        assert((*wordList)[*wordListSize-1]);
        strcpy((*wordList)[*wordListSize-1], line);
    }
    fclose(fp);
}

int sample_isValid(char **list, int list_size, char *sample, int ploidy, int *valid_count)
{
    if(!list_size)
    {
        (*valid_count)+=ploidy;
        return 1;
    }
    for(int i=0; i < list_size; i++)
    {
        if((strcmp(list[i],sample)) == 0)
        {
            (*valid_count)+=ploidy;
            return 1;
        }
    }
    return 0;
}

void readHeaderFile(pre_t *preData, header_t *headerData)
{
    int eol=0, eof=0, status, headerFound=0, index, validList_index=0;
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

    if((dp=opendir(preData->input)) == NULL)
    {
        fprintf(stderr,"cannot open directory: %s\n", preData->input);
        assert(0);
    }
    while((entry=readdir(dp)) != NULL)
    {
        strcpy(tmp_buff,preData->input);
        strcat(tmp_buff,entry->d_name);
        stat(tmp_buff,&statbuf);
        if(S_ISREG(statbuf.st_mode))
        {
            if(strstr(entry->d_name, "header") != NULL)
            {
                if(headerFound == 0)
                {
                    strcpy(headerFile,tmp_buff);
                    sscanf(entry->d_name,"%[^'_']",headerData->alignmentID);
                    sprintf(headerFile,"%s%s",preData->input,entry->d_name);
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

    if(headerFound == 0)
    {
        fprintf(stderr,"ERROR: Could not find Header file in %s\n",preData->input);
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
    printf("alignment %s found\n",headerData->alignmentID);
#endif

    inFileType fpHeader=FOPEN(headerFile,"r");
    status=getNextLine(fpHeader, &line, &eol, &eof, &lineLength);
    assert(line);
    headerData->headerLine1=(char *)realloc_buff((void *)headerData->headerLine1,
                                                 lineLength, "char");
    if(status == 1 && eol == 1)
    {
        strcpy(headerData->headerLine1, line);
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
                        headerData->headerLine2=(char *)realloc_buff(
                                                        (void *)headerData->headerLine2,
                                                        lineLength, "char");
                        strcpy(headerData->headerLine2,line);
                        index=0;
                        while(getWordFromString(line,&word,&eol,&wordLength,&index)==1)
                        {
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

                        headerData->valid_mask_size=1;

                        while(getWordFromString(line,&word,&eol,&wordLength,&index) == 1)
                        {
                            if(validList_index == headerData->valid_mask_size-1)
                            {
                                headerData->valid_mask_size+=BUFFERSIZE;
                                headerData->valid_mask=(int *)realloc_buff((void *)headerData->valid_mask,
                                                            (headerData->valid_mask_size), "int");
                            }
                            headerData->valid_mask[validList_index++]=sample_isValid(preData->sampleList,
                                                                                      preData->sampleListSize,
                                                                                      word,
                                                                                      preData->ploidy,
                                                                                      (&headerData->valid_count));
                            VCFsamples++;
                            if(eol == 1)
                                break;
                        }
                        assert(VCFsamples!=0);
                        status=getNextLine(fpHeader, &line, &eol, &eof, &lineLength);
                        if(status == 1 && eol==1)
                        {
                            sscanf(line,"%s %d %d %d", alignmentId_temp,
                                                       &(headerData->snipsPerFile),
                                                       &(headerData->snipSize),
                                                       &(headerData->totalSnips));
                            if((strcmp(alignmentId_temp,headerData->alignmentID))!=0)
                            {
#ifdef VERBOSE
                                fprintf(stderr, "Different internal and \
external ID: %s != %s\n", alignmentId_temp, headerData->alignmentID);
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
    char command[8+INFILENAMESIZE*2]; // Allocate space for command
    sprintf(command,"cp -r %s %s", headerFile, preData->output);
    system(command);
    free(line);
    free(word);
}

//this function checks the file names and keeps the ones that have snips
//included in the windows
void findFiles(pre_t *preData, header_t *headerData)
{
    DIR* dp;
    char file_found[INFILENAMESIZE];
    struct dirent *entry;
    struct stat statbuf;
    char tmpId[STRINGLENGTH];

    if((dp=opendir(preData->input)) == NULL)
    {
        fprintf(stderr,"cannot open directory: %s\n", preData->input);
        assert(0);
    }
    while((entry=readdir(dp)) != NULL)
    {
        strcpy(file_found,preData->input);
        strcat(file_found,entry->d_name);
        stat(file_found,&statbuf);
        if(S_ISREG(statbuf.st_mode))
        {
            if(strstr(entry->d_name, headerData->alignmentID) != NULL &&
               strstr(entry->d_name, "header") == NULL &&
               strstr(entry->d_name, ".vcf.gz") != NULL)
            {
                sscanf(entry->d_name,"%[^_]",tmpId);
                if(!strcmp(tmpId,headerData->alignmentID))
                {
                        headerData->filesListNum++;
                        headerData->filesList=(char**)realloc_buff((void *)headerData->filesList,
                                                   headerData->filesListNum,"char*");
                        headerData->filesList[headerData->filesListNum-1]=(char *)malloc(
                                                 sizeof(char)*INFILENAMESIZE);
                        sprintf(headerData->filesList[headerData->filesListNum-1],"%s%s",preData->input,
                                entry->d_name);
                }
            }
        }
    }
    closedir(dp);
    sort_files(headerData->filesList,headerData->filesListNum);
}
