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
#include "header.h"

gzFile createHeaderFile(gzFile fpIn,
        char* outputPathName,
        char ** headerLine1,
        char ** headerLine2,
        char **line,
        int* lineLength,
        char **word,
        int* wordLength)
{
    int eol=0, eof=0, status;
    char headerFile[INFILENAMESIZE+30];
    sprintf(headerFile,"%sheader.vcf.gz",outputPathName);

    gzFile fpOut=NULL;
    printf("Creating Header: %s\n\n",headerFile);
    fpOut = gzopen(headerFile,"w");

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        *headerLine1 = (char *) malloc(sizeof(char)*(strlen(*line)));
        strcpy(*headerLine1,*line);
        gzprintf(fpOut,"%s\n",*line);
        while(1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                writeLine(fpOut,line);
                if(strlen(*line) >=6)
                {
                    if((*line)[0] == '#' &&
                            (*line)[1] == 'C' &&
                            (*line)[2] == 'H' &&
                            (*line)[3] == 'R' &&
                            (*line)[4] == 'O' &&
                            (*line)[5] == 'M')
                    {
                        *headerLine2 = (char *) malloc(sizeof(char)*(strlen(*line)));
                        strcpy(*headerLine2,*line);
                        return fpOut;
                    }
                }
            }
            else if(eof==1)
                assert(0);
        }
    }
    else if(eof==1)
        assert(0);
    return NULL;
}

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
        int* posMax)
{
    char outputFile[INFILENAMESIZE+30];
    int status=-1,
        eof=0,
        eol=0,
        index = 0,
        pos=0,
        lineCounter=0,
        startingValue=0,
        endingValue=0,
        firstPos=0,
        newFile=0,
        printDelay=0;

    gzFile fpOut;
    char * allID = (char*)malloc(256*sizeof(char));
    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        strcpy(allID,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                (*snipSize)++;
            index++;
        }
        lineCounter=1;
        (*maxcount) = (size*1024*1024)/(double)((*snipSize)/8);

        if((*maxcount) < 0 || (*maxcount) > lines)
            (*maxcount) = lines;

        printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d\n\n\n",
                lines,(*snipSize),(*maxcount));

        (*counter)++;
        sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,alignmentId,1,
                (*maxcount),pos);
        startingValue = 1;
        endingValue = (*maxcount);
        //firstPos = pos;
        //(*posMin) = pos;
        //(*posMax) = pos;

        printf(" %d%% -- Creating File: %s\n", ((*counter)*100)/lines,outputFile);

        fpOut = gzopen(outputFile,"w");
        if(fpOut == NULL)
        {
            fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
            free(allID);
            return 0;
        }
        writeLine(fpOut,headerLine1);
        writeLine(fpOut,headerLine2);
        if(!strcmp(alignmentId,allID))
            writeLine(fpOut,line);
        newFile = 0;

        if( lineCounter >= lines)
        {
            gzclose(fpOut);

            char outputFileNew[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("100%% -- Creating File: %s\n",outputFile);
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    startingValue,(*counter),
                    firstPos,(*posMax));
            printf("100%% -- Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                free(allID);
                return 0;
            }
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        free(allID);
        return 0;
    }

    while(1)
    {
        index = 0;
        if(newFile == 1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                //char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;
                if(!strcmp(allID,alignmentId))
                {
                    (*counter)++;
                    (*posMax)=pos;
                    sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,
                            alignmentId,(*counter),
                            (*counter)+(*maxcount)-1,
                            pos);
                    startingValue = (*counter);
                    endingValue = (*counter)+(*maxcount)-1;
                    firstPos = pos;
                    if((*posMin) == 0)
                        (*posMin) = pos;

                    //printf("\033[A\33[2K");
                    printf("\033[A\33[2K");
                    printf(" %d%% -- Creating File: %s\n",((*counter)*100)/lines,
                            outputFile);

                    fpOut = gzopen(outputFile,"w");
                    if(fpOut == NULL)
                    {
                        fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
                        free(allID);
                        return 0;
                    }
                    writeLine(fpOut,headerLine1);
                    writeLine(fpOut,headerLine2);
                    writeLine(fpOut,line);
                    newFile = 0;

                    if(lineCounter >= lines)
                    {
                        gzclose(fpOut);

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        printf("100%% -- Creating File: %s\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                (*posMin),(*posMax));
                        printf("100%% -- Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            free(allID);
                            return 0;
                        }
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }
            else if(eof ==1)
            {
                printf("%d Snips written\n",(*counter));
                free(allID);
                return 1;
            }
        }
        else
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                //char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;

                if(((*counter)*100/lines)%5 == 0 &&
                        printDelay != ((*counter)*100/lines))
                {
                    printf("\033[A\33[2K");
                    printf(" %d%% -- Creating File: %s\n",((*counter)*100)/lines,
                            outputFile);
                    printDelay = ((*counter)*100/lines);
                }
                
                if(!strcmp(allID,alignmentId))
                {     
                    writeLine(fpOut,line);
                    (*counter)++;
                    (*posMax)=pos;
                    if((*posMin) == 0)
                        (*posMin) = pos;
                    if(lineCounter >= lines)
                    {
                        gzclose(fpOut);

                        printf("\033[A\33[2K");
                        char outputFileNew[INFILENAMESIZE+30];
                        printf("100%% -- Creating File: %s\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                (*posMin),(*posMax));
                        printf("100%% -- Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            free(allID);
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        free(allID);
                        return 1;
                    }
                    else if((*counter) == endingValue)
                    {
                        newFile = 1;
                        gzclose(fpOut);

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                (*posMin),(*posMax));
                        printf(" %d%% -- Renaming %s to %s\n",((*counter)*100)/lines,
                                outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            free(allID);
                            return 0;
                        }
                    }
                }
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                //free(allID);
            }
            else if(eof ==1)
            {
                gzclose(fpOut);

                char outputFileNew[INFILENAMESIZE+30];
                printf("100%% -- Creating File: %s\n",outputFile);
                sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                        alignmentId,
                        startingValue,
                        (*counter),(*posMin),(*posMax));
                printf("100%% -- Renaming %s to %s\n",outputFile,outputFileNew);
                int ret = rename(outputFile, outputFileNew);
                if(ret != 0)
                {
                    fprintf(stderr,"\n ERROR: unable to rename the file\n");
                    free(allID);
                    return 0;
                }
                printf("%d Snips written\n",(*counter));
                free(allID);
                return 1;
            }
        }
    }
    printf("%d Snips written\n",(*counter));
    free(allID);
    return 1;
}

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
        int Wmax)
{
    char outputFile[INFILENAMESIZE+30];
    int status=-1,
        eof=0,
        eol=0,
        index = 0,
        pos=0,
        lineCounter=0,
        startingValue=0,
        endingValue=0,
        firstPos=0,
        newFile=0,
        printDelay=0;

    gzFile fpOut;

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                (*snipSize)++;
            index++;
        }
        lineCounter++;
        *maxcount = (size*1024*1024)/((double)((*snipSize)/8));

        if(*maxcount < 0 || *maxcount > lines)
            *maxcount = lines;

        printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, Wmin: %d, \
                Wmax: %d\n\n\n",lines,*snipSize,*maxcount, Wmin, Wmax);

        if(lineCounter >= Wmin)
        {
            (*counter)++;
            sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,alignmentId,1,
                    (*maxcount),pos);
            startingValue = 1;
            endingValue = (*maxcount);
            firstPos = pos;
            *posMin = pos;

            printf("Creating File: %s, %d%%\n",outputFile,((*counter)*100)/(Wmax-Wmin+1));

            fpOut = gzopen(outputFile,"w");
            if(fpOut == NULL)
            {
                fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
                return 0;
            }
            writeLine(fpOut,headerLine1);
            writeLine(fpOut,headerLine2);
            writeLine(fpOut,line);
            newFile = 0;

            if(lineCounter >= Wmax)
            {
                *posMax = pos;
                gzclose(fpOut);

                char outputFileNew[INFILENAMESIZE+30];
                printf("\033[A\33[2K");
                printf("Creating File: %s, 100%%\n",outputFile);
                sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                        alignmentId,
                        startingValue,(*counter),
                        firstPos,pos);
                printf("Renaming %s to %s\n",outputFile,outputFileNew);
                int ret = rename(outputFile, outputFileNew);
                if(ret != 0)
                {
                    fprintf(stderr,"\n ERROR: unable to rename the file\n");
                    return 0;
                }
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
        else
            newFile = 1;


    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }

    while(1)
    {
        index = 0;
        if(newFile == 1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;
                if(!strcmp(allID,alignmentId))
                {
                    if(lineCounter >= Wmin)
                    {
                        (*counter)++;
                        sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,
                                alignmentId,
                                (*counter),
                                (*counter)+(*maxcount)-1,
                                pos);
                        startingValue = (*counter);
                        endingValue = (*counter)+(*maxcount)-1;
                        firstPos = pos;
                        if(*posMin == 0)
                        {
                            *posMin = pos;
                        }

                        printf("\033[A\33[2K");
                        printf("Creating File: %s, %d%%\n",outputFile,
                                ((*counter)*100)/(Wmax-Wmin+1));

                        fpOut = gzopen(outputFile,"w");
                        if(fpOut == NULL)
                        {
                            fprintf(stderr,"\n ERROR: failed to open file %s",
                                    outputFile);
                            return 0;
                        }
                        writeLine(fpOut,headerLine1);
                        writeLine(fpOut,headerLine2);
                        writeLine(fpOut,line);
                        newFile = 0;
                    }
                    else
                        newFile = 1;

                    if(lineCounter >= Wmax)
                    {
                        *posMax = pos;
                        gzclose(fpOut);
                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: No snips found inside window\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        printf("Creating File: %s, 100%%\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }else if(eof ==1)
            {
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
        else {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;

                if(((*counter)*100/(Wmax-Wmin+1) )%5 == 0 &&
                        printDelay != (*counter)*100/(Wmax-Wmin+1))
                {
                    printf("\033[A\33[2K");
                    printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/(Wmax-Wmin+1));
                    printDelay = (((*counter)*100)/(Wmax-Wmin+1));
                }

                if(!strcmp(allID,alignmentId))
                {
                    writeLine(fpOut,line);
                    (*counter)++;

                    if(lineCounter >= Wmax)
                    {
                        *posMax = pos;
                        if(fpOut != NULL)
                        {
                            gzclose(fpOut);
                        }
                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: No snips found inside window\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("Creating File: %s, 100%%\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                    else if((*counter) == endingValue)
                    {
                        newFile = 1;
                        gzclose(fpOut);

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }
            else if(eof ==1)
            {
                gzclose(fpOut);

                if((*counter) == 0)
                {
                    fprintf(stderr,"\n ERROR: No snips found inside window\n");
                    return 0;
                }
                char outputFileNew[INFILENAMESIZE+30];
                printf("Creating File: %s, 100%%\n",outputFile);
                sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                        alignmentId,
                        startingValue,(*counter),
                        firstPos,pos);
                printf("Renaming %s to %s\n",outputFile,outputFileNew);
                int ret = rename(outputFile, outputFileNew);
                if(ret != 0)
                {
                    fprintf(stderr,"\n ERROR: unable to rename the file\n");
                    return 0;
                }
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
    }
    printf("%d Snips written\n",(*counter));
    return 1;
}

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
        int posWmax)
{
    char outputFile[INFILENAMESIZE+30];
    int status=-1,
        eof=0,
        eol=0,
        index = 0,
        pos=0,
        lineCounter=0,
        startingValue=0,
        endingValue=0,
        firstPos=0,
        newFile=0;

    gzFile fpOut;

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                (*snipSize)++;
            index++;
        }
        lineCounter++;
        *maxcount = (size*1024*1024)/((double)((*snipSize)/8));

        if(*maxcount < 0 || *maxcount > lines)
            *maxcount = lines;

        printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, posWmin: %d, \
                posWmax: %d\n\n\n",lines,*snipSize,*maxcount, posWmin, posWmax);

        if(pos >= posWmin)
        {
            (*counter)++;
            sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,alignmentId,1,
                    (*maxcount),pos);
            startingValue = 1;
            endingValue = (*maxcount);
            firstPos = pos;
            *posMin = pos;

            printf("Creating File: %s\n",outputFile);

            fpOut = gzopen(outputFile,"w");
            if(fpOut == NULL)
            {
                fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
                return 0;
            }
            writeLine(fpOut,headerLine1);
            writeLine(fpOut,headerLine2);
            writeLine(fpOut,line);
            (*counter)++;
            newFile = 0;
        }
        else
            newFile = 1;

        if(pos >= posWmax)
        {
            *posMax = pos;
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }

            char outputFileNew[INFILENAMESIZE+30];
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    startingValue,(*counter),
                    firstPos,pos);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }

    while(1)
    {
        index = 0;
        if(newFile == 1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;
                if(!strcmp(allID,alignmentId))
                {
                    if(pos >= posWmin)
                    {
                        (*counter)++;
                        sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,
                                alignmentId,
                                (*counter),
                                (*counter)+(*maxcount)-1,
                                pos);
                        startingValue = (*counter);
                        endingValue = (*counter)+(*maxcount)-1;
                        firstPos = pos;
                        if(*posMin == 0)
                        {
                            *posMin = pos;
                        }

                        printf("\033[A\33[2K");
                        printf("Creating File: %s\n",outputFile);

                        fpOut = gzopen(outputFile,"w");
                        if(fpOut == NULL)
                        {
                            fprintf(stderr,"\n ERROR: failed to open file %s",outputFile);
                            return 0;
                        }
                        writeLine(fpOut,headerLine1);
                        writeLine(fpOut,headerLine2);
                        writeLine(fpOut,line);
                        newFile = 0;
                    }
                    else
                        newFile = 1;

                    if(pos >= posWmax)
                    {
                        *posMax = pos;
                        gzclose(fpOut);
                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: No snips found inside window\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }else if(eof ==1)
            {
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
        else
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;

                if(!strcmp(allID,alignmentId))
                {
                    writeLine(fpOut,line);
                    (*counter)++;

                    if(pos >= posWmax)
                    {
                        *posMax = pos;
                        gzclose(fpOut);

                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: No snips found inside window\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,startingValue,(*counter),firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                    else if((*counter) == endingValue)
                    {
                        newFile = 1;
                        gzclose(fpOut);

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }
            else if(eof ==1)
            {
                gzclose(fpOut);

                if((*counter) == 0)
                {
                    fprintf(stderr,"\n ERROR: No snips found inside window\n");
                    return 0;
                }
                char outputFileNew[INFILENAMESIZE+30];
                sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                        alignmentId,
                        startingValue,(*counter),
                        firstPos,pos);
                printf("Renaming %s to %s\n",outputFile,outputFileNew);
                int ret = rename(outputFile, outputFileNew);
                if(ret != 0)
                {
                    fprintf(stderr,"\n ERROR: unable to rename the file\n");
                    return 0;
                }
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
    }
}

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
        int listSize)
{
    char outputFile[INFILENAMESIZE+30];
    int status=-1,
        eof=0,
        eol=0,
        index = 0,
        pos=0,
        lineCounter=0,
        startingValue=0,
        endingValue=0,
        firstPos=0,
        newFile=0,
        inListPos=0,
        printDelay=0;

    gzFile fpOut;

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                (*snipSize)++;
            index++;
        }
        lineCounter++;
        *maxcount = (size*1024*1024)/((double)((*snipSize)/8));

        if(*maxcount < 0 || *maxcount > lines)
            *maxcount = lines;

        printf("Total snips: %d, Snip Size: %d bits, Snips per file: %d, using input \
                List\n\n\n",lines,*snipSize,*maxcount);

        if(pos == (*inList)[inListPos])
        {
            (*counter)++;
            sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,alignmentId,1,
                    (*maxcount),pos);
            startingValue = 1;
            endingValue = (*maxcount);
            firstPos = pos;
            *posMin = pos;

            printf("Creating File: %s, %d%%\n",outputFile, ((*counter)*100)/listSize);

            fpOut = gzopen(outputFile,"w");
            if(fpOut == NULL)
            {
                fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
                return 0;
            }
            writeLine(fpOut,headerLine1);
            writeLine(fpOut,headerLine2);
            writeLine(fpOut,line);
            (*inList)[inListPos] = 0;
            inListPos++;
            newFile = 0;
        }
        else
        {
            newFile = 1;
        }

        if(pos > (*inList)[inListPos])
        {
            inListPos++;
        }

        if(inListPos >= listSize)
        {
            *posMax = pos;
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: no position found\n");
                return 0;
            }

            char outputFileNew[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("Creating File: %s, 100%%\n",outputFile);
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    startingValue,(*counter),
                    firstPos,pos);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }


    while(1)
    {
        index = 0;
        if(newFile == 1)
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;
                if(!strcmp(allID,alignmentId))
                {
                    if(pos == (*inList)[inListPos])
                    {
                        (*counter)++;
                        sprintf(outputFile,"%s%s_%d_%d_%d_.vcf.gz",outputPathName,
                                alignmentId,
                                (*counter),
                                (*counter)+(*maxcount)-1,
                                pos);
                        startingValue = (*counter);
                        endingValue = (*counter)+(*maxcount)-1;
                        firstPos = pos;
                        if(*posMin == 0)
                        {
                            *posMin = pos;
                        }

                        printf("\033[A\33[2K");
                        printf("Creating File: %s, %d%%\n",outputFile,
                                ((*counter)*100)/listSize);

                        fpOut = gzopen(outputFile,"w");
                        if(fpOut == NULL)
                        {
                            fprintf(stderr,"\n ERROR: failed to open file %s",outputFile);
                            return 0;
                        }
                        writeLine(fpOut,headerLine1);
                        writeLine(fpOut,headerLine2);
                        writeLine(fpOut,line);
                        (*inList)[inListPos] = 0;
                        inListPos++;

                        newFile = 0;
                    }
                    else
                    {
                        newFile = 1;
                    }

                    if(pos > (*inList)[inListPos])
                    {
                        inListPos++;
                    }

                    if(inListPos >= listSize)
                    {
                        *posMax = pos;
                        gzclose(fpOut);
                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: no position found\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        printf("Creating File: %s, 100%%\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }else if(eof ==1)
            {
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
        else
        {
            status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
            if(status == 1 && eol==1)
            {
                status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
                assert(status);
                char * allID = (char*)malloc(strlen(*word)*sizeof(char));
                strcpy(allID,*word);
                status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
                assert(status);
                pos = atoi(*word);
                lineCounter++;

                if((((*counter)*100)/listSize)%5 == 0 &&
                        printDelay != (((*counter)*100)/listSize))
                {
                    printf("\033[A\33[2K");
                    printf("Creating File: %s, %d%%\n",outputFile,
                            ((*counter)*100)/listSize);
                    printDelay = (((*counter)*100)/listSize);
                }

                if(!strcmp(allID,alignmentId))
                {
                    if (pos == (*inList)[inListPos])
                    {
                        writeLine(fpOut,line);
                        (*counter)++;
                        (*inList)[inListPos] = 0;
                        inListPos++;
                    }
                    else if(pos > (*inList)[inListPos])
                        inListPos++;

                    if(inListPos >= listSize)
                    {
                        *posMax = pos;
                        gzclose(fpOut);
                        if((*counter) == 0)
                        {
                            fprintf(stderr,"\n ERROR: no position found\n");
                            return 0;
                        }

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("Creating File: %s, 100%%\n",outputFile);
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                        printf("%d Snips written\n",(*counter));
                        return 1;
                    }
                    else if((*counter) == endingValue)
                    {
                        newFile = 1;
                        gzclose(fpOut);

                        char outputFileNew[INFILENAMESIZE+30];
                        printf("\033[A\33[2K");
                        sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                                alignmentId,
                                startingValue,
                                (*counter),
                                firstPos,pos);
                        printf("Renaming %s to %s\n",outputFile,outputFileNew);
                        int ret = rename(outputFile, outputFileNew);
                        if(ret != 0)
                        {
                            fprintf(stderr,"\n ERROR: unable to rename the file\n");
                            return 0;
                        }
                    }
                }
                else
                //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
                //            "skipped\n",lineCounter,allID);
                free(allID);
            }
            else if(eof ==1)
            {
                gzclose(fpOut);

                if((*counter) == 0)
                {
                    fprintf(stderr,"\n ERROR: no position found\n");
                    return 0;
                }
                char outputFileNew[INFILENAMESIZE+30];
                sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                        alignmentId,
                        startingValue,
                        (*counter),firstPos,pos);
                printf("Renaming %s to %s\n",outputFile,outputFileNew);
                int ret = rename(outputFile, outputFileNew);
                if(ret != 0)
                {
                    fprintf(stderr,"\n ERROR: unable to rename the file\n");
                    return 0;
                }
                printf("%d Snips written\n",(*counter));
                return 1;
            }
        }
    }
    printf("%d Snips written\n",(*counter));
    return 1;
}

int createSingleFileIndexW(gzFile fpIn,
        char* outputPathName,
        int *counter,
        char **line,
        int* lineLength,
        char **word,
        int* wordLength,
        int lines,
        int Wmin,
        int Wmax)
{
    int eol=0,
        eof=0,
        status,
        index = 0,
        snipSize = 0,
        lineCounter = 0,
        printDelay = 0,
        pos,
        firstPos=0,
        firstIndex=0;

    char alignmentId[INFILENAMESIZE];
    char outputFile[INFILENAMESIZE+30];
    gzFile fpOut;

    sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
    printf("Creating Single File: %s\n",outputFile);
    fpOut = gzopen(outputFile,"w");
    if(fpOut == NULL)
    {
        fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
        return 0;
    }

    while(1)
    {
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            writeLine(fpOut,line);
            if(strlen(*line) >=6)
            {
                if((*line)[0] == '#' &&
                        (*line)[1] == 'C' &&
                        (*line)[2] == 'H' &&
                        (*line)[3] == 'R' &&
                        (*line)[4] == 'O' &&
                        (*line)[5] == 'M')
                    break;
            }
        }
        else if(eof==1)
        {
            fprintf(stderr,"\n ERROR: Header malformed\n");
            return 0;
        }
    }
    printf("Header copied\n");

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                snipSize++;
            index++;
        }
        lineCounter++;
        printf("Total snips: %d, Snip Size: %d bits, Wmin: %d, Wmax: %d\n\n\n",lines,
                snipSize,
                Wmin,
                Wmax);

        if(lineCounter >= Wmin)
        {
            writeLine(fpOut,line);
            (*counter)++;
            if(firstPos == 0)
            {
                firstPos = pos;
                firstIndex = lineCounter;
            }
        }

        printf("Processing: %d%%\n",((*counter)*100)/(Wmax-Wmin+1));

        if(lineCounter >= Wmax)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }

            char outputFileNew[INFILENAMESIZE+30];
            char outputFile[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    firstIndex,
                    firstIndex+(*counter)-1,
                    firstPos,pos);
            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }

    while(1)
    {
        index = 0;
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
            assert(status);
            char * allID = (char*)malloc(strlen(*word)*sizeof(char));
            strcpy(allID,*word);
            status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
            assert(status);
            pos = atoi(*word);
            lineCounter++;

            if((((*counter)*100)/(Wmax-Wmin+1))%5 == 0 &&
                    printDelay != ((*counter)*100)/(Wmax-Wmin+1))
            {
                printf("\033[A\33[2K");
                printf("Processing: %d%%\n",((*counter)*100)/(Wmax-Wmin+1));
                printDelay = (((*counter)*100)/(Wmax-Wmin+1));
            }

            if(!strcmp(allID,alignmentId))
            {
                if(lineCounter >= Wmin)
                {
                    writeLine(fpOut,line);
                    (*counter)++;
                    if(firstPos == 0)
                    {
                        firstPos = pos;
                        firstIndex = lineCounter;
                    }
                }

                if(lineCounter >= Wmax)
                {
                    gzclose(fpOut);
                    if((*counter) == 0)
                    {
                        fprintf(stderr,"\n ERROR: No snips found inside window\n");
                        return 0;
                    }
                    char outputFileNew[INFILENAMESIZE+30];
                    char outputFile[INFILENAMESIZE+30];
                    printf("\033[A\33[2K");
                    printf("Processing: 100%%\n");
                    sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                            alignmentId,
                            firstIndex,
                            firstIndex+(*counter)-1,
                            firstPos,pos);
                    sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
                    printf("Renaming %s to %s\n",outputFile,outputFileNew);
                    int ret = rename(outputFile, outputFileNew);
                    if(ret != 0)
                    {
                        fprintf(stderr,"\n ERROR: unable to rename the file\n");
                        return 0;
                    }
                    printf("%d Snips written\n",(*counter));
                    return 1;
                }
            }
            else
            //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be "
            //            "skipped\n",lineCounter,allID);
            free(allID);
        }
        else if(eof ==1)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }
            char outputFileNew[INFILENAMESIZE+30];
            char outputFile[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    firstIndex,
                    firstIndex+(*counter)-1,
                    firstPos,pos);
            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
}

int createSingleFilePosW(gzFile fpIn,
        char* outputPathName,
        int *counter,
        char **line,
        int* lineLength,
        char **word,
        int* wordLength,
        int lines,
        int posWmin,
        int posWmax)
{
    int eol=0,
        eof=0,
        status,
        index = 0,
        snipSize = 0,
        lineCounter = 0,
        printDelay = 0,
        printCount = 0,
        pos,
        firstPos=0,
        firstIndex=0,
        i;

    char alignmentId[INFILENAMESIZE];
    char outputFile[INFILENAMESIZE+30];
    gzFile fpOut;

    sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
    printf("Creating Single File: %s\n",outputFile);
    fpOut = gzopen(outputFile,"w");
    if(fpOut == NULL)
    {
        fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
        return 0;
    }

    while(1)
    {
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            writeLine(fpOut,line);
            if(strlen(*line) >=6)
            {
                if((*line)[0] == '#' &&
                        (*line)[1] == 'C' &&
                        (*line)[2] == 'H' &&
                        (*line)[3] == 'R' &&
                        (*line)[4] == 'O' &&
                        (*line)[5] == 'M')
                    break;
            }
        }
        else if(eof==1)
        {
            fprintf(stderr,"\n ERROR: Header malformed\n");
            exit(0);
        }
    }
    printf("Header copied\n");

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                snipSize++;
            index++;
        }
        lineCounter++;
        printf("Total snips: %d, Snip Size: %d bits, posWmin: %d, posWmax: %d\n\n\n",
                lines,snipSize, posWmin, posWmax);

        if(pos >= posWmin)
        {
            if(firstPos == 0)
            {
                firstPos = pos;
                firstIndex = lineCounter;
            }
            writeLine(fpOut,line);
            (*counter)++;
        }

        printf("Processing");
        for(i=0;i<=printCount;i++)
            printf(".");
        printf("\n");
        printCount++;
        if(printCount == 5)
            printCount = 0;

        if(pos >= posWmax)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }

            char outputFileNew[INFILENAMESIZE+30];
            char outputFile[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    firstIndex,
                    firstIndex+(*counter)-1,
                    firstPos,pos);
            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
            return 1;
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }

    while(1)
    {
        index = 0;
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
            assert(status);
            char * allID = (char*)malloc(strlen(*word)*sizeof(char));
            strcpy(allID,*word);
            status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
            assert(status);
            pos = atoi(*word);
            lineCounter++;

            if((((*counter)*100)/lines)%5 == 0 && printDelay != (((*counter)*100)/lines))
            {
                printf("\033[A\33[2K");
                printf("Processing");
                for(i=0;i<=printCount;i++)
                    printf(".");
                printf("\n");
                printCount++;
                if(printCount == 5)
                    printCount = 0;
                printDelay = (((*counter)*100)/lines);
            }

            if(!strcmp(allID,alignmentId))
            {
                if(pos >= posWmin)
                {
                    if(firstPos == 0)
                    {
                        firstPos = pos;
                        firstIndex = lineCounter;
                    }
                    writeLine(fpOut,line);
                    (*counter)++;
                }

                if(pos >= posWmax)
                {
                    gzclose(fpOut);
                    if((*counter) == 0)
                    {
                        fprintf(stderr,"\n ERROR: No snips found inside window\n");
                        return 0;
                    }
                    // hit overflow warning, temporarily grew the buffer
                    // MPAMPIS
                    char outputFileNew[4*INFILENAMESIZE+30];
                    char outputFile[INFILENAMESIZE+30];
                    printf("\033[A\33[2K");
                    printf("Processing: 100%%\n");
                    sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,
                            alignmentId,
                            firstIndex,
                            firstIndex+(*counter)-1,
                            firstPos,pos);
                    sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
                    printf("Renaming %s to %s\n",outputFile,outputFileNew);
                    int ret = rename(outputFile, outputFileNew);
                    if(ret != 0)
                    {
                        fprintf(stderr,"\n ERROR: unable to rename the file\n");
                        return 0;
                    }
                    return 1;
                }
            }
            else
            //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be skipped\n",
            //            lineCounter,allID);
            free(allID);
        }
        else if(eof ==1)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }
            // hit overflow warning, temporarily grew the buffer
            // MPAMPIS
            char outputFileNew[4*INFILENAMESIZE+30];
            char outputFile[INFILENAMESIZE+30];
            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            sprintf(outputFileNew,"%s%s_%d_%d_%d_%d.vcf.gz",outputPathName,alignmentId,
                    firstIndex,
                    firstIndex+(*counter)-1,
                    firstPos,pos);
            sprintf(outputFile,"%stemp.vcf.gz",outputPathName);
            printf("Renaming %s to %s\n",outputFile,outputFileNew);
            int ret = rename(outputFile, outputFileNew);
            if(ret != 0)
            {
                fprintf(stderr,"\n ERROR: unable to rename the file\n");
                return 0;
            }
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
}

int createSingleFilePosL(gzFile fpIn,
        char* outputPathName,
        int *counter,
        char **line,
        int* lineLength,
        char **word,
        int* wordLength,
        int lines,
        int ** inList,
        int listSize)
{
    int eol=0,
        eof=0,
        status,
        index = 0,
        snipSize = 0,
        lineCounter = 0,
        printDelay = 0,
        pos,
        firstPos=0,
        inListPos=0;

    char alignmentId[INFILENAMESIZE];
    char outputFile[INFILENAMESIZE+30];
    gzFile fpOut;

    sprintf(outputFile,"%spositions_from_list.vcf.gz",outputPathName);
    printf("Creating Single File: %s\n",outputFile);
    fpOut = gzopen(outputFile,"w");
    if(fpOut == NULL)
    {
        fprintf(stderr,"\n ERROR: failed to open file %s", outputFile);
        return 0;
    }

    while(1)
    {
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            writeLine(fpOut,line);
            if(strlen(*line) >=6)
            {
                if((*line)[0] == '#' &&
                        (*line)[1] == 'C' &&
                        (*line)[2] == 'H' &&
                        (*line)[3] == 'R' &&
                        (*line)[4] == 'O' &&
                        (*line)[5] == 'M')
                    break;
            }
        }
        else if(eof==1)
        {
            fprintf(stderr,"\n ERROR: Header malformed\n");
            exit(0);
        }
    }
    printf("Header copied\n");

    status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
    if(status == 1 && eol==1)
    {
        status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
        assert(status);
        //strcpy(alignmentId,*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
        assert(status);
        pos = atoi(*word);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//id
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//ref
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//alt
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//qual
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//filter
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//info
        assert(status);
        status = getWordFromString(*line,word, &eol, wordLength, &index);//format
        assert(status);
        while(index < strlen(*line))
        {
            if ((*line)[index] == '1' || (*line)[index] == '0' || (*line)[index] == '.')
                snipSize++;
            index++;
        }
        lineCounter++;
        printf("Total snips: %d, Snip Size: %d bits, using input List\n\n\n",lines,
                snipSize);

        if(pos == (*inList)[inListPos])
        {
            if(firstPos == 0)
                firstPos = pos;
            writeLine(fpOut,line);
            (*counter)++;
            (*inList)[inListPos] = 0;
            inListPos++;
        }
        else if(pos > (*inList)[inListPos])
            inListPos++;

        printf("Processing: %d%%\n",((*counter)*100)/listSize);

        if(inListPos >= listSize)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: no position found\n");
                return 0;
            }

            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
    else
    {
        fprintf(stderr,"\n ERROR: no snips found\n");
        return 0;
    }

    while(1)
    {
        index = 0;
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1)
        {
            status = getWordFromString(*line,word, &eol, wordLength, &index);//chrom
            assert(status);
            char * allID = (char*)malloc(strlen(*word)*sizeof(char));
            strcpy(allID,*word);
            status = getWordFromString(*line,word, &eol, wordLength, &index);//pos
            assert(status);
            pos = atoi(*word);
            lineCounter++;

            if((((*counter)*100)/listSize)%5 == 0 &&
                    printDelay != (((*counter)*100)/listSize))
            {
                printf("\033[A\33[2K");
                printf("Processing: %d%%\n",((*counter)*100)/listSize);
                printDelay = (((*counter)*100)/listSize);
            }

            if(!strcmp(allID,alignmentId))
            {
                if(pos == (*inList)[inListPos])
                {
                    if(firstPos == 0)
                        firstPos = pos;
                    writeLine(fpOut,line);
                    (*counter)++;
                    (*inList)[inListPos] = 0;
                    inListPos++;
                }
                else if(pos > (*inList)[inListPos])
                    inListPos++;

                if(inListPos >= listSize || lineCounter >= lines)
                {
                    gzclose(fpOut);
                    if((*counter) == 0)
                    {
                        fprintf(stderr,"\n ERROR: no position found\n");
                        return 0;
                    }
                    printf("\033[A\33[2K");
                    printf("Processing: 100%%\n");
                    printf("%d Snips written\n",(*counter));
                    return 1;
                }
            }
            else
            //    fprintf(stderr,"WARNING: snip %d is alignment %s and will be skipped\n",
            //            lineCounter,allID);
            free(allID);
        }
        else if(eof ==1)
        {
            gzclose(fpOut);
            if((*counter) == 0)
            {
                fprintf(stderr,"\n ERROR: No snips found inside window\n");
                return 0;
            }
            printf("\033[A\33[2K");
            printf("Processing: 100%%\n");
            printf("%d Snips written\n",(*counter));
            return 1;
        }
    }
}
