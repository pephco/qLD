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

int getFileFormat (gzFile fp)
{
    
    char tmp;
    tmp=gzgetc(fp);
    char macsFLAG[] = "COMMAND";
    int flaglength = 7,j;
    char vcfFLAG[] = "##fileformat=VCF";
    int flaglengthVCF = 16;
    int format = OTHER_FORMAT;
    while(tmp!=EOF)
    {
        if(tmp=='/')
        {
            tmp = gzgetc(fp);
            
            if(tmp=='/')
            {
                format = MS_FORMAT;
                break;
            }
            else
                tmp = gzgetc(fp);				
        }
        else
        {
            if(tmp=='>')
            {
                format = FASTA_FORMAT;
                break;
            }
            else
            {
                int counter = 0;
                
                while(counter < flaglength)
                {
                    if(tmp != macsFLAG[counter])
                        break;
                    
                    tmp = gzgetc(fp);
                    
                    ++counter;
                }
                
                j = counter;
                
                
                if(j == flaglength)
                {
                    format = MACS_FORMAT;
                    break;
                }
                
                else
                {
                    
                    gzseek(fp, -j - 1, SEEK_CUR);
                    
                    tmp = gzgetc( fp);
                    
                    int counter = 0;
                    
                    while(counter < flaglengthVCF)
                    {
                        if(tmp != vcfFLAG[counter])
                            break;
                        
                        tmp = gzgetc(fp);
                        
                        ++counter;
                    }
                    
                    j = counter;
                    
                    if(j == flaglengthVCF)
                    {
                        format = VCF_FORMAT;
                        break;
                    }
                    
                    else
                        tmp = gzgetc(fp);
                }
            }
        }		
    }
    
    return format;
}

void printHelp (FILE *fp)
{
    fprintf(fp,
    "VCF_parser manual\n"
    "-----------------\n"
    "\t-input       inputFile\n"
    "\t-output      outputDir\n"
    "\t-size        parts_size\n"
    "\t-Wmin        snip index\n"
    "\t-Wmax        snip index\n"
    "\t-posWmin     snip pos\n"
    "\t-posWmax     snip pos\n"
    "\t-inputList   inputFile\n"
    "\t-chrom       chromosome\n"
    "\t-toSingleOutput\n"
    "\nDescription:\n"
    "\t-input     <STRING>  Specifies the name of the input alignment file.\n"
    "\t-output    <STRING>  Specifies the path of the output alignment files.\n"
    "\t-size      <INT>     Specifies the size of the memory footprint of the output\n"
    "\t                     alignment files in MB, if toSingleOutput is set,\n"
    "\t                     size is not needed. Supported file formats: VCF.\n\n"
    "\t-Wmin      <INT>     index of the minimum snip to be included, minimum 1 (default)\n"
    "\t-Wmax      <INT>     index of the maximum snip to be included,\n"
    "\t                     maximum total-Snips (default)\n"
    "\t-posWmin   <INT>     pos of the minimum snip to be included, must be valid\n"
    "\t-posWmax   <INT>     pos of the maximum snip to be included, must be valid\n"
    "\t-inputList <STRING>  input text file with the pos to keep\n"
    "\t-chrom     <STRING>  Specifies the chromosome to be extracted from the original VCF\n"
    "\t-toSingleOutput      Used to generate a new VCF that is part of the input file,\n"
    "\t                     -Wmin and -Wmax mandatory with this command\n"
    "\n\n");
}

void commandLineParser(int argc, char** argv, 
                       char * infile, 
                       char * outPath, 
                       int * fileFormat,
                       int * size,
                       int * Wmin,
                       int * Wmax,
                       int *Wset,
                       int * posWmin,
                       int * posWmax,
                       int *posWset, 
                       char * inList, 
                       int * inListSet,
                       int * toSingleOutput,
                       char * alignmentId)
{
    int i, pathSet = 0, fileSet=0, sizeSet=0, CHset=0;
    gzFile fp;
    
    for(i=1; i<argc; ++i)
    {
        if(!strcmp(argv[i], "-input")) 
        { 
            if (i!=argc-1)
            {			
                strcpy(infile,argv[++i]);
                
                fp=gzopen(infile,"r");
                
                if (fp==NULL)
                {
                    fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",infile);
                    exit(0);
                }
                else
                {
                    *fileFormat = getFileFormat (fp);
                    if (*fileFormat!=VCF_FORMAT)
                    {
                        fprintf(stderr, "\n ERROR: File %s is not VCF format.\n\n",
                                infile);
                        exit(0);
                    }
                    
                    gzclose(fp);
                    
                    fileSet=1;
                }
            }
            continue;
        }
        if(!strcmp(argv[i], "-inputList")) 
        { 
            if (i!=argc-1)
            {			
                strcpy(inList,argv[++i]);
                
                fp=gzopen(inList,"r");
                
                if (fp==NULL)
                {
                    fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",inList);
                    exit(0);
                }
                else
                {
                    gzclose(fp);
                    *inListSet=1;
                }
            }
            continue;
        }
        if(!strcmp(argv[i], "-output")) 
        {
            if (i!=argc-1)
            {			
    			char outputPathName[INFILENAMESIZE];
                strcpy(outputPathName,argv[++i]);
                struct stat st = {0};
                
                if (stat(outputPathName, &st) == -1) 
                {
                    int status = mkdir(outputPathName, 0700);
                    if (status == -1) 
                    {
                        fprintf(stderr, "\n ERROR: Directory %s could not be created.\n",
                                outputPathName);
                        exit(0);
                    }
                    else
                        fprintf(stdout, "\n Directory %s created.\n",outputPathName);
                }
                else 
                {
                    fprintf(stdout, "\n Directory %s exists, renaming.\n",outputPathName);
                    int pos = strlen(outputPathName)-1;
                    if(outputPathName[pos] == '/')
                        outputPathName[pos] = '\0';
                    // hit overflow warnings @sprintf, changed the size 
                    // to mute them temporarily
                    // MPAMPIS
                    int i=1;
    				char outputPathNameNew[INFILENAMESIZE+4*sizeof(int)];
                    do
                    {
                    	sprintf(outputPathNameNew,"%s_%d",outputPathName,i);
                    	i++;
                    }while(stat(outputPathNameNew, &st) != -1);

                	strcpy(outputPathName,outputPathNameNew);
                    int status = mkdir(outputPathName, 0700);
                    if (status == -1) 
                    {
                        fprintf(stderr, "\n ERROR: Directory %s could not be created.\n",
                                outputPathName);
                        exit(0);
                    }
                    else
                        fprintf(stdout, "\n Directory %s created.\n",outputPathName);
                }
                int pos = strlen(outputPathName)-1;
                if(outputPathName[pos] != '/') {
                    outputPathName[pos+1] = '/';
                    outputPathName[pos+2] = '\0';
                }
                strcpy(outPath,outputPathName);
                pathSet=1;
            }
            continue;
        }
        if(!strcmp(argv[i], "-size"))
        {
            if (i!=argc-1)
            {
                *size = atoi(argv[++i]);
                fprintf(stdout, "\n Part size: %d MB.\n",*size);
                
                if(*size!=0)
                    sizeSet=1;
            } 
            continue;
        }

        if(!strcmp(argv[i], "-Wmin")) 
        {
            if (i!=argc-1)
            {           
                *Wmin = atoi(argv[++i]);
                *Wset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-Wmax")) 
        {
            if (i!=argc-1)
            {           
                *Wmax = atoi(argv[++i]);
                *Wset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-posWmin")) 
        {
            if (i!=argc-1)
            {           
                *posWmin = atoi(argv[++i]);
                *posWset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-posWmax")) 
        {
            if (i!=argc-1)
            {           
                *posWmax = atoi(argv[++i]);
                *posWset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-chrom")) 
        {
            if (i!=argc-1)
            {           
                strcpy(alignmentId,argv[++i]);
                CHset = 1;
            }
            continue;
        }

        if(!strcmp(argv[i], "-toSingleOutput")) 
        {
        	*toSingleOutput = 1;
        	continue;
        }
        
        if(!strcmp(argv[i], "-help")||!strcmp(argv[i], "-h")||!strcmp(argv[i], "--help"))
        {             
            printHelp (stdout);
            
            exit(0);
        }
        
        fprintf(stderr, "\n ERROR: %s is not a valid command line parameter\n",argv[i]);
        exit(0);
    }
    
    if (pathSet==0)
    {
        fprintf(stderr, "\n ERROR: Please specify a path for the output with -output\n");
        exit(0);
    }
    
    if (CHset==0)
    {
        fprintf(stderr, "\n ERROR: Please specify a chromosome with -chrom\n");
        exit(0);
    }

    if (fileSet==0)
    {
        fprintf(stderr, "\n ERROR: Please specify an alignment with -input\n");
        exit(0);
    }
    
    if (( sizeSet==0) && *toSingleOutput == 0)
    {
        fprintf(stderr, "\n ERROR: Please specify a size in MB for the output with \
-size\n");
        exit(0);
    }

    if(*toSingleOutput == 1 && (*Wset)==0 && (*posWset)==0 && (*inListSet)==0)
    {
        fprintf(stderr, "\n ERROR: Please specify an index window with -Wmin and/or Wmax,\
or a pos window with -posWmin and/or -posWmax, or an input List with -inputList\n");
        exit(0);
    }

    if(((*Wset)==1 && (*posWset)==1) || 
       ((*Wset)==1 && (*inListSet)==1) || 
       ((*inListSet)==1 && (*posWset)==1))
    {
        fprintf(stderr, "\n ERROR: Please specify only an index window, or a pos window, \
or an input List\n");
        exit(0);
    }
}

int main(int argc, char** argv)
{
    int fileFormat=OTHER_FORMAT,
    	size=0,
    	Wmin=0,
    	Wmax=0,
    	Wset = 0,
    	posWmin=0,
    	posWmax=0,
    	posWset = 0,
    	inListSet = 0,
    	toSingleOutput=0,
    	lineLength = STRINGLENGTH,
    	wordLength = STRINGLENGTH,
    	lines = 0,
    	listSize = 0,
    	posMin = 0,
    	posMax = 0,
    	status=0,
    	maxcount=0,
    	counter=0,
    	snipSize=0;

    double time1,time2,time3;

    char inputFileName[INFILENAMESIZE];
    char inputListName[INFILENAMESIZE];
    char outputPathName[INFILENAMESIZE];
    char alignmentId[INFILENAMESIZE];
    char headerFileOld[INFILENAMESIZE+30];
    // hit overflow warnings @sprintf, changed the size 
    // to mute them temporarily
    // MPAMPIS
    char headerFileNew[2*INFILENAMESIZE+30];

    gzFile fpIn=NULL, fpHeader=NULL;
    char ** headerLine1, ** headerLine2;

    int **inList = (int**)malloc(sizeof(int*));
    *inList = (int*)malloc(sizeof(int));

    commandLineParser(argc, 
                      argv, 
                      inputFileName, 
                      outputPathName, 
                      &fileFormat, 
                      &size, 
                      &Wmin, 
                      &Wmax, 
                      &Wset, 
                      &posWmin, 
                      &posWmax, 
                      &posWset,
                      inputListName,
                      &inListSet, 
                      &toSingleOutput,
                      alignmentId);

    char **line = (char**)malloc(sizeof(char*));
    *line = (char*)malloc(sizeof(char)*STRINGLENGTH);
    char **word = (char**)malloc(sizeof(char*));
    *word = (char*)malloc(sizeof(char)*STRINGLENGTH);

    fpIn = gzopen(inputFileName,"r");
    time1 = gettime();
    lines = countLines(fpIn, line, &lineLength);
    gzclose(fpIn);
    if(Wset == 1) 
    {
	    if(Wmin < 1 || Wmin > lines) 
        {
	    	fprintf(stderr,"\n WARNING: Wmin must be between 1 and totalSnips, setting \
Wmin to 1 (default)\n");
	    	Wmin = 1;
	    }

	    if(Wmax < 1 || Wmax > lines) 
        {
	    	fprintf(stderr,"\n WARNING: Wmax must be between 1 and totalSnips, setting \
Wmax to total Snips (default)\n");
	    	Wmax = lines;
	    }

	    if(Wmin > Wmax) 
        {
	    	fprintf(stderr,"\n WARNING: Wmin must be <= Wmax, setting Wmin = Wmax\n");
	    	Wmin = Wmax;
	    }
	}
	else if(posWset == 1) 
    {

	    if(posWmin < 1) 
        {
	    	fprintf(stderr,"\n WARNING: posWmin must be between 1 and 2147483646 (max \
positive int), setting posWmin to 1 (default)\n");
	    	posWmin = 1;
	    }

	    if(posWmax < 1) 
        {
	    	fprintf(stderr,"\n WARNING: posWmax must be between 1 and 2147483646 (max \
positive int), setting posWmax to 2147483646 (default)\n");
	    	posWmax = 2147483646;
	    }

	    if(posWmin > posWmax) 
        {
	    	fprintf(stderr,"\n WARNING: posWmin must be <= posWmax, setting posWmin = \
posWmax\n");
	    	posWmin = posWmax;
	    }
	}
	else if(inListSet == 1) 
    {
		listSize = sortList(inputListName,inList, line, &lineLength);
	}


    fpIn = gzopen(inputFileName,"r");
    time2 = gettime();
    if(toSingleOutput == 0) 
    {
    	printf("Creating Multiple Files\n");
    	headerLine1 = (char **) malloc(sizeof(char*));
    	headerLine2 = (char **) malloc(sizeof(char*));
	    fpHeader = createHeaderFile(fpIn, 
                                    outputPathName, 
                                    headerLine1, 
                                    headerLine2, 
                                    line, 
                                    &lineLength, 
                                    word, 
                                    &wordLength);
	    
	    if(Wset == 1)
	    	status = createPartIndexW(fpIn, 
                                      outputPathName, 
                                      headerLine1, 
                                      headerLine2, 
                                      line, 
                                      &lineLength, 
                                      word, 
                                      &wordLength, 
                                      lines,
                                      size, 
                                      &snipSize, 
                                      alignmentId, 
                                      &maxcount,
                                      &counter, 
                                      &posMin, 
                                      &posMax,
                                      Wmin,
                                      Wmax);
	    else if(posWset == 1)
			status = createPartPosW(fpIn, 
                                    outputPathName, 
                                    headerLine1, 
                                    headerLine2, 
                                    line, 
                                    &lineLength, 
                                    word, 
                                    &wordLength, 
                                    lines,
                                    size, 
                                    &snipSize, 
                                    alignmentId, 
                                    &maxcount, //-Wint-conversion MPAMPIS
                                    &counter, 
                                    &posMin, 
                                    &posMax,posWmin,posWmax);
	    else if(inListSet == 1) 
        {
			status = createPartPosL(fpIn, 
                                    outputPathName, 
                                    headerLine1, 
                                    headerLine2, 
                                    line, 
                                    &lineLength, 
                                    word, 
                                    &wordLength, 
                                    lines,
                                    size, 
                                    &snipSize, 
                                    alignmentId, 
                                    &maxcount,
                                    &counter, 
                                    &posMin, 
                                    &posMax, 
                                    inList,
                                    listSize);
	    	printf("\n");
			for(inListSet=0;inListSet<listSize;inListSet++)
				if((*inList)[inListSet] != 0)
	        		fprintf(stderr,"\tWARNING: list pos %d not found\n",
                            (*inList)[inListSet]);
			printf("\n\n");
	    }
	    else
	    	status = createPart(fpIn, 
                                outputPathName, 
                                headerLine1, 
                                headerLine2, 
                                line, 
                                &lineLength, 
                                word, 
                                &wordLength, 
                                lines,
                                size, 
                                &snipSize, 
                                alignmentId, 
                                &maxcount,
                                &counter, 
                                &posMin, 
                                &posMax);

	    if(status == 0) 
        {
	    	fprintf(stderr,"\n ERROR encountered, abnormal termination\n");
	    	return 0;
	    }

	    if(maxcount <0 || maxcount>counter)
	    	maxcount = counter;
	    gzprintf(fpHeader,"%s %d %d %d %d %d\n",alignmentId, maxcount, snipSize, 
                                                counter, posMin, posMax);
	    printf("%s %d %d %d %d %d\n",alignmentId, maxcount, snipSize, counter, posMin, 
                                     posMax);
	    gzclose(fpHeader);
	    sprintf(headerFileNew,"%s%s_header.vcf.gz",outputPathName,alignmentId);

	    sprintf(headerFileOld,"%sheader.vcf.gz",outputPathName);
		printf("Renaming %s to %s\n",headerFileOld, headerFileNew);
	    int ret = rename(headerFileOld, headerFileNew);
		if(ret != 0) 
        {
			fprintf(stderr,"\n ERROR: unable to rename the file\n");
	        return 0;
		}
    	free(*headerLine1);
    	free(*headerLine2);
    	free(headerLine1);
    	free(headerLine2);
	}
	else if(toSingleOutput==1 && Wset ==1)
        status = createSingleFileIndexW(fpIn, 
                                        outputPathName, 
                                        &counter, 
                                        line, 
                                        &lineLength, 
                                        word, 
                                        &wordLength, 
                                        lines, 
                                        Wmin, 
                                        Wmax);
	else if(toSingleOutput==1 && posWset ==1)
        status = createSingleFilePosW(fpIn, 
                                      outputPathName, 
                                      &counter, 
                                      line, 
                                      &lineLength, 
                                      word, 
                                      &wordLength, 
                                      lines, 
                                      posWmin, 
                                      posWmax);
	else if(toSingleOutput==1 && inListSet ==1) 
    {
        status = createSingleFilePosL(fpIn, 
                                      outputPathName, 
                                      &counter, 
                                      line, 
                                      &lineLength, 
                                      word, 
                                      &wordLength, 
                                      lines, 
                                      inList, 
                                      listSize);
		printf("\n");
		for(inListSet=0;inListSet<listSize;inListSet++)
			if((*inList)[inListSet] != 0)
        		fprintf(stderr,"\tWARNING: list pos %d not found\n",(*inList)[inListSet]);
		printf("\n\n");
	}
		
	if(status == 0) 
    {
    	fprintf(stderr,"\n ERROR encountered, abnormal termination\n");
    	return 0;
    }
    gzclose(fpIn);
    free(*line);
    free(line);
    free(*word);
    free(word);
    free(*inList);
    free(inList);
    time3 = gettime();
    printf("get lines: %fs, create Parts: %fs, Total Time: %fs \n",time2 - time1, 
                                                                   time3-time2, 
                                                                   time3-time1);
    return 1;
}
