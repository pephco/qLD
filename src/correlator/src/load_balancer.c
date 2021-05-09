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
#include "../include/load_balancer.h"
#include "../include/correlate_IO.h"

// Takes two lists sorted in decreasing order, and merge their nodes
// together to make one big sorted list which is returned
t_node* SortedMerge(t_node* a, t_node* b)
{
    // Base cases
    if (a == NULL)
        return b;

    else if (b == NULL)
        return a;

    t_node* result = NULL;

    // Pick either a or b, and recur
    if((a->filesListNum + a->filesListNum2) > (b->filesListNum + b->filesListNum2))
    {
        result = a;
        result->next = SortedMerge(a->next, b);
    }
    else if((a->filesListNum + a->filesListNum2) < (b->filesListNum + b->filesListNum2))
    {
        result = b;
        result->next = SortedMerge(a, b->next);
    }
    else
    {
        if((a->posWmax - a->posWmin) >= (b->posWmax - b->posWmin))
        {
            result = a;
            result->next = SortedMerge(a->next, b);
        }
        else
        {
            result = b;
            result->next = SortedMerge(a, b->next);
        }
    }

    return result;
}

/*
   Split the nodes of the given list into front and back halves,
   and return the two lists using the reference parameters.
   If the length is odd, the extra node should go in the front list.
   It uses the fast/slow pointer strategy
*/
void FrontBackSplit(t_node* source, t_node** frontRef,
        t_node** backRef)
{
    // if length is less than 2, handle separately
    if (source == NULL || source->next == NULL)
    {
        *frontRef = source;
        *backRef = NULL;
        return;
    }

    t_node* slow = source;
    t_node* fast = source->next;

    // Advance 'fast' two nodes, and advance 'slow' one node
    while (fast != NULL)
    {
        fast = fast->next;
        if (fast != NULL)
        {
            slow = slow->next;
            fast = fast->next;
        }
    }

    // 'slow' is before the midpoint in the list, so split it in two
    // at that point.
    *frontRef = source;
    *backRef = slow->next;
    slow->next = NULL;
}

// Sort given linked list using Merge sort algorithm
void MergeSort(t_node** head)
{
    // Base case -- length 0 or 1
    if (*head == NULL || (*head)->next == NULL)
        return;

    t_node* a;
    t_node* b;

    // Split head into 'a' and 'b' sublists
    FrontBackSplit(*head, &a, &b);

    // Recursively sort the sublists
    MergeSort(&a);
    MergeSort(&b);

    // answer = merge the two sorted lists together
    *head = SortedMerge(a, b);
}

void preprocess_data(FILE *fpInRep,
        char *inputPathName,
        char *inputPath2Name,
        char *outputFileName,
        sample_t *sampleList,
        sample_t *sampleList2,
        int inPath2set,
        int posWset1,
        int posWset2,
        int posWmin1,
        int posWmax1,
        int posWmin2,
        int posWmax2,
        int mdf,
        int *task_count)
{
    int totalSnips=0, totalSnips2=0, snipsPerFile=0;
    int snipsPerFile2=0, posMin=0, posMax=0, posMin2=0, posMax2=0;

    char alignmentId[INFILENAMESIZE];
    char alignmentId2[INFILENAMESIZE];

    int filesListNum=0, filesListNum2=0;
    int snipSize=0, snipSize2=0;

    char *headerLine1=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerLine1);
    char *headerLine2=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerLine2);
    char *headerLine1_2=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerLine1_2);
    char *headerLine2_2=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(headerLine2_2);

    char **filesList=(char **)malloc(sizeof(char *));
    assert(filesList);
    *filesList=NULL;
    char **filesList2=(char **)malloc(sizeof(char *));
    assert(filesList2);
    *filesList2=NULL;

    //we read the header file to get the two important header lines and
    //(last line of header) the snips per file, the snip size in bits, the total snips,
    //and the minimum and maximum positions inside the file
    readHeaderFile(inputPathName,
            &headerLine1,
            &headerLine2,
            alignmentId,
            &snipsPerFile,
            &snipSize,
            &totalSnips,
            &posMin,
            &posMax);
    if(inPath2set == 1)
    {
        readHeaderFile(inputPath2Name,
                &headerLine1_2,
                &headerLine2_2,
                alignmentId2,
                &snipsPerFile2,
                &snipSize2,
                &totalSnips2,
                &posMin2,
                &posMax2);
        if(snipSize != snipSize2)
        {
            fprintf(stderr, "\n ERROR: The two inputs have different snip sizes.\n\t%d !=\
                    %d\n",snipSize,snipSize2);
            exit(1);
        }
    }
    else
    {
        readHeaderFile(inputPathName,
                &headerLine1_2,
                &headerLine2_2,
                alignmentId2,
                &snipsPerFile2,
                &snipSize2,
                &totalSnips2,
                &posMin2,
                &posMax2);
    }

    if(posWset1 == 1)
    {
        if(posWmin1 < posMin || posWmin1 > posMax)
            posWmin1 = posMin;

        if(posWmax1 < posMin || posWmax1 > posMax)
            posWmax1 = posMax;

        if(posWset2 == 1 && inPath2set == 1)
        {
            if(posWmin2 < posMin2 || posWmin2 > posMax2)
                posWmin2 = posMin2;

            if(posWmax2 < posMin2 || posWmax2 > posMax2)
                posWmax2 = posMax2;
        }
        else if(posWset2 == 1)
        {
            if(posWmin2 < posMin || posWmin2 > posMax)
                posWmin2 = posMin;

            if(posWmax2 < posMin || posWmax2 > posMax)
                posWmax2 = posMax;
        }
        else if(inPath2set == 1)
        {
            posWset2 = 1;
            if(posWmin1 < posMin2 || posWmin1 > posMax2)
                posWmin2 = posMin2;
            else
                posWmin2 = posWmin1;

            if(posWmax1 < posMin2 || posWmax1 > posMax2)
                posWmax2 = posMax2;
            else
                posWmax2 = posWmax1;
        }
        else
        {
            if(posWmin2 < posMin || posWmin2 > posMax)
                posWmin2 = posMin;

            if(posWmax2 < posMin || posWmax2 > posMax)
                posWmax2 = posMax;

        }
    }
    else
    {
        posWset1=1;
        posWmin1=posMin;
        posWmax1=posMax;
        posWmin2=posMin;
        posWmax2=posMax;
    }

    //little hack for correct output format in split-file inputLists.
    if(inPath2set == 1 && 
       posWset1 == 1 &&
       posWset2 == 1 && 
       posWmin1 == posWmin2 && 
       posWmax1 == posWmax2)
        posWset2=0;

    //this check is needed if the parser is given a large enough file size that the
    //total snips are less than the calculated snip size
    if(snipsPerFile < 0 || snipsPerFile > totalSnips)
    {
        snipsPerFile = totalSnips;
    }

    if(inPath2set == 1 && (snipsPerFile2 < 0 || snipsPerFile2 > totalSnips2))
    {
        snipsPerFile2 = totalSnips2;
    }

    //with the given window(s) or input list we find the files we need to
    //open to read the snips we need
    findFiles(inputPathName,
              alignmentId,
              snipsPerFile,
              totalSnips,
              posWmin1,
              posWmax1,
              &filesList,
              &filesListNum,
              mdf);

    if(inPath2set == 1)
    {
        findFiles(inputPath2Name,
                  alignmentId2,
                  snipsPerFile2,
                  totalSnips2,
                  posWmin2,
                  posWmax2,
                  &filesList2,
                  &filesListNum2,
                  mdf);
    }
    else
    {
        findFiles(inputPathName,
                  alignmentId2,
                  snipsPerFile2,
                  totalSnips2,
                  posWmin2,
                  posWmax2,
                  &filesList2,
                  &filesListNum2,
                  mdf);
    }
    
    // Count number of tasks, to use on task allocation to threads 
    // (except in competing mode)
    (*task_count)++;

#ifdef VERBOSE
    printf("Task[%d] created\n",(*task_count));
#endif

    enqueue_task(fpInRep,
                 inputPathName,
                 posWmin1,
                 posWmax1,
                 inputPath2Name,
                 posWmin2,
                 posWmax2,
                 outputFileName,
                 inPath2set,
                 posWset1,
                 posWset2,
                 filesList,
                 filesListNum,
                 filesList2,
                 filesListNum2,
                 headerLine1,
                 headerLine2,
                 headerLine1_2,
                 headerLine2_2,
                 snipSize,
                 snipSize2,
                 sampleList,
                 sampleList2);

    for(int i=0; i < filesListNum; i++)
            free(filesList[i]);
    free(filesList);
    for(int i=0; i < filesListNum2; i++)
            free(filesList2[i]);
    free(filesList2);
}
