/******************************************************************************
 * EMTV3DMain.cpp: EMTV for 3D image reconstruction
 *   
 ******************************************************************************
 * Copyright (C) 2010~2011 EMTV 3D Reconstruction project
 * Authors: "Jianwen Chen" <jwchen@ee.ucla.edu>
 *          "Ming Yan" <basca.yan@gmail.com>
 *          "Yi Zou" <zouyi@cs.ucla.edu>
 *
 * Version : 
 ******************************************************************************/


#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>


#ifdef WIN32
#include "extras/getopt.h"
#else
#include <sys/time.h>
#include <getopt.h>
#endif

#ifdef _MSC_VER
  typedef __int64 int64_t;
#else
  #include <stdint.h>		// Use the C99 official header
#endif

#define SINGLEPOINT 
#ifdef SINGLEPOINT
#define DataType float 
#else
#define DataType double 
#endif

char inputfilename0[100];
char inputfilename1[100];
int datasize;

void readFileToVector(char* in_file, DataType *data, int size)
{
   FILE* inputfile;	
   int i;

   inputfile = fopen(in_file, "r"); 
   if( inputfile == NULL ) 
   {
     printf("Cannot open the input file!\n");
     exit(1);
   }

   for (i = 0; i < size; i++)
         {
#ifdef SINGLEPOINT
           fscanf(inputfile,"%f",&data[i]);
#else           
           fscanf(inputfile,"%lf",&data[i]);
#endif           
         }
   fclose(inputfile);
}

#define ABS(i)  (i>0)?i:(-i)
void CompareVector(DataType *data0, DataType *data1, int size)
{
    int i;
    DataType diff=0.0;
    DataType diffall=0.0;

    for(i=0;i<size;i++)
    {
        diff = data0[i]-data1[i]; 

        if(diff > 0.1) 
#ifdef SINGLEPOINT
          printf("%d, data0[%d] = %f data1[%d]= %f \n", i, i, data0[i], i, data1[i]);
#else           
          printf("%d, data0[%d] = %lf data1[%d]= %lf \n", i, i, data0[i], i, data1[i]);
#endif          

        diffall += ABS(diff);
    }

#ifdef SINGLEPOINT
    printf("diffall = %f\n", diffall);
#else           
    printf("diffall = %lf\n", diffall);
#endif          
}

int help()
{
  printf("./filecompare -i inputfilename0 -e inputfilename1 -s size\n");
  return 0;
}


int configure(int argc, char*argv[])
{
  char option_def[] = "?vh:i:e:s:";
  int ch;

  for(;;){
    ch = getopt(argc, argv, option_def);
    if (ch == -1)
        break;

    switch (ch)
    {
      case '?':
      case 'v':
      case 'h':
          help();
          return -1;
      case 'i':
          strcpy(inputfilename0, optarg); 
          printf("input filename0 : %s \n", inputfilename0);  
          break;
      case 'e':
          strcpy(inputfilename1, optarg); 
          printf("input filename1 : %s \n", inputfilename1);  
          break;
      case 's':
          datasize = atoi(optarg); 
          printf("vector size:  %d\n",datasize);  
          break;
      default:
          printf("Unknow option detection\n");  
          return 0; 
    }
  }

  return 0;
}

int main(int argc, char* argv[])
{
        DataType *data0;
        DataType *data1;
	if(configure(argc, argv)==-1)
            return -1;

        if(argc<7)
        {
            help();
            return -1;
        }

#ifdef SINGLEPOINT
           printf("Using single point data-type\n");
#else             
           printf("Using double point data-type\n");
#endif
    
        //memory allocation
        data0 = (DataType*)malloc(datasize * sizeof(DataType));
        data1 = (DataType*)malloc(datasize * sizeof(DataType));

        //read input data 
	readFileToVector(inputfilename0, data0, datasize);
	readFileToVector(inputfilename1, data1, datasize);

        //compare the data 
        CompareVector(data0, data1, datasize);

        //free the memory 
        free(data0);
        free(data1);
 
        return 0;
}



