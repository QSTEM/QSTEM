/*************************************************************
 * file: readparams.c
 *
 * contains functions for reading parameters from a data file
 * These parameters are specified by a title string
 * The title string will be case sensitive
 * Several data files can be kept open (up to STACK_SIZE)
 * parameter files can be pushed on/pulled off the stack
 *
 *************************************************************/

#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include "readparams.h"

#define COMMENT '%'
#define PAR_BUF_LEN 1024
#define STACK_SIZE 5

FILE *fp=NULL;
FILE **fpStack = NULL;
char parBuf[PAR_BUF_LEN];
char commentChar=COMMENT;

/********************************************************
 * open the parameter file and return 1 for success,
 * 0 for failure
 * If a file is already open, close it first
 *******************************************************/
int parOpen(const char *fileName) {
  int i;

  if (fpStack == NULL) {
    fpStack = (FILE **)malloc(STACK_SIZE*sizeof(FILE *));
    for (i=0;i<STACK_SIZE;i++)
      fpStack[i] = NULL;
    fpStack[0] = fp;
  }

  if (fp!=NULL)
    fclose(fp);

  fp=fopen(fileName,"r");
  return (fp != NULL);
}

/******************************************************
 * push the current FILE pointer on the stack, in order 
 * open a second (or up to 5th) parameter file
 */
void parFpPush() {
  int i;

  if (fpStack == NULL) {
    fpStack = (FILE **)malloc(STACK_SIZE*sizeof(FILE *));
    for (i=0;i<STACK_SIZE;i++)
      fpStack[i] = NULL;
  }
  else {
    for (i=1;i<STACK_SIZE;i++)
      fpStack[i] = fpStack[i-1];
    fpStack[0] = fp; 
  }
  fp = NULL;
}

/******************************************************
 * discard the current file pointer (close, if not zero)
 * and make the previous one on the stack active
 ****************************************************/
void parFpPull() {
  int i;

  parClose();
  fp = fpStack[0];
  for (i=1;i<STACK_SIZE;i++)
    fpStack[i-1] = fpStack[i];
  fpStack[STACK_SIZE-1] = NULL;
}


/********************************************************
 * Close the parameter file, if it is open
 *******************************************************/
void parClose() {
  if (fp!=NULL)
    fclose(fp);
  fp = NULL;
}

/*******************************************************
 * if for some reason a function needs to see the file 
 * pointer, it can obtain it by calling this function 
 ******************************************************/
FILE *getFp() {
  return fp;
}


/***********************************************************
 * setComment(char newComment) will set a new character as 
 * the comment character
 **********************************************************/
void setComment(char newComment) {
  commentChar = newComment;
}

/************************************************************
 * function: int readparam(char *title, char *parString)
 * 
 * returns 1 for sucess, 0 for failure
 * fp: file pointer to (open) parameter file
 * title: string defining the parameter
 * parString: string containing the parameter
 * wrapFlag: decides whether we can wrap around the end of the file
 *
 * The function will start at the current file pointer but 
 * start from the beginning if it could not find the parameter
 * It will therefore work faster if the parameters are called 
 * in order.
 ************************************************************/
int readparam(const char *title, char *parString, int wrapFlag) {
  char *result;
  char *str=NULL,*comment;
  // int iresult;

  if (fp == NULL)
    return 0;

  do {
    result = fgets(parBuf,(int)PAR_BUF_LEN,fp);
    /* cut off at the comment */
    comment = strchr(parBuf,COMMENT);
    if (comment!=NULL)
      *comment = '\0';
    if (result !=NULL) 
      str = strstr(parBuf,title);
  } while((result != NULL) && (str == NULL));

  if ((result == NULL) && (wrapFlag)) {
    /* reset the file position pointer, if needed */
    fseek(fp, 0L, SEEK_SET);
    do {
      result = fgets(parBuf,PAR_BUF_LEN,fp);
      comment = strchr(parBuf,COMMENT);
      if (comment!=NULL)
	*comment = '\0';
      if (result !=NULL) {
	str = strstr(parBuf,title);
      }
    } while((result != NULL) && (str == NULL));
  }
  if (result == NULL) {
    // printf("Could not find parameter %s\n",title);
    return 0;
  }

  strcpy(parString,str+strlen(title));
  
  return 1;
} 

/*********************************************************
 * reset the parameter file pointer to zero, 
 * in order to look for the first occurence of something.
 *********************************************************/
void resetParamFile() {
  if (fp !=NULL)
    fseek(fp, 0L, SEEK_SET);
}

/* This function returns a pointer to the next word in the string str
 * It does not alter str.  Gaps between words are defined by any of
 * the characters in delim. (e.g. delim=" \t")
 * returns NULL, if end of string has been found
 */
char *strnext(char *str,char *delim) {
  int found=0;
  char *str2;


  for(str2 = str;(*str2 != '\0');str2++) {
    if (strchr(delim,*str2)) found=1;
    if ((found) && (strchr(delim,*str2) == NULL))
      break;
  }
  if (*str2 == '\0')
    return NULL;
  if (*str2 == '\n')
    return NULL;
  return str2;
}


/************************************************************
 * function: int readparam(char *title, char *parString)
 * 
 * returns 1 for sucess, 0 for failure
 * title: string defining the parameter
 * parString: string containing the parameter
 *
 * The function will start at the current file pointer and
 * write the title and value (parString) of the next 
 * (legal, i.e. at least 2 words, at least one after the colon)
 * parameter in the corresponding strings.
 ************************************************************/
int readNextParam(char *title, char *parString) {
  char *result;
  char *str,*comment;
  // int iresult;

  if (fp == NULL)
    return 0;

  do {
    result = fgets(parBuf,(int)PAR_BUF_LEN,fp);
    if (result == NULL)
      return 0;

    /* cut off at the comment */
    comment = strchr(parBuf,COMMENT);
    if (comment!=NULL)
      *comment = '\0';
    /* find the colon */
    str = strchr(parBuf,':');
    if (str != NULL) {
      str = strnext(str," \t");
      // printf("%s\n",str);
    }
  } while(str == NULL);
  /* Now we found a legal parameter in $(parBuf): str
   */
  strcpy(parString,str);
  strncpy(title,parBuf,str-parBuf);
  title[str-parBuf] = '\0';
  return 1;  /* success */
} 



/******************************************************************
 * This function reads the next line of the file and copies it into
 * buf (at most bufLen characters).
 * Content of buf not altered, if unsuccessful
 ******************************************************************/ 
int readNextLine(char *buf,int bufLen) {
  char *result;
 
  if (fp == NULL)
    return 0;

  result = fgets(parBuf,(int)PAR_BUF_LEN,fp);
  if (result == NULL)
      return 0;
  result = strchr(parBuf,'%');
  if (result != NULL)
    *result = '\0';
  strncpy(buf,parBuf,bufLen);

  return 1;  /* success */
} 




