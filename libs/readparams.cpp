/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#include "readparams.hpp"

#define COMMENT '%'
#define PAR_BUF_LEN 1024
//#define STACK_SIZE 5

//FILE *fp=NULL;
//FILE **fpStack = NULL;
static char parBuf[PAR_BUF_LEN];
static char commentChar=COMMENT;

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
int readparam(FILE *fp, const char *title, char *parString, int wrapFlag) {
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
void resetParamFile(FILE *fp) {
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
int readNextParam(FILE *fp, char *title, char *parString) {
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
int readNextLine(FILE *fp, char *buf,int bufLen) {
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

