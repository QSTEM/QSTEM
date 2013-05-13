/*************************************************************
 * file: readparams.c
 *
 * contains functions for reading parameters from a data file
 * These parameters are specified by a title string
 * The title string will be case sensitive
 *
 *************************************************************/
#ifndef READPARAMS_H
#define READPARAMS_H


/********************************************************
 * open the parameter file and return 1 for success,
 * 0 for failure
 * If a file is already open, close it first
 *******************************************************/
int parOpen(const char *fileName);

/********************************************************
 * Close the parameter file, if it is open
 *******************************************************/
void parClose();

/******************************************************
 * push the current FILE pointer on the stack, in order 
 * open a second (or up to 5th) parameter file
 */
void parFpPush();

/******************************************************
 * discard the current file pointer (close, if not zero)
 * and make the previous one on the stack active
 ****************************************************/
void parFpPull();

/*******************************************************
 * if for some reason a function needs to see the file 
 * pointer, it can obtain it by calling this function 
 ******************************************************/
FILE *getFp();

/***********************************************************
 * setComment(char newComment) will set a new character as 
 * the comment character
 **********************************************************/
void setComment(char newComment);

/************************************************************
 * function: int readparam(char *title, char *parString)
 * 
 * returns 1 for sucess, 0 for failure
 * fp: file pointer to (open) parameter file
 * title: string defining the parameter
 * parString: string containing the parameter
 *
 * The function will start at the current file pointer but 
 * start from the beginning if it could not find the parameter
 * It will therefore work faster if the parameters are called 
 * in order.
 ************************************************************/
int readparam(char *title, char *parString,int wrapFlag);


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
int readNextParam(char *title, char *parString);

/******************************************************************
 * This function reads the next line of the file and copies it into
 * buf (at most bufLen characters).
 * Content of buf not altered, if unsuccessful
 ******************************************************************/ 
int readNextLine(char *buf,int bufLen);


void resetParamFile();

char *strnext(char *str,char *delim);


#endif /* READPARAMS_H */
