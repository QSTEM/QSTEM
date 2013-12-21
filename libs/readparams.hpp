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
 *
 *************************************************************/
#ifndef READPARAMS_H
#define READPARAMS_H

#define	NCMAX	132	/* characters per line to read */

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
int readparam(FILE *fp, const char *title, char *parString,int wrapFlag);


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
int readNextParam(FILE *fp, char *title, char *parString);

/******************************************************************
 * This function reads the next line of the file and copies it into
 * buf (at most bufLen characters).
 * Content of buf not altered, if unsuccessful
 ******************************************************************/ 
int readNextLine(FILE *fp, char *buf,int bufLen);


void resetParamFile(FILE *fp);

char *strnext(char *str,char *delim);


#endif /* READPARAMS_H */
