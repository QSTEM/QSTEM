
#include <stdlib.h>
#include <stdio.h>
#include "fftw3.h"
#include "memory_fftw3.h"
 
#ifndef WIN32
#include <stdint.h>
#endif
/*
#define PRINT_MESSAGE
*/

/*---------------------------- any2D() -------------------------------*/
/*
	2D array allocator for any type of size 'size'
	make space for m[0...(nx-1)][0..(ny-1)]

*/
void **any2D( int nx, int ny,int size, const char *message )
{	void **m;
	int i;

	m = (void **)fftw_malloc( nx * sizeof(void *)); 
	if( m == NULL ) {
		printf("any2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] =  fftw_malloc( ny *nx* size );
	if( m[0] == NULL ){
	  printf("any2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (void *)((intptr_t)(m[0])+(i*ny*size));
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif

	return m;

}  /* end any2D() */


/*---------------------------- any3D() -------------------------------*/
/*
	3D array allocator for any type 
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
void ***any3D( int nx, int ny,int nz,int size, const char *message )
{	
  void ***m;
  int i,j;
  
  m = (void***) fftw_malloc( nx * sizeof(void **)); 
  if( m == NULL ) {
    printf("any3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (void **)fftw_malloc(ny*sizeof(void *));
    if (m[i] == NULL) {
      printf("any3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = fftw_malloc(nz*ny*nx* size );
  if( m[0][0] == NULL ){
	printf("any3Df cannot allocate consecutive memory %d MB (for array %s)\n",
	    nz*ny*nx*size/(1024*1024),message);
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] =(void *)((intptr_t)(m[0][0])+size*nz*(i*ny+j));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end any3D() */


