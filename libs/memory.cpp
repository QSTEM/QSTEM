
#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "fftw3.h"
 
/*
#define PRINT_MESSAGE
*/
typedef float float_t;
// typedef float[2] fftw_complex;

/*---------------------------- char1D() -------------------------------*/
/*
	1D array allocator for type char
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
char* char1D( int n, const char *message )
{
	char *m;
	
	m = (char*) malloc( n * sizeof( char ) );
	if( m == NULL ) {
		printf("char1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return( m );

}  /* end char1D() */

/*---------------------------- int1D() -------------------------------*/
/*
	1D array allocator for type int
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
int* int1D( int n, const char *message )
{
	int *m;
	
	m = (int*) malloc( n * sizeof( int ) );
	if( m == NULL ) {
		printf("int1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif
	return( m );

}  /* end int1D() */

/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
float_t *float1D( int n, const char *message )
{
	float_t *m;
	
	m = (float_t*) malloc( n * sizeof( float_t) );
	if( m == NULL ) {
		printf("float1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return( m );

}  /* end float1D() */

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
double* double1D( int n, const char *message )
{
	double *m;
	
	m = (double*) malloc( n * sizeof( double ) );
	if( m == NULL ) {
		printf("double1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return( m );

} /* end double1D() */

/*---------------------------- char2D() -------------------------------*/
/*
	2D array allocator for type char
	make space for m[0...(nx-1)][0..(ny-1)]
	
*/
char **char2D( int nx, int ny, const char *message )
{	char **m;
	int i;

	m = (char**) malloc( nx * sizeof( char* ) ); 
	if( m == NULL ) {
		printf("char2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (char*) malloc( ny * sizeof( char ) );
		if( m[i] == NULL ){
			printf("char2D cannot allocate arrays, size=%d: %s\n",
			       ny, message );
			exit(0);
		}
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;
}  /* end char2D() */

/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
short **short2D( int nx, int ny, const char *message )
{	short **m;
	int i;

	m = (short**) malloc( nx * sizeof( short* ) ); 
	if( m == NULL ) {
		printf("short2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}
	m[0] = (short *) malloc( ny *nx* sizeof(short) );
	if( m[0] == NULL ){
	  printf("long2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;
}  /* end short2d() */

/*---------------------------- int2D() -------------------------------*/
/*
	2D array allocator for type int
	make space for m[0...(nx-1)][0..(ny-1)]

*/
int **int2D( int nx, int ny, const char *message )
{	int **m;
	int i;

	m = (int**) malloc( nx * sizeof(int* ) ); 
	if( m == NULL ) {
		printf("int2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (int *) malloc( ny *nx* sizeof(int) );
	if( m[0] == NULL ){
	  printf("int2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (int *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (int) = %d\n",message,(int)m);
#endif

	return m;

}  /* end int2D() */

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message )
{	long **m;
	int i;

	m = (long**) malloc( nx * sizeof( long* ) ); 
	if( m == NULL ) {
		printf("long2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (long *) malloc( ny *nx* sizeof(long) );
	if( m[0] == NULL ){
	  printf("long2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;

}  /* end long2d() */

/*---------------------------- float32_2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/

float **float32_2D( int nx, int ny, const char *message )
{	float **m;
	int i;

	m = (float**) malloc( nx * sizeof( float* ) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (float *) malloc( ny *nx* sizeof( float ) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (float *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (float_t) = %d\n",message,(int)m);
#endif

	return m;

}  /* end float2D() */


/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_t **float2D( int nx, int ny, const char *message )
{	float_t **m;
	int i;

	m = (float_t**) malloc( nx * sizeof( float_t* ) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (float_t *) malloc( ny *nx* sizeof( float_t ) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (float_t *)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (float_t) = %d\n",message,(int)m);
#endif

	return m;

}  /* end float2D() */
/*---------------------------- float3D() -------------------------------*/
/*
	3D array allocator for type float_t
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
float_t ***float3D( int nx, int ny,int nz, const char *message )
{	
  float_t ***m;
  int i,j;
  
  m = (float_t***) malloc( nx * sizeof(float_t**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float_t**)malloc(ny*sizeof(float_t*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float_t*) malloc(nz*ny*nx* sizeof(float_t) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (float_t*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (float_t) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end float3D() */


/*---------------------------- float32_3D() -------------------------------*/
/*
	3D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
float ***float32_3D( int nx, int ny,int nz, const char *message )
{	
  float ***m;
  int i,j;
  
  m = (float***) malloc( nx * sizeof(float**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float **)malloc(ny*sizeof(float*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float*) malloc(nz*ny*nx* sizeof(float) );
  if( m[0] == NULL ){
    printf("float32_3D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (float*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (float32) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end float32_3D() */




/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type doubel
	make space for m[0...(nx-1)][0..(ny-1)]

*/
double **double2D( int nx, int ny, const char *message )
{	double **m;
	int i;

	m = (double**) malloc( nx * sizeof(double* ) ); 
	if( m == NULL ) {
		printf("double2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (double *) malloc( ny *nx* sizeof(double) );
	if( m[0] == NULL ){
	  printf("double2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = &(m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s\n",message);
#endif

	return m;

}  /* end double2D() */

/*---------------------------- complex2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
fftw_complex **complex2D( int nx, int ny, const char *message )
{	fftw_complex **m;
	int i;

	m = (fftw_complex**) malloc( nx * sizeof(fftw_complex*) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (fftw_complex*) malloc( ny *nx* sizeof(fftw_complex) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (fftw_complex*)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif

	return m;

}  /* end complex2D() */

/*---------------------------- complex3D() -------------------------------*/
/*
	3D array allocator for type fftw_complex
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
fftw_complex ***complex3D( int nx, int ny,int nz, const char *message )
{	
  fftw_complex ***m;
  int i,j;
  
  m = (fftw_complex***) malloc( nx * sizeof(fftw_complex**) ); 
  if( m == NULL ) {
    printf("complex3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (fftw_complex**)malloc(ny*sizeof(fftw_complex*));
    if (m[i] == NULL) {
      printf("complex3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (fftw_complex*) malloc(nz*ny*nx* sizeof(fftw_complex) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (fftw_complex*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end complex3D() */

/*---------------------------- complex3Df() -------------------------------*/
/*
	3D array allocator for type fftwf_complex
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
fftw_complex ***complex3Df( int nx, int ny,int nz, const char *message)
{	
  return complex3D(nx,ny,nz,message );
}  /* end complex3Df() */


/*---------------------------- any2D() -------------------------------*/
/*
	2D array allocator for any type of size 'size'
	make space for m[0...(nx-1)][0..(ny-1)]

*/
void **any2D( int nx, int ny,int size, const char *message )
{	void **m;
	int i;

	m =  (void **)malloc( nx * sizeof(void *)); 
	if( m == NULL ) {
		printf("any2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] =  malloc( ny *nx* size );
	if( m[0] == NULL ){
	  printf("any2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (void *)((int)(m[0])+(i*ny*size));
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
  
  m = (void***)malloc( nx * sizeof(void **)); 
  if( m == NULL ) {
    printf("any3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (void **)malloc(ny*sizeof(void **));
    if (m[i] == NULL) {
      printf("any3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = malloc(nz*ny*nx* size );
  if( m[0] == NULL ){
    printf("any3D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] =(void *)((int)(m[0][0])+size*nz*(i*ny+j));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end any3D() */


