
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

// template design inspired by:
// https://svn.ssec.wisc.edu/repos/lib_ifts/trunk/lib_ifts/src/fftw.hxx

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif





/*************  1D arrays ****************/

boost::shared_ptr<float_1D_type> float1D(int size, std::string message)
{
	float1D_type arr(boost::extents[size]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<float1D_type>(arr);
}

/*************   2D arrays  *************/

boost::shared_ptr<float2D_type> float2D(int nx, int ny, std::string message)
{
	float2D_type arr(boost::extents[nx][ny]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<float2D_type>(arr);
}

boost::shared_ptr<double2D_type> double2D(int nx, int ny, std::string message)
{
	double2D_type arr(boost::extents[nx][ny]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<double2D_type>(arr);
}

boost::shared_ptr<short2D_type> int2D(int nx, int ny, std::string message)
{
	short2D_type arr(boost::extents[nx][ny]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<short2D_type>(arr);
}

/*************    3D arrays *************/


// complex double matrix aligned for SIMD instructions
typedef boost::multi_array<std::complex<double>, 2, fftw_allocator< std::complex<double> > > fftw_cplx_doubleses;


/*---------------------------- float3D() -------------------------------*/
/*
	3D array allocator for type float_tt
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
float_tt ***float3D( int nx, int ny,int nz, const char *message )
{	
  float_tt ***m;
  int i,j;
  
  m = (float_tt***) fftw_malloc( nx * sizeof(float_tt**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float_tt**)fftw_malloc(ny*sizeof(float_tt*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float_tt*) fftw_malloc(nz*ny*nx* sizeof(float_tt) );
  if( m[0] == NULL ){
    printf("float2D cannot allocate arrays, size=%d: %s\n",
	   ny*nx, message );
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (float_tt*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (float_tt) = %d\n",message,(int)m);
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
  
  m = (float***) fftw_malloc( nx * sizeof(float**) ); 
  if( m == NULL ) {
    printf("float3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (float **)fftw_malloc(ny*sizeof(float*));
    if (m[i] == NULL) {
      printf("float3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (float*) fftw_malloc(nz*ny*nx* sizeof(float) );
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




/*---------------------------- complex2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
fftw_complex **complex2D( int nx, int ny, const char *message)
{
  fftw_complex **m;
  int i;
  
  m = (fftw_complex**) fftw_malloc( nx * sizeof(fftw_complex*) ); 
  if( m == NULL ) {
    printf("float2D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  
  m[0] = (fftw_complex*) fftw_malloc( ny *nx* sizeof(fftw_complex) );
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


/*---------------------------- complex2Df() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
fftwf_complex **complex2Df( int nx, int ny, const char *message)
{	fftwf_complex **m;
	int i;

	m = (fftwf_complex**) fftw_malloc( nx * sizeof(fftwf_complex*) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	m[0] = (fftwf_complex*) fftw_malloc( ny *nx* sizeof(fftwf_complex) );
	if( m[0] == NULL ){
	  printf("float2D cannot allocate arrays, size=%d: %s\n",
		 ny*nx, message );
	  exit(0);
	}
	for (i=1; i<nx; i++){
	  m[i] = (fftwf_complex*)(&m[0][i*ny]);
	}
#ifdef PRINT_MESSAGE
	printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif

	return m;

}  /* end complex2Df() */

/*---------------------------- complex3D() -------------------------------*/
/*
	3D array allocator for type fftw_complex
	make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

*/
fftw_complex ***complex3D( int nx, int ny,int nz, const char *message)
{	
  fftw_complex ***m;
  int i,j;
  
  m = (fftw_complex***)fftw_malloc( nx * sizeof(fftw_complex**) ); 
  if( m == NULL ) {
    printf("complex3D cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (fftw_complex**)fftw_malloc(ny*sizeof(fftw_complex*));
    if (m[i] == NULL) {
      printf("complex3D cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (fftw_complex*) fftw_malloc(nz*ny*nx* sizeof(fftw_complex) );
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
fftwf_complex ***complex3Df( int nx, int ny,int nz, const char *message)
{	
  fftwf_complex ***m;
  int i,j;
  
  m = (fftwf_complex***)fftwf_malloc( nx * sizeof(fftwf_complex**) ); 
  if( m == NULL ) {
    printf("complex3Df cannot allocate pointers, size=%d: %s\n",
	   nx, message );
    exit(0);
  }
  for (i=0;i<nx;i++) {
    m[i] = (fftwf_complex**)fftwf_malloc(ny*sizeof(fftwf_complex*));
    if (m[i] == NULL) {
      printf("complex3Df cannot allocate pointers (stage2), "
	     "size=%d: %s\n",ny, message );
      exit(0);
    }
  }
  
  m[0][0] = (fftwf_complex*) fftwf_malloc(nz*ny*nx* sizeof(fftwf_complex) );
  if( m[0][0] == NULL ){
    printf("complex3Df cannot allocate consecutive memory %d MB (for array %s)\n",
	    nz*ny*nx* sizeof(fftwf_complex)/(1024*1024),message);
    exit(0);
  }
  for (i=0; i<nx; i++) for (j=0;j<ny;j++) {
    m[i][j] = (fftwf_complex*)(&(m[0][0][nz*(i*ny+j)]));
  }
#ifdef PRINT_MESSAGE
  printf("allocated memory for %s (fftw_complex) = %d\n",message,(int)m);
#endif
  
  return m;
    
}  /* end complex3D() */



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


