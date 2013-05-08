
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

boost::shared_ptr<float1D_type> float1D(int size, std::string message)
{
	float1D_type arr(boost::extents[size]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<float1D_type>(arr);
}

boost::shared_ptr<double1D_type> double1D(int size, std::string message)
{
	double1D_type arr(boost::extents[size]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<double1D_type>(arr);
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
	int2D_type arr(boost::extents[nx][ny]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<int2D_type>(arr);
}

boost::shared_ptr<complex2D_type> complex2D(int nx, int ny, std::string message)
{
	complex2D_type arr(boost::extents[nx][ny]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<complex2D_type>(arr);
}


/*************    3D arrays *************/

boost::shared_ptr<float3D_type> float3D(int nx, int ny, int nz, std::string message)
{
	float3D_type arr(boost::extents[nx][ny][nz]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<float3D_type>(arr);
}

boost::shared_ptr<complex3D_type> complex3D(int nx, int ny,int nz, const char *message)
{
	complex3D_type arr(boost::extents[nx][ny][nz]);
	printf("allocated memory for %s\n",message.c_str());
	#ifdef PRINT_MESSAGE
		printf("allocated memory for %s\n",message.c_str());
	#endif
	return boost::shared_ptr<complex3D_type>(arr);
}



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


