#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <wchar.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif

float *rotate(int Xrot_in, int Yrot_in, int Zrot_in, int NumOfCoords, float *Coords)
{
int i;
float cx,cy,cz,sx,sy,sz;
static float *tempCoords;

tempCoords=malloc(sizeof(float[3])*NumOfCoords);

cx=cos(((float)Xrot_in)/180*M_PI);cy=cos(((float)Yrot_in)/180*M_PI);cz=cos(((float)Zrot_in)/180*M_PI);
sx=sin(((float)Xrot_in)/180*M_PI);sy=sin(((float)Yrot_in)/180*M_PI);sz=sin(((float)Zrot_in)/180*M_PI);

for(i=0;i<NumOfCoords;i++){
	/*
	*(tempCoords+3*i)	=(cy*cz)*(*(Coords+3*i))		+(cy*sz)*(*(Coords+3*i+1))		-(sy)*(*(Coords+3*i+2));
	*(tempCoords+3*i+1)	=(cz*sy*sx-cx*sz)*(*(Coords+3*i)) 	+(cx*cz+sx*sy*sz)*(*(Coords+3*i+1)) 	+(sx*cy)*(*(Coords+3*i+2));
	*(tempCoords+3*i+2)	=(sx*sz+cx*sy*cz)*(*(Coords+3*i)) 	+(cx*sy*sz-sx*cz)*(*(Coords+3*i+1))	+(cx*cy)*(*(Coords+3*i+2));
	*/
	*(tempCoords+3*i)	=(cy*cz)*(*(Coords+3*i))	+(cz*sy*sx-cx*sz)*(*(Coords+3*i+1))	+(sx*sz+cx*sy*cz)*(*(Coords+3*i+2));
	*(tempCoords+3*i+1)	=(cy*sz)*(*(Coords+3*i))	+(cx*cz+sx*sy*sz)*(*(Coords+3*i+1))	+(cx*sy*sz-sx*cz)*(*(Coords+3*i+2));
	*(tempCoords+3*i+2)	=-(sy)*(*(Coords+3*i))		+(sx*cy)*(*(Coords+3*i+1))		+(cx*cy)*(*(Coords+3*i+2));

	}

return tempCoords;
}

float *bubsort_Coord( float *Coord, int NumOfCoords)
{
static float extremum[6]={1,2,3,4,5,6};
int i;
/*Initialize*/
extremum[0]=*(Coord);
extremum[2]=*(Coord+1);
extremum[4]=*(Coord+2);
extremum[1]=*(Coord);
extremum[3]=*(Coord+1);
extremum[5]=*(Coord+2);

/*Biggests*/
for(i=0;i<NumOfCoords;i++){
	extremum[0]=*(Coord+3*i)>extremum[0]?*(Coord+3*i):extremum[0];
	extremum[2]=*(Coord+3*i+1)>extremum[2]?*(Coord+3*i+1):extremum[2];
	extremum[4]=*(Coord+3*i+2)>extremum[4]?*(Coord+3*i+2):extremum[4];

/*Smallests*/
	extremum[1]=*(Coord+3*i)<extremum[1]?*(Coord+3*i):extremum[1];
	extremum[3]=*(Coord+3*i+1)<extremum[3]?*(Coord+3*i+1):extremum[3];
	extremum[5]=*(Coord+3*i+2)<extremum[5]?*(Coord+3*i+2):extremum[5];
	}

return extremum;
}

int main(int argc,char *argv[])
{
int NumOfAtoms;
int i,j,dummy=0;

float dimensions[3];
float *oBox;
float *rBox;
float oBox_ar[][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
float *oCoord,*rCoord;
float *extremum_acceptor;
float extremum_box[6]={0,0,0,0,0,0};
float extremum_main[6]={0,0,0,0,0,0};
float Spacing;
int Xrot,Yrot,Zrot;
char nullchar[255];

FILE *in_stream;

/* File stream opening*/
if ((in_stream=fopen(argv[1],"r"))==NULL) {
	perror("Error when opening file");
	exit(EXIT_FAILURE);
	}

/* Coordicates and Rotation Parameters Load in */
fscanf(in_stream,"Number of particles = %d\n",&NumOfAtoms);

oCoord=malloc(sizeof(float[3])*NumOfAtoms);
rCoord=malloc(sizeof(float[3])*NumOfAtoms);
extremum_acceptor=malloc(sizeof(float)*6);
rBox=malloc(sizeof(float[3])*8);
oBox=oBox_ar;

	for(j=0;j<2;j++) fgets(nullchar,255,in_stream);
fscanf(in_stream,"H0(1,1) = %f A\n",&dimensions[0]);
	for(j=0;j<3;j++) fgets(nullchar,255,in_stream);
fscanf(in_stream,"H0(2,2) = %f A\n",&dimensions[1]);
	for(j=0;j<3;j++) fgets(nullchar,255,in_stream);
fscanf(in_stream,"H0(3,3) = %f A\n",&dimensions[2]);
	for(j=0;j<4;j++) fgets(nullchar,255,in_stream);

	dummy=0;
while(1){
        if(fscanf(in_stream,"%f %f %f \n",(oCoord+3*dummy),(oCoord+3*dummy+1),(oCoord+3*dummy+2))==EOF)
		break;
        else{
/*		printf("%f %f %f\n",*(oCoord+3*dummy),*(oCoord+3*dummy),*(oCoord+3*dummy+2));*/
		dummy++;
        }
	}


Xrot=atoi(argv[2]);
Yrot=atoi(argv[3]);
Zrot=atoi(argv[4]);
Spacing=atof(argv[5]);

/* Preload Box and Atom coordicates with dimensions */
for(i=0;i<8;i++){
	*(oBox+3*i)=*(oBox+3*i)*dimensions[0];
	*(oBox+3*i+1)=*(oBox+3*i+1)*dimensions[1];
	*(oBox+3*i+2)=*(oBox+3*i+2)*dimensions[2];
	}

for(i=0;i<NumOfAtoms;i++){
	*(oCoord+3*i)=*(oCoord+3*i)*dimensions[0];
	*(oCoord+3*i+1)=*(oCoord+3*i+1)*dimensions[1];
	*(oCoord+3*i+2)=*(oCoord+3*i+2)*dimensions[2];
	}

/* Bubble Sorting for the extremums */

rBox=rotate(Xrot,Yrot,Zrot,(int)8,oBox);
extremum_acceptor=bubsort_Coord(rBox,(int)8);
for(i=0;i<6;i++) extremum_box[i]=*(extremum_acceptor+i);

rCoord=rotate(Xrot,Yrot,Zrot,NumOfAtoms,oCoord);
extremum_acceptor=bubsort_Coord(rCoord,NumOfAtoms);
for(i=0;i<6;i++) extremum_main[i]=*(extremum_acceptor+i);

/* Testing Scripts */
/*
printf("--- %s ---\n",argv[1]);
printf("Num of Atoms = %d\n",NumOfAtoms);
printf("Dimensions X-> %f; Y-> %f;  Z-> %f\n",dimensions[0],dimensions[1],dimensions[2]);
printf("\n");
printf("Rotation Parameters X-> %d; Y-> %d; Z-> %d\n",Xrot,Yrot,Zrot);
printf("\n");
for(i=0;i<8;i++)
        printf("Orig Box Num %d = X-> %f; Y-> %f;  Z-> %f\n",i,*(oBox+3*i),*(oBox+3*i+1),*(oBox+3*i+2));
for(i=0;i<8;i++)
	printf("Rota Box Num %d = X-> %f; Y-> %f;  Z-> %f\n",i,*(rBox+3*i),*(rBox+3*i+1),*(rBox+3*i+2));
printf("\n");
for(i=0;i<NumOfAtoms;i++)
	printf("Orig atom Num %d = X-> %f; Y-> %f;  Z-> %f \n",i,*(oCoord+3*i),*(oCoord+3*i+1),*(oCoord+3*i+2));
for(i=0;i<NumOfAtoms;i++)
        printf("Rota atom Num %d = X-> %f; Y-> %f;  Z-> %f \n",i,*(rCoord+3*i),*(rCoord+3*i+1),*(rCoord+3*i+2));
printf("Extremums Main %f %f %f %f %f %f\n",extremum_main[0],extremum_main[1],extremum_main[2],extremum_main[3],extremum_main[4],extremum_main[5]);
printf("Extremums Box %f %f %f %f %f %f\n",extremum_box[0],extremum_box[1],extremum_box[2],extremum_box[3],extremum_box[4],extremum_box[5]);
*/

/* Final outputs of extremem values in X and Y direction */
printf("Xmax-> %f Xmin-> %f\n",(extremum_main[0]-extremum_box[1]),(extremum_main[1]-extremum_box[1]));
printf("Ymax-> %f Ymin-> %f\n",(extremum_main[2]-extremum_box[3]),(extremum_main[3]-extremum_box[3]));
printf("Spaced Xmax-> %f Xmin-> %f\n",(extremum_main[0]-extremum_box[1]+Spacing),(extremum_main[1]-extremum_box[1]-Spacing));
printf("Spaced Ymax-> %f Ymin-> %f\n",(extremum_main[2]-extremum_box[3]+Spacing),(extremum_main[3]-extremum_box[3]-Spacing));
printf("%f %f %f %f\n",(extremum_main[0]-extremum_box[1]+Spacing),(extremum_main[1]-extremum_box[1]-Spacing),(extremum_main[2]-extremum_box[3]+Spacing),(extremum_main[3]-extremum_box[3]-Spacing));

/* Ending Part */
free(oCoord);
free(rCoord);
free(rBox);
fclose(in_stream);
return 0;
}
