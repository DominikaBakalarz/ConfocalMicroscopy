#include "lodepng.h"
#include <bits/stdc++.h>

#define MAXCELLS 1048576

using namespace std;
vector <unsigned char> image;
float cells[MAXCELLS*3];
int cStatus[MAXCELLS];
int cSize[MAXCELLS];
int cTypes[MAXCELLS];
int numberOfCells;
int numberOfTrials;
long imageLen;
double matchPercent;
double interPercent;
double cellSize;
unsigned width,height;
vector <int> pixelsInCell;

void decodeOneStep(const char* filename)
{
	unsigned error = lodepng::decode(image, width, height,filename);
	if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
        imageLen=(long)width*(long)height;
}

bool black(int r, int g, int b)
{
	if(r<=20 && g<=20 && b<=20)
		return true;
	else
		return false;

}

bool red(int r, int g, int b)
{
	if(r>g+5)
		return true;
	else
		return false;
}

bool green(int r, int g, int b)
{
	if (g-r>-20 && g>30)
		return true;
	else
		return false;
}
bool white(int r, int g, int b)
{
	if(r==255 && g==255 && b==255)
		return true;
	else
		return false;

}

long selectPixel(long b)
{
	int m=rand();
	return m%b;
}

int DrawPixel(int x, int y)
{
	int a=4*(1024*y+x);
	if (a>=0 && a<=4*1024*1024)
		pixelsInCell.push_back(a);

}

void drawTrialCell(int x0, int y0, int r)
{
	for(int i=-r; i< r; i++)
	{
		for(int j=-r; j< r; j++)
		{
			if( j*j + i*i< r*r)
			DrawPixel(x0-j, y0-i);
		}
	}
}

bool enoughGreen( vector <int> a)
{
	int i=0;
	for(int j=0; j<a.size(); j++)
	{
		if(green(image[a[j]], image[a[j]+1], image[a[j]+2]))
			i++;
	}
	if ((float)i/a.size()> matchPercent)
		return true;
	else
		return false;
}
bool enoughRed( vector <int> a)
{
	int i=0;
	for(int j=0; j<a.size(); j++)
	{
		if(red(image[a[j]], image[a[j]+1], image[a[j]+2]))
			i++;
	}
	if ((float)i/a.size()> matchPercent)
		return true;
	else
		return false;
}

void drawCell(vector <int> a)
{
	for (int i=0; i<a.size(); i++)
	{
		image[a[i]]=255;
		image[a[i]+1]=255;
		image[a[i]+2]=255;
		image[a[i]+3]=255;

	}
}
void drawDeathCell(vector <int> a)
{
	for (int i=0; i<a.size(); i++)
	{
		image[a[i]]=0;
		image[a[i]+1]=0;
		image[a[i]+2]=255;
		image[a[i]+3]=255;

	}
}

bool smallIntersection( vector <int> a)
{
	int i=0;
	for(int j=0; j<(int)a.size(); j++)
	{
		if(white(image[a[j]], image[a[j]+1], image[a[j]+2]))
			i++;
	}
	if ((float)i/a.size()<interPercent)
		return true;
	else
		return false;
}

void saveCell(int x,int y,int z,int status,int size) {
   cells[3*numberOfCells]=(float)x;
   cells[3*numberOfCells+1]=(float)y;
   cells[3*numberOfCells+2]=(float)z;
   cStatus[numberOfCells]=status;
   cSize[numberOfCells]=size;
   cTypes[numberOfCells]=1;
   numberOfCells++;
}

void swap_Nbyte(char *data, int n, int m)
{
  int i, j;
  char old_data[16];

  for(j = 0; j < n; j++) {
    memcpy(&old_data[0], &data[j * m], m);
    for(i = 0; i < m; i++)
      data[j * m + i] = old_data[m - i - 1];
  }
}


void writeVTK() {
  FILE *fhandle;
   char header[256];
  int i;

  fhandle=fopen("result.vtk","wb");

  sprintf(header,"# vtk DataFile Version 2.0\nStemscan output\nBINARY\nDATASET UNSTRUCTURED_GRID\n");
  fwrite(header,sizeof(char),strlen(header),fhandle);
  memset(header,0,256);
  sprintf(header,"\nPOINTS %d float\n",numberOfCells);
  fwrite(header,sizeof(char),strlen(header),fhandle);
  swap_Nbyte((char*)cells,numberOfCells*3,sizeof(float));
  fwrite(cells,sizeof(float),numberOfCells*3,fhandle);
  memset(header,0,256);
  sprintf(header,"\nCELL_TYPES %d\n",numberOfCells); 
  fwrite(header,sizeof(char),strlen(header),fhandle); 
  swap_Nbyte((char*)cTypes,numberOfCells,sizeof(int));
  fwrite(cTypes,sizeof(int),numberOfCells,fhandle);
  memset(header,0,256);
  sprintf(header,"\nPOINT_DATA %d",numberOfCells);
  fwrite(header,sizeof(char),strlen(header),fhandle);
  memset(header,0,256);
  sprintf(header,"\nSCALARS status int 1\nLOOKUP_TABLE default\n"); 
  fwrite(header,sizeof(char),strlen(header),fhandle);
  swap_Nbyte((char*)cStatus,numberOfCells,sizeof(int));
  fwrite(cStatus,sizeof(int),numberOfCells,fhandle);
  memset(header,0,256);
  sprintf(header,"\nSCALARS size int 1\nLOOKUP_TABLE default\n");
  fwrite(header,sizeof(char),strlen(header),fhandle);
  swap_Nbyte((char*)cSize,numberOfCells,sizeof(int));
  fwrite(cSize,sizeof(int),numberOfCells,fhandle);

  fclose(fhandle);
}

int main(int argc,char **argv)
{
    vector <int> vtk;

    numberOfCells=0;

    if(argc!=5) {
      printf("Usage: ./stemscan matchPercent interPercent cellSize numberOfTrials\n");
      exit(1);
    }

    matchPercent=(double)atof(argv[1]);
    interPercent=(double)atof(argv[2]);
    cellSize=(double)atof(argv[3]);
    numberOfTrials=atoi(argv[4]);

    for (int i=0; i<1024*1024;i++)
        vtk.push_back(0);

    decodeOneStep("03_rgb_z19.png");

    for(int i=0;i<numberOfTrials;i++) {
      long pixNum;
      int x,y,z;
      z=0;
      pixNum=selectPixel(imageLen);
      x=pixNum%width;
      y=pixNum/width;
      pixNum*=4;
      if(green(image[pixNum],image[pixNum+1],image[pixNum+2])) {
        drawTrialCell(x,y,cellSize);
        if(smallIntersection(pixelsInCell) && enoughGreen(pixelsInCell)) {
          drawCell(pixelsInCell);
          saveCell(x,y,z,1,cellSize);
        }  
        pixelsInCell.clear();
      }
      if(red(image[pixNum],image[pixNum+1],image[pixNum+2])) {
        drawTrialCell(x,y,cellSize);
        if(smallIntersection(pixelsInCell) && enoughRed(pixelsInCell)) {
          drawDeathCell(pixelsInCell); 
          saveCell(x,y,z,0,cellSize);
        }
        pixelsInCell.clear();
      }
    }     

  writeVTK();

  unsigned error = lodepng::encode("result.png", image, 1024, 1024);

  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

}




