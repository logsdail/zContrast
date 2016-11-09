#include "CContrastImage.h"

using namespace std;

CContrastImage::CContrastImage(string sz, vector<string> *o)
// Constructor
// Input: string(sz) - File to be read containing intensity details
{
	output_content = o;
	readImage(sz);
}

CContrastImage::CContrastImage(int nRows, int nColumns, float gridSize, vector<string> *o)
// Constructor
// Input: int(nRows) - Y Total
// 	  int(ncolumns) - X Total
//	  float(gridSize) - Point spacing
{
	output_content = o;
	
	m_iRows = nRows;
	m_iColumns = nColumns;
	m_fGridSize = gridSize;
	
	int total = nRows*nColumns;	
	m_dArray.resize(total);
}

void CContrastImage::fourier(float nScale)
// Perform fourier transform on data
// This is lifted from an example file supplied with fftw3
// Though cout's have been removed to make viewing easier
// Data is centralised afterwards
// Inputs: float(nScale) - This is our scaling factor for the images,
//			   Otherwise the spread is too great.
{
	int i = 0;
	int j = 0;
  	double *in;
  	int nx = getWidth();
  	int ny = getHeight();
	
 	fftw_complex *out;
  	fftw_plan plan_forward;
	
	// INPUT ARRAY
  	in = (double*) malloc ( sizeof ( double ) * nx * ny );
	
  	for ( i = 0; i < nx; i++ )
		for ( j = 0; j < ny; j++ )
		{
			in[i*nx+j] = getIntensity(j,i);
			// Addition to clear array ready for out put Fourier
			setIntensity(j,i,0);
		}
	
	// OUTPUT ARRAY
  	int nyh = ( ny / 2 ) + 1;	//half+1 for Fourier transform
 	int ny_half = ( ny / 2 );	//true halves
  	int nx_half = ( nx / 2 );	//true halves
	
	// PERFORM FFT
  	out = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * nx * nyh );
  	plan_forward = fftw_plan_dft_r2c_2d ( nx, ny, in, out, FFTW_ESTIMATE );
  	fftw_execute ( plan_forward );
	
	// TRANSFORM COMPLETE  //
	// GET MAX FOR SCALING //
	
	double MAX = 0;
	
  	for ( i = 0; i < nx; i++ )
		for ( j = 0; j < ny_half; j++ )
		{
			if (out[i*nyh+j][0] > 0 && out[i*nyh+j][0] > MAX)
			{
				MAX = out[i*nyh+j][0];
			}
			else if (0-out[i*nyh+j][0] > MAX)
			{
				MAX = 0-out[i*nyh+j][0];
			}
		}
	
	// IMAGED CENTRALISING
	// SCALE is 3000 / log (1+INPUT)
	// Puts everything in close proximity for imaging
	double SCALE_VALUE = 3000; // COULD SOFTCODE THIS
	double SCALAR = SCALE_VALUE / log(1 + MAX*nScale);
	
  	for ( i = 0; i < nx; i++ )
	{
		for ( j = 0; j < (ny_half); j++)
		{
			if ( i < nx_half)
			{
				if (out[i*nyh+j][0] > 0)
				{
					setIntensity(ny_half-1+j,nx_half-i,SCALAR*log(1+out[i*nyh+j][0]*nScale));
					setIntensity(ny_half-j,nx_half-1+i,SCALAR*log(1+out[i*nyh+j][0]*nScale));
				}
				else
				{
					setIntensity(ny_half-1+j,nx_half-i,SCALAR*log(1-out[i*nyh+j][0]*nScale));
					setIntensity(ny_half-j,nx_half-1+i,SCALAR*log(1-out[i*nyh+j][0]*nScale));
				}
			}
			else
			{
				if (out[i*nyh+j][0] > 0)
				{
					setIntensity(ny_half-1+j,nx_half-1+(nx-i),SCALAR*log(1+out[i*nyh+j][0]*nScale));
					setIntensity(ny_half-j,nx_half-(nx-i),SCALAR*log(1+out[i*nyh+j][0]*nScale));
				}
				else
				{
					setIntensity(ny_half-1+j,nx_half-1+(nx-i),SCALAR*log(1-out[i*nyh+j][0]*nScale));
					setIntensity(ny_half-j,nx_half-(nx-i),SCALAR*log(1-out[i*nyh+j][0]*nScale));
				}
				
			}
			
		}
	}
	
	// CLEAR OUT MEMORY //
  	fftw_destroy_plan(plan_forward);
  	free(in);
  	fftw_free(out);
}

void CContrastImage::readImage(string sz)
// Reads in data from file
// Input: string(sz) - Filename
{
	Information new_file;
	new_file.readFile(sz);
	// SET DATA //
	setHeight((int)new_file.getN());
	setWidth((int)new_file.getM());
	setStep(new_file.getStep());
	setDataArray(new_file.getDataArray());
}

void CContrastImage::setIntensity(int y, int x, double value)
// Sets intensity at point
// Inputs: int(y) - Y COORDINATE
//	   int(x) - X COORDINATE
//	   double(value) - VALUE TO INSERT
{
	if(y >= 0 && y < m_iRows && x >= 0 && x < m_iColumns)
	{
		int position = (y*m_iRows) + x;
		m_dArray[position] = value;
	}
	else
	{
		string xS;
		string yS;
		NumberToString(x,xS);
		NumberToString(y,yS);
		
#pragma omp critical
		{
		output_content->push_back("Failed to write 2D array out of bounds: x =" + xS + ", y = " + yS);
		}
	}
}

double CContrastImage::getIntensity(int y, int x)
// Returns intensity at point
// Inputs: int(y) - Y coordinate
//	   int(x) - X coordinate
// Output: double - Value at point
{
	if(y >= 0 && y < m_iRows && x >= 0 && x < m_iColumns)
	{
		int position = (y*m_iRows) + x;
		return m_dArray[position];
	}
	else
	{
#pragma omp critical
		{
		output_content->push_back("Out of 2D Array Bounds returning 0.0");
		}
		return 0.0;
	}
}

void CContrastImage::saveImageMatrix(string sz)
// Saves image matrix (akin to .txt file)
// Input: string(sz) - Filename
{
	Information new_file;
	new_file.saveImageMatrix(sz,m_dArray,m_iColumns,m_iRows,m_fGridSize);
}

void CContrastImage::saveImage(string sz)
// Saves image (akin to .zcon file)
// Input: string(sz) - Filename
{
	Information new_file;
	new_file.saveImage(sz,m_dArray,m_iColumns,m_iRows,m_fGridSize);
}

void CContrastImage::saveCrossSection(const string &sz,const int &n)	
{
	vector<float> crossSection;
	for (int i = 0; i < m_iColumns; i++) 
		crossSection.push_back(getIntensity(n,i));
	
	Information new_file;
	new_file.saveImage(sz,crossSection,m_iColumns,1,m_fGridSize);
}

Information CContrastImage::getInformation(const string &filename)
// Returns data from file
// Output - Information class of file data
{
	Information new_file;
	new_file.setDataArray(m_dArray,m_iColumns,m_iRows,m_fGridSize);
	new_file.setFilename(filename);
	
	return new_file;
}
