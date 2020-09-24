///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>

using namespace std;

// constants
const int           RED = 0;                // red channel
const int           GREEN = 1;                // green channel
const int           BLUE = 2;                // blue channel
const unsigned char BACKGROUND[3] = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
	double        res;

	res = 1;
	for (int i = 1; i <= s; i++)
		res = (n - i + 1) * res / i;

	return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
	data = new unsigned char[width * height * 4];
	ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
	int i;

	width = w;
	height = h;
	data = new unsigned char[width * height * 4];

	for (i = 0; i < width * height * 4; i++)
		data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image)
{
	width = image.width;
	height = image.height;
	data = NULL;
	if (image.data != NULL) {
		data = new unsigned char[width * height * 4];
		memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
	}
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
	if (data)
		delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
	unsigned char   *rgb = new unsigned char[width * height * 3];
	int		    i, j;

	if (!data)
		return NULL;

	// Divide out the alpha
	for (i = 0; i < height; i++)
	{
		int in_offset = i * width * 4;
		int out_offset = i * width * 3;

		for (j = 0; j < width; j++)
		{
			RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
		}
	}

	return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
	TargaImage	*out_image = Reverse_Rows();

	if (!out_image)
		return false;

	if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
	{
		cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
		return false;
	}

	delete out_image;

	return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
	unsigned char   *temp_data;
	TargaImage	    *temp_image;
	TargaImage	    *result;
	int		        width, height;

	if (!filename)
	{
		cout << "No filename given." << endl;
		return NULL;
	}// if

	temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
	if (!temp_data)
	{
		cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
		width = height = 0;
		return NULL;
	}
	temp_image = new TargaImage(width, height, temp_data);
	free(temp_data);

	result = temp_image->Reverse_Rows();

	delete temp_image;

	return result;
}// Load_Image


//By hahaaga



///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	//Whole image
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{

			//Every pixel to gray;
			int result = (int)floor((this->getColor(j, i, R) * 0.299f) + (this->getColor(j, i, G) * 0.587f) + (this->getColor(j, i, B) * 0.114f));
			this->getColor(j, i, R) = result;
			this->getColor(j, i, G) = result;
			this->getColor(j, i, B) = result;
		}
	}
	return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	//Whole image
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{

			//every pixel R->3bit G->3bit B->2bit

			this->getColor(j, i, R) &= 0xE0;

			this->getColor(j, i, G) &= 0xE0;

			this->getColor(j, i, B) &= 0xC0;
		}
	}
	return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	std::map<uint16_t, uint64_t> list;//14-10R 5-9G  0-4B 
	//Whole image to 5-bit color space and record in color set
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{

			//every pixel to 5-bit

			this->getColor(j, i, R) &= 0xF8;

			this->getColor(j, i, G) &= 0xF8;

			this->getColor(j, i, B) &= 0xF8;
			//build 15-bit color
			uint16_t color = this->getColor(j, i, R) << 7;
			color |= this->getColor(j, i, G) << 2;
			color |= this->getColor(j, i, B) >> 3;
			//find in list  if not exist init to 1
			auto iterFindInlist = list.find(color);
			if (iterFindInlist == list.end())
			{
				list[color] = 1;
			}
			//count++
			else
			{
				(*iterFindInlist).second++;
			}

		}
	}

	//256 popular colur list
	vector<pair<uint16_t, uint64_t>> sorted_256_list;
	//copy map to vector
	for (auto v : list)
	{
		sorted_256_list.push_back(v);
	}
	//sort by bigger front
	sort(sorted_256_list.begin(), sorted_256_list.end(),
		[](const pair<uint16_t, uint64_t>& p1, const pair<uint16_t, uint64_t>& p2)
	{
		return p1.second > p2.second;
	}
	);
	// >256 colors remove
	if (sorted_256_list.size() > 256)
	{
		sorted_256_list.erase(sorted_256_list.begin() + 256, sorted_256_list.end());
	}

	//Whole image to closest color
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//build 15-bit color
			uint8_t red = this->getColor(j, i, R) >> 3;
			uint8_t green = this->getColor(j, i, G) >> 3;
			uint8_t blue = this->getColor(j, i, B) >> 3;

			//find min Distance color
			int index = 0;
			double minDistance = INFINITY;

			int k = 0;
			for (const auto& c : sorted_256_list)
			{
				//bit mask out individual RGB
				double dist = pow((red - ((c.first & 0x7C00) >> 10)), 2) + pow(green - ((c.first & 0x03E0) >> 5), 2) + pow(blue - (c.first & 0x1F), 2);
				if (dist < minDistance)
				{
					minDistance = dist;
					index = k;
					if (minDistance == 0.0)
					{
						break;
					}
				}
				k++;
			}
			//save color
			this->getColor(j, i, R) = (sorted_256_list[index].first & 0x7C00) >> 7; //R
			this->getColor(j, i, G) = (sorted_256_list[index].first & 0x03E0) >> 2; //G
			this->getColor(j, i, B) = (sorted_256_list[index].first & 0x001F) << 3; //B
		}
	}

	return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	if (this->To_Grayscale())
	{
		//Whole image to black or white threshold 127(half)
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				//if pixel >=128  convert to 255
				if (this->getColor(j, i, R) > 127)
				{
					this->getColor(j, i, R) = 255;
					this->getColor(j, i, G) = 255;
					this->getColor(j, i, B) = 255;
				}
				//or to 0
				else
				{
					this->getColor(j, i, R) = 0;
					this->getColor(j, i, G) = 0;
					this->getColor(j, i, B) = 0;
				}
			}
		}
		return true;
	}
	else
	{
		ClearToBlack();
		return false;
	}
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	srand(time(NULL));
	if (this->To_Grayscale())
	{
		//Whole image add random value(-51,51) first, then threshold by 0.5
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{

				int v = this->getColor(j, i, R) + ((rand() % 103) - 51);
				if (v < 0)
				{
					v = 0;
				}
				if (v > 255)
				{
					v = 255;
				}
				//if pixel >=128  convert to 255
				if (v > 127)
				{
					this->getColor(j, i, R) = 255;
					this->getColor(j, i, G) = 255;
					this->getColor(j, i, B) = 255;
				}
				//or to 0
				else
				{
					this->getColor(j, i, R) = 0;
					this->getColor(j, i, G) = 0;
					this->getColor(j, i, B) = 0;
				}
			}
		}
		return true;
	}
	else
	{
		ClearToBlack();
		return false;
	}
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	if (this->To_Grayscale())
	{
		//make a buffer
		double* temp = new double[this->width*this->height];

		//Copy display to Buffer and convert to float point with [0,1]
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				temp[i*width + j] = (double)this->getColor(j, i, R) / 255.0;
			}
		}


		//Floyd-Steinberg
		for (int y = 0; y < this->height; y++)
		{
			//even line, left to right
			if (y % 2 == 0)
			{
				for (int x = 0; x < this->width; x++)
				{
					double oldPixel = temp[y*width + x];
					double newPixel = (oldPixel > 0.5) ? 1 : 0;
					double error = oldPixel - newPixel;

					temp[y*width + x] = newPixel;

					if (x < this->width - 1)
					{
						temp[y*width + (x + 1)] += error * (7.0 / 16.0);
					}
					if (y < this->height - 1)
					{
						temp[(y + 1)*width + x] += error * (5.0 / 16.0);
						if (x > 0)
						{
							temp[(y + 1)*width + (x - 1)] += error * (3.0 / 16.0);
						}
						if (x < this->width - 1)
						{
							temp[(y + 1)*width + (x + 1)] += error * (1.0 / 16.0);
						}
					}
				}
			}
			//odd line, right to left
			else
			{
				for (int x = this->width - 1; x >= 0; x--)
				{
					double newPixel = (temp[y*width + x] > 0.5) ? 1 : 0;
					double error = temp[y*width + x] - newPixel;

					temp[y*width + x] = newPixel;

					if (x < this->width - 1)
					{
						temp[y*width + (x + 1)] += error * (7.0 / 16.0);
					}
					if (y < this->height - 1)
					{
						temp[(y + 1)*width + x] += error * (5.0 / 16.0);
						if (x > 0)
						{
							temp[(y + 1)*width + (x - 1)] += error * (3.0 / 16.0);
						}
						if (x < this->width - 1)
						{
							temp[(y + 1)*width + (x + 1)] += error * (1.0 / 16.0);
						}
					}
				}
			}

			//Copy back to display
			for (int i = 0; i < this->height; i++)
			{
				for (int j = 0; j < this->width; j++)
				{
					
					this->getColor(j, i, R) = temp[i*this->width + j] * 255;
					this->getColor(j, i, G) = temp[i*this->width + j] * 255;
					this->getColor(j, i, B) = temp[i*this->width + j] * 255;
				}
			}

		}
		//recycle memory
		delete[] temp;
		return true;
	}
	else
	{
		return false;
	}
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	if (this->To_Grayscale())
	{
		//measure brightness
		uint64_t sum = 0;
		//Whole image
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				//add to sum;
				sum += this->getColor(j, i, R);
			}
		}

		//binary serach
		int upBound = 255;
		int downBound = 0;
		uint8_t* tempData = new uint8_t[this->width*this->height * 4];
		do
		{

			memcpy_s(tempData, this->width*this->height * 4, this->data, this->width*this->height * 4);
			uint64_t  newSum = 0;

			//Whole image to mid
			for (int i = 0; i < this->height; i++)
			{
				for (int j = 0; j < this->width; j++)
				{
					//if pixel >=128  convert to 255
					if (tempData[(i*this->width + j) * 4] > (upBound + downBound) / 2)
					{
						tempData[(i*this->width + j) * 4] = 255;
						tempData[(i*this->width + j) * 4 + 1] = 255;
						tempData[(i*this->width + j) * 4 + 2] = 255;
					}
					//or to 0
					else
					{
						tempData[(i*this->width + j) * 4] = 0;
						tempData[(i*this->width + j) * 4 + 1] = 0;
						tempData[(i*this->width + j) * 4 + 2] = 0;
					}
				}
			}

			//cal Brightness
			for (int i = 0; i < this->height; i++)
			{
				for (int j = 0; j < this->width; j++)
				{
					newSum += tempData[(i*this->width + j) * 4];
				}
			}

			if (newSum < sum)
			{
				upBound = (upBound + downBound) / 2 - 1;
			}
			else if (newSum > sum)
			{
				downBound = (upBound + downBound) / 2 + 1;
			}
			else
			{
				break;
			}


		} while (downBound <= upBound);


		memcpy_s(this->data, this->width*this->height * 4, tempData, this->width*this->height * 4);
		delete[] tempData;


		return true;
	}
	else
	{
		ClearToBlack();
		return false;
	}
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	if (this->To_Grayscale())
	{
		//The partten matrix
		double thr_mask[4][4] =
		{
			{
				0.7059, 0.0588, 0.4706, 0.1765
			},

			{
				0.3529, 0.9412, 0.7647, 0.5294
			},

			{
				0.5882, 0.8235, 0.8824, 0.2941
			},

			{
				0.2353, 0.411833335655526322494068836022, 0.1176, 0.6471
			}
		};
		//Whole image
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				//normalize this pixel's value
				double normalVal = this->getColor(j, i, R);
				normalVal /= 255.0;

				//bright than threshold
				if (normalVal >= thr_mask[j % 4][i % 4] * 1.0664889125)
				{
					this->getColor(j, i, R) = 255;
					this->getColor(j, i, G) = 255;
					this->getColor(j, i, B) = 255;
				}
				else
				{
					this->getColor(j, i, R) = 0;
					this->getColor(j, i, G) = 0;
					this->getColor(j, i, B) = 0;
				}

			}
		}

		return true;
	}
	else
	{
		return false;
	}
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	ClearToBlack();
	return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Over: Images not the same size\n";
		return false;
	}

	ClearToBlack();
	return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_In: Images not the same size\n";
		return false;
	}

	ClearToBlack();
	return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Out: Images not the same size\n";
		return false;
	}

	ClearToBlack();
	return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Atop: Images not the same size\n";
		return false;
	}

	ClearToBlack();
	return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
	if (width != pImage->width || height != pImage->height)
	{
		cout << "Comp_Xor: Images not the same size\n";
		return false;
	}

	ClearToBlack();
	return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
	if (!pImage)
		return false;

	if (width != pImage->width || height != pImage->height)
	{
		cout << "Difference: Images not the same size\n";
		return false;
	}// if

	for (int i = 0; i < width * height * 4; i += 4)
	{
		unsigned char        rgb1[3];
		unsigned char        rgb2[3];

		RGBA_To_RGB(data + i, rgb1);
		RGBA_To_RGB(pImage->data + i, rgb2);

		data[i] = abs(rgb1[0] - rgb2[0]);
		data[i + 1] = abs(rgb1[1] - rgb2[1]);
		data[i + 2] = abs(rgb1[2] - rgb2[2]);
		data[i + 3] = 255;
	}

	return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	ClearToBlack();
	return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	ClearToBlack();
	return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	ClearToBlack();
	return false;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
	ClearToBlack();
	return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	ClearToBlack();
	return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	ClearToBlack();
	return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
	ClearToBlack();
	return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
	ClearToBlack();
	return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	ClearToBlack();
	return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	ClearToBlack();
	return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
	ClearToBlack();
	return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0; i < 3; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
	unsigned char   *dest = new unsigned char[width * height * 4];
	TargaImage	    *result;
	int 	        i, j;

	if (!data)
		return NULL;

	for (i = 0; i < height; i++)
	{
		int in_offset = (height - i - 1) * width * 4;
		int out_offset = i * width * 4;

		for (j = 0; j < width; j++)
		{
			dest[out_offset + j * 4] = data[in_offset + j * 4];
			dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
			dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
			dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
		}
	}

	result = new TargaImage(width, height, dest);
	delete[] dest;
	return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
	memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;
			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
				int dist_squared = x_off * x_off + y_off * y_off;
				if (dist_squared <= radius_squared) {
					data[(y_loc * width + x_loc) * 4 + 0] = s.r;
					data[(y_loc * width + x_loc) * 4 + 1] = s.g;
					data[(y_loc * width + x_loc) * 4 + 2] = s.b;
					data[(y_loc * width + x_loc) * 4 + 3] = s.a;
				}
				else if (dist_squared == radius_squared + 1) {
					data[(y_loc * width + x_loc) * 4 + 0] =
						(data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
					data[(y_loc * width + x_loc) * 4 + 1] =
						(data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
					data[(y_loc * width + x_loc) * 4 + 2] =
						(data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
					data[(y_loc * width + x_loc) * 4 + 3] =
						(data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
				}
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
	radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}

