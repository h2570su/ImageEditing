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
#include <ppl.h>
#include <iomanip>

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
			int result = (int)floor((this->getColor(j, i, R) * 0.299) + (this->getColor(j, i, G) * 0.587) + (this->getColor(j, i, B) * 0.114));
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
			//build 15-bit color
			uint16_t color = (this->getColor(j, i, R) & 0xF8) << 7;
			color |= (this->getColor(j, i, G) & 0xF8) << 2;
			color |= (this->getColor(j, i, B) & 0xF8) >> 3;
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

	//Whole image to closest color, using parallel for
	concurrency::parallel_for(0, this->width*this->height, [this, sorted_256_list](int idx)
	{
		//convert xy
		int i = idx / this->width;
		int j = idx % this->width;
		//build 15-bit color

		uint8_t red = this->getColor(j, i, R);
		uint8_t green = this->getColor(j, i, G);
		uint8_t blue = this->getColor(j, i, B);

		//find min Distance color
		int index = 0;
		double minDistance = INFINITY;

		int k = 0;
		uint16_t findColor = 0;
		for (const auto& c : sorted_256_list)
		{
			//bit mask out individual RGB
			int r = red - ((c.first & 0x7C00) >> 7);
			int g = green - ((c.first & 0x03E0) >> 2);
			int b = blue - ((c.first & 0x1F) << 3);
			double dist = r * r + g * g + b * b;
			if (dist < minDistance)
			{
				minDistance = dist;
				index = k;
				findColor = c.first;
				if (minDistance == 0.0)
				{
					break;
				}
			}
			k++;
		}
		//save color
		this->getColor(j, i, R) = (findColor & 0x7C00) >> 7; //R
		this->getColor(j, i, G) = (findColor & 0x03E0) >> 2; //G
		this->getColor(j, i, B) = (findColor & 0x001F) << 3; //B

	});

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
	const float threshold = 0.5;
	if (this->To_Grayscale())
	{
		//make a buffer
		float* temp = new float[this->width*this->height];

		//Copy display to Buffer and convert to float point with [0,1]
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				temp[i*width + j] = (float)this->getColor(j, i, R) / 255.0;
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
					float oldPixel = temp[y*this->width + x];
					float newPixel = (oldPixel > threshold) ? 1.0 : 0.0;
					float error = oldPixel - newPixel;

					temp[y*this->width + x] = newPixel;

					if (x + 1 < this->width)
					{
						temp[y*this->width + (x + 1)] += (error * (7.0 / 16.0));
					}

					if (y + 1 < this->height)
					{
						if (x - 1 >= 0)
						{
							temp[(y + 1)*this->width + (x - 1)] += (error * (3.0 / 16.0));
						}

						temp[(y + 1)*this->width + x] += (error * (5.0 / 16.0));
						if (x + 1 < this->width)
						{
							temp[(y + 1)*this->width + (x + 1)] += (error * (1.0 / 16.0));
						}
					}
				}

			}
			//odd line, right to left
			else
			{
				for (int x = this->width - 1; x >= 0; x--)
				{
					float oldPixel = temp[y*width + x];
					float newPixel = (oldPixel > threshold) ? 1.0 : 0.0;
					float error = oldPixel - newPixel;

					temp[y*this->width + x] = newPixel;

					if (x - 1 >= 0)
					{
						temp[y*this->width + (x - 1)] += (error * (7.0 / 16.0));
					}

					if (y + 1 < this->height)
					{
						if (x + 1 < this->width)
						{
							temp[(y + 1)*this->width + (x + 1)] += (error * (3.0 / 16.0));
						}
						temp[(y + 1)*this->width + x] += (error * (5.0 / 16.0));
						if (x - 1 >= 0)
						{
							temp[(y + 1)*this->width + (x - 1)] += (error * (1.0 / 16.0));
						}
					}

				}
			}
		}

		//Copy back to display
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{

				if (temp[i*this->width + j] > threshold)
				{
					temp[i*this->width + j] = 255;
				}
				else
				{
					temp[i*this->width + j] = 0;
				}

				this->getColor(j, i, R) = temp[i*this->width + j];
				this->getColor(j, i, G) = temp[i*this->width + j];
				this->getColor(j, i, B) = temp[i*this->width + j];
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
		uint64_t counter[256] = { 0 };
		uint64_t sum = 0;
		//Whole image
		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				//add to sum;
				sum += this->getColor(j, i, R);

				//add to indiv counter
				counter[this->getColor(j, i, R)]++;
			}
		}

		sum /= 255;


		int threshold = 255;
		uint64_t whiteCount = 0;
		for (threshold = 255; threshold >= 0; threshold--)
		{
			whiteCount += counter[threshold];
			if (whiteCount >= sum)
			{
				break;
			}
		}

		for (int i = 0; i < this->height; i++)
		{
			for (int j = 0; j < this->width; j++)
			{
				this->getColor(j, i, R) = (this->getColor(j, i, R) >= threshold) ? 255 : 0;
				this->getColor(j, i, G) = this->getColor(j, i, R);
				this->getColor(j, i, B) = this->getColor(j, i, R);
			}
		}

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
	const float RGshades[]{ 0 / 255.0, 36 / 255.0, 73 / 255.0, 109 / 255.0, 146 / 255.0, 182 / 255.0, 219 / 255.0, 255 / 255.0 };
	const float Bshades[]{ 0 / 255.0, 85 / 255.0, 170 / 255.0, 255 / 255.0 };

	float thresholdOffset = 0 / 255.0;
	float errorRate = 1;

	//make a buffer
	float* temp = new float[this->width*this->height * 3];

	//Copy display to Buffer and convert to float point with [0,1]
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			temp[(i*width + j) * 3] = (float)this->getColor(j, i, R) / 255.0;
			temp[(i*width + j) * 3 + 1] = (float)this->getColor(j, i, G) / 255.0;
			temp[(i*width + j) * 3 + 2] = (float)this->getColor(j, i, B) / 255.0;
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
				//record old color
				float oldPixelR = temp[(y*this->width + x) * 3];
				float oldPixelG = temp[(y*this->width + x) * 3 + 1];
				float oldPixelB = temp[(y*this->width + x) * 3 + 2];
				float newPixelR = 0;
				float newPixelG = 0;
				float newPixelB = 0;
				//pick new color
				//R
				for (int i = 1; i < 8; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelR >= RGshades[i - 1] && oldPixelR <= RGshades[i])
					{
						if (oldPixelR >= ((RGshades[i - 1] + RGshades[i]) / 2) + thresholdOffset)
						{
							newPixelR = RGshades[i];
						}
						else
						{
							newPixelR = RGshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelR > 1.0)
				{
					newPixelR = 1.0;
				}
				if (oldPixelR < 0)
				{
					newPixelR = 0;
				}

				//G
				for (int i = 1; i < 8; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelG >= RGshades[i - 1] && oldPixelG <= RGshades[i])
					{
						if (oldPixelG >= ((RGshades[i - 1] + RGshades[i]) / 2) + thresholdOffset)
						{
							newPixelG = RGshades[i];
						}
						else
						{
							newPixelG = RGshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelG > 1.0)
				{
					newPixelG = 1.0;
				}
				if (oldPixelG < 0)
				{
					newPixelG = 0;
				}

				//B
				for (int i = 1; i < 4; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelB >= Bshades[i - 1] && oldPixelB <= Bshades[i])
					{
						if (oldPixelB >= ((Bshades[i - 1] + Bshades[i]) / 2) + thresholdOffset)
						{
							newPixelB = Bshades[i];
						}
						else
						{
							newPixelB = Bshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelB > 1.0)
				{
					newPixelB = 1.0;
				}
				if (oldPixelB < 0)
				{
					newPixelB = 0;
				}

				float errorR = (oldPixelR - newPixelR)*errorRate;
				float errorG = (oldPixelG - newPixelG)*errorRate;
				float errorB = (oldPixelB - newPixelB)*errorRate;

				temp[(y*this->width + x) * 3] = newPixelR;
				temp[(y*this->width + x) * 3 + 1] = newPixelG;
				temp[(y*this->width + x) * 3 + 2] = newPixelB;

				if (x + 1 < this->width)
				{
					temp[(y*this->width + (x + 1)) * 3] += (errorR * (7.0 / 16.0));
					temp[(y*this->width + (x + 1)) * 3 + 1] += (errorG * (7.0 / 16.0));
					temp[(y*this->width + (x + 1)) * 3 + 2] += (errorB * (7.0 / 16.0));
				}

				if (y + 1 < this->height)
				{
					if (x - 1 >= 0)
					{
						temp[((y + 1)*this->width + (x - 1)) * 3] += (errorR * (3.0 / 16.0));
						temp[((y + 1)*this->width + (x - 1)) * 3 + 1] += (errorG * (3.0 / 16.0));
						temp[((y + 1)*this->width + (x - 1)) * 3 + 2] += (errorB * (3.0 / 16.0));
					}

					temp[((y + 1)*this->width + x) * 3] += (errorR * (5.0 / 16.0));
					temp[((y + 1)*this->width + x) * 3 + 1] += (errorG * (5.0 / 16.0));
					temp[((y + 1)*this->width + x) * 3 + 2] += (errorB * (5.0 / 16.0));

					if (x + 1 < this->width)
					{
						temp[((y + 1)*this->width + (x + 1)) * 3] += (errorR * (1.0 / 16.0));
						temp[((y + 1)*this->width + (x + 1)) * 3 + 1] += (errorG * (1.0 / 16.0));
						temp[((y + 1)*this->width + (x + 1)) * 3 + 2] += (errorB * (1.0 / 16.0));
					}
				}
			}

		}
		//odd line, right to left
		else
		{
			for (int x = this->width - 1; x >= 0; x--)
			{
				//record old color
				float oldPixelR = temp[(y*this->width + x) * 3];
				float oldPixelG = temp[(y*this->width + x) * 3 + 1];
				float oldPixelB = temp[(y*this->width + x) * 3 + 2];
				float newPixelR = 0;
				float newPixelG = 0;
				float newPixelB = 0;
				//pick new color
				//R
				for (int i = 1; i < 8; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelR >= RGshades[i - 1] && oldPixelR <= RGshades[i])
					{
						if (oldPixelR >= ((RGshades[i - 1] + RGshades[i]) / 2) + thresholdOffset)
						{
							newPixelR = RGshades[i];
						}
						else
						{
							newPixelR = RGshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelR > 1.0)
				{
					newPixelR = 1.0;
				}
				if (oldPixelR < 0)
				{
					newPixelR = 0;
				}

				//G
				for (int i = 1; i < 8; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelG >= RGshades[i - 1] && oldPixelG <= RGshades[i])
					{
						if (oldPixelG >= ((RGshades[i - 1] + RGshades[i]) / 2) + thresholdOffset)
						{
							newPixelG = RGshades[i];
						}
						else
						{
							newPixelG = RGshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelG > 1.0)
				{
					newPixelG = 1.0;
				}
				if (oldPixelG < 0)
				{
					newPixelG = 0;
				}

				//B
				for (int i = 1; i < 4; i++)
				{
					//if pixel in Range [i-1]to[i]
					if (oldPixelB >= Bshades[i - 1] && oldPixelB <= Bshades[i])
					{
						if (oldPixelB >= ((Bshades[i - 1] + Bshades[i]) / 2) + thresholdOffset)
						{
							newPixelB = Bshades[i];
						}
						else
						{
							newPixelB = Bshades[i - 1];
						}
						break;
					}
				}
				if (oldPixelB > 1.0)
				{
					newPixelB = 1.0;
				}
				if (oldPixelB < 0)
				{
					newPixelB = 0;
				}

				float errorR = (oldPixelR - newPixelR)*errorRate;
				float errorG = (oldPixelG - newPixelG)*errorRate;
				float errorB = (oldPixelB - newPixelB)*errorRate;

				temp[(y*this->width + x) * 3] = newPixelR;
				temp[(y*this->width + x) * 3 + 1] = newPixelG;
				temp[(y*this->width + x) * 3 + 2] = newPixelB;

				if (x - 1 >= 0)
				{
					temp[(y*this->width + (x - 1)) * 3] += (errorR * (7.0 / 16.0));
					temp[(y*this->width + (x - 1)) * 3 + 1] += (errorG * (7.0 / 16.0));
					temp[(y*this->width + (x - 1)) * 3 + 2] += (errorB * (7.0 / 16.0));
				}

				if (y + 1 < this->height)
				{
					if (x + 1 < this->width)
					{
						temp[((y + 1)*this->width + (x + 1)) * 3] += (errorR * (3.0 / 16.0));
						temp[((y + 1)*this->width + (x + 1)) * 3 + 1] += (errorG * (3.0 / 16.0));
						temp[((y + 1)*this->width + (x + 1)) * 3 + 2] += (errorB * (3.0 / 16.0));
					}

					temp[((y + 1)*this->width + x) * 3] += (errorR * (5.0 / 16.0));
					temp[((y + 1)*this->width + x) * 3 + 1] += (errorG * (5.0 / 16.0));
					temp[((y + 1)*this->width + x) * 3 + 2] += (errorB * (5.0 / 16.0));
					if (x - 1 >= 0)
					{
						temp[((y + 1)*this->width + (x - 1)) * 3] += (errorR * (1.0 / 16.0));
						temp[((y + 1)*this->width + (x - 1)) * 3 + 1] += (errorG * (1.0 / 16.0));
						temp[((y + 1)*this->width + (x - 1)) * 3 + 2] += (errorB * (1.0 / 16.0));
					}
				}

			}
		}
	}

	//Copy back to display
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			this->getColor(j, i, R) = temp[(i*this->width + j) * 3] * 255.0;
			this->getColor(j, i, G) = temp[(i*this->width + j) * 3 + 1] * 255.0;
			this->getColor(j, i, B) = temp[(i*this->width + j) * 3 + 2] * 255.0;
		}
	}

	//recycle memory
	delete[] temp;
	return true;

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
	uint64_t diff = 0;

	for (int i = 0; i < width * height * 4; i += 4)
	{
		unsigned char        rgb1[3];
		unsigned char        rgb2[3];

		RGBA_To_RGB(data + i, rgb1);
		RGBA_To_RGB(pImage->data + i, rgb2);

		data[i] = abs(rgb1[0] - rgb2[0]);
		data[i + 1] = abs(rgb1[1] - rgb2[1]);
		data[i + 2] = abs(rgb1[2] - rgb2[2]);
		if (data[i] != 0 || data[i + 1] != 0 || data[i + 2] != 0)
		{
			diff += data[i] + data[i + 1] + data[i + 2];
			/*data[i] = 255;
			data[i + 1] = 255;
			data[i + 2] = 255;*/
		}
		data[i + 3] = 255;

	}
	cout << "Diffence: " << diff << " val*Pixel, in: " << (double)diff / (this->width*this->height * 3 * 255) << " %" << endl;
	return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	//new image canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	//prevent over pixel be wrong color;
	this->invaildPixel = 0;
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//color sum in RGB
			int sum[3] = { 0 };
			//5x5 matrix
			for (int k = -2; k <= 2; k++)
			{
				for (int l = -2; l <= 2; l++)
				{
					//offseting xy
					int x = j + l;
					int y = i + k;
					//out of bound, reflect
					if (x < 0)
					{
						x *= -1;
					}
					if (x >= this->width)
					{
						int ref = x - this->width + 1;
						x -= 2 * ref;
					}

					if (y < 0)
					{
						y *= -1;
					}
					if (y >= this->height)
					{
						int ref = y - this->height + 1;
						y -= 2 * ref;
					}
					//summing
					sum[0] += this->getColor(x, y, R);
					sum[1] += this->getColor(x, y, G);
					sum[2] += this->getColor(x, y, B);

				}
			}
			//div by 5x5
			sum[0] /= 25;
			sum[1] /= 25;
			sum[2] /= 25;

			//to canvas
			canvas[(i*this->width + j) * 4] = sum[0];
			canvas[(i*this->width + j) * 4 + 1] = sum[1];
			canvas[(i*this->width + j) * 4 + 2] = sum[2];
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//copy canvas to display
	memcpy_s(this->data, this->width*this->height * 4, canvas, this->width*this->height * 4);
	//recycle memory
	delete[] canvas;
	return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	//Bartlett filter matrix
	double weightMat[5][5] =
	{
		{0.027,0.111,0.194,0.111,0.027},
		{0.111,0.361,0.611,0.361,0.111},
		{0.194,0.611,1.000,0.611,0.194},
		{0.111,0.361,0.611,0.361,0.111},
		{0.027,0.111,0.194,0.111,0.027}
	};
	/*
		{0.0277,0.1111,0.1944,0.1111,0.0277},
		{0.1111,0.3611,0.6111,0.3611,0.1111},
		{0.1944,0.6111,1.0000,0.6111,0.1944},
		{0.1111,0.3611,0.6111,0.3611,0.1111},
		{0.0277,0.1111,0.1944,0.1111,0.0277}

	*/
	double divider = 0;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			divider += weightMat[i][j];
		}
	}

	//new image canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	//prevent over pixel be wrong color;
	this->invaildPixel = 0;

	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//color sum in RGB
			int sum[3] = { 0 };
			//5x5 matrix
			for (int k = -2; k <= 2; k++)
			{
				for (int l = -2; l <= 2; l++)
				{

					//offseting xy
					int x = j + l;
					int y = i + k;
					//out of bound, reflect
					if (x < 0)
					{
						x *= -1;
					}
					if (x >= this->width)
					{
						int ref = x - this->width + 1;
						x -= 2 * ref;
					}

					if (y < 0)
					{
						y *= -1;
					}
					if (y >= this->height)
					{
						int ref = y - this->height + 1;
						y -= 2 * ref;
					}
					//summing
					sum[0] += this->getColor(x, y, R)*weightMat[l + 2][k + 2];
					sum[1] += this->getColor(x, y, G)*weightMat[l + 2][k + 2];
					sum[2] += this->getColor(x, y, B)*weightMat[l + 2][k + 2];

				}
			}
			//div by 5x5
			sum[0] /= divider;
			sum[1] /= divider;
			sum[2] /= divider;

			//to canvas
			canvas[(i*this->width + j) * 4] = sum[0];
			canvas[(i*this->width + j) * 4 + 1] = sum[1];
			canvas[(i*this->width + j) * 4 + 2] = sum[2];
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//copy canvas to display
	memcpy_s(this->data, this->width*this->height * 4, canvas, this->width*this->height * 4);
	//recycle memory
	delete[] canvas;
	return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	return Filter_Gaussian_N(5);
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////


bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
	//generate Gaussian filter matrix
	double** weightMat = new double*[N];
	for (int i = 0; i < N; i++)
	{
		weightMat[i] = new double[N];
	}

	//generate 1D binomial table
	uint64_t* bin = new uint64_t[N];
	for (int i = 0; i < N; i++)
	{
		bin[i] = Binomial(N - 1, i);
	}

	//calculate 2D Gaussian martix
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			weightMat[i][j] = bin[i] * bin[j];
		}
	}

	//recycle memory
	delete[] bin;



	//calculate divider
	double divider = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			divider += weightMat[i][j];
		}
	}
	//pre divide
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			weightMat[i][j] /= divider;
		}
	}


	//new image canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	//prevent over pixel be wrong color;
	this->invaildPixel = 0;

	int halfN = (int)(N / 2);
	Concurrency::parallel_for(0, this->width*this->height, [&](int idx)
	{
		int Xx = idx % this->width;
		int Yy = idx / this->width;
			//color sum in RGB
			double sum[3] = { 0 };
			//NxN matrix

			for (int k = -halfN; k <= halfN; k++)
			{
				for (int l = -halfN; l <= halfN; l++)
				{

					//offseting xy
					int x = Xx + l;
					int y = Yy + k;
					//out of bound, reflect
					if (x < 0)
					{
						x *= -1;
					}
					if (x >= this->width)
					{
						int ref = x - this->width + 1;
						x -= 2 * ref;
					}

					if (y < 0)
					{
						y *= -1;
					}
					if (y >= this->height)
					{
						int ref = y - this->height + 1;
						y -= 2 * ref;
					}
					//summing
					sum[0] += this->getColor(x, y, R)*weightMat[l + halfN][k + halfN];
					sum[1] += this->getColor(x, y, G)*weightMat[l + halfN][k + halfN];
					sum[2] += this->getColor(x, y, B)*weightMat[l + halfN][k + halfN];

				}
			}
			//div by 5x5
			//sum[0] /= divider;
			//sum[1] /= divider;
			//sum[2] /= divider;

			//to canvas
			canvas[(Yy*this->width + Xx) * 4] = sum[0];
			canvas[(Yy*this->width + Xx) * 4 + 1] = sum[1];
			canvas[(Yy*this->width + Xx) * 4 + 2] = sum[2];
			canvas[(Yy*this->width + Xx) * 4 + 3] = 255;
		
	});

	//copy canvas to display
	memcpy_s(this->data, this->width*this->height * 4, canvas, this->width*this->height * 4);
	//recycle memory
	delete[] canvas;
	for (int i = 0; i < N; i++)
	{
		delete[] weightMat[i];
	}
	delete[] weightMat;

	return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	//Bartlett filter matrix
	double weightMat[5][5] =
	{
		{0.027,0.111,0.194,0.111,0.027},
		{0.111,0.361,0.611,0.361,0.111},
		{0.194,0.611,1.000,0.611,0.194},
		{0.111,0.361,0.611,0.361,0.111},
		{0.027,0.111,0.194,0.111,0.027}
	};

	double divider = 0;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			divider += weightMat[i][j];
		}
	}

	//new image canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	//prevent over pixel be wrong color;
	this->invaildPixel = 0;

	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//color sum in RGB
			int sum[3] = { 0 };
			//5x5 matrix
			for (int k = -2; k <= 2; k++)
			{
				for (int l = -2; l <= 2; l++)
				{

					//offseting xy
					int x = j + l;
					int y = i + k;
					//out of bound, reflect
					if (x < 0)
					{
						x *= -1;
					}
					if (x >= this->width)
					{
						int ref = x - this->width + 1;
						x -= 2 * ref;
					}

					if (y < 0)
					{
						y *= -1;
					}
					if (y >= this->height)
					{
						int ref = y - this->height + 1;
						y -= 2 * ref;
					}
					//summing
					sum[0] += this->getColor(x, y, R)*weightMat[l + 2][k + 2];
					sum[1] += this->getColor(x, y, G)*weightMat[l + 2][k + 2];
					sum[2] += this->getColor(x, y, B)*weightMat[l + 2][k + 2];

				}
			}
			//div by 5x5
			sum[0] /= divider;
			sum[1] /= divider;
			sum[2] /= divider;

			//to canvas
			canvas[(i*this->width + j) * 4] = sum[0];
			canvas[(i*this->width + j) * 4 + 1] = sum[1];
			canvas[(i*this->width + j) * 4 + 2] = sum[2];
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//subtract source
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//to canvas
			int Rr = this->getColor(j, i, R) - canvas[(i*this->width + j) * 4];
			int Gg = this->getColor(j, i, G) - canvas[(i*this->width + j) * 4 + 1];
			int Bb = this->getColor(j, i, B) - canvas[(i*this->width + j) * 4 + 2];
			//if<0 to 0
			canvas[(i*this->width + j) * 4] = (Rr < 0) ? 0 : Rr;
			canvas[(i*this->width + j) * 4 + 1] = (Gg < 0) ? 0 : Gg;
			canvas[(i*this->width + j) * 4 + 2] = (Bb < 0) ? 0 : Bb;
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//copy canvas to display
	memcpy_s(this->data, this->width*this->height * 4, canvas, this->width*this->height * 4);
	//recycle memory
	delete[] canvas;
	return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	//Bartlett filter matrix
	double weightMat[5][5] =
	{
		{0.027,0.111,0.194,0.111,0.027},
		{0.111,0.361,0.611,0.361,0.111},
		{0.194,0.611,1.000,0.611,0.194},
		{0.111,0.361,0.611,0.361,0.111},
		{0.027,0.111,0.194,0.111,0.027}
	};

	double divider = 0;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			divider += weightMat[i][j];
		}
	}

	//new image canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	//prevent over pixel be wrong color;
	this->invaildPixel = 0;

	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//color sum in RGB
			int sum[3] = { 0 };
			//5x5 matrix
			for (int k = -2; k <= 2; k++)
			{
				for (int l = -2; l <= 2; l++)
				{

					//offseting xy
					int x = j + l;
					int y = i + k;
					//out of bound, reflect
					if (x < 0)
					{
						x *= -1;
					}
					if (x >= this->width)
					{
						int ref = x - this->width + 1;
						x -= 2 * ref;
					}

					if (y < 0)
					{
						y *= -1;
					}
					if (y >= this->height)
					{
						int ref = y - this->height + 1;
						y -= 2 * ref;
					}
					//summing
					sum[0] += this->getColor(x, y, R)*weightMat[l + 2][k + 2];
					sum[1] += this->getColor(x, y, G)*weightMat[l + 2][k + 2];
					sum[2] += this->getColor(x, y, B)*weightMat[l + 2][k + 2];

				}
			}
			//div by 5x5
			sum[0] /= divider;
			sum[1] /= divider;
			sum[2] /= divider;

			//to canvas
			canvas[(i*this->width + j) * 4] = sum[0];
			canvas[(i*this->width + j) * 4 + 1] = sum[1];
			canvas[(i*this->width + j) * 4 + 2] = sum[2];
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//subtract source then add it on source
	for (int i = 0; i < this->height; i++)
	{
		for (int j = 0; j < this->width; j++)
		{
			//to canvas
			int Rr = this->getColor(j, i, R) * 2 - canvas[(i*this->width + j) * 4];
			int Gg = this->getColor(j, i, G) * 2 - canvas[(i*this->width + j) * 4 + 1];
			int Bb = this->getColor(j, i, B) * 2 - canvas[(i*this->width + j) * 4 + 2];

			//ceiling

			Rr = (Rr > 255) ? 255 : Rr;
			Gg = (Gg > 255) ? 255 : Gg;
			Bb = (Bb > 255) ? 255 : Bb;

			//if<0 to 0
			canvas[(i*this->width + j) * 4] = (Rr < 0) ? 0 : Rr;
			canvas[(i*this->width + j) * 4 + 1] = (Gg < 0) ? 0 : Gg;
			canvas[(i*this->width + j) * 4 + 2] = (Bb < 0) ? 0 : Bb;
			canvas[(i*this->width + j) * 4 + 3] = 255;
		}
	}

	//copy canvas to display
	memcpy_s(this->data, this->width*this->height * 4, canvas, this->width*this->height * 4);
	//recycle memory
	delete[] canvas;
	return true;
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
	const float filter33[3][3] =
	{
		{1 / 16.0, 1 / 8.0, 1 / 16.0},
		{1 / 8.0, 1 / 4.0, 1 / 8.0},
		{1 / 16.0, 1 / 8.0, 1 / 16.0},
	};

	//cal new size
	int newWidth = (int)(this->width / 2);
	int newHeight = (int)(this->height / 2);

	//create new canvas
	uint8_t* canvas = new uint8_t[newWidth*newHeight * 4];
	memset(canvas, 0, sizeof(uint8_t)*newWidth*newHeight * 4);


	//process
	for (int y = 0; y < newHeight; y++)
	{
		for (int x = 0; x < newWidth; x++)
		{
			//a new pixel
			float sum[3] = { 0.0 };
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dx = -1; dx <= 1; dx++)
				{
					int srcX = (x * 2) + dx;
					int srcY = (y * 2) + dy;

					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}


					//summing
					sum[0] += this->getColor(srcX, srcY, R)*filter33[dx + 1][dy + 1];
					sum[1] += this->getColor(srcX, srcY, G)*filter33[dx + 1][dy + 1];
					sum[2] += this->getColor(srcX, srcY, B)*filter33[dx + 1][dy + 1];
				}

			}
			canvas[(y*newWidth + x) * 4] = sum[0];
			canvas[(y*newWidth + x) * 4 + 1] = sum[1];
			canvas[(y*newWidth + x) * 4 + 2] = sum[2];
			canvas[(y*newWidth + x) * 4 + 3] = 255;
		}
	}



	//delete original image
	delete[] this->data;
	//point data to new image
	this->data = canvas;
	this->width = newWidth;
	this->height = newHeight;

	return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	const float filter33[3][3] =
	{
		{1 / 16.0, 1 / 8.0, 1 / 16.0},
		{1 / 8.0, 1 / 4.0, 1 / 8.0},
		{1 / 16.0, 1 / 8.0, 1 / 16.0},
	};

	const float filter34[3][4] =
	{
		{1 / 32.0, 3 / 32.0, 3 / 32.0, 1 / 32.0},
		{2 / 32.0, 6 / 32.0, 6 / 32.0, 2 / 32.0},
		{1 / 32.0, 3 / 32.0, 3 / 32.0, 1 / 32.0},
	};

	const float filter44[4][4] =
	{
		{1 / 64.0, 3 / 64.0, 3 / 64.0, 1 / 64.0},
		{3 / 64.0, 9 / 64.0, 9 / 64.0, 3 / 64.0},
		{3 / 64.0, 9 / 64.0, 9 / 64.0, 3 / 64.0},
		{1 / 64.0, 3 / 64.0, 3 / 64.0, 1 / 64.0}
	};

	//cal new size
	int newWidth = (int)(this->width * 2);
	int newHeight = (int)(this->height * 2);

	//create new canvas
	uint8_t* canvas = new uint8_t[newWidth*newHeight * 4];
	memset(canvas, 0, sizeof(uint8_t)*newWidth*newHeight * 4);


	//process
	Concurrency::parallel_for(0, newWidth*newHeight, [&](int idx)
	{
		int x = idx % newWidth;
		int y = idx / newWidth;
		//a new pixel
		float sum[3] = { 0.0 };
		if (x % 2 == 0 && y % 2 == 0)
		{
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dx = -1; dx <= 1; dx++)
				{
					int srcX = (x / 2) + dx;
					int srcY = (y / 2) + dy;

					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}


					//summing
					sum[0] += this->getColor(srcX, srcY, R)*filter33[dx + 1][dy + 1];
					sum[1] += this->getColor(srcX, srcY, G)*filter33[dx + 1][dy + 1];
					sum[2] += this->getColor(srcX, srcY, B)*filter33[dx + 1][dy + 1];
				}
			}
		}
		else if (x % 2 == 1 && y % 2 == 0)
		{
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dx = -1; dx <= 2; dx++)
				{
					int srcX = (x / 2) + dx;
					int srcY = (y / 2) + dy;

					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}


					//summing
					sum[0] += this->getColor(srcX, srcY, R)*filter34[dy + 1][dx + 1];
					sum[1] += this->getColor(srcX, srcY, G)*filter34[dy + 1][dx + 1];
					sum[2] += this->getColor(srcX, srcY, B)*filter34[dy + 1][dx + 1];
				}
			}
		}
		else if (x % 2 == 0 && y % 2 == 1)
		{
			for (int dy = -1; dy <= 2; dy++)
			{
				for (int dx = -1; dx <= 1; dx++)
				{
					int srcX = (x / 2) + dx;
					int srcY = (y / 2) + dy;

					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}


					//summing
					sum[0] += this->getColor(srcX, srcY, R)*filter34[dx + 1][dy + 1];
					sum[1] += this->getColor(srcX, srcY, G)*filter34[dx + 1][dy + 1];
					sum[2] += this->getColor(srcX, srcY, B)*filter34[dx + 1][dy + 1];
				}
			}
		}
		else
		{
			for (int dy = -1; dy <= 2; dy++)
			{
				for (int dx = -1; dx <= 2; dx++)
				{
					int srcX = (x / 2) + dx;
					int srcY = (y / 2) + dy;

					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}


					//summing
					sum[0] += this->getColor(srcX, srcY, R)*filter44[dx + 1][dy + 1];
					sum[1] += this->getColor(srcX, srcY, G)*filter44[dx + 1][dy + 1];
					sum[2] += this->getColor(srcX, srcY, B)*filter44[dx + 1][dy + 1];
				}
			}
		}
		canvas[(y*newWidth + x) * 4] = sum[0];
		canvas[(y*newWidth + x) * 4 + 1] = sum[1];
		canvas[(y*newWidth + x) * 4 + 2] = sum[2];
		canvas[(y*newWidth + x) * 4 + 3] = 255;
	});




	//delete original image
	delete[] this->data;
	//point data to new image
	this->data = canvas;
	this->width = newWidth;
	this->height = newHeight;

	return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	float** filter = NULL;
	int filterWidthHalf = (scale * 2 - 1);
	if (filterWidthHalf < 1)
	{
		filterWidthHalf = 1;
	}

	//new matrix
	filter = new float*[filterWidthHalf * 2 + 1];
	for (int i = 0; i < filterWidthHalf * 2 + 1; i++)
	{
		filter[i] = new float[filterWidthHalf * 2 + 1];
	}

	//center
	filter[filterWidthHalf][filterWidthHalf] = 1 / 4.0;
	//top half
	for (int i = 0; i < filterWidthHalf; i++)
	{
		filter[filterWidthHalf][i] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*(i + 1.0));
		for (int j = 0; j < filterWidthHalf; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
		}
		for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
		}
	}
	//center line
	for (int j = 0; j < filterWidthHalf; j++)
	{
		filter[j][filterWidthHalf] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
	}
	for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
	{
		filter[j][filterWidthHalf] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
	}
	//bottom half
	for (int i = filterWidthHalf + 1; i < filterWidthHalf * 2 + 1; i++)
	{
		filter[filterWidthHalf][i] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - i));
		for (int j = 0; j < filterWidthHalf; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
		}
		for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
		}
	}

	//calculate matrix sum for divding usage
	float matrixSum = 0;
	for (int i = 0; i < 2 * filterWidthHalf + 1; i++)
	{
		for (int j = 0; j < 2 * filterWidthHalf + 1; j++)
		{
			matrixSum += filter[i][j];
		}
	}

	//DEBUG Output matrix
	/*cout << fixed << setprecision(7);
	for (int i = 0; i < 2 * filterWidthHalf + 1; i++)
	{
		for (int j = 0; j < 2 * filterWidthHalf + 1; j++)
		{
			cout << filter[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << "Sum: " << matrixSum << endl;*/
	//DEBUG




	//cal new size
	int newWidth = (int)(this->width*scale);
	int newHeight = (int)(this->height*scale);

	//create new canvas
	uint8_t* canvas = new uint8_t[newWidth*newHeight * 4];
	memset(canvas, 0, sizeof(uint8_t)*newWidth*newHeight * 4);


	//copy to bigger canvas
	Concurrency::parallel_for(0, newWidth*newHeight, [&](int idx)
	{
		int x = idx % newWidth;
		int y = idx / newWidth;
		int srcX = (x / scale);
		int srcY = (y / scale);
		/*if (ceil(srcX) == floor(srcX) && floor(srcY) == ceil(srcY))*/
		{
			canvas[(y*newWidth + x) * 4] = this->data[(int)(srcY*this->width + srcX) * 4];
			canvas[(y*newWidth + x) * 4 + 1] = this->data[(int)(srcY*this->width + srcX) * 4 + 1];
			canvas[(y*newWidth + x) * 4 + 2] = this->data[(int)(srcY*this->width + srcX) * 4 + 2];
			canvas[(y*newWidth + x) * 4 + 3] = 255;
		}

	});
	//create new canvas for filter
	uint8_t* canvas2 = new uint8_t[newWidth*newHeight * 4];
	memset(canvas2, 0, sizeof(uint8_t)*newWidth*newHeight * 4);

	//bartlett
	Concurrency::parallel_for(0, newWidth*newHeight, [&](int idx)
	{
		int x = idx % newWidth;
		int y = idx / newWidth;
		//a new pixel
		float sum[3] = { 0.0 };

		for (int dy = -filterWidthHalf; dy <= filterWidthHalf; dy++)
		{
			for (int dx = -filterWidthHalf; dx <= filterWidthHalf; dx++)
			{
				int srcX = x + dx;
				int srcY = y + dy;

				//out of bound, reflect
				if (srcX < 0)
				{
					srcX *= -1;
				}
				if (srcX >= newWidth)
				{
					int ref = srcX - newWidth + 1;
					srcX -= 2 * ref;
				}

				if (srcY < 0)
				{
					srcY *= -1;
				}
				if (srcY >= newHeight)
				{
					int ref = srcY - newHeight + 1;
					srcY -= 2 * ref;
				}


				//summing
				sum[0] += canvas[(srcX + newWidth * srcY) * 4] * filter[dx + filterWidthHalf][dy + filterWidthHalf];
				sum[1] += canvas[(srcX + newWidth * srcY) * 4 + 1] * filter[dx + filterWidthHalf][dy + filterWidthHalf];
				sum[2] += canvas[(srcX + newWidth * srcY) * 4 + 2] * filter[dx + filterWidthHalf][dy + filterWidthHalf];
			}

		}
		canvas2[(y*newWidth + x) * 4] = sum[0] / matrixSum;
		canvas2[(y*newWidth + x) * 4 + 1] = sum[1] / matrixSum;
		canvas2[(y*newWidth + x) * 4 + 2] = sum[2] / matrixSum;
		canvas2[(y*newWidth + x) * 4 + 3] = 255;

	});

	//delete matrix

	for (int i = 0; i < filterWidthHalf * 2 + 1; i++)
	{
		delete[] filter[i];
	}
	delete[] filter;

	//delete original image
	delete[] this->data;
	delete[] canvas;
	//point data to new image
	this->data = canvas2;
	this->width = newWidth;
	this->height = newHeight;

	return true;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
	float** filter = NULL;
	int filterWidthHalf = 1;

	//new matrix
	filter = new float*[filterWidthHalf * 2 + 1];
	for (int i = 0; i < filterWidthHalf * 2 + 1; i++)
	{
		filter[i] = new float[filterWidthHalf * 2 + 1];
	}

	//center
	filter[filterWidthHalf][filterWidthHalf] = 1 / 4.0;
	//top half
	for (int i = 0; i < filterWidthHalf; i++)
	{
		filter[filterWidthHalf][i] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*(i + 1.0));
		for (int j = 0; j < filterWidthHalf; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
		}
		for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
		}
	}
	//center line
	for (int j = 0; j < filterWidthHalf; j++)
	{
		filter[j][filterWidthHalf] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
	}
	for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
	{
		filter[j][filterWidthHalf] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
	}
	//bottom half
	for (int i = filterWidthHalf + 1; i < filterWidthHalf * 2 + 1; i++)
	{
		filter[filterWidthHalf][i] = filter[filterWidthHalf][filterWidthHalf] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - i));
		for (int j = 0; j < filterWidthHalf; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*(j + 1.0));
		}
		for (int j = filterWidthHalf + 1; j < filterWidthHalf * 2 + 1; j++)
		{
			filter[j][i] = filter[filterWidthHalf][i] * (1.0 / (filterWidthHalf + 1)*((2 * filterWidthHalf + 1) - j));
		}
	}

	//calculate matrix sum for divding usage
	float matrixSum = 0;
	for (int i = 0; i < 2 * filterWidthHalf + 1; i++)
	{
		for (int j = 0; j < 2 * filterWidthHalf + 1; j++)
		{
			matrixSum += filter[i][j];
		}
	}

	//DEBUG Output matrix
	/*cout << fixed << setprecision(7);
	for (int i = 0; i < 2 * filterWidthHalf + 1; i++)
	{
		for (int j = 0; j < 2 * filterWidthHalf + 1; j++)
		{
			cout << filter[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << "Sum: " << matrixSum << endl;*/
	//DEBUG




	float theta = angleDegrees * 2 * acos(-1) / 360.0;

	//create new canvas
	uint8_t* canvas2 = new uint8_t[this->width*this->height * 4];
	memset(canvas2, 0, sizeof(uint8_t)*this->width*this->height * 4);


	//process
	for (int y = 0; y < this->height; y++)
	{
		for (int x = 0; x < this->width; x++)
		{
			int srcX = ((x - (this->width / 2))*cos(-theta) - (y - (this->height / 2)) * sin(-theta)) + (this->width / 2);
			int srcY = ((x - (this->width / 2))*sin(-theta) + (y - (this->height / 2)) * cos(-theta)) + (this->height / 2);

			if (srcX >= 0 && srcX < this->width&&srcY >= 0 && srcY < this->height)
			{
				canvas2[(y*this->width + x) * 4] = this->getColor(srcX, srcY, R);
				canvas2[(y*this->width + x) * 4 + 1] = this->getColor(srcX, srcY, G);
				canvas2[(y*this->width + x) * 4 + 2] = this->getColor(srcX, srcY, B);
				canvas2[(y*this->width + x) * 4 + 3] = 255;
			}
		}
	}

	//create new canvas
	uint8_t* canvas = new uint8_t[this->width*this->height * 4];
	memset(canvas, 0, sizeof(uint8_t)*this->width*this->height * 4);

	for (int y = 0; y < this->height; y++)
	{
		for (int x = 0; x < this->width; x++)
		{
			//a new pixel
			float sum[3] = { 0.0 };
			for (int dy = -filterWidthHalf; dy <= filterWidthHalf; dy++)
			{
				for (int dx = -filterWidthHalf; dx <= filterWidthHalf; dx++)
				{
					int srcX = x + dx;
					int srcY = y + dy;
					//out of bound, reflect
					if (srcX < 0)
					{
						srcX *= -1;
					}
					if (srcX >= this->width)
					{
						int ref = srcX - this->width + 1;
						srcX -= 2 * ref;
					}

					if (srcY < 0)
					{
						srcY *= -1;
					}
					if (srcY >= this->height)
					{
						int ref = srcY - this->height + 1;
						srcY -= 2 * ref;
					}

					//summing
					sum[0] += canvas2[(srcX + this->width*srcY) * 4] * filter[dx + filterWidthHalf][dy + filterWidthHalf];
					sum[1] += canvas2[(srcX + this->width*srcY) * 4 + 1] * filter[dx + filterWidthHalf][dy + filterWidthHalf];
					sum[2] += canvas2[(srcX + this->width*srcY) * 4 + 2] * filter[dx + filterWidthHalf][dy + filterWidthHalf];


				}

			}
			canvas[(y*this->width + x) * 4] = sum[0] / matrixSum;
			canvas[(y*this->width + x) * 4 + 1] = sum[1] / matrixSum;
			canvas[(y*this->width + x) * 4 + 2] = sum[2] / matrixSum;
			canvas[(y*this->width + x) * 4 + 3] = 255;
		}
	}

	//delete matrix

	for (int i = 0; i < filterWidthHalf * 2 + 1; i++)
	{
		delete[] filter[i];
	}
	delete[] filter;

	//delete original image
	delete[] this->data;
	delete[] canvas2;
	//point data to new image
	this->data = canvas;

	//this->Filter_Bartlett();
	return true;
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

