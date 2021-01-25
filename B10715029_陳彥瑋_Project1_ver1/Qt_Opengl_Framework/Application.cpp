#include "Application.h"
#include "qt_opengl_framework.h"
#include <vector>
#include <algorithm>
#include <functional>

Application::Application()
{

}
Application::~Application()
{

}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene( void )
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	
}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
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

		for (i = 0 ; i < 3 ; i++)
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
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];

			for (int k=0; k<3; k++)
				img_data[offset_rgba+k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			//unsigned char uni = rgb[offset_rgb + rr] / 32.0 + rgb[offset_rgb + gg] / 32.0 + rgb[offset_rgb + bb] / 64.0;

			//for (int k = 0; k < 3; k++)
				//img_data[offset_rgba + k] = uni;
			img_data[offset_rgba + rr] = rgb[offset_rgb + rr] >> 5 << 5;
			img_data[offset_rgba + gg] = rgb[offset_rgb + gg] >> 5 << 5;
			img_data[offset_rgba + bb] = rgb[offset_rgb + bb] >> 6 << 6;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	//*add

	//make histogram & storage 5bit data
	long long int histogram[32768] = { 0 };
	unsigned char *rgb_5bit = new unsigned char[img_height * img_width * 3];
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			
			rgb_5bit[offset_rgb + rr] = rgb[offset_rgb + rr] >> 3;
			rgb_5bit[offset_rgb + gg] = rgb[offset_rgb + gg] >> 3;
			rgb_5bit[offset_rgb + bb] = rgb[offset_rgb + bb] >> 3;

			histogram[rgb_5bit[offset_rgb + rr] * 32 * 32 + rgb_5bit[offset_rgb + gg] * 32 + rgb_5bit[offset_rgb + bb]] ++;
		}
	}

	//find most popular 256 color & convert index to 5bit rgb
	unsigned int popular256[3 * 256];

	for (int i = 0; i < 256; i++)
	{
		long long int highestVote = -1;
		int index = 0;
		for (int j = 0; j < 32768; j++)
		{
			if (histogram[j] > highestVote)
			{
				highestVote = histogram[j];
				index = j;
			}
		}

		histogram[index] = -2; //mark as selected

		popular256[i * 3 + rr] = index / 32 / 32; // red
		popular256[i * 3 + gg] = (index / 32) % 32; // green
		popular256[i * 3 + bb] = index % 32; //blue
		
	}

	//find nearest color & assign
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			
			double nearestDistance = 99999; //sqrt( (32-0)^2 + (32-0)^2 + (32-0)^2 ) = 55.425...
			short index = 0;
			for (int k = 0; k < 256; k++)
			{
				double currDist = sqrt(pow((rgb_5bit[offset_rgb + rr] - popular256[k * 3 + rr]), 2)
									 + pow((rgb_5bit[offset_rgb + gg] - popular256[k * 3 + gg]), 2)
									 + pow((rgb_5bit[offset_rgb + bb] - popular256[k * 3 + bb]), 2));
				if (currDist <= nearestDistance)
				{
					nearestDistance = currDist;
					index = k;
				}
			}

			img_data[offset_rgba + rr] = popular256[index * 3 + rr] << 3;
			img_data[offset_rgba + gg] = popular256[index * 3 + gg] << 3;
			img_data[offset_rgba + bb] = popular256[index * 3 + bb] << 3;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb_5bit;
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			double intensity = gray / (double)256;

			if (intensity >= 0.5)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = WHITE;
					
			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = BLACK;
			}
			
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	srand(time(NULL));

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			double intensity = gray / (double)256;
			intensity += ((rand() % 5) / (double)10) - 0.2;

			if (intensity >= 0.5)
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = WHITE;

			}
			else
			{
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = BLACK;
			}

			img_data[offset_rgba + aa] = WHITE;
		}
	}
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	double *error = new double[img_height * img_width];
	for (int i = 0; i < (img_height * img_width); i++)
	{
		error[i] = 0;
	}

	bool toRight = true; // for zig-zag. to right=true, to left = false
	int j;


	for (int i = 0; i < img_height; i++)
	{
		if (toRight)
		{
			j = 0;
		}
		else
		{
			j = img_width - 1;
		}

		while (j < img_width && j >= 0)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			double intensity = gray / (double)256 + error[i * img_width + j];

			//fill color
			bool newColor;
			if (intensity >= 0.5)
			{
				newColor = 1;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = WHITE;
			}
			else
			{
				newColor = 0;
				for (int k = 0; k < 3; k++)
					img_data[offset_rgba + k] = BLACK;
			}
			img_data[offset_rgba + aa] = WHITE;

			//evaluate new error
			double newError = intensity - newColor;

			bool right, down, left; // check position can move to or not
			if (j + 1 >= img_width)
			{
				right = false;
			}
			else
			{
				right = true;
			}
			if (i + 1 >= img_height)
			{
				down = false;
			}
			else
			{
				down = true;
			}
			if (j - 1 < 0)
			{
				left = false;
			}
			else
			{
				left = true;
			}

			// fill error
			if (down == false)
			{
				//at the buttom of the image

				if ( (right == false && toRight == true) || (left == false && toRight == false) )
				{
					// at the end of the image
					
				}
				else
				{
					if (toRight == true)
					{
						//next is going to right
						error[i * img_width + (j + 1)] += newError; //other road is closed, get all the error
					}
					else
					{
						//next is going to left
						error[i * img_width + (j - 1)] += newError; //other road is closed, get all the error
					}
				}
			}
			else
			{
				//can move to down

				if (right == false && toRight == true)
				{
					// at the right border
					
					error[(i + 1) * img_width + j] += newError * 10 / 16; //down get 10/16 error
					error[(i + 1) * img_width + (j - 1)] += newError * 6 / 16; //down-left get 6/16 error
				}
				else if (left == false && toRight == false)
				{
					// at the left border

					error[(i + 1) * img_width + j] += newError * 10 / 16; //down get 10/16 error
					error[(i + 1) * img_width + (j + 1)] += newError * 6 / 16; //down-right get 6/16 error
				}
				else
				{
					// at the middle

					if (toRight == true)
					{
						//next move to right
						error[i * img_width + (j + 1)] += newError * 7 / 16; //right get 7/16 error
						error[(i + 1) * img_width + (j + 1)] += newError * 1 / 16; //down-right get 1/16 error
						error[(i + 1) * img_width + j] += newError * 5 / 16; //down get 5/16 error
						error[(i + 1) * img_width + (j - 1)] += newError * 3 / 16; //down-left get 3/16 error
					}
					else
					{
						//next move to left
						error[i * img_width + (j - 1)] += newError * 7 / 16; //left get 7/16 error
						error[(i + 1) * img_width + (j - 1)] += newError * 1 / 16; //down-left get 1/16 error
						error[(i + 1) * img_width + j] += newError * 5 / 16; //down get 5/16 error
						error[(i + 1) * img_width + (j + 1)] += newError * 3 / 16; //down-right get 3/16 error
					}
				}
			}

			//update
			if (toRight)
			{
				j++;
			}
			else
			{
				j--;
			}
		}

		if (toRight)
		{
			toRight = false;
		}
		else
		{
			toRight = true;
		}
	}

	delete[] error;
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	double *intensity_copy = new double[img_height * img_width];
	double *intensity_sortcopy = new double[img_height * img_width];
	double intensity_sum = 0;

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			double intensity = gray / (double)256;
			intensity_sum += intensity;
			intensity_copy[i * img_width + j] = intensity;
			intensity_sortcopy[i * img_width + j] = intensity;
		}
	}

	double intensity_avg = intensity_sum / (double)(img_width * img_height);
	int medianIndex = img_height * img_width * intensity_avg;

	std::sort(intensity_sortcopy, intensity_sortcopy + (img_height * img_width), std::greater<double>());

	double threshold = intensity_sortcopy[medianIndex];

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			
			if (intensity_copy[i * img_width + j] > threshold)
			{
				img_data[offset_rgba + rr] = WHITE;
				img_data[offset_rgba + gg] = WHITE;
				img_data[offset_rgba + bb] = WHITE;
			}
			else
			{
				img_data[offset_rgba + rr] = BLACK;
				img_data[offset_rgba + gg] = BLACK;
				img_data[offset_rgba + bb] = BLACK;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] intensity_sortcopy;
	delete[] intensity_copy;
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	double *intensity_copy = new double[img_height * img_width];
	const double mask[4][4] = { 0.7059, 0.3529, 0.5882, 0.2353,
								0.0588, 0.9412, 0.8235, 0.4118,
								0.4706, 0.7647, 0.8824, 0.1176,
								0.1765, 0.5294, 0.2941, 0.6471 };

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			unsigned char gray = 0.299 * rgb[offset_rgb + rr] + 0.587 * rgb[offset_rgb + gg] + 0.114 * rgb[offset_rgb + bb];
			double intensity = gray / (double)256;

			intensity_copy[i * img_width + j] = intensity;
		}
	}

	int chunk_i = img_height / 4;
	if (img_height % 4 != 0)
	{
		chunk_i++;
	}
	int chunk_j = img_width / 4;
	if (img_width % 4 != 0)
	{
		chunk_j++;
	}

	for (int i = 0; i < chunk_i; i++)
	{
		for (int j = 0; j < chunk_j; j++)
		{
			for (int m = 0; m < 4; m++)
			{
				for (int n = 0; n < 4; n++)
				{
					if (((i * 4 + m) < img_height) && ((j * 4 + n) < img_width))
					{
						int offset_rgba = (i * 4 + m) * img_width * 4 + (j * 4 + n) * 4;

						if (intensity_copy[(i * 4 + m) * img_width + (j * 4 + n)] > mask[m][n])
						{
							img_data[offset_rgba + rr] = WHITE;
							img_data[offset_rgba + gg] = WHITE;
							img_data[offset_rgba + bb] = WHITE;
						}
						else
						{
							img_data[offset_rgba + rr] = BLACK;
							img_data[offset_rgba + gg] = BLACK;
							img_data[offset_rgba + bb] = BLACK;
						}

						img_data[offset_rgba + aa] = WHITE;
					}
				}
			}
		}
	}

	delete[] intensity_copy;
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();

	//*add
	unsigned char redTable[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
	unsigned char greenTable[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
	unsigned char blueTable[4] = { 0, 85, 170, 255 };

	double *redError = new double[img_height * img_width];
	double *greenError = new double[img_height * img_width];
	double *blueError = new double[img_height * img_width];
	for (int i = 0; i < (img_height * img_width); i++)
	{
		redError[i] = 0;
		greenError[i] = 0;
		blueError[i] = 0;
	}

	bool toRight = true; // for zig-zag. to right=true, to left = false
	int j;


	for (int i = 0; i < img_height; i++)
	{
		if (toRight)
		{
			j = 0;
		}
		else
		{
			j = img_width - 1;
		}

		while (j < img_width && j >= 0)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			unsigned char red = rgb[offset_rgb + rr] + (char)redError[i * img_width + j];
			unsigned char green = rgb[offset_rgb + gg] + (char)greenError[i * img_width + j];
			unsigned char blue = rgb[offset_rgb + bb] + (char)blueError[i * img_width + j];

			//get colsed color
			unsigned char newRed = redTable[red >> 5];
			unsigned char newGreen = greenTable[green >> 5];
			unsigned char newBlue = blueTable[blue >> 6];

			//fill color
			img_data[offset_rgba + rr] = newRed;
			img_data[offset_rgba + gg] = newGreen;
			img_data[offset_rgba + bb] = newBlue;
			img_data[offset_rgba + aa] = WHITE;

			//evaluate error
			double newRedError = red - newRed;
			double newGreenError = green - newGreen;
			double newBlueError = blue - newBlue;

			bool right, down, left; // check position can move to or not
			if (j + 1 >= img_width)
			{
				right = false;
			}
			else
			{
				right = true;
			}
			if (i + 1 >= img_height)
			{
				down = false;
			}
			else
			{
				down = true;
			}
			if (j - 1 < 0)
			{
				left = false;
			}
			else
			{
				left = true;
			}

			// fill error
			if (down == false)
			{
				//at the buttom of the image

				if ((right == false && toRight == true) || (left == false && toRight == false))
				{
					// at the end of the image

				}
				else
				{
					if (toRight == true)
					{
						//next is going to right

						//other road is closed, get all the error
						redError[i * img_width + (j + 1)] += newRedError;
						greenError[i * img_width + (j + 1)] += newGreenError;
						blueError[i * img_width + (j + 1)] += newBlueError;
					}
					else
					{
						//next is going to left

						//other road is closed, get all the error
						redError[i * img_width + (j - 1)] += newRedError;
						greenError[i * img_width + (j - 1)] += newGreenError;
						blueError[i * img_width + (j - 1)] += newBlueError;
					}
				}
			}
			else
			{
				//can move to down

				if (right == false && toRight == true)
				{
					// at the right border

					//down get 10/16 error
					redError[(i + 1) * img_width + j] += newRedError * 10 / 16;
					greenError[(i + 1) * img_width + j] += newGreenError * 10 / 16;
					blueError[(i + 1) * img_width + j] += newBlueError * 10 / 16;

					//down-left get 6/16 error
					redError[(i + 1) * img_width + (j - 1)] += newRedError * 6 / 16;
					greenError[(i + 1) * img_width + (j - 1)] += newGreenError * 6 / 16;
					blueError[(i + 1) * img_width + (j - 1)] += newBlueError * 6 / 16;
				}
				else if (left == false && toRight == false)
				{
					// at the left border

					//down get 10/16 error
					redError[(i + 1) * img_width + j] += newRedError * 10 / 16;
					greenError[(i + 1) * img_width + j] += newGreenError * 10 / 16;
					blueError[(i + 1) * img_width + j] += newBlueError * 10 / 16;

					//down-right get 6/16 error
					redError[(i + 1) * img_width + (j + 1)] += newRedError * 6 / 16;
					greenError[(i + 1) * img_width + (j + 1)] += newGreenError * 6 / 16;
					blueError[(i + 1) * img_width + (j + 1)] += newBlueError * 6 / 16;
				}
				else
				{
					// at the middle

					if (toRight == true)
					{
						//next move to right

						//right get 7/16 error
						redError[i * img_width + (j + 1)] += newRedError * 7 / 16;
						greenError[i * img_width + (j + 1)] += newGreenError * 7 / 16;
						blueError[i * img_width + (j + 1)] += newBlueError * 7 / 16;

						//down-right get 1/16 error
						redError[(i + 1) * img_width + (j + 1)] += newRedError * 1 / 16;
						greenError[(i + 1) * img_width + (j + 1)] += newGreenError * 1 / 16;
						blueError[(i + 1) * img_width + (j + 1)] += newBlueError * 1 / 16;

						//down get 5/16 error
						redError[(i + 1) * img_width + j] += newRedError * 5 / 16;
						greenError[(i + 1) * img_width + j] += newGreenError * 5 / 16;
						blueError[(i + 1) * img_width + j] += newBlueError * 5 / 16;

						//down-left get 3/16 error
						redError[(i + 1) * img_width + (j - 1)] += newRedError * 3 / 16;
						greenError[(i + 1) * img_width + (j - 1)] += newGreenError * 3 / 16;
						blueError[(i + 1) * img_width + (j - 1)] += newBlueError * 3 / 16;
					}
					else
					{
						//next move to left

						//left get 7/16 error
						redError[i * img_width + (j - 1)] += newRedError * 7 / 16;
						greenError[i * img_width + (j - 1)] += newGreenError * 7 / 16;
						blueError[i * img_width + (j - 1)] += newBlueError * 7 / 16;

						//down-left get 1/16 error
						redError[(i + 1) * img_width + (j - 1)] += newRedError * 1 / 16;
						greenError[(i + 1) * img_width + (j - 1)] += newGreenError * 1 / 16;
						blueError[(i + 1) * img_width + (j - 1)] += newBlueError * 1 / 16;

						//down get 5/16 error
						redError[(i + 1) * img_width + j] += newRedError * 5 / 16;
						greenError[(i + 1) * img_width + j] += newGreenError * 5 / 16;
						blueError[(i + 1) * img_width + j] += newBlueError * 5 / 16;

						//down-right get 1/16 error
						redError[(i + 1) * img_width + (j + 1)] += newRedError * 3 / 16;
						greenError[(i + 1) * img_width + (j + 1)] += newGreenError * 3 / 16;
						blueError[(i + 1) * img_width + (j + 1)] += newBlueError * 3 / 16;
					}
				}
			}//*/
			
			//fill error
			/*if (j + 1 < img_width)
			{
				//right get 7/16 error
				redError[i * img_width + (j + 1)] += newRedError * 7 / 16;
				greenError[i * img_width + (j + 1)] += newGreenError * 7 / 16;
				blueError[i * img_width + (j + 1)] += newBlueError * 7 / 16;
			}
			if ((j + 1 < img_width) && (i + 1 < img_height))
			{
				//down-right get 1/16 error
				redError[(i + 1) * img_width + (j + 1)] += newRedError * 1 / 16;
				greenError[(i + 1) * img_width + (j + 1)] += newGreenError * 1 / 16;
				blueError[(i + 1) * img_width + (j + 1)] += newBlueError * 1 / 16;
			}
			if (i + 1 < img_height)
			{
				//down get 5/16 error
				redError[(i + 1) * img_width + j] += newRedError * 5 / 16;
				greenError[(i + 1) * img_width + j] += newGreenError * 5 / 16;
				blueError[(i + 1) * img_width + j] += newBlueError * 5 / 16;
			}
			if ((j - 1 >= 0) && (i + 1 < img_height))
			{
				//down-left get 3/16 error
				redError[(i + 1) * img_width + (j - 1)] += newRedError * 3 / 16;
				greenError[(i + 1) * img_width + (j - 1)] += newGreenError * 3 / 16;
				blueError[(i + 1) * img_width + (j - 1)] += newBlueError * 3 / 16;
			}//*/

			//update
			if (toRight)
			{
				j++;
			}
			else
			{
				j--;
			}
		}

		if (toRight)
		{
			toRight = false;
		}
		else
		{
			toRight = true;
		}
	}

	delete[] redError;
	delete[] greenError;
	delete[] blueError;
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering( double filter[][5] )
{
	unsigned char *rgb = this->To_RGB();

	//*add
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			//looking neighberhood
			double newRed = 0;
			double newGreen = 0;
			double newBlue = 0;
			for (int m = 0; m < 5; m++)
			{
				for (int n = 0; n < 5; n++)
				{
					if ((i - 2 + m) >= 0 && (i - 2 + m) < img_height && (j - 2 + n) >= 0 && (j - 2 + n) < img_width)
					{
						int offset_rgb = (i - 2 + m) * img_width * 3 + (j - 2 + n) * 3;
						newRed += rgb[offset_rgb + rr] * filter[m][n];
						newGreen += rgb[offset_rgb + gg] * filter[m][n];
						newBlue += rgb[offset_rgb + bb] * filter[m][n];
					}
					else
					{
						newRed += WHITE * filter[m][n];
						newGreen += WHITE * filter[m][n];
						newBlue += WHITE * filter[m][n];
					}
				}
			}

			//fill color
			int offset_rgba = i * img_width * 4 + j * 4;

			if (newRed < 0.0)
			{
				newRed = 0;
			}
			else if (newRed > 255.0)
			{
				newRed = 255.0;
			}
			if (newGreen < 0.0)
			{
				newGreen = 0;
			}
			else if (newGreen > 255.0)
			{
				newGreen = 255.0;
			}
			if (newBlue < 0.0)
			{
				newBlue = 0;
			}
			else if (newBlue > 255.0)
			{
				newBlue = 255.0;
			}

			img_data[offset_rgba + rr] = (unsigned char)newRed;
			img_data[offset_rgba + gg] = (unsigned char)newGreen;
			img_data[offset_rgba + bb] = (unsigned char)newBlue;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::filtering( double **filter, int n )
{
	unsigned char *rgb = this->To_RGB();

	//*add
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			//looking neighberhood
			double newRed = 0;
			double newGreen = 0;
			double newBlue = 0;
			for (int a = 0; a < n; a++)
			{
				for (int b = 0; b < n; b++)
				{
					if ((i - (n / 2) + a) >= 0 && (i - (n / 2) + a) < img_height && (j - (n / 2) + b) >= 0 && (j - (n / 2) + b) < img_width)
					{
						int offset_rgb = (i - (n / 2) + a) * img_width * 3 + (j - (n / 2) + b) * 3;
						newRed += rgb[offset_rgb + rr] * filter[a][b];
						newGreen += rgb[offset_rgb + gg] * filter[a][b];
						newBlue += rgb[offset_rgb + bb] * filter[a][b];
					}
					else
					{
						newRed += WHITE * filter[a][b];
						newGreen += WHITE * filter[a][b];
						newBlue += WHITE * filter[a][b];
					}
				}
			}

			//fill color
			int offset_rgba = i * img_width * 4 + j * 4;
			
			if (newRed < 0.0)
			{
				newRed = 0;
			}
			else if (newRed > 255.0)
			{
				newRed = 255.0;
			}
			if (newGreen < 0.0)
			{
				newGreen = 0;
			}
			else if (newGreen > 255.0)
			{
				newGreen = 255.0;
			}
			if (newBlue < 0.0)
			{
				newBlue = 0;
			}
			else if (newBlue > 255.0)
			{
				newBlue = 255.0;
			}

			img_data[offset_rgba + rr] = (unsigned char)newRed;
			img_data[offset_rgba + gg] = (unsigned char)newGreen;
			img_data[offset_rgba + bb] = (unsigned char)newBlue;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	//add*

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	//*add
	double filter[5][5] = { 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25,
							1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25,
							1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25,
							1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25,
							1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25};
	Application::filtering(filter);
	//add*
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	//*add
	double filter[5][5] = { 1.0 / 81, 2.0 / 81, 3.0 / 81, 2.0 / 81, 1.0 / 81,
							2.0 / 81, 4.0 / 81, 6.0 / 81, 4.0 / 81, 2.0 / 81,
							3.0 / 81, 6.0 / 81, 9.0 / 81, 6.0 / 81, 3.0 / 81,
							2.0 / 81, 4.0 / 81, 6.0 / 81, 4.0 / 81, 2.0 / 81,
							1.0 / 81, 2.0 / 81, 3.0 / 81, 2.0 / 81, 1.0 / 81 };
	Application::filtering(filter);
	//add*
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	//*add
	double filter[5][5] = { 1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256,
							4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
							6.0 / 256, 24.0 / 256, 36.0 / 256, 24.0 / 256, 6.0 / 256,
							4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
							1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256 };
	Application::filtering(filter);
	//add*
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N( unsigned int N )
{
	//*add
	unsigned long long int weight = pow(2, (N-1) * 2); // N must < 33

	int *filter_1d;

	double **filter_2d = new double*[N];
	for (int i = 0; i < N; i++)
	{
		filter_2d[i] = new double[N];
	}

	//calculate pascal triangle
	int *currTriangle = new int[2];
	currTriangle[0] = currTriangle[1] = 1;

	int *nextTriangle;
	for (int i = 2; i < N; i++)
	{
		nextTriangle = new int[i + 1];
		nextTriangle[0] = nextTriangle[i] = 1; //head and tail = 1

		for (int j = 0; j < i - 1; j++)
		{
			nextTriangle[j + 1] = currTriangle[j] + currTriangle[j + 1];
		}

		//stroage
		delete[] currTriangle;
		currTriangle = nextTriangle;
	}
	filter_1d = currTriangle;

	//evaluate 2d filter
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			filter_2d[i][j] = filter_1d[i] * (double)filter_1d[j] / weight;
		}
	}

	Application::filtering(filter_2d, N);

	delete[] filter_1d;
	delete[] filter_2d;
	//add*
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	//*add
	double lowPassFilter[5][5] = { 1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256,
								   4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
								   6.0 / 256, 24.0 / 256, 36.0 / 256, 24.0 / 256, 6.0 / 256,
								   4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
								   1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256 };
	double original[5][5] = { 0, 0, 0, 0, 0,
							  0, 0, 0, 0, 0,
							  0, 0, 1, 0, 0,
							  0, 0, 0, 0, 0,
							  0, 0, 0, 0, 0 };

	double filter[5][5];
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			filter[i][j] = original[i][j] - lowPassFilter[i][j];
		}
	}

	Application::filtering(filter);
	//add*
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	//unsigned char *rgb = this->To_RGB();

	//*add
	double lowPassFilter[5][5] = { 1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256,
								   4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
								   6.0 / 256, 24.0 / 256, 36.0 / 256, 24.0 / 256, 6.0 / 256,
								   4.0 / 256, 16.0 / 256, 24.0 / 256, 16.0 / 256, 4.0 / 256,
								   1.0 / 256, 4.0 / 256, 6.0 / 256, 4.0 / 256, 1.0 / 256 };
	double original[5][5] = { 0, 0, 0, 0, 0,
							  0, 0, 0, 0, 0,
							  0, 0, 1, 0, 0,
							  0, 0, 0, 0, 0,
							  0, 0, 0, 0, 0 };

	double filter[5][5];
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			filter[i][j] = original[i][j] - lowPassFilter[i][j];
			filter[i][j] += original[i][j];
		}
	}

	Application::filtering(filter);
	//add*

	//delete[] rgb;
	//mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	//renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	//*add
	unsigned char *rgb = this->To_RGB();

	int newHeight = (img_height + 1) / 2.0;
	int newWidth  = (img_width  + 1) / 2.0;
	img_data = new unsigned char[newWidth * 4 * newHeight * 4];

	double filter[3][3] = { 1.0 / 16, 1.0 / 8, 1.0 / 16,
							1.0 / 8, 1.0 / 4, 1.0 / 8,
							1.0 / 16, 1.0 / 8, 1.0 / 16 };

	for (int i = 0; i < img_height; i+=2)
	{
		for (int j = 0; j < img_width; j+=2)
		{
			//looking neighberhood
			double newRed = 0;
			double newGreen = 0;
			double newBlue = 0;
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					if ((i - 1 + m) >= 0 && (i - 1 + m) < img_height && (j - 1 + n) >= 0 && (j - 1 + n) < img_width)
					{
						int offset_rgb = (i - 1 + m) * img_width * 3 + (j - 1 + n) * 3;
						newRed += rgb[offset_rgb + rr] * filter[m][n];
						newGreen += rgb[offset_rgb + gg] * filter[m][n];
						newBlue += rgb[offset_rgb + bb] * filter[m][n];
					}
					else
					{
						newRed += BLACK * filter[m][n];
						newGreen += BLACK * filter[m][n];
						newBlue += BLACK * filter[m][n];
					}
				}
			}

			//fill color
			int offset_rgba = (i / 2) * newWidth * 4 + (j / 2) * 4;
			img_data[offset_rgba + rr] = (unsigned char)newRed;
			img_data[offset_rgba + gg] = (unsigned char)newGreen;
			img_data[offset_rgba + bb] = (unsigned char)newBlue;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	mImageDst = QImage(img_data, newWidth, newHeight, QImage::Format_ARGB32);
	this->img_height = newHeight;
	this->img_width = newWidth;

	delete[] rgb;
	//add*

	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	//*add
	unsigned char *rgb = this->To_RGB();

	const int brushNum = 3;
	unsigned char radius[brushNum] = { 7, 3, 1 }; //biggest to smallest
	unsigned char *canvas = new unsigned char[img_width * 4 * img_height * 4];
	unsigned char *refrenceImg = new unsigned char[img_width * 4 * img_height * 4];
	unsigned char *copy_imgData = new unsigned char[img_width * 4 * img_height * 4];
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			canvas[offset_rgba + rr] = 0;
			canvas[offset_rgba + gg] = 0;
			canvas[offset_rgba + bb] = 0;
			canvas[offset_rgba + aa] = WHITE;

			copy_imgData[offset_rgba + rr] = img_data[offset_rgba + rr];
			copy_imgData[offset_rgba + gg] = img_data[offset_rgba + gg];
			copy_imgData[offset_rgba + bb] = img_data[offset_rgba + bb];
			copy_imgData[offset_rgba + aa] = img_data[offset_rgba + aa];
		}
	}

	for (int i = 0; i < brushNum; i++)
	{
		this->Gaussion_Filter(copy_imgData, refrenceImg, img_height, img_width, 2 * radius[i] + 1);
		this->NPR_Paint_Layer(canvas, refrenceImg, radius[i]);
		
	}
	
	delete[] refrenceImg;
	delete[] copy_imgData;
	delete[] canvas;
	delete[] rgb;
	//add*
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{
	//*add
	const char fg = 1;
	const char T = 25;

	std::vector<Stroke> strokeList;
	srand(time(NULL));

	double *diff = new double[img_height * img_width];
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			diff[i*img_width + j] = sqrt(pow((tCanvas[offset_rgba + rr] - tReferenceImage[offset_rgba + rr]), 2)
				+ pow((tCanvas[offset_rgba + gg] - tReferenceImage[offset_rgba + gg]), 2)
				+ pow((tCanvas[offset_rgba + bb] - tReferenceImage[offset_rgba + bb]), 2));
		}
	}

	int grid = fg * tBrushSize;

	for (int y = 0; y < img_height; y += grid)
	{
		for (int x = 0; x < img_width; x += grid)
		{
			//sum errors and get max error
			int maxX, maxY;
			double maxError = -1;
			double sumError = 0.0;

			for (int m = 0; m < grid; m++)
			{
				for (int n = 0; n < grid; n++)
				{
					int i = y - grid / 2 + m;
					int j = x - grid / 2 + n;
					if (i >= 0 && i < img_height && j >= 0 && j < img_width)
					{
						double err = diff[i*img_width + j];
						sumError += err;
						if (err > maxError)
						{
							maxX = j;
							maxY = i;
						}
					}
				}
			}

			sumError /= grid * grid;

			if (sumError > T)
			{
				int offset_rgba = maxY * img_width * 4 + maxX * 4;
				Stroke s(tBrushSize, maxX, maxY, tReferenceImage[offset_rgba + rr], tReferenceImage[offset_rgba + gg],
					tReferenceImage[offset_rgba + bb], tReferenceImage[offset_rgba + aa]);
				strokeList.push_back(s);
			}
		}
	}

	//random change psition
	for (int i = 0; i < strokeList.size(); i++)
	{
		int randomIndex = rand() % strokeList.size();
		Stroke s = strokeList[i];
		strokeList[i] = strokeList[randomIndex];
		strokeList[randomIndex] = s;
	}

	//print stroke
	for (int i = 0; i < strokeList.size(); i++)
	{
		this->Paint_Stroke(strokeList[i]);
	}

	//refresh canvas
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgba = i * img_width * 4 + j * 4;
			tCanvas[offset_rgba + rr] = img_data[offset_rgba + rr];
			tCanvas[offset_rgba + gg] = img_data[offset_rgba + gg];
			tCanvas[offset_rgba + bb] = img_data[offset_rgba + bb];
			tCanvas[offset_rgba + aa] = img_data[offset_rgba + aa];
		}
	}

	delete[] diff;
	//add*
}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}

//gaussion filter
void Application::Gaussion_Filter(unsigned char *srcImage, unsigned char *dstImage, int height, int width, int n)
{
	unsigned long long int weight = pow(2, (n - 1) * 2); // N must < 33

	int *filter_1d;
	double **filter_2d = new double*[n];
	for (int i = 0; i < n; i++)
	{
		filter_2d[i] = new double[n];
	}

	//calculate pascal triangle
	int *currTriangle = new int[2];
	currTriangle[0] = currTriangle[1] = 1;

	int *nextTriangle;
	for (int i = 2; i < n; i++)
	{
		nextTriangle = new int[i + 1];
		nextTriangle[0] = nextTriangle[i] = 1; //head and tail = 1

		for (int j = 0; j < i - 1; j++)
		{
			nextTriangle[j + 1] = currTriangle[j] + currTriangle[j + 1];
		}

		//stroage
		delete[] currTriangle;
		currTriangle = nextTriangle;
	}
	filter_1d = currTriangle;

	//evaluate 2d filter
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			filter_2d[i][j] = filter_1d[i] * (double)filter_1d[j] / weight;
		}
	}

	//-------------------------------
	//main filter
	//-------------------------------
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			//looking neighberhood
			double newRed = 0;
			double newGreen = 0;
			double newBlue = 0;
			for (int a = 0; a < n; a++)
			{
				for (int b = 0; b < n; b++)
				{
					if ((i - (n / 2) + a) >= 0 && (i - (n / 2) + a) < height && (j - (n / 2) + b) >= 0 && (j - (n / 2) + b) < img_width)
					{
						int offset_rgba = (i - (n / 2) + a) * width * 4 + (j - (n / 2) + b) * 4;
						newRed += srcImage[offset_rgba + rr] * filter_2d[a][b];
						newGreen += srcImage[offset_rgba + gg] * filter_2d[a][b];
						newBlue += srcImage[offset_rgba + bb] * filter_2d[a][b];
					}
					else
					{
						newRed += WHITE * filter_2d[a][b];
						newGreen += WHITE * filter_2d[a][b];
						newBlue += WHITE * filter_2d[a][b];
					}
				}
			}

			//fill color
			int offset_rgba = i * width * 4 + j * 4;

			if (newRed < 0.0)
			{
				newRed = 0;
			}
			else if (newRed > 255.0)
			{
				newRed = 255.0;
			}
			if (newGreen < 0.0)
			{
				newGreen = 0;
			}
			else if (newGreen > 255.0)
			{
				newGreen = 255.0;
			}
			if (newBlue < 0.0)
			{
				newBlue = 0;
			}
			else if (newBlue > 255.0)
			{
				newBlue = 255.0;
			}

			dstImage[offset_rgba + rr] = (unsigned char)newRed;
			dstImage[offset_rgba + gg] = (unsigned char)newGreen;
			dstImage[offset_rgba + bb] = (unsigned char)newBlue;
			dstImage[offset_rgba + aa] = WHITE;
		}
	}
	//-----------------------------------------------------

	delete[] filter_1d;
	delete[] filter_2d;
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
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}



