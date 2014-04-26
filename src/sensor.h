
#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "color.h"

#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfChromaticities.h>
#include <OpenEXR/ImfStandardAttributes.h>
using namespace Imf;
using namespace Imath;

inline int min(int a, int b) {
	return a < b ? a : b;
}

inline int max(int a, int b) {
	return a < b ? b : a;
}

struct Sensor {
	
	struct Pixel {
		rgba sum;
		float weight;
		unsigned int samples;

		Pixel() : sum(0,0,0,1), weight(0), samples(0) {}
	};


	int xres, yres;
	float half;
	float inv_size;
	std::vector<Pixel> array;
	int samples;
	float scale;
	
	Sensor(int xres, int yres, float filtersize) : xres(xres), yres(yres), half(filtersize*0.5f), inv_size(1.f/filtersize), array(xres*yres, Pixel()), samples(0), scale(sqrt(xres*yres)) {}

	void detect(float x, float y, float wavelength, float radiance) {
		samples++;
		rgba c = Spectrum::torgb(wavelength, radiance);
		
		x = x*scale + xres/2;  
		y = y*scale + xres/2;

		int left = max(0,ceil(x-half)), right = min(xres-1,floor(x+half));
		int bottom = max(0,ceil(y-half)), top = min(yres-1,floor(y+half));
		
		for(int iy = bottom; iy <= top; iy++) {
			for(int ix = left; ix <= right; ix++) {
				float weight = mitchell((x-ix)*inv_size)*mitchell((y-iy)*inv_size)*inv_size*inv_size;
				pixel(ix, iy).sum += c*weight;
				pixel(ix, iy).weight += weight;
				pixel(ix, iy).samples++;
			}
		}
	}

	float mitchell(float d) {
		static const float B = 1.f/3;
		static const float C = 1.f/3;
		float x = 4*abs(d);
		if(x >= 2) return 0;
		else if(x >= 1) return ((-B-6*C)*x*x*x + (6*B + 30*C)*x*x + (-12*B-48*C)*x + (8*B+24*C)) * (1.f/6); 
		else return ((12-9*B-6*C)*x*x*x + (-18+12*B+6*C)*x*x + (6-2*B)) * (1.f/6);
	}

	void outputEXR(std::string file_name) {
		std::vector<Rgba> exr_pixels(xres*yres);
		float s = ((float)xres*yres)/samples;
		for(int i = 0; i < xres*yres; i++) {
			exr_pixels[i].r = array[i].sum.r * s;
			exr_pixels[i].g = array[i].sum.g * s;
			exr_pixels[i].b = array[i].sum.b * s;
			exr_pixels[i].a = 1;
		}
		Header header(xres, yres);
		Chromaticities chroma(V2f(1,0), V2f(0,1), V2f(0,0), V2f(1./3, 1./3));
		addChromaticities(header, chroma);
		RgbaOutputFile file(file_name.c_str(), header, WRITE_RGBA);
		file.setFrameBuffer(&exr_pixels[0], 1, xres);
		file.writePixels(yres);
	}

	Pixel& pixel(int x, int y) {
		return array[y*xres + x];
	}

/*
	bool outputPFM(std::string file_name) {
		rgba c;

		std::ofstream out(file_name.c_str(), std::ios::out | std::ios::binary);

		if( !out.is_open() ) return false;

		out << "PF\n";

		out << xres << " " << yres << "\n";    
		out << -1.0 << "\n";

		for ( int i = 0; i < yres; i++ ) {

			for ( int j = 0; j < xres; j++ ) {

				c = pixel(j, i).sum;

				out.write( (char*)&c.r, sizeof( float ) );
				out.write( (char*)&c.g, sizeof( float ) );
				out.write( (char*)&c.b, sizeof( float ) );
			}
		}
		
		out.close();

		return true;
	}
*/
};

#endif

