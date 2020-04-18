/***************************************************************************
 *   Copyright (C) 2019 by Nuno Cardoso                                    *
 *   nmrcardoso@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef SAVE_IMAGE_H
#define SAVE_IMAGE_H


#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <map>
#include <vector>
#include <cmath>
#include <float.h>




#define LE16(n) (uchar)((n) & 0xFF), (uchar)(((n)>>8) & 0xFF)
#define LE32(n) (uchar)((n) & 0xFF), (uchar)(((n)>>8) & 0xFF), (uchar)(((n)>>16) & 0xFF), (uchar)(((n)>>24) & 0xFF)

#define INCHES_PER_METER 39.3701

class Image {
public:
	enum Format { PNG, BMP, TGA, PPM };
	static const char *FILE_CHOOSER_FILTER;
private:
	FILE *_file;
	uchar *_buffer;
	size_t _width, _height;
	int _error;
public:
	Image(const char *f, size_t w, size_t h, Fl_Widget *wgt = NULL);
	~Image();
	size_t write(Format f);
	inline int error(void) { return ferror(_file) | _error; }
	inline void close(void) { fclose(_file); _file = NULL; }
private:
	Image(const Image &image); // Unimplemented copy constructor
	Image &operator=(const Image &image); // Unimplemented assignment operator
	size_t write_png(void);
	size_t write_bmp(void);
	size_t write_tga(void);
	size_t write_ppm(void);
};




using namespace std;



#endif
