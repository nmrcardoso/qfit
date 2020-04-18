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

#ifndef QFIT_H
#define QFIT_H


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



#include "structs.h"






using namespace std;


class Program;


class MyGlWindow : public Fl_Widget {
	friend class Program;
	DataLine data;
	ConstFit& fit;
	int width;
	int height;
	bool plotVr;
	double Tmin, Tmax;
	double border;
	double w_x0, w_x1, w_y0, w_y1, w_w, w_h;
	int fontboxx, fontboxy;
	double dx, dy;
	int ldxy;
	FitRange rangeVr;
	double Vx_min, Vx_max;
	GSLfitRes res;
	bool fitVr;
	double xmin, xmax;	
	double ymin, ymax;

	public:
	MyGlWindow(int X, int Y, int W, int H, ConstFit& _fit) ;
  private:
	void draw() ;

	void DrawTicksXY(int fAxisLength, int div1, int div2, double fWmin, double fWmax, int yaxis);

	void createdraw() ;	

	int convertLeftButton(int pos, double x0, double x1, double tmin, double tmax) const;
	int convertRightButton(int pos, double x0, double x1, double tmin, double tmax) const;
  	void Mouse_Handler_GLWiN(int pos, int whichButtonCliked);
	int handle(int event);

};


#endif
