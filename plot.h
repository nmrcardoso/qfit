/***************************************************************************
 *   Copyright (C) 2020 by Nuno Cardoso                                    *
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

#ifndef PLOT_QFIT_H
#define PLOT_QFIT_H

#include <FL/fl_draw.H>

#include "structs.h"
#include "ticks.h"

using namespace std;

class DrawErrorPlot : public DrawTicks{	
	double min;
	double max;
	int x, y, w, h;
	int fAxisLength_x, fAxisLength_y;
	Fl_Color IntPointsColor, ExtPointsColor;
	Fl_Color bcolor, bground, fitRcolor;
	int Errorlinewidth;
	int dxy;
	public:
	int w_x0, w_x1, w_y0, w_y1;
	double xmin, xmax, ymin, ymax;
	DrawErrorPlot(int x, int y, int w, int h);
	void setErrorLineWidth(int size);
	void setPointSize(int size);
	void setBackColor(Fl_Color color);
	void setBackgroundColor(Fl_Color color);
	void setFRangeColor(Fl_Color color);
	void setPointsColor(Fl_Color IntPoints, Fl_Color ExtPoints);
	void draw(DataLine data, double t1, double t2);
	
};


#endif

