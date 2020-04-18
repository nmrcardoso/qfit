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


//#define USE_FL_WIDGET_DRAW

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

#include "plot.h"

#include "log.h"

#ifndef uint
#define uint unsigned int
#endif


using namespace std;






DrawErrorPlot::DrawErrorPlot(int x, int y, int w, int h): x(x), y(y), w(w), h(h){
	bground = FL_LIGHT2;
	bcolor = FL_LIGHT3;
	fitRcolor = FL_LIGHT1;
	IntPointsColor = FL_RED;
	ExtPointsColor = FL_BLUE;
	Errorlinewidth = 2;
	dxy = 4;
}


void DrawErrorPlot::setErrorLineWidth(int size){
	Errorlinewidth = size;
}
void DrawErrorPlot::setPointSize(int size){
	dxy = size;
}

void DrawErrorPlot::setPointsColor(Fl_Color IntPoints, Fl_Color ExtPoints){
	IntPointsColor = IntPoints;
	ExtPointsColor = ExtPoints;
}
void DrawErrorPlot::setBackgroundColor(Fl_Color color){
	bground = color;
}
void DrawErrorPlot::setBackColor(Fl_Color color){
	bcolor = color;
}

void DrawErrorPlot::setFRangeColor(Fl_Color color){
	fitRcolor = color;
}




template<typename type>
void fontscreenbox(type val, Fl_Font fontname, int fontsize, int &wi, int &hi){
	wi=0, hi=0;
	stringstream s; 
	s << val;
	fl_font(fontname, fontsize);
	fl_measure(s.str().c_str(), wi, hi);
}

void DrawErrorPlot::draw(DataLine data, double t1, double t2){


	xmin = data[0].x;
	xmax = data[data.size()-1].x;
	ymin = data.begin()->y, ymax = data.begin()->y;
	for( DataLine::iterator it = data.begin(); it != data.end(); ++it ) {
		ymin = std::min( it->y - it->error, ymin );
		ymax = std::max( it->y + it->error, ymax );
	}
	//xmin -= .25; xmax += .25;
	int newbinsX = data.size();
	
	OptimizeTicksLimits(data.size(), newbinsX, xmin, xmax, false);
	
	//qlog << xmin << ":::.:::" << xmax << ":::.:::" << data.size() << endl;
	//for(int i = 0; i < data.size(); i++) qlog << i << "::" << data[i].x << "::..." << data[i].y << endl;
	
	int newbinsY = 5;
	//needs this before draw area and call plot for the best number of ticks
	OptimizeTicksLimits(8, newbinsY, ymin, ymax, false);
	setTicks(newbinsX, 5, newbinsY, 5);

	// Calculate the best space around plot
	int border = 10;
	int w0=0, h0=0, w1=0, h1=0, x0 = 0, y0 = 0;
	if(xaxis){
		fontscreenbox(xmin, getTickFont(), getTickFontSize(), w0, h0);
		fontscreenbox(xmax, getTickFont(), getTickFontSize(), w1, h1);
		y0 = std::max( h0, h1);
		w0 = std::max( w0, w1)/2;
		border = std::max( border, w0);
		border = std::max( border, h1/2);
	}
	if(yaxis){
		fontscreenbox(ymin, getTickFont(), getTickFontSize(), w0, h0);
		fontscreenbox(ymax, getTickFont(), getTickFontSize(), w1, h1);
		w0 = std::max( w0, w1);
		double ytemp = (ymax-ymin)/newbinsY;
		fontscreenbox(ytemp, getTickFont(), getTickFontSize(), w1, h1);
		x0 = std::max( w0, w1);
		border = std::max( border, h1/2);
	}
	
	w_x0 = x+border+x0;
	w_x1 = x+w-border;
	w_y0 = y+border;
	w_y1 = y+h-border-y0;
	int w_w = w_x1-w_x0;
	int w_h = w_y1-w_y0;
	InitTicks(w_x0, w_x1, w_y0, w_y1);

	double xrange = (xmax-xmin), yrange = (ymax-ymin);
	fAxisLength_x = w_x1-w_x0;
	fAxisLength_y = w_y1-w_y0;
	double factorx = fAxisLength_x/xrange;
	double factory = fAxisLength_y/yrange;

	//background color around the plot
	fl_color(bcolor);
	fl_rectf(x, y, w, h);

	//background color inside the plot
	fl_line_style(0, 1);
	fl_color(bground);
	fl_rectf(w_x0, w_y0, w_w, w_h);

	//highlight zone of the fit range
	if( t2 > t1  && t1 >= data[0].x && t2 <= data[data.size()-1].x){
		double recxmin =  (t1-xmin - 0.5)*factorx;
		double recxmax = (t2-xmin + 0.5)*factorx - recxmin ;
		if(t1 == data[0].x){
			recxmin = 0.;
			recxmax = (t2-xmin + 0.5)*factorx;
		}
		if(t2 == data[data.size()-1].x){
			recxmin = (t1-xmin - 0.5)*factorx;
			recxmax = w_w-recxmin;
		}
		if(t1 == data[0].x && t2== data[data.size()-1].x){
			recxmin = 0.;
			recxmax = w_w;
		}
		fl_color(fitRcolor);
		fl_rectf(w_x0 + recxmin, w_y0, recxmax, w_h );
	}

	//Draw axis plot
	drawAxis(xmin, xmax, ymin, ymax);
	//DRAW POINTS
	for( uint i = 0; i < data.size(); ++i ) {
		//qlog << data[i].x << " :: " << data[i].y << endl;
		double xx = data[i].x;
		int x = w_x0+round((xx-xmin)*factorx);
		double yy = data[i].y;
		int y = w_y1-round((yy-ymin)*factory);
		double error = data[i].error;
		if (t2 > t1 && data[i].x >= t1 && data[i].x <= t2 )
			fl_color(IntPointsColor);
		else
			fl_color(ExtPointsColor);
		//draw square point
		fl_line_style(0, Errorlinewidth);
		fl_rectf(x - dxy, y-dxy, 2*dxy, 2*dxy );
		//draw error bar
		fl_line(x, y + error * factory, x, y - error * factory);
		////draw error bar top bar			
		fl_line(x - dxy, y + error * factory, x + dxy, y + error * factory);
		////draw error bar bottom bar
		fl_line(x - dxy, y - error * factory, x + dxy, y - error * factory);
	}
}



