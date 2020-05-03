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

#include "ticks.h"


#include "log.h"

using namespace std;


void DrawTicks::setDefaultValues(){
	fontcolor = FL_BLACK;
	axislinecolor = FL_BLACK;
	fontaxis = FL_HELVETICA_BOLD;
	fontsize = 12;
	NxDivLvl1 = 10;
	NxDivLvl2 = 5;
	NyDivLvl1 = 10;
	NyDivLvl2 = 5;
	ticks_color = FL_BLACK;
	ticks_length = 20;
	ticks_width = 1;
	axislinewidth = 1;
	gridlines = false;
	gridlinewidth = 1;
	xaxis = false;
	yaxis = false;
}


DrawTicks::DrawTicks() {
	w_x0 = 0;
	w_x1 = 0;
	w_y0 = 0;
	w_y1 = 0;
	ScreenAxisLength_x = w_x1-w_x0;
	ScreenAxisLength_y = w_y1-w_y0;
	setDefaultValues();
}

void DrawTicks::InitTicks(int w_x0_, int w_x1_, int w_y0_, int w_y1_){
	w_x0 = w_x0_;
	w_x1 = w_x1_;
	w_y0 = w_y0_;
	w_y1 = w_y1_;
	ScreenAxisLength_x = w_x1-w_x0;
	ScreenAxisLength_y = w_y1-w_y0;
}


void DrawTicks::setTickLength(int size){
	ticks_length = size;
}
void DrawTicks::setTickWidth(int width){
	ticks_width = width;
}
void DrawTicks::setTickColor(Fl_Color color){
	ticks_color = color;
}
void DrawTicks::setAxisLineWidth(int size){
	axislinewidth = size;
}
void DrawTicks::setAxisLineColor(Fl_Color axislinecolor_){
	axislinecolor = axislinecolor_;
}
void DrawTicks::setAxisFontSize(int size){
	fontsize = size;
}
void DrawTicks::setAxisFontColor(Fl_Color fontcolor_){
	fontcolor = fontcolor_;
}
void DrawTicks::setGridLines(int size){
	gridlines = true;
	gridlinewidth = size;
}
void DrawTicks::unsetGridLines(){
	gridlines = false;
}
	
void DrawTicks::setAxisX() { xaxis = true; }
void DrawTicks::setAxisY() { yaxis = true; }
void DrawTicks::setTicks(int NxDivLvl1_, int NxDivLvl2_, int NyDivLvl1_, int NyDivLvl2_) { 
		NxDivLvl1 = NxDivLvl1_, NxDivLvl2 = NxDivLvl2_; 
		NyDivLvl1 = NyDivLvl1_, NyDivLvl2 = NyDivLvl2_;
}

void DrawTicks::drawAxis(double xmin, double xmax, double ymin, double ymax){
	fl_font(fontaxis, fontsize);
	if( xaxis) {  // Draw x axis
		//draw top and bottom lines
		fl_color(axislinecolor);
		fl_line_style(0, axislinewidth);
		fl_line(w_x0, w_y0, w_x1, w_y0);
		fl_line(w_x0, w_y1, w_x1, w_y1);


		//draw tick lines
		if(NxDivLvl1){
			Ticks xticks(ScreenAxisLength_x, NxDivLvl1, NxDivLvl2, xmin, xmax);
			xticks.CalculateTicks();
			for (int ii=0; ii<xticks.NumTicksLvl1; ii++){
				int posx = round(double(w_x0)+xticks.TicksPosLvl1[ii]);
				if(gridlines){
					if(posx > w_x0 && posx < w_x1){
						fl_color(48);
						fl_line_style(1,gridlinewidth);
						fl_line(posx, w_y0, posx, w_y1);
					}
				}
				fl_color(ticks_color);
				fl_line_style(0,ticks_width);
				fl_line(posx, w_y1-ticks_length, posx, w_y1);
				fl_line(posx, w_y0, posx, w_y0+ticks_length);
				//draw values in x ticks
				int wi=0, hi=0;
				stringstream s; 
				s << xticks.TicksValuesLvl1[ii];
				fl_measure(s.str().c_str(), wi, hi);
				fl_color(fontcolor);
			 	fl_draw(s.str().c_str(), posx-wi/2, w_y1+ticks_length/4, wi, hi, FL_ALIGN_CENTER);
			}
			// draw minor ticks
			fl_color(ticks_color);
			fl_line_style(0,ticks_width);
			for (int ii=0; ii<xticks.NumTicksLvl2; ii++){
				int iniv = round(double(w_x0)+xticks.TicksPosLvl2[ii]);
				fl_line(iniv, w_y1-ticks_length/2, iniv, w_y1);
				fl_line(iniv, w_y0, iniv, w_y0+ticks_length/2);
			}
		}
	}
	if( yaxis ){ // Draw y axis
		//set left and right lines	
		fl_color(axislinecolor);
		fl_line_style(0, axislinewidth);
		fl_line(w_x0, w_y0, w_x0, w_y1);
		fl_line(w_x1, w_y0, w_x1, w_y1);
		//draw tick lines
		if(NyDivLvl1){
			Ticks yticks(ScreenAxisLength_y, NyDivLvl1, NyDivLvl2, ymin, ymax);
			yticks.CalculateTicks();
			fl_line_style(0,5);
			for (int ii=0; ii<yticks.NumTicksLvl1; ii++){
				int posy = round(double(w_y1)-yticks.TicksPosLvl1[ii]);
				if(gridlines){
					if(posy < w_y1 && posy > w_y0) {
						fl_color(48);
						fl_line_style(1,gridlinewidth);
						fl_line(w_x0, posy, w_x1, posy);
					}
				}
				//if(posy < w_y0) qlog << posy << ":::" << w_y0 << endl;
				fl_color(ticks_color);
				fl_line_style(0,ticks_width);
				//right tick lines
				fl_line(w_x0, posy, w_x0+ticks_length, posy);
				//left tick lines
				fl_line(w_x1, posy, w_x1-ticks_length, posy);
				//draw values in x ticks
				int wi=0, hi=0;
				stringstream s; 
				s << yticks.TicksValuesLvl1[ii];
				//const char *s = ToString(yticks.TicksValuesLvl1[ii]).c_str();
				fl_measure(s.str().c_str(), wi, hi);
				fl_color(fontcolor);
			 	fl_draw(s.str().c_str(), w_x0-wi-ticks_length/4, posy-hi/2, wi, hi, FL_ALIGN_RIGHT);
			}
			// draw minor ticks
			fl_color(ticks_color);
			fl_line_style(0,ticks_width);
			for (int ii=0; ii<yticks.NumTicksLvl2; ii++){
				int iniy = round(double(w_y1)-yticks.TicksPosLvl2[ii]);
				fl_line(w_x0, iniy, w_x0+ticks_length/2, iniy);
				fl_line(w_x1, iniy, w_x1-ticks_length/2, iniy);
			}
		}
	}
}




double Log10(double val){ return std::log(val)/std::log(10.); }


Ticks::Ticks(){
	TicksPosLvl1 = 0;
	TicksValuesLvl1 = 0;
	TicksPosLvl2 = 0;
	NumTicksLvl1 = 0;
	NumTicksLvl2 = 0;
}
Ticks::Ticks(int ScreenAxisLength, int NDivLvl1, int NDivLvl2, double funcMin, double funcMax): ScreenAxisLength(ScreenAxisLength), NDivLvl1(NDivLvl1), NDivLvl2(NDivLvl2), funcMin(funcMin), funcMax(funcMax){
	TicksPosLvl1 = 0;
	TicksValuesLvl1 = 0;
	TicksPosLvl2 = 0;
	NumTicksLvl1 = 0;
	NumTicksLvl2 = 0;
}
Ticks::~Ticks(){
	if(TicksPosLvl1) delete[] TicksPosLvl1;
	if(TicksValuesLvl1) delete[] TicksValuesLvl1;
	if(TicksPosLvl2) delete[] TicksPosLvl2;
}

void Ticks::CalculateTicks(int ScreenAxisLength_, int NDivLvl1_, int NDivLvl2_, double funcMin_, double funcMax_){
	ScreenAxisLength = ScreenAxisLength_;
	NDivLvl1 = NDivLvl1_;
	NDivLvl2 = NDivLvl2_;
	funcMin = funcMin_;
	funcMax = funcMax_;
	CalculateTicks();
}

void Ticks::CalculateTicks(){
	if(TicksPosLvl1){ delete[] TicksPosLvl1; TicksPosLvl1 = 0;}
	if(TicksValuesLvl1){ delete[] TicksValuesLvl1; TicksValuesLvl1 = 0;}
	if(TicksPosLvl2){ delete[] TicksPosLvl2; TicksPosLvl2 = 0;}

	int nDivOpt;
	double step1=0, step2=0, fmin2=0, fmax2=0;
	double fmin = funcMin;
	double fmax = funcMax;
	int fNDiv = NDivLvl1;
	int fNDivLvl1 = fNDiv%100;
	int fNDivLvl2 = NDivLvl2;

	// Level 1 tick marks.
	OptimizeTicks(fmin,  fmax, fNDivLvl1, funcMin, funcMax, nDivOpt, step1);
	fNDivLvl1   = nDivOpt;
	NumTicksLvl1 = fNDivLvl1+1;
	TicksPosLvl1  = new double[NumTicksLvl1];
	TicksValuesLvl1  = new double[NumTicksLvl1];
	//qlog << NumTicksLvl1 << endl;
	double r = ScreenAxisLength/(fmax-fmin);
	int step11 = round(step1);
	for(int i=0; i < NumTicksLvl1; ++i){
		TicksValuesLvl1[i] = funcMin + i*step1;
		TicksPosLvl1[i] = r*(TicksValuesLvl1[i]-fmin);
	}	
	// Level 2 tick marks.
	if (fNDivLvl2) {
		double t2;
		OptimizeTicks(funcMin, funcMin+step1, fNDivLvl2, fmin2, fmax2, nDivOpt, step2);
		fNDivLvl2       = nDivOpt;
		step2        = abs((TicksPosLvl1[1]-TicksPosLvl1[0])/fNDivLvl2);
		int nTickl = (int)(TicksPosLvl1[0]/step2);
		int nTickr = (int)((ScreenAxisLength-TicksPosLvl1[NumTicksLvl1-1])/step2);
		NumTicksLvl2     = fNDivLvl1*(fNDivLvl2-1)+nTickl+nTickr;
		TicksPosLvl2      = new double[NumTicksLvl2];
		int k = 0;
		for (int i=0; i<NumTicksLvl1-1; i++) {
			t2 = TicksPosLvl1[i]+step2;
			for (int j=0; j<fNDivLvl2-1; j++) {
				TicksPosLvl2[k] = round(int(t2));
				k++;
				t2 = t2+step2;
			}
		}
		if (nTickl) {
			t2 = TicksPosLvl1[0]-step2;
			for (int i=0; i<nTickl; i++) {
			TicksPosLvl2[k] = round(int(t2));
			k++;
			t2 = t2-step2;
			}
		}
		if (nTickr) {
			t2 = TicksPosLvl1[NumTicksLvl1-1]+step2;
			for (int i=0; i<nTickr; i++) {
				TicksPosLvl2[k] = round(int(t2));
				k++;
				t2 = t2+step2;
			}
		}
		/*qlog << "________________________________" << endl;
		for (int i=0; i<NumTicksLvl2-1; i++) {
			qlog << TicksPosLvl2[i] << ":::" << TicksPosLvl2[i+1]-TicksPosLvl2[i] << endl;
		}*/
		
	}
}


void CallL20(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, int &nold, double &BinWidth, int &ntemp, int optionTime);
void CallL90(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, int &nold, double &BinWidth, int &ntemp, int optionTime);
void CallLOK(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, double &BinWidth, int optionTime);


void OptimizeTicks(double A1, double A2, int nold, double &BinLow, double &BinHigh, int &nbins, double &BinWidth){
   nbins = nold;
   int ntemp = 0;
   int optionTime = 0;
   double al = min(A1,A2);
   double ah = max(A1,A2);
   if (al == ah) ah = al + 1;
   // if nold  ==  -1 , program uses binwidth input from calling routine
   if (nold == -1 && BinWidth > 0 ){
   		CallL90(al, ah, BinLow, BinHigh, nbins, nold, BinWidth, ntemp, optionTime);
		return;
   }
   ntemp = max(nold, 2);
   if (ntemp < 1) ntemp = 1;
   CallL20(al, ah, BinLow, BinHigh, nbins, nold, BinWidth, ntemp, optionTime);
}



void CallL20(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, int &nold, double &BinWidth, int &ntemp, int optionTime){
   double awidth = (ah-al)/double(ntemp);
   double timemulti = 1;
   if (awidth >= FLT_MAX) {   //in float.h
		CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
		return;
   }
   if (awidth <= 0) { 
		CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
		return;
   }
// Get nominal bin width in exponential form
   int jlog   = int(Log10(awidth));
   if (jlog <-200 || jlog > 200) {
	  BinLow   = 0;
	  BinHigh  = 1;
	  BinWidth = 0.01;
	  nbins    = 100;
	  CallL90(al, ah, BinLow, BinHigh, nbins, nold, BinWidth, ntemp, optionTime);
	  return;
   }
   if (awidth <= 1 && (!optionTime || timemulti==1) ) jlog--;
   double sigfig = awidth*pow(10,-jlog) -1e-10;
   // in the above statement, it is important to subtract 1e-10
   // to avoid precision problems if the tests below

//  Round mantissa
   double siground = 0;
	if      (sigfig <= 1)    siground = 1;
	else if (sigfig <= 2)    siground = 2;
	else if (sigfig <= 5 && (!optionTime || jlog<1))  siground = 5;
	else if (sigfig <= 6 && optionTime && jlog==1)    siground = 6;
	else { siground = 1;   jlog++; }

   BinWidth = siground*pow(10,jlog);
   if (optionTime) BinWidth *= timemulti;

   CallL90(al, ah, BinLow, BinHigh, nbins, nold, BinWidth, ntemp, optionTime);
}

void CallL90(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, int &nold, double &BinWidth, int &ntemp, int optionTime){
   double alb  = al/BinWidth;
   if (abs(alb) > 1e9) {
	  BinLow  = al;
	  BinHigh = ah;
	  if (nbins > 10 * nold && nbins > 10000) nbins = nold;
	  CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
	  return;
   }
   int lwid   = int(alb);
   if (alb < 0) lwid--;
   BinLow     = BinWidth*double(lwid);
   alb        = ah/BinWidth + 1.00001;
   int kwid = int(alb);
   if (alb < 0) kwid--;
   BinHigh = BinWidth*double(kwid);
   nbins = kwid - lwid;
   if (nold == -1){
		CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
		return;
   }
   if (nold <= 5) {          //    Request for one bin is difficult case
	  if (nold > 1 || nbins == 1){
		CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
		return;
	  }
	  BinWidth = BinWidth * 2;
	  nbins    = 1;
	  CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
	  return;
   }
   if (2*nbins == nold && !optionTime) {
		ntemp++; 
		CallL20(al, ah, BinLow, BinHigh, nbins, nold, BinWidth, ntemp, optionTime);
		return;
	}
	CallLOK(al, ah, BinLow, BinHigh, nbins, BinWidth, optionTime);
}


void CallLOK(double &al, double &ah, double &BinLow, double &BinHigh, int &nbins, double &BinWidth, int optionTime){
   double oldBinLow = BinLow;
   double oldBinHigh = BinHigh;
   int oldnbins = nbins;

   double atest = BinWidth * 0.0001;
   if (al-BinLow  >= atest) { BinLow  += BinWidth;  nbins--; }
   if (BinHigh-ah >= atest) { BinHigh -= BinWidth;  nbins--; }
   if (!optionTime && BinLow >= BinHigh) {
	  //this case may happen when nbins <=5
	  BinLow = oldBinLow;
	  BinHigh = oldBinHigh;
	  nbins = oldnbins;
   }
   else if (optionTime && BinLow>=BinHigh) {
	  nbins = 2 * oldnbins;
	  BinHigh = oldBinHigh;
	  BinLow = oldBinLow;
	  BinWidth = (oldBinHigh - oldBinLow)/nbins;
	  atest = BinWidth * 0.0001;
	  if (al-BinLow  >= atest) { BinLow  += BinWidth;  nbins--; }
	  if (BinHigh-ah >= atest) { BinHigh -= BinWidth;  nbins--; }
   }
}





void OptimizeTicksLimits(int nbins, int &newbins, double &xmin, double &xmax, bool isInteger){
   double binlow = 0,binhigh = 0,binwidth=0;
   int n=0;
   double dx = 0.1 * (xmax - xmin);
   if (isInteger) dx = 5 * (xmax - xmin)/nbins;
   double fmin = xmin - dx;
   double fmax = xmax + dx;
   if (fmin < 0 && xmin >= 0) fmin = 0;
   if (fmax > 0 && xmax <= 0) fmax = 0;

   OptimizeTicks(fmin, fmax, nbins, binlow, binhigh, n, binwidth);

	nbins = n;
	if (binwidth <= 0 || binwidth > 1.e+39) {
		xmin = -1;
		xmax = 1;
	}
	else if(xmin < binlow || xmax > binhigh){
		fmin = binlow;
		while(xmin < fmin) {
			fmin -= binwidth;
			nbins++;
		}
		fmax = binhigh;
		while(xmax > fmax) {
			fmax += binwidth;
			nbins++;
		}
		xmin = fmin;
		xmax = fmax;
	} 
	else {
		xmin    = binlow;
		xmax    = binhigh;
	}
	newbins = nbins;
}















