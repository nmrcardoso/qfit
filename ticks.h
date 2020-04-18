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

#ifndef TICKS_QFIT_H
#define TICKS_QFIT_H






#include <FL/fl_draw.H>


using namespace std;





// Color is the type we use to represent color. We can use Color like this:
//    grid.set_color(Color::red);
struct Color {
    enum Color_type {
        red=FL_RED,
        blue=FL_BLUE,
        green=FL_GREEN,
        yellow=FL_YELLOW,
        white=FL_WHITE,
        black=FL_BLACK,
        magenta=FL_MAGENTA,
        cyan=FL_CYAN,
        dark_red=FL_DARK_RED,
        dark_green=FL_DARK_GREEN,
        dark_yellow=FL_DARK_YELLOW,
        dark_blue=FL_DARK_BLUE,
        dark_magenta=FL_DARK_MAGENTA,
        dark_cyan=FL_DARK_CYAN
    };

    enum Transparency { invisible = 0, visible=255 };

    Color(Color_type cc) :c(Fl_Color(cc)), v(visible) { }
    Color(Color_type cc, Transparency vv) :c(Fl_Color(cc)), v(vv) { }
    Color(int cc) :c(Fl_Color(cc)), v(visible) { }
    Color(Transparency vv) :c(Fl_Color()), v(vv) { }    // default color

    int as_int() const { return c; }

    char visibility() const { return v; } 
    void set_visibility(Transparency vv) { v=vv; }
private:
    char v;    // invisible and visible for now
    Fl_Color c;
};


/*
    enum Font_type {
        helvetica=FL_HELVETICA,
        helvetica_bold=FL_HELVETICA_BOLD,
        helvetica_italic=FL_HELVETICA_ITALIC,
        helvetica_bold_italic=FL_HELVETICA_BOLD_ITALIC,
        courier=FL_COURIER,
        courier_bold=FL_COURIER_BOLD,
        courier_italic=FL_COURIER_ITALIC,
        courier_bold_italic=FL_COURIER_BOLD_ITALIC,
        times=FL_TIMES,
        times_bold=FL_TIMES_BOLD,
        times_italic=FL_TIMES_ITALIC,
        times_bold_italic=FL_TIMES_BOLD_ITALIC,
        symbol=FL_SYMBOL,
        screen=FL_SCREEN,
        screen_bold=FL_SCREEN_BOLD,
        zapf_dingbats=FL_ZAPF_DINGBATS
 */



class DrawTicks {
	//Ticks properties
	int NxDivLvl1, NxDivLvl2;
	int NyDivLvl1, NyDivLvl2;
	int ticks_length;
	int ticks_width;
	Fl_Font fontaxis;
	int fontsize;
	Fl_Color fontcolor;
	Fl_Color ticks_color;

	//Gridlines properties
	bool gridlines;
	int gridlinewidth;

	//Axis line properties
	int axislinewidth;
	Fl_Color axislinecolor;

	//Axis properties
	int w_x0, w_x1, w_y0, w_y1;
	int ScreenAxisLength_x, ScreenAxisLength_y;

	void setDefaultValues();
	protected:
	bool xaxis;
	bool yaxis;

	int getTickFontSize() const { return fontsize;}
	Fl_Font getTickFont() const { return fontaxis;}

	virtual void setTicks(int NxDivLvl1_, int NxDivLvl2_, int NyDivLvl1_, int NyDivLvl2_);
	public:

	DrawTicks();
	virtual ~DrawTicks() { };
	virtual void InitTicks(int w_x0_, int w_x1_, int w_y0_, int w_y1_);

	virtual void setAxisFontSize(int size);
	virtual void setAxisFontColor(Fl_Color fontcolor_);
	virtual void setAxisX();
	virtual void setAxisY();
	virtual void setAxisLineWidth(int size);
	virtual void setAxisLineColor(Fl_Color axislinecolor_);

	virtual void setTickLength(int size);
	virtual void setTickWidth(int width);
	virtual void setTickColor(Fl_Color color);

	virtual void setGridLines(int size);
	virtual void unsetGridLines();

	virtual void drawAxis(double xmin, double xmax, double ymin, double ymax);
};










void OptimizeTicksLimits(int nbins, int &newbins, double &xmin, double &xmax, bool isInteger);
void OptimizeTicks(double A1, double A2, int nold, double &BinLow, double &BinHigh, int &nbins, double &BinWidth);


class Ticks{
  private:
	int NDivLvl1; 		  //number of major ticks if possible
	int NDivLvl2; 		  //number of minor ticks if possible, if 0 no minor ticks
	int ScreenAxisLength;  //plot size in screen coordinates
	double funcMin;	  //minimum value of the function to plot
	double funcMax;     //maximum value of the function to plot

  public:
	int NumTicksLvl1;	  //final number of major ticks
	double *TicksPosLvl1;  //major ticks in screen coordinates
	double *TicksValuesLvl1; //major ticks in plot coordinates for labels
	int NumTicksLvl2;	  //final number of minor ticks
	double *TicksPosLvl2;  //minor ticks in screen coordinates

	Ticks();
	Ticks(int ScreenAxisLength, int NDivLvl1, int NDivLvl2, double funcMin, double funcMax);
	~Ticks();
	void CalculateTicks(int ScreenAxisLength_, int NDivLvl1_, int NDivLvl2_, double funcMin_, double funcMax_);
	void CalculateTicks();
};























#endif

