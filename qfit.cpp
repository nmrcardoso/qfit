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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cassert>
#include <map>
#include <vector>
#include <cmath>
#include <float.h>


#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Hor_Value_Slider.H>
#include <FL/Fl_Input.H>
//#include <FL/Fl_run.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Select_Browser.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/fl_ask.H>
#include <FL/Fl_RGB_Image.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Printer.H>
#include <FL/Fl_Check_Button.H>


#include "ticks.h"
#include "fit.h"
#include "plot.h"
#include "qfit.h"

#include "log.h"

#include "save_image.h"

#include "jackerr.h"


#ifndef uint
#define uint unsigned int
#endif


using namespace std;

mylog qlog;





MyGlWindow::MyGlWindow(int X, int Y, int W, int H, ConstFit& _fit) : Fl_Widget(X,Y,W,H), fit(_fit), width(W), height(H) { 
	plotVr = false; 
	fitVr = false;
	
	dx = 4;
	dy = 4;
	ldxy = 6;
}
void MyGlWindow::draw() {
    fl_push_clip( x(), y(), w(), h() );
	createdraw();	
    fl_pop_clip();	
}


void MyGlWindow::createdraw() {

	//Draw plot
	double t1, t2;
	t1 = fit.range.t1;
	t2 = fit.range.t2;
	if(plotVr){ t1 = rangeVr.t1; t2 = rangeVr.t2; }
	//Plot points
	DrawErrorPlot plot(x(), y(), w(), h());
	plot.setBackColor(FL_LIGHT3);
	plot.setBackgroundColor(FL_LIGHT2);
	plot.setFRangeColor(FL_LIGHT1);
	plot.setErrorLineWidth(2);
	plot.setPointSize(3);
	plot.setAxisFontSize(16);
	plot.setAxisFontColor(FL_BLACK);
	plot.setAxisLineColor(FL_BLACK);
	plot.setAxisLineWidth(1);
	plot.setAxisX();
	plot.setAxisY();
	plot.setGridLines(1);
	plot.setTickLength(20);
	plot.setTickWidth(1);
	plot.setTickColor(FL_BLACK);
	plot.draw(data, t1, t2);

	xmin = plot.xmin;
	xmax = plot.xmax;
	ymin = plot.ymin;
	ymax = plot.ymax;
	w_x0 = plot.w_x0;
	w_x1 = plot.w_x1;
	w_y0 = plot.w_y0;
	w_y1 = plot.w_y1;
	w_w = w_x1-w_x0;
	w_h = w_y1-w_y0;

	double xrange = (xmax-xmin), yrange = (ymax-ymin);
	double factorx = w_w/xrange;
	double factory = w_h/yrange;


	if(!plotVr){
		//Plot fit line
		fl_line_style(0,2);
		fl_color(FL_RED);
		fl_line(w_x0, w_y1-(fit.mean-ymin)*factory, w_x1, w_y1-(fit.mean-ymin)*factory);
	}	
	if(plotVr && fitVr){
		//Fit and draw results and fit line for V(r) final
		DataLine data1;
		for( uint i = 0; i < data.size(); ++i ) {
			if ( data[i].x >= rangeVr.t1 && data[i].x <= rangeVr.t2 )
				data1.push_back(data[i]);
		}
		res = GSLfit(data1);

		fl_line_style(0,2);
		fl_color(FL_RED);
		double stepxx = xrange/1000.;
		for( double i = xmin; i < xmax; i+=stepxx ){
			double y = res.val[0] + res.val[1] / i + res.val[2] * i; 
			double yy = res.val[0] + res.val[1] / (i+stepxx) + res.val[2] * (i+stepxx); 
			//fl_circle(w_x0+(i-xmin)*factorx, w_y1-(y-ymin)*factory, .4);
			fl_line(w_x0+(i-xmin)*factorx, w_y1-(y-ymin)*factory, w_x0+(i+stepxx-xmin)*factorx, w_y1-(yy-ymin)*factory);
		}
		fl_font(FL_BOLD, 16);
		string beg[3] = {"A = ", "B = ", "C = "};
		string leg = "Fit: A + B / r + C * r\n";
		leg += "Range: [" + ToString(data1[0].x) + "," + ToString(data1[data1.size()-1].x) + "]\n";
		for(int i = 0; i < 3; ++i){
			leg += beg[i] + ToString(res.val[i]) + " ± " + ToString(res.errorJ[i]) + "\n";
		}
		leg += "chi2/dof = " + ToString(res.chi2_dofJ);
		qlog << leg << endl;
		int wi=0, hi=0;
		fl_measure(leg.c_str(), wi, hi);
		int startx = w_x1-wi-4*dx;
		int starty = w_y1-hi-4*dy;
		fl_color(52);
		fl_rectf(startx, starty, wi, hi );
		fl_color(FL_BLACK);
	 	fl_draw(leg.c_str(), startx, starty, wi, hi, FL_ALIGN_LEFT);

	}
}


int MyGlWindow::convertLeftButton(int pos, double x0, double x1, double tmin, double tmax) const  {
	return (int)(0.5 + tmin + (pos - x0) / x1 * ( 1 + tmax - tmin ));
}
int MyGlWindow::convertRightButton(int pos, double x0, double x1, double tmin, double tmax) const  {
	return (int)(-0.5 + tmin + double(pos - x0) / x1 * ( 1 + tmax - tmin ));
}

void MyGlWindow::Mouse_Handler_GLWiN(int pos, int whichButtonCliked){
	if(pos>=x() && pos<=x()+w()){// && !plotVr){
		switch(whichButtonCliked){
			case 1: 
				if(!plotVr){
					int val = convertLeftButton(pos, w_x0, w_w, xmin, xmax);
					if( val < fit.range.t2) {
						if(val < Tmin) val = Tmin;
						fit.range.t1 = val;
					}
				}
				else if(fitVr){ 
					int val = convertLeftButton(pos, w_x0, w_w, Vx_min, Vx_max);
					int n=0;
					for(int i = 0; i < data.size(); i++)
						if(data[i].x <= rangeVr.t2 && data[i].x >= val) n++;
					if( val < rangeVr.t2 && n > 3) {
						if(val < Vx_min) val = Vx_min;
						rangeVr.t1 = val;
					}
				}
				break;
			case 3:
				if(!plotVr){
					int val = convertRightButton(pos, w_x0, w_w, xmin, xmax);
					if( val > fit.range.t1) {
						if(val > Tmax) val = Tmax;
						fit.range.t2 = val;
					}
				}
				else if(fitVr){
					int val = convertRightButton(pos, w_x0, w_w, Vx_min, Vx_max);
					int n=0;
					for(int i = 0; i < data.size(); i++)
						if(data[i].x <= val && data[i].x >= rangeVr.t1) n++;
					if( val > rangeVr.t1 && n > 3) {
						if(val > Vx_max) val = Vx_max;
						rangeVr.t2 = val;
					}
				} 
				break;
			default: 
				qlog << "NONE"<< endl; 
				break;
		}
	}
}

int MyGlWindow::handle(int event) {
	switch (event) {
        case FL_ENTER:
            //window()->cursor(FL_CURSOR_CROSS);
            window()->cursor(FL_CURSOR_HAND);
            return 1;
		case FL_LEAVE:
            window()->cursor(FL_CURSOR_DEFAULT);
			return 1;
		case FL_FOCUS :
			Fl::focus(this);
			redraw();
			return 1;
		case FL_UNFOCUS :
			redraw();
			return 1;
		case FL_PUSH :
			switch(Fl::event_button()){
	        	case FL_LEFT_MOUSE:
					Mouse_Handler_GLWiN((int)Fl::event_x(), 1);
					redraw();
  					do_callback();
					return 1;
		        break;
			    case FL_RIGHT_MOUSE:
					Mouse_Handler_GLWiN((int)Fl::event_x(), 3);
					redraw();
  					do_callback();
					return 1;
				break;
			}
	  	break;
	}
	return Fl_Widget::handle(event);
}






















class Program{
	Fl_Window window;
	MyGlWindow *glwin;
	Fl_Slider *minslider, *maxslider;
	Fl_Button *nextR, *prevR;
	Fl_Output* boxR1R2;
	Fl_Output* fitres;
	Fl_Input* fileinp;
	Fl_Button *save_res, *save_data, *exit_prog, *plotVr, *printData, *save_png;
	Fl_Select_Browser *browser;
	Fl_Check_Button *fitVr;

	FileDataCols fcols;
	Data data;
	Data::iterator dataIter;
	ConstFit fit;
	map<RPoint,FitRange,LessRPoint> ranges;
	int Tmax = 1000000;
	int Tmin = 0;
	int x0, x1;
	string rawname;

	RPoint R1R2() {
		return dataIter->first;
	}


	static void changeR1R2( Fl_Widget* w, void* data ) {
		Fl_Button* button = (Fl_Button*)w;
		Program* program = (Program*)data;
		assert( button && program );
		if ( button == program->nextR )
			program->nextR1R2();
		else if ( button == program->prevR )
			program->prevR1R2();
		else
			assert(0);
		program->glwin->data = program->dataIter->second;
		program->updateAll();
		program->glwin->redraw();
		Fl::redraw();
		Fl::flush();
	}
	
	void prevR1R2() {
		if (dataIter != data.begin()){
			--dataIter;
			int loc = browser->value();
			browser->select(loc,0);
			browser->value(loc-1);
		}
	}
	void nextR1R2() {
		Data::iterator new_it = dataIter;
		if (++new_it != data.end()){
			dataIter = new_it;
			int loc = browser->value();
			browser->select(loc,0);
			browser->value(loc+1);		
		}
	}

	static void checkSliderRange( Fl_Widget* w, void* data ) {
  		Fl_Slider* slider = (Fl_Slider*)w;
		Program* program = (Program*)data;
		if ( slider == program->minslider ) {
			if ( slider->value() < program->Tmin )
				slider->value( program->Tmin );
			else if ( slider->value() >= program->maxslider->value() - 1 )
				slider->value( program->maxslider->value() - 1 );
		} else if ( slider == program->maxslider )	{
			if ( slider->value() <= program->minslider->value() + 1 )
				slider->value( program->minslider->value() + 1 );
			else if ( slider->value() > program->Tmax )
				slider->value( program->Tmax );
		}
		program->ranges[program->R1R2()].t1 = int( program->minslider->value() );
		program->ranges[program->R1R2()].t2 = int( program->maxslider->value() );
		program->updateAll();
	}


	static void glwinDraw( Fl_Widget* w, void* data ) {
  		MyGlWindow* glmywin = (MyGlWindow*)w;
		Program* program = (Program*)data;
		program->ranges[program->R1R2()].t1 = glmywin->fit.range.t1;
		program->minslider->value(glmywin->fit.range.t1);
		program->ranges[program->R1R2()].t2 = glmywin->fit.range.t2;
		program->maxslider->value(glmywin->fit.range.t2);
		program->updateAll();
	}


	static void PlotResults( Fl_Widget* w, void* data ) {
		Fl_Button* button = (Fl_Button*)w;
		Program* program = (Program*)data;
		if( !program->glwin->plotVr ){
			DataLine dataVr;
			for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
				RPoint rpoint = it->first;
				FitRange range = program->ranges[rpoint];
				ConstFit fit = constFit( it->second, range.t1, range.t2 );
				DataPoint d; 
				d.x = rpoint.r2;
				d.y = fit.mean;
				d.error = fit.error;
				dataVr.push_back(d);
			}
			program->glwin->data = dataVr;
			program->glwin->plotVr = true;
			program->glwin->Vx_min = dataVr[0].x;
			program->glwin->Vx_max = dataVr[dataVr.size()-1].x;

			if(program->glwin->rangeVr.t1 > dataVr[dataVr.size()-1].x-3 || program->glwin->rangeVr.t1 < dataVr[0].x) 
				program->glwin->rangeVr.t1 = dataVr[0].x;
			if(program->glwin->rangeVr.t2 > dataVr[dataVr.size()-1].x || program->glwin->rangeVr.t2 < dataVr[0].x+3) 
				program->glwin->rangeVr.t2 = dataVr[dataVr.size()-1].x;

			button->label( "Back" );
			program->disableEnableButtons();
			program->updateAll();		
		}
		else{
			program->glwin->data = program->dataIter->second;
			program->glwin->plotVr = false;
			button->label( "Plot V(r)" );
			program->disableEnableButtons();
			program->updateAll();
		}
	}


	static void PrintResultsPNG( Fl_Widget* w, void* data ) {
		Fl_Button* button = (Fl_Button*)w;
		Program* program = (Program*)data;
		// save plot in png format
		string filename = program->file_output_name() + ".png";
		int ww=program->glwin->w(), hh=program->glwin->h();
		Image im(filename.c_str(), ww, hh, program->glwin);
		qlog << "Saving actual plot in: " << filename << endl;
		im.write(Image::PNG);		
		// save plot in tga format
		/*filename = program->file_output_name() + ".tga";
		qlog << "Saving actual plot in: " << filename << endl;
		Image im1(filename.c_str(), ww, hh, program->glwin);
		im1.write(Image::TGA);		
		// save plot in bmp format
		filename = program->file_output_name() + ".bmp";
		qlog << "Saving actual plot in: " << filename << endl;
		Image im1(filename.c_str(), ww, hh, program->glwin);
		im1.write(Image::BMP);
		// save plot ppm format
		filename = program->file_output_name() + ".ppm";
		qlog << "Saving actual plot in: " << filename << endl;
		Image im1(filename.c_str(), ww, hh, program->glwin);
		im1.write(Image::PPM);*/
	}

	static void PrintResults( Fl_Widget* w, void* data ) {
		Fl_Button* button = (Fl_Button*)w;
		Program* program = (Program*)data;

		Fl_Printer *printer = new Fl_Printer();
		if (printer->start_job(1) == 0) {
			int width, height;
			printer->start_page();
			printer->printable_rect(&width, &height);
			
			int fontsize = 17;
			fl_font(FL_COURIER_BOLD, fontsize);
			time_t now; time(&now);
			int wi=0, hi=0;
			fl_measure(ctime(&now), wi, hi);

			float sw = float(width)/float(program->glwin->w()*2);
			float sh = float(height)/float((program->glwin->h()+3*hi)*(program->data.size()+2)/2);
			sw = float(width)/float(program->glwin->w());
			sh = float(height)/float((program->glwin->h()+2*hi));
			
			
			int numg = program->data.size() + 1;
			if( program->fcols.nc != 1 ) numg = program->data.size();
		    float scale;			
			int nx, ny, newnx,newny;
			float plotsize = 0;
			float extrasp = program->glwin->w()+width;
			for(int i = 1; i < numg; ++i){
				newnx = i;
				newny = (numg+newnx - 1) / newnx;
				sw = float(width)/float(program->glwin->w()*newnx);
				sh = float(height)/float((program->glwin->h()+3*hi)*newny);
				scale = sw;
				if (sw > sh) scale = sh;	
				float plotsize0 = program->glwin->w()*scale;	
				float extrasp0 = width - program->glwin->w()*newnx*scale;
				if(plotsize0 >= plotsize && extrasp0 <= extrasp) { 
					nx = newnx; 
					ny = newny; 
					plotsize = plotsize0; 
					extrasp = extrasp0;
				}
			}
			sw = float(width)/float(program->glwin->w()*nx);
			sh = float(height)/float((program->glwin->h()+3*hi)*ny);
			if (sw > sh) sw = sh;
			printer->scale(sw);

			if(width > height){
				printer->rotate(180);
				printer->translate(0,0);
			}
			printer->origin(0, 0);
			int hh=0;
			int ww=0, iw=0;
			Data::iterator dataIter1 = program->dataIter;
			bool plotVrr = program->glwin->plotVr;
			bool fitVrr = program->glwin->fitVr;
			if( program->fcols.nc == 1 ) {
				program->glwin->plotVr = false;
				program->glwin->fitVr = true;
			}
			for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
				program->dataIter = it;
				program->glwin->data = program->dataIter->second;
				program->updateAll();
				program->glwin->redraw();
				printer->origin(ww, hh+hi);
				fl_color(FL_BLACK);
				fl_font(FL_COURIER_BOLD, fontsize);
				stringstream text;
				text << "V(";
				if(program->fcols.nc != 1) text << program->dataIter->first.r1 << ",";
				text << program->dataIter->first.r2 << ") = " << program->fit.mean << " ± " << program->fit.error << "   chi²/dof = " << program->fit.chi2dof << "\n[t1,t2] = [" << program->ranges[program->dataIter->first].t1 << "," << program->ranges[program->dataIter->first].t2 << "]";
				
				fl_draw(text.str().c_str(), 0, 0);
				printer->origin(ww, hh+2*hi);
				printer->print_widget(program->glwin, 0, 0);

				if(iw < nx-1){
					iw++;
					ww+= width/sw/nx;
				}
				else{ 
					iw=0;
					hh+=program->glwin->h()+3*hi;
					ww=0;
				}
			}
			printer->origin(ww, hh+2*hi);

			if( program->fcols.nc == 1 ) {
				DataLine dataVr;
				for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
					RPoint rpoint = it->first;
					FitRange range = program->ranges[rpoint];
					ConstFit fit = constFit( it->second, range.t1, range.t2 );
					DataPoint d; 
					d.x = rpoint.r2;
					d.y = fit.mean;
					d.error = fit.error;
					dataVr.push_back(d);
				}
				program->glwin->data = dataVr;
				program->glwin->plotVr = true;			
				program->glwin->fitVr = fitVrr;
				program->glwin->Vx_min = dataVr[0].x;
				program->glwin->Vx_max = dataVr[dataVr.size()-1].x;
				if(program->glwin->rangeVr.t1 > dataVr[dataVr.size()-1].x-3 || program->glwin->rangeVr.t1 < dataVr[0].x) program->glwin->rangeVr.t1 = dataVr[0].x;
				if(program->glwin->rangeVr.t2 > dataVr[dataVr.size()-1].x || program->glwin->rangeVr.t2 < dataVr[0].x+3) program->glwin->rangeVr.t2 = dataVr[dataVr.size()-1].x;
				program->updateAll();

				printer->print_widget(program->glwin, 0, 0);

				program->glwin->plotVr = plotVrr;
				program->glwin->fitVr = fitVrr;
				program->dataIter = dataIter1;
				if(!program->glwin->plotVr) program->glwin->data = program->dataIter->second;
				program->disableEnableButtons();
				}
			program->updateAll();
			if(width > height){
					printer->untranslate();
			}
			printer->end_page();
			printer->end_job();
		}		
		delete printer;
	}




	static ConstFit constFit( const DataLine& dataline, int t1, int t2 ) { //implementar erros jackniffe
		ConstFit fit;
		fit.range.t1 = t1;
		fit.range.t2 = t2;
		double s0 = 0, s1 = 0;
		int num = 0;
		for( DataLine::const_iterator it = dataline.begin(); it != dataline.end(); ++it )
			if ( it->x >= t1 && it->x <= t2 ) {
				double err2 = it->error * it->error;
				s0 += 1.0 / err2;
				s1 += it->y / err2;
				++num;
			}
		fit.mean = s1 / s0;
		vector<double> trials;
		double s2 = 0;
		for( auto p : dataline ) {
			if (p.x >= t1 && p.x <= t2) {
				double err2 = p.error * p.error;
				trials.push_back( (s1 - p.y/err2) / (s0 - 1./err2) );
				s2 += pow( (p.y - fit.mean) / p.error, 2 );
			}
		}
		fit.error = jackerr(trials);
		//fit.error = jackerr(trials, fit.mean);
		fit.chi2dof = s2 / (num - 1);
		return fit;
	}



	static void exitprog( Fl_Widget* w, void* data ) {
		switch ( fl_choice("Are you sure you want to quit?", "Yes", "No", 0) ) {
		  case 0: exit(0); // Yes
		  case 1:  return;// No (default)
		}
	}



	static void saveResults( Fl_Widget* w, void* data ) {
		Program* program = (Program*)data;
		string filename = program->file_output_name() + ".res";
		qlog << "Saving V to: " << filename << endl; 
		ofstream out(filename.c_str());
		if( !out.is_open()){
			qlog << "Error opening file...\nExiting..." << endl;
			exit(1);
		}
		out.precision(12);
		//out.setf(ios_base::scientific);
		for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
			RPoint rpoint = it->first;
			FitRange range = program->ranges[rpoint];
			ConstFit fit = constFit( it->second, range.t1, range.t2 );
			if(program->fcols.nc == 1)
				out << rpoint.r2 << '\t' << fit.mean << '\t' << fit.error << '\t' << fit.chi2dof << '\t' << range.t1 << '\t' << range.t2 << '\n';
			else
				out << rpoint.r1 << '\t' << rpoint.r2 << '\t' << fit.mean << '\t' << fit.error << '\t' << fit.chi2dof << '\t' << range.t1 << '\t' << range.t2 << '\n';
		}
	}
	static void saveData( Fl_Widget* w, void* data ) {
		Program* program = (Program*)data;
		string filename = program->file_output_name() + ".qfit";
		qlog << "Saving results to: " << filename << endl;
		ofstream out(filename.c_str());
		if( !out.is_open()){
			qlog << "Error opening file...\nExiting..." << endl;
			exit(1);
		}
		out.precision(12);
		//out.setf(ios_base::scientific);
		out << program->data.size() << '\n';

		for( map<RPoint,FitRange,LessRPoint>::iterator it = program->ranges.begin(); it != program->ranges.end(); ++it ) {
			RPoint rpoint = it->first;
			FitRange range = it->second;
			if(program->fcols.nc != 1) out << rpoint.r1 << '\t';
			out << rpoint.r2 << '\t' << range.t1 << '\t' << range.t2 << '\n'; 
			}

		for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
			RPoint rpoint = it->first;
			const DataLine& line = it->second;
			for( size_t i = 0; i < line.size(); ++i ){
				if(program->fcols.nc != 1) out << rpoint.r1 << '\t';
				out << rpoint.r2 << '\t' << line[i].x << '\t' << line[i].y << '\t' << line[i].error << '\n';
			}
		}
	}
	void updateAll() {
		const DataLine& dataline = data[R1R2()];
		fit = constFit( dataline, ranges[R1R2()].t1, ranges[R1R2()].t2 );
		updateResult();
		updateSliders();
		redrawGraph();
	}

	void updateSliders() {
		minslider->value( ranges[R1R2()].t1 );
		maxslider->value( ranges[R1R2()].t2 );
	}

	void updateResult() {
		stringstream text;
		text << "V(";
		if(fcols.nc != 1) text << R1R2().r1 << ",";
		text << R1R2().r2 << ") = " << fit.mean << " ± " << fit.error << "   chi²/dof = " << fit.chi2dof << "    [t1,t2] = [" << ranges[R1R2()].t1 << "," << ranges[R1R2()].t2 << "]";
		fitres->value(text.str().c_str());
		text.str("");
		text.clear();
		if(fcols.nc == 1) text << R1R2().r2;
		else text << R1R2().r1 << "," << R1R2().r2;
		boxR1R2->value( text.str().c_str() );
	}


 	void CreateBrowserContents() {
		for( Data::iterator it = data.begin(); it != data.end(); ++it ) {
			RPoint rpoint = it->first;
			stringstream text;
			if( fcols.nc == 1) text << "@c@b@ " << rpoint.r2;
			else text << "@c@b@ " << rpoint.r1 << "," << rpoint.r2;
			browser->add(text.str().c_str());
		}
		browser->value(1);
	}

	static void changefitVr(Fl_Widget* w, void* data ) {
  		Fl_Check_Button* but = (Fl_Check_Button*)w;
		Program* program = (Program*)data;
		if(but->value()) {
			but->set();
			program->glwin->fitVr = true;
		}
		else{
			but->clear();
			program->glwin->fitVr = false;
		}
		program->redrawGraph();
	}


 	static void select_brow_cb(Fl_Widget* w, void* data ) {
  		Fl_Select_Browser* browse = (Fl_Select_Browser*)w;
		Program* program = (Program*)data;
		assert( browse && program );
		int location = ((int)program->browser->value());
		int pos = 1;
		program->browser->deselect();
		program->browser->select(location);
		for( Data::iterator it = program->data.begin(); it != program->data.end(); ++it ) {
			if(pos == location){
				program->dataIter = it;				
				break;
			}
			pos+=1;
		}
		program->glwin->data = program->dataIter->second;
		program->updateAll();
		program->glwin->redraw();
		Fl::redraw();
		Fl::flush();
	}

 	string file_output_name() {
		string name = fileinp->value();
		if(name.size() < 1 ) name = rawname;
		return name;
	}

	double convertvalue(string str){
		stringstream ss(str);
		double val;
		ss >> val;
		//No need to proper read "inf" and "nan" values for now		
		if(false){
			if(str == "Inf" || str == "inf"){
				 val = std::numeric_limits<double>::infinity();
			}
			if(str == "-Inf" || str == "-inf"){
				 val = -std::numeric_limits<double>::infinity();
			}
			if(str == "Nan" || str == "nan"){
				 val = std::numeric_limits<double>::quiet_NaN();
			}
			if(str == "-Nan" || str == "-nan"){
				 val = -std::numeric_limits<double>::quiet_NaN();
			}
		}
		return val;
	}
	
	
	
	void loadvalues(ifstream &in){
		RPoint rpoint;
		DataPoint datapoint;
    	std::string str;
		double *data_read = new double[fcols.tc];
		if(fcols.nc == 1) qlog << "r \tt \tV(r,t) \tV_error" << endl;
		else qlog << "r1 \tr2 \tt \tV(r1,r2,t) \tV_error" << endl;
		while( ! in.eof() ) {
			for( int i = 0; i < fcols.tc; ++i){
				in >> str; 
				data_read[i] = convertvalue(str);
			}
			if(in.fail()) break;
			if(fcols.nc == 1){
				rpoint.r1 = 0;
				rpoint.r2 = data_read[0];
				datapoint.x = data_read[1];
				datapoint.y = data_read[fcols.d];
				datapoint.error = data_read[fcols.d+1];
			}
			else{
				rpoint.r1 = data_read[0];
				rpoint.r2 = data_read[1];
				datapoint.x = data_read[2];
				datapoint.y = data_read[fcols.d];
				datapoint.error = data_read[fcols.d+1];
			}
			if(fcols.nc != 1)
				qlog << rpoint.r1 << '\t';
			qlog << rpoint.r2 << '\t' << datapoint.x << '\t' << datapoint.y << '\t' << datapoint.error << endl;
			if ( data.find(rpoint) == data.end() )
				data[rpoint] = DataLine();
			data[rpoint].push_back( datapoint );
		}
		qlog << "------------------------------------------------------" << endl;
		delete[] data_read;
		qlog << flush;	
	}
	void setTlimits(bool setrange = false){
		Tmax = 0;
		Tmin = 100;
		double tmaxx = 0, tmin=100;
		for( Data::iterator it = data.begin(); it != data.end(); ++it ) {
			FitRange range;
			tmaxx = it->second[it->second.size()-1].x;
			if(Tmax < tmaxx) Tmax = tmaxx;
			tmin = it->second[0].x;
			if(Tmin > tmin) Tmin = tmin;
			if(setrange){
				range.t1 = Tmin;
				range.t2 = Tmax;
				ranges[it->first] = range;
			}
		}
		glwin->Tmin=Tmin;
		glwin->Tmax=Tmax;	
	}

	public:
	void load( string filename ) {
		qlog << "Reading data from: " << filename << endl;
		ifstream input( filename.c_str() );
		if (!input.is_open()){
            qlog <<"Error saving configuration: " << filename << endl;
            exit(1);
        }
		data.clear();
		loadvalues( input);
		input.close();	
		dataIter = data.begin();
		glwin->data = dataIter->second;
		setTlimits(true);
	}
	void loadmfit( string filename ) {
		qlog << "Reading data from: " << filename << endl;
		ifstream input(filename.c_str());
		if (!input.is_open()){
            qlog << "Error saving configuration: " << filename << endl;
            exit(1);
        }
		size_t sz = 0;
		input >> sz;
		ranges.clear();
		data.clear();
		if(fcols.nc == 1) qlog << "r \t";
		else qlog << "r1 \tr2 \t" << endl;
		qlog << "tmin \ttmax" << endl;		
		for( size_t i = 0; i < sz; ++i ) {
			RPoint rpoint;
			FitRange range;
			if(fcols.nc==1){
				rpoint.r1 = 0.;
				input >> rpoint.r2 >> range.t1 >> range.t2;
			}
			else input >> rpoint.r1 >> rpoint.r2 >> range.t1 >> range.t2;
			ranges[rpoint] = range;
			if(fcols.nc != 1)
				qlog << rpoint.r1 <<'\t';
			qlog << rpoint.r2 << '\t' << range.t1 << '\t' << range.t2 << endl;
		}
		qlog << "------------------------------------------------------" << endl;		
		loadvalues( input);
		input.close();	
		
		assert( data.size() == sz );
		for( Data::iterator it = data.begin(); it != data.end(); ++it ) {
			RPoint rpoint = it->first;
			assert( ranges.find(rpoint) != ranges.end() );
		}
		dataIter = data.begin();
		glwin->data = dataIter->second;
		setTlimits(false);
	}

	void disableEnableButtons(){
		if(glwin->plotVr){
			plotVr->label( "Back" );
			minslider->deactivate();
			maxslider->deactivate();
			browser->deactivate();
			prevR->deactivate();
			nextR->deactivate();
			minslider->hide();
			maxslider->hide();
			browser->hide();
			prevR->hide();
			nextR->hide();
			fitres->hide();
			boxR1R2->hide();
			fitres->deactivate();
			boxR1R2->deactivate();
			fitVr->show();
			fitVr->activate();
		}
		else{
			plotVr->label( "Plot V(r)" );
			minslider->show();
			maxslider->show();
			browser->show();
			prevR->show();
			nextR->show();
			minslider->activate();
			maxslider->activate();
			browser->activate();
			prevR->activate();
			nextR->activate();
			fitVr->deactivate();
			fitVr->hide();
			fitres->activate();
			fitres->show();
			boxR1R2->show();
			boxR1R2->activate();
		}
	}

	Program( int argc,  const char **argv, int x = 100, int y = 100, int w = 800, int h = 600, const char* name = "qfit" ) : window( x, y, w, h, name  ){
		// Defining input file reading parameters
		// File input columns can be defined as:
		// 1- r1 r2 t V(r1,r2,t) error  -> nc = 2
		//   ./program inputfile 2
		// 2- r t V(r,t) error -> no nc or nc = 1
		//   ./program inputfile
		//   ./program inputfile 1
		// 3- r t V1(r,t) error_V1 V2(r,t) error_V2 ...
		//   ./program inputfile nc tc d 
		//   tc total number of columns, d column with V to read, error should be d+1, d column id starts at 0.
		string filename = argv[1];
		//qlog << "number of parameters: " << argc << endl;
		if( argc == 2){
			fcols.nc=1;
			fcols.tc = 4;
			fcols.d = 2;
		}
		else if(argc==3){
			fcols.nc=atoi(argv[2]);
			if( fcols.nc == 1 ){
				fcols.tc = 4;
				fcols.d = 2;
			}
			else if( fcols.nc == 2 ){
				fcols.tc = 5;
				fcols.d = 3;
			}
			else{
				qlog << "Not implemented for nc = " << fcols.nc << endl;
				qlog << "Only implemented for nc = 1 or 2 (columns of r)" << endl;
				exit(1);
			}
		}
		else if( argc == 5){
			fcols.nc=atoi(argv[2]);
			fcols.tc = atoi(argv[3]);
			fcols.d = atoi(argv[4]);
		}
		else{
			qlog << "Not implemented..." << endl;
			exit(1);
		}
		window.begin();
			//Window Icon
			uchar buffer[32*32*3];
			Fl_RGB_Image icon(buffer, 32, 32, 3);
			//icon.color_average(FL_BLUE, 0.0);
			window.icon(&icon);

			//Browse r's
			browser = new Fl_Select_Browser(0,10,120,h-210,0);
			browser->type(FL_MULTI_BROWSER);
			browser->selection_color(FL_BLUE);
			//browser->color(239);
			browser->color(52);
			browser->textsize(18);
			browser->callback(select_brow_cb, this);
			
			//Print fit result
			fitres = new Fl_Output( 10, h-195, w-20, 40 );
			//fitres->textfont(FL_COURIER);
			fitres->textsize(18);
			fitres->textfont(FL_BOLD);
			//fitres->color(FL_BLUE);
			fitres->color(52);
	  		fitres->tooltip("Display fit result...");

			//Graphical part of the fit
			x0 = 120;
			x1 = w-x0;
			glwin = new MyGlWindow(  x0, 10, x1, h - 200-10, fit);
			glwin->callback(glwinDraw, this);
			//glwin->tooltip("Display data...");

			//Slider to choose lower fit range
			minslider = new Fl_Slider(  220, h - 120, w-100-220, 20 );
			minslider->type(FL_HOR_NICE_SLIDER);
			minslider->box(FL_FLAT_BOX);
			minslider->slider_size(.05);
			minslider->step( 1 );
			minslider->callback( checkSliderRange, this );

			//Slider to choose upper fit range
			maxslider = new Fl_Slider( 220, h - 140, w-100-220, 20 );
			maxslider->type(FL_HOR_NICE_SLIDER);
			maxslider->box(FL_FLAT_BOX);		
			maxslider->slider_size(.05);
			maxslider->step( 1 );
			maxslider->callback( checkSliderRange, this );

			//Plot V(r)
			if( fcols.nc == 1 ){
				plotVr = new Fl_Button( w-90, h - 140, 80, 40, "Plot V(r)" );
				plotVr->callback( PlotResults, this );
			}

			fitVr = new Fl_Check_Button(w-190, h - 140, 80, 40, "Fit V(r)");
			fitVr->callback( changefitVr, this );
			fitVr->deactivate();
			fitVr->hide();

			//Button to go back in r
			prevR = new Fl_Button( 10, h - 140, 40, 40, "@<-" );
			prevR->callback( changeR1R2, this );

			//Button to go foward in r
			nextR = new Fl_Button(160, h - 140, 40, 40, "@->" );
			nextR->callback( changeR1R2, this );

			//Print actual r in fit
			boxR1R2 = new Fl_Output( 55, h - 140, 100, 40 );
			boxR1R2->align(FL_ALIGN_RIGHT);
			boxR1R2->textsize(18);
			boxR1R2->textfont(FL_BOLD);
			//boxR1R2->color(239);
			boxR1R2->color(52);

			//Button to save V result
			save_res = new Fl_Button( 240, h - 60, 100, 40, "Save V" );
			save_res->callback( saveResults, this );

			//Button to save all results, and fit ranges, allow to read and retune fits
			save_data = new Fl_Button( 360, h - 60, 100, 40, "Save Results" );
			save_data->callback( saveData, this );

			//Button to save png
			save_png = new Fl_Button( w-320, h - 60, 90, 40, "Save in PNG" );
			save_png->callback( PrintResultsPNG, this );

			//Button to print
			printData = new Fl_Button( w-220, h - 60, 90, 40, "Print" );
			printData->callback( PrintResults, this );

			//Button to exit
			exit_prog = new Fl_Button( w-100, h - 60, 90, 40, "Exit" );
			exit_prog->callback( exitprog, this );

			//Field to change output file name to save results, witout extension
			fileinp = new Fl_Input( 20, h - 60, 200, 40, "Output filename" );
			fileinp->align(FL_ALIGN_TOP);
			//fileinp->color(222);
			fileinp->color(52);			
		window.end();

		size_t lastindex = filename.find_last_of("."); 
		rawname = filename.substr(0, lastindex);
		string extension = filename.substr(lastindex, filename.size());

		if( argc == 5) rawname += "_" + ToString(fcols.d/2);

		//load data
		if ( extension == ".qfit" )
			loadmfit( filename );
		else
			load( filename );

		//set output filename
		fileinp->value(rawname.c_str());

		//Set slider range
		minslider->range( Tmin, Tmax );
		minslider->value( Tmin );
		maxslider->range( Tmin, Tmax );
		maxslider->value( Tmax );
		
		//Add r's do browser
		CreateBrowserContents();

		window.show();
		updateAll();
	}

	void redrawGraph() {
		glwin->redraw();
		//Fl::redraw();
	}
};





int main(int argc,  const char **argv){
	qlog.set("qfit.log");
	if ( argc > 1 ){
		Program* program = new Program( argc, argv );
	}
	else{
		qlog << "Error: missing arguments...\n Exiting..." << endl;
		exit(0);
	}
	Fl::run();
}
