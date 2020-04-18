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

#ifndef STRUCTS_H
#define STRUCTS_H


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



#ifndef uint
//typedef unsigned int uint;
#define uint unsigned int
#endif


template<class T>
inline std::string ToString(const T number){
    std::stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}
template<>
inline std::string ToString<double>(const double number){
    std::stringstream ss;//create a stringstream
	//ss.precision(2);
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}
template<>
inline std::string ToString<float>(const float number){
    std::stringstream ss;//create a stringstream
	//ss.precision(2);
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}


double jackerr(const std::vector<double>& trials);

struct RPoint {
	RPoint() {}
	RPoint( double r1, double r2 ) : r1(r1), r2(r2) {}
	double r1, r2;
	RPoint inverse() const {
		RPoint res;
		res.r1 = r2;
		res.r2 = r1;
		return res;
	}
};

struct LessRPoint {
	bool operator()( const RPoint& x, const RPoint& y ) const {
		return ( x.r1 == y.r1 ? x.r2 < y.r2 : x.r1 < y.r1 );
	}
};

struct FitRange {
	int t1, t2;
};

struct DataPoint {
	double x, y, error;
};

struct ConstFit {
	FitRange range;
	double mean, error;
	double chi2dof;
};

struct FileDataCols{
	int nc;
	int tc;
	int d;
};

typedef std::vector<DataPoint> DataLine;
typedef std::map<RPoint,DataLine,LessRPoint> Data;



struct GSLfitRes{
	double val[3];
	double error[3];
	double errorJ[3];
	double chi2_dof;
	double chi2_dofJ;
};
	




#endif
