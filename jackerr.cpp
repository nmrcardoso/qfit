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


#include <numeric>
#include <vector>
#include <cmath>



double jackerr(std::vector<double> const & trials){
    double mean = std::accumulate(trials.begin(), trials.end(), 0.0) / trials.size();
    double sq_sum = std::inner_product(trials.begin(), trials.end(), trials.begin(), 0.0,
        [](double const & x, double const & y) { return x + y; },
        [mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
    return std::sqrt((trials.size() - 1) * sq_sum/trials.size());
}
double jackerr(std::vector<double> const & trials, const double mean){
    double sq_sum = std::inner_product(trials.begin(), trials.end(), trials.begin(), 0.0,
        [](double const & x, double const & y) { return x + y; },
        [mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
    return std::sqrt((trials.size() - 1) * sq_sum/trials.size());
}


/*
double jackerr(const vector<double>& trials){
    auto len = trials.size();
	double s1 = 0, s2 = 0;
	for( double x : trials ) {
		s1 += x;
		s2 += x*x;
	}
	double sumres = s2 - s1 * s1 / len;
    //return (len - 1) * std::sqrt(sumres) / len;
    return std::sqrt((len - 1) * sumres/len);
}
*/
