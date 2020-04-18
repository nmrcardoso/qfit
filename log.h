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
#ifndef LOG_QFIT_H
#define LOG_QFIT_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>


class mylog {
	void create_logfile(std::string filename){
    	file.open(filename.c_str(), std::ios::out);
        if (!file.is_open()){
            std::cout << "Error opening file: " << filename << std::endl;
            exit(1);
        }
        std::cout << "Setting log file to: " << filename << std::endl;
    	writefile = true;
    }   
  public:
	std::ofstream file;
	bool writefile = false;
	
	// this is the type of std::cout
	typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
	// this is the function signature of std::endl
	typedef CoutType& (*StandardEndLine)(CoutType&);
	// define an operator<< to take in std::endl
	mylog& operator<<(StandardEndLine manip){
		// call the function, but we cannot return it's value
		manip(std::cout);
		//if( writefile ) file << std::endl;
		if( writefile ){
			file << "\n";
			file.flush();
		}
		return *this;
	}	
    mylog(std::string filename){ create_logfile( filename ); }  
    mylog(){ writefile = false; }
    ~mylog(){
    	if( writefile ) file.close();
    }
    void set(std::string filename){
    	if( writefile ) file.close();
    	create_logfile( filename );
    }      
};

template <class T>
mylog& operator<< (mylog& st, T val) {
  if( st.writefile ) st.file << val;
  std::cout << val;
  return st;
}

extern mylog qlog;


	
#endif

