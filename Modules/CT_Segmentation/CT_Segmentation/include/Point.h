/*
 * This file is part of libTIM.
 *
 * Copyright (©) 2005-2013  Benoit Naegel
 * Copyright (©) 2013 Theo de Carpentier
 *
 * libTIM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libTIM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/gpl>.
 */

#ifndef Point_h
#define Point_h
// ----------------------------------

#include <iostream>
#include "Types.h"

namespace LibTIM {

/** \defgroup Point Point 
	\ingroup DataStructures 
**/

/*@{*/

/// Point Structure
/** Basic structure to manipulate 3D points of type T
**/

template <class T>
struct Point
{
	T x, y, z;
	
	Point(TCoord x=0, TCoord y=0, TCoord z=0) {this->x=x; this->y=y; this->z=z;}
	Point& operator+=(Point <T> q) {x+=q.x; y+=q.y; z+=q.z; return *this;}
	Point& operator-=(Point <T> q) {x-=q.x; y-=q.y; z-=q.z; return *this;}
	bool operator==(Point <T> q) {return (x==q.x && y== q.y && z==q.z);}
	
	void operator=(const Point <T> &q) {
		this->x=q.x;
		this->y=q.y;
		this->z=q.z;
		}
	
	Point(const Point &q) {operator=(q);}
	void operator()(T x, T y, T z) {this->x=x; this->y=y; this->z=z; }
	
	void print() {std::cout << x << " " << y << " " << z <<"\n";}
};

template <class T>
Point<T> operator+(Point <T> p, Point <T> q)
{
	Point <T> r=p;
	return r+=q;
}

template <class T>
Point<T> operator-(Point <T> p, Point <T> q)
{
	Point <T> r=p;
	return r-=q;
}

/*@}*/

}

#endif
