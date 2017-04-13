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

#ifndef FlatSE_h
#define FlatSE_h

#include <iostream>
#include <vector>
#include <limits>

#include "Image.h"
#include "Point.h"


///LibTIM library
/**
Philosophy is:
	\li reduced set of data structures (principal are: Image() FlatSE() NonFlatSE() )
	\li user oriented, easy to use
	\li large choice of mathematical morphology functions (recent algorithms, research code,...)
	\li efficient implementation 
**/

namespace LibTIM {

/** \defgroup FlatSE Flat Structuring Elements
	\ingroup DataStructures
**/
/*@{*/

///Container base class for flat structuring elements (or binary masks)
/**
	\par Example:
	\verbatim
	FlatSE se;
	se.make2DN9();
	\endverbatim
	creates a 2D structuring element containing a 3x3 square. The origin is at the center.
	\verbatim
	FlatSE se;
	se.makeBallEuclidian2D(r, im);
	\endverbatim
	creates 2D structuring element containing a ball of radius r according
	to the voxels spacing of im. The origin is at the center.
	
	WARNING: some algorithms require a  \e connexity rather than a structuring element in parameters.
	To this end, use for example make2DN8() to compute a 8-neighborhood (now the center is \e not included
	in the structuring element)
	Adding a \e connexity structure is in project.
**/

class FlatSE
{
	
	protected:
		//The flat SE is stored as a vector of points and offsets
		std::vector<Point<TCoord> > points;
		std::vector<TOffset> offsets;
		
		//Negative and positive offsets of SE 
		//Useful to add borders
		TCoord negativeOffsets[3];
		TCoord positiveOffsets[3];
		
	public:
	
	FlatSE(){}
	FlatSE(const Image<U8> &im);
	
	FlatSE &operator=(const FlatSE &se);
	FlatSE(const FlatSE &se)
		{
		operator=(se);
		}
	
	/// returns the number of points contained in the structuring element (cardinal of the set)
	int getNbPoints() const;

	/// computes the offset of each point, according to "size"
	void setContext(const TSize *size);

	void setNegPosOffsets();
	
	Point<TCoord> getPoint(int i) const {return points[i];}
	void addPoint(Point <TCoord> p) {points.push_back(p);}
	
	/// returns the offset of "point"
	TOffset getOffset(int point) {return offsets[point];}

	const TCoord* getNegativeOffsets() const;
	const TCoord* getPositiveOffsets() const;

	void makeSymmetric();

	Image<U8> toImage();
	
	FlatSE &operator+=(FlatSE &b){for (unsigned int i=0; i<b.points.size(); i++) points.push_back(b.points[i]); return *this;}
	
	//Iterators
	
	typedef std::vector<Point<TCoord> >::iterator iterator_point;
	typedef std::vector<TOffset >::iterator iterator_offset;
	typedef iterator_offset iterator;
	
	iterator begin() {return offsets.begin();}
	iterator end() {return offsets.end(); }
	
	iterator_point begin_point() {return points.begin();}
	iterator_point end_point() {return points.end(); }
	
	//Various neighborhoods
	
	/// Basic neighborhoods in 2D N4 and N8
	void make2DN4();
	void make2DN5();
	void make2DN8();
	void make2DN9();
	
	void make2DEuclidianBall(int r);
	
	///In 3D N6,18,26
	void make3DN6();
	void make3DN18();
	void make3DN26();
	void make3DN27();
	
	void make3DAxialSegment(int l);
	
	//Methods to create various structuring elements	
	
	template <class VoxelType> void makeBallEuclidian2D ( Image <VoxelType> &img, double r);
	template <class VoxelType> void makeBallChessboard2D(Image <VoxelType> &img, double rx, double ry);
	
	template <class VoxelType> void makeBallEuclidian3D(Image <VoxelType> &img, double r);
	
	template <class VoxelType> void makeCircle2D(Image <VoxelType> &img, double r, double t);
		
	void print() {
		std::cout << "FlatSE \n"; 
		for(unsigned int i=0; i< points.size(); i++) points[i].print(); 
		std::cout <<" "; std::cout << "\n";}

	void reserve(size_t size) {points.reserve(size); offsets.reserve(size);}
	void clear(){points.clear(); offsets.clear();}
	
};

/*@}*/

}

#include "FlatSE.hxx"

#endif
