/*
 * Levelset2d.h
 *
 *  Created on: May 31, 2014
 *      Author: Reid
 		Modified by: Konstantinos
 */

#ifndef LEVELSET2D_H_
#define LEVELSET2D_H_

#include "definitions.h"

class Levelset2d {

public:

	// default constructor
	Levelset2d() {
		_xdim = 0; _ydim = 0;
	}
	// actual constructor
	Levelset2d(const vector<float> & levelset, const size_t & xdim, const size_t & ydim):
		_levelset(levelset), _xdim(xdim), _ydim(ydim) {
		if (_xdim*_ydim != _levelset.size()) {
			cout << "ERROR: levelset size not consistent with dimensions" << endl;
		}
	}
	
	// checks if there is penetration, if there is, finds the penetration amount and the normalized gradient
	// and stores them in input fields value and gradient.  If not, returns false.
	bool isPenetration(const Vector2d & point, double & value, Vector2d & gradient) const {
		double x = point(0);
		double y = point(1);

		// check if the point exists in the level set, if not, return false
		if (x+1 > double(_xdim) || y+1 > double(_ydim) || x < 0 || y < 0){
			return false;
		}
		size_t xr = (size_t)round(x);
		size_t yr = (size_t)round(y);

		// check if the point is close to the surface, if not, return false
		if (getGridValue(xr,yr) > 1) {
			return false;
		}

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(xf, yf);
		double b2 	= getGridValue(xc, yf) - b1;
		double b3 	= getGridValue(xf, yc) - b1;
		double b4 	= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);

		value = b1 + b2*dx + b3*dy + b4*dx*dy;
		if (value > 0) {
			return false;
		}
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return true;
	}

	// Methods to deform the level sets
	void shearVertical(const double & gamma) {
		for (size_t j = 0; j < _ydim; j++) {
			for (size_t i = 0; i < _xdim; i++) {
//				_levelset[j*_xdim + i] += ((double)j - (double)_ydim/2.)*gamma;
				_levelset[j*_xdim + i] += ((double)j - 45.)*gamma;
			}
		}
	}
	
	void shearHorizontal(const double & gamma) {
		for (size_t j = 0; j < _ydim; j++) {
			for (size_t i = 0; i < _xdim; i++) {
//				_levelset[j*_xdim + i] += ((double)i - (double)_xdim/2.)*gamma;
				_levelset[j*_xdim + i] += ((double)i - 45.)*gamma;
			}
		}
	}
	
//  public get methods for debugging
//	double getXdim() const {
//		return _xdim;
//	}
//	double getYdim() const {
//		return _ydim;
//	}
//	vector<double> getLevelset() const {
//		return _levelset;
//	}
	
private:

	double getGridValue(size_t & x, size_t & y) const {
		return _levelset[y*_xdim + x];
	}
	
	vector<float> _levelset;
	size_t _xdim;
	size_t _ydim;
};

#endif /* LEVELSET2D_H_ */
