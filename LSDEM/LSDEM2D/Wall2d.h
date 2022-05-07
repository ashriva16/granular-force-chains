// Wall2d.h

// Implementation of a flat wall (line) in 2D

// Author: Konstantinos Karapiperis
// Date: October 23, 2016

#ifndef WALL2D_H_
#define WALL2D_H_

#include "definitions.h"
#include "Grain2d.h"

class Wall2d {
	
public:
	Wall2d() {
		_d = 0; _kn = 0; _velocity = 0; _mu = 0; _id = INT_MAX-1;
	}
	Wall2d(const Vector2d & normal, const Vector2d & position, const double & kn, const double & mu, const size_t & id): 
		_normal(normal), _position(position), _kn(kn), _mu(mu), _id(id)  {
		_normal = _normal/_normal.norm();
		_d = -(normal.dot(position) );
		_velocity = 0;
	}
	
	bool bCircleContact(const Grain2d & grain) const {
		if	 ( _normal.dot(grain.getPosition()) + _d > grain.getRadius() ) {
			return false;
		}
		return true;
	}
	
	// frictionless wall
	bool findWallForceMoment(const Grain2d & grain, Vector2d & force, double & moment, size_t & ncontacts, Vector3d & stress) const {
		// zero the input vars
		force << 0., 0.;
		moment = 0.;
		stress << 0,0,0;
		ncontacts = 0;
		bool 	 checkflag = false;			// changes to true if penetration exists
		double 	 penetration; 				// penetration amount
		Vector2d df;						// force increment on grain
		Vector2d v;							// relative velocity
		Vector2d sdot;						// projection of relative velocity into tangential direction
		Vector2d Fs;						// shear force
		Vector2d ptcm;						// point wrt the center of mass of the grain in real space

		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				ncontacts++;
				ptcm = grain.getPointList()[ptidx] - grain.getPosition();
				ptcm = ptcm + penetration*ptcm/ptcm.norm();
				checkflag = true;
				df = -penetration*_normal*_kn;
				force += df;
				// moment -= ptcm.cross(df);
				moment += ptcm(0)*df(1) - ptcm(1)*df(0);
				
				stress(0) += df(0)*ptcm(0);
				stress(1) += df(1)*ptcm(1);
				stress(2) += 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
			}
		}
		return checkflag;
	}
	
	// wall with static friction
	bool findWallForceMomentFriction(Grain2d & grain, Vector2d & force, double & moment, size_t & ncontacts, Vector3d & stress, const double & dt) const {
		// zero the input vars
		force << 0., 0.;
		moment = 0.;
		stress << 0,0,0;
		ncontacts = 0;
		bool 	 checkflag = false;			// changes to true if penetration exists
		double 	 penetration; 				// penetration amount
		Vector2d df;						// force increment on grain
		Vector2d v;							// relative velocity
		Vector2d tangent; 		            // surface tangent
		double   ds;						// relative displacement of a point of contact
		double sdot;						// projection of relative velocity into tangential direction
		Vector2d Fs;						// shear force
		Vector2d ptcm;						// point wrt the center of mass of the grain in real space
		double Fsmag;						// magnitude of shear force

		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				ncontacts++;
				ptcm = grain.getPointList()[ptidx] - grain.getPosition();
				ptcm = ptcm + penetration*ptcm/ptcm.norm();
				checkflag = true;
				df = -penetration*_normal*_kn;
				force += df;
				// moment -= ptcm.cross(df);
				moment += ptcm(0)*df(1) - ptcm(1)*df(0);

				stress(0) += df(0)*ptcm(0);
				stress(1) += df(1)*ptcm(1);
				stress(2) += 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));

				// Friction contribution
				// Find relative tangential kinematics
				v << grain.getVelocity()(0) - grain.getOmega()*ptcm(1) - _velocity*_normal(0),
					 grain.getVelocity()(1) + grain.getOmega()*ptcm(0) - _velocity*_normal(1);
				tangent << -_normal(1), _normal(0);
				sdot = tangent.dot(-v);
				ds = sdot*dt;
				// TODO: Compare the above with the following to make sure they are equivalent
				// ds = (v - v.dot(_normal)*_normal)*dt; 
				
				grain.getNodeShearsNonConst()[ptidx] -= ds*_kn; // technically should be ks but walls don't have ks
				grain.getNodeContactNonConst()[ptidx] = _id;

				if (grain.getNodeShearsNonConst()[ptidx] > 0) {
					Fsmag = std::min(grain.getNodeShearsNonConst()[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(grain.getNodeShearsNonConst()[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				grain.getNodeShearsNonConst()[ptidx] = Fsmag;
				force += Fs;
			    moment += ptcm(0)*Fs(1) - ptcm(1)*Fs(0);
			}
			else if (grain.getNodeContactNonConst()[ptidx] == _id ){
				grain.getNodeContactNonConst()[ptidx] = INT_MAX;
				grain.getNodeShearsNonConst()[ptidx] = 0.0;
			}
		}
		return checkflag;
	}

	bool isContact(const Grain2d & grain) const {
		double 	penetration;
		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				return true;
			}
		}
		return false;
	}
	
	// methods to move/rotate the wall
	void moveWall(const Vector2d & amount) {
		_position += amount;
		_d = -(_normal.dot(_position) );
	}

	void rotateWall(const Matrix2d & R) {
		_normal = R*_normal;
		_d = -(_normal.dot(_position) );
	}
	
	// newton-raphson timestep in the normal direction
	void takeTimestep(const double & netForce, const size_t & ncontacts) {
		static const double alpha = 0.7;
		//cout <<"Wall Check " << netForce << "\n" << _normal << "amount" << alpha*netForce/_kn/5.*_normal << endl << endl; 
		if (ncontacts < 5)
			moveWall( alpha*netForce/_kn/5.*_normal );
		else
			moveWall( alpha*netForce/_kn/double(ncontacts)*_normal );
	}
	
	// change methods
	void changeKn(const double & newkn) {
		_kn = newkn;
	}
	
	void changeMu(const double & newmu) {
		_mu = newmu;
	}
	
	void changeVelocity(const double & velocity) {
		_velocity = velocity;
	}

	// get methods
	const Vector2d & getPosition() const {
		return _position;
	}
	const Vector2d & getNormal() const {
		return _normal;
	}
	
private:
	Vector2d _normal;		// normal direction, points away from wall
	Vector2d _position;		// 'center' of wall
	double 	_kn;			// wall stiffness
	double	_mu;			// wall friction
	double 	_velocity;		// velocity in the normal direction
	double 	_d;				// such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal
	size_t 	_id;
};

# endif
