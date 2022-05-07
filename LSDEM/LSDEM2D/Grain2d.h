/*
 * Grain2d.h
 *
 *  Created on: October 25, 2016
 *      Author: Konstantinos Karapiperis
 */

#ifndef GRAIN2D_H_
#define GRAIN2D_H_

#include "definitions.h"
#include "Levelset2d.h"
#include "WorldStates.h"

class Grain2d {
public:
	// constructors
	Grain2d() {
		_radius = 0; _mass = 0; _momentInertia = 0; _omega = 0; _theta = 0; _id = 0; _density = 1;
		_morphologyID = 0; _kn = 0; _ks = 0; _mu = 0; _ncontacts = 0;
	}

	Grain2d(const double & mass, const Vector2d & position, const Vector2d & velocity,
			  const double & momentInertia, const double & theta, const double & omega,
			  const Vector2d & cmLset, const vector<Vector2d> & pointList,
			  const double & radius, const Levelset2d & lset, const int & id,
			  const int & morphologyID, const double & kn, const double & ks,
			  const double & mu, const double & cresN, const double & cresS):
			  _mass(mass), _position(position), _velocity(velocity),  _momentInertia(momentInertia),
			  _theta(theta), _omega(omega), _cmLset(cmLset), _radius(radius), _lset(lset), _id(id),
			  _morphologyID(morphologyID), _kn(kn), _ks(ks), _mu(mu), _cresN(cresN), _cresS(cresS) {
		
		_pointList = pointList;
		_ncontacts = 0;
		_density = 1;
		_nodeShears.resize(_pointList.size());
		_nodeContact.resize(_pointList.size());
		// _nodeNormals.resize(_pointList.size());

		// Compute rotation matrix required for updating the pointlist from ref. to actual config. below
		Matrix2d rotMatrix;
		rotMatrix << cos(_theta), -sin(_theta), sin(_theta), cos(_theta);

		for (size_t i = 0; i < _pointList.size(); i++) {
			_nodeShears[i] = 0;
			_nodeContact[i] = 0;
			_pointList[i] = rotMatrix*_pointList[i] + _position;
			// _nodeNormals[i] << 0,0;
		}
	}

	// Checks if bounding circles intersect between this and other (1st level contact check)
	bool bcircleGrainIntersection(const Grain2d & other) const {
		if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
			 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		return false;
	}
	
	// Actual (2nd level) contact check between *this and other. If there is no contact, returns false.
	// Compares points of *this to the level set of other.
	// If there is contact, returns true and updates force, which is the force on *this,
	// thisMoment, which is the moment on *this, and otherMoment, the moment on other.
	bool findInterparticleForceMoment(const Grain2d & other, const double & dt, Vector2d & force, double & thisMoment, double & otherMoment) {

		// Declare temporary variables
	    force << 0., 0.;	    								 // force at a contact point (shared by the two grains in contact)
	    thisMoment = 0.; 										 // moment on this particle due to a contact
	    otherMoment = 0.;										 // moment on the other particle due to a contact
		Vector2d		   ptOtherCM; 		                     // point wrt the center of mass of other in real space
		Vector2d		   ptOtherLset; 	                     // point in the reference config of other's level set
		double		       penetration;	                         // penetration amount (is negative by convention of the level set)
		Vector2d		   normal; 			                     // surface normal pointing from other to *this
		const double       cos2 = cos(other.getTheta());
		const double       sin2 = sin(other.getTheta());
		Vector2d		   ptThisCM; 		                     // point wrt the center of mass of this in real space
		Vector2d		   df; 				                     // force increment from a single point
		Vector2d		   tangent; 		                     // surface tangent
		double		       sdot; 			                     // relative velocity of a point of contact
		double		       ds;									 // relative displacement of a point of contact
		Vector2d		   Fs; 				                     // vector of frictional/shear force
		Vector2d           relativeVel;							 // vector of relative velocity between the two grains
        
        double Kn = (_kn+other.getKn())/2.;						 // average normal contact stiffness between the particles in contact 			 
		double Ks = (_ks+other.getKs())/2.;						 // average shear contact stiffness between the particles in contact
		double Fsmag;											 // magnitude of shear force
		bool isContact = false;									 // return (assigned to true if contact exists)
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(_cresN)/sqrt(M_PI*M_PI + log(_cresN)*log(_cresN));
		double GamaS = -2*sqrt(_ks*_mass*other.getMass()/(_mass+other.getMass()))*log(_cresS)/sqrt(M_PI*M_PI + log(_cresS)*log(_cresS));

		// Do we need to include GamaS in the contact law?

		// Iterate through all of the points of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
		
			ptThisCM = _pointList[ptidx] - _position;
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			//	if there is penetration, finds forces and moments due to each contact
			if ( other.getLset().isPenetration(ptOtherLset, penetration, normal) ) {
				//cout << "penetration: " << penetration << endl; 
				isContact = true;
                _ncontacts++;
				// rotate the normal from the reference config of 2's level set to real space
				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				// find the tangent
				tangent << -normal(1), normal(0);
				// update force: normal force contribution
				relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
				               other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
				// viscoelastic contact law
				df = penetration*normal*Kn - GamaN*normal.dot(relativeVel)*normal;
				force -= df;
				// update moments: eccentric loading contribution
				otherMoment -= df(0)*ptOtherCM(1) - df(1)*ptOtherCM(0);
				thisMoment += df(0)*ptThisCM(1) - df(1)*ptThisCM(0);
				// force/moment calculations based on friction
				//sdot = tangent.dot(_velocity - other.getVelocity()) - (_omega*ptThisCM.norm() + other.getOmega()*ptOtherCM.norm());
				sdot = tangent.dot(-relativeVel);
				ds = sdot*dt;
				_nodeContact[ptidx] = other.getId();
				// Static friction law
				_nodeShears[ptidx] += _ks*ds + GamaS*sdot;	// elastic predictor
				if (_nodeShears[ptidx] > 0) {
					Fsmag = std::min(_nodeShears[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(_nodeShears[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				_nodeShears[ptidx] = Fsmag;
				force += Fs;
				thisMoment -= Fs(0)*ptThisCM(1) - Fs(1)*ptThisCM(0);
				otherMoment += Fs(0)*ptOtherCM(1) - Fs(1)*ptOtherCM(0);
			}
			// if there is no contact between the point and other, reset the shear force if other was the last contact
			else if (_nodeContact[ptidx] == other.getId() ){
				_nodeContact[ptidx] = 0;
				_nodeShears[ptidx] = 0;
				// _nodeNormals[ptidx] << 0,0;
			}
		} // end for loop iterating through contact points
		return isContact;
	} // end findInterparticleForceMoment


	// Finds contact data between *this and other - for output purposes
	CData findContactData(const Grain2d & other) const {

		// Declare temp variables
		CData 			   cData;								 // output
		Vector2d 		   force;								 // total force
		Vector2d		   ptOtherCM; 		                     // point wrt the center of mass of other in real space
		Vector2d		   ptOtherLset; 	                     // point in the reference config of other's level set
		double		       penetration;	                         // penetration amount (is negative by convention of the level set)
		Vector2d		   normal; 			                     // surface normal pointing from other to *this
		const double       cos2 = cos(other.getTheta());
		const double       sin2 = sin(other.getTheta());
		Vector2d		   ptThisCM; 		                     // point wrt the center of mass of this in real space
		Vector2d		   df; 				                     // force increment from a single point
		Vector2d		   tangent; 		                     // surface tangent
		double		       sdot; 			                     // relative velocity of a point of contact
		double		       ds;									 // relative displacement of a point of contact
		Vector2d		   Fs; 				                     // vector of frictional/shear force
		Vector2d           relativeVel;							 // vector of relative velocity between the two grains
        double Kn = (_kn+other.getKn())/2.;						 // average normal contact stiffness between the particles in contact 			 
		double Ks = (_ks+other.getKs())/2.;						 // average shear contact stiffness between the particles in contact
		double Fsmag;											 // magnitude of shear force
		//double cres = .6;										 // auxilliary variable for computing contact damping coeff below
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(_cresN)/sqrt(M_PI*M_PI + log(_cresN)*log(_cresN));
		double GamaS = -2*sqrt(_ks*_mass*other.getMass()/(_mass+other.getMass()))*log(_cresS)/sqrt(M_PI*M_PI + log(_cresS)*log(_cresS));
		// Initialize!
		force.fill(0.0);
		normal.fill(0.0);

		// iterate through all of the points of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {

			ptThisCM = _pointList[ptidx] - _position;
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();

			//	if there is penetration, get contact info into the cData struct 
			if ( other.getLset().isPenetration(ptOtherLset, penetration, normal) ) {

				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
				tangent << -normal(1), normal(0);
				relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
				               other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
				df = penetration*normal*Kn - GamaN*normal.dot(relativeVel)*normal;
				force -= df;
				if (_nodeShears[ptidx] > 0) {
					Fsmag = std::min(_nodeShears[ptidx], df.norm()*_mu );
				}
				else {
					Fsmag = std::max(_nodeShears[ptidx], -df.norm()*_mu );
				}
				Fs = -tangent*Fsmag;
				force += Fs;

				cData._cpairs.push_back(Vector2i(_id, other.getId()));
				cData._forces.push_back(force);
				cData._normals.push_back(normal);
				cData._clocs.push_back(_pointList[ptidx]);
			}
		} // end for loop iterating through contact points
		return cData;
	} // findContactData


	// Compute kinetic energy of the grain
	double computeKineticEnergy() const {
		double ke = 0;
		ke += .5*_mass*_velocity.squaredNorm();
		ke += .5*_momentInertia*_omega*_omega;
		return ke;
	}


	// Explicit time integration
	void takeTimestep(const Vector2d & force, const double & moment, const double & gDamping, const double & dt) {
		_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*force/_mass   );
		_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*moment/_momentInertia);
		double cosd = cos(_omega*dt);
		double sind = sin(_omega*dt);
		// must update the points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind, 
								(_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position + _velocity*dt;
		}
		_position = _position + dt*_velocity;
		_theta = _theta + dt*_omega;
	}
	
	void moveGrain(const Vector2d & amount) {
		_position = _position + amount;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += amount;
		}
	}
	
	// Change methods
	void changeMu(const double & newmu) {
		_mu = newmu;
	}

	void changePos(const Vector2d & pos) {
		Vector2d disp = pos - _position;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += disp;
		}
		_position = pos;	
	}

	void changeRot(const double & rot) {
		double dtheta = rot - _theta;
		_theta = rot;
		double cosd = cos(dtheta);
		double sind = sin(dtheta);
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind, 
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position;
		}
	}

	void changeKn(const double & kn) {
		_kn = kn;
	}

	void changeKs(const double & ks) {
		_ks = ks;
	}

	void changeCres(const double & cresN,const double & cresS){
		_cresN = cresN;
		_cresS = cresS; 
	}

	void changeDensity(const double & density) {
		_mass *= density/_density;
		_momentInertia *= density/_density;
		_density = density;
	}
    
	void changeId(const size_t & id) {
		_id = id;
	}

	void changeVelocity(const Vector2d & velocity){
		_velocity = velocity;
	}

    void changeShearHist(const size_t & nodeidx, const size_t & contacting, const double & shear) {
		_nodeShears[nodeidx] = shear;
		_nodeContact[nodeidx] = contacting;
	}
	//generate random velocity between -maxVel and maxVel in both coordinates 

	void setRandomVel(const double & maxVel) {
		double v1,v2;  
		v1 = 2*maxVel*( ((double)rand() / RAND_MAX) - 0.5);
		v2 = 2*maxVel*( ((double)rand() / RAND_MAX) - 0.5); 
		Vector2d newVel(v1,v2);
		_velocity = newVel;
	}

	// erase friction history
	void clearFriction() {
		for (size_t i = 0; i < _nodeShears.size(); i++) {
			_nodeShears[i] = 0.;
			// _nodeNormals[i] << 0,0,0;
			_nodeContact[i] = 0;
		}
	}

	// Helper methods
	const double & getMass() const {
		return _mass;
	}
	const Vector2d & getPosition() const {
		return _position;
	}
	const Vector2d & getVelocity() const {
		return _velocity;
	}
	const double & getTheta() const {
		return _theta;
	}
	const double & getOmega() const {
		return _omega;
	}
	const Vector2d & getCmLset() const {
		return _cmLset;
	}	

	const double & getRadius() const {
		return _radius;
	}
	const Levelset2d & getLset() const {
		return _lset;
	}	
	const double & getKn() const {
		return _kn;
	}
	const double & getKs() const {
		return _ks;
	}
	const double & getMu() const {
		return _mu;
	}
	const double & getCresN() const{
		return _cresN;
	}
	const double & getCresS() const{
		return _cresS;
	}
	const size_t & getId() const {
		return _id;
	}
    const size_t & getmorphologyID()const {
        return _morphologyID;
    }
    const vector<Vector2d> getPointList() const {
    	return _pointList;
    }
	const vector<size_t> & getNodeContact() const {
		return _nodeContact;
	}
	const vector<double> & getNodeShears() const {
		return _nodeShears;
	}
	// const vector<Vector2d> & getNodeNormals() const {
	// 	return _nodeNormals;
	// }
	 
	// non const
	vector<double> & getNodeShearsNonConst() {
		return _nodeShears;
	}
	vector<size_t> & getNodeContactNonConst() {
		return _nodeContact;
	}

private:

	double 		_mass;				// particle mass
	Vector2d 	_position; 			// location of center of mass in real space
	Vector2d 	_velocity;			// velocity of center of mass
	double 		_momentInertia;		// moment of inertia in principal frame (purely diagonal terms)
	double 		_theta;				// particle orientation
	double 		_omega;				// particle angular velocity
	Vector2d	_cmLset; 			// center of mass wrt the level set reference configuration (cm: at (0,0), I: diagonal)
	vector<Vector2d>  _pointList; 	// list of points comprising the grain in real space (translated and rotated)
	double 		_radius;			// radius of bounding sphere
	Levelset2d 	_lset;				// level set of grain
	double		_kn;				// normal stiffness
	double		_ks;				// shear stiffness
	double		_mu;				// interparticle friction
	double		_density;			// particle density (default 1?)
	double 		_cresN; 			//coefficient of restitution normal
	double 		_cresS;				//'' shear
	size_t      _morphologyID;		// ID representing the morphology type of the particle
	size_t 		_ncontacts;			// number of contacts of the grain (wals + other grains)
	size_t 		_id;				// ID (numbering) of the grain
	vector<double>	_nodeShears;	// shears at each node
	vector<size_t>	_nodeContact;  	// index of grain the node is contacting
	// vector<Vector2d> _nodeNormals;
};

#endif /* Grain2D_H_ */
