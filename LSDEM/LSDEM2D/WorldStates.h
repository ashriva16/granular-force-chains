/*
 * WorldStates.h
 *
 *  Created on: October 25, 2016
 *  Author: Konstantinos Karapiperis
 */

#ifndef WORLDSTATES_H_
#define WORLDSTATES_H_

struct GrainState2d {

	GrainState2d() {}
	GrainState2d(const Vector3d & svtemp, const vector<Vector2d> & gftemp, const vector<double> & gmtemp):
		_stressVoigt(svtemp), _grainForces(gftemp), _grainMoments(gmtemp) {}

	void resize(const int & ngrains) {
		_grainForces.resize(ngrains);
		_grainMoments.resize(ngrains);
	}

	void reset() {
		_stressVoigt << 0., 0., 0.;
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainForces[i] << 0., 0.;
			_grainMoments[i] = 0.;
		}
	}

	void operator+=(const GrainState2d & w) {
		_stressVoigt += w._stressVoigt; 
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainForces[i] += w._grainForces[i];
			_grainMoments[i] += w._grainMoments[i];
		}
	}

	Vector3d 			_stressVoigt;			// macroscopic stress of assembly
	vector<Vector2d> 	_grainForces;			// forces on grain
	vector<double>   	_grainMoments;			// moments on grain
};


struct WallState2d {

	WallState2d() {}
	WallState2d(const vector<Vector2d> & wtemp):
		_wallForces(wtemp) {}
	
	void reset() {
		for (size_t i = 0; i < _wallForces.size(); i++) {
			_wallForces[i] << 0., 0.;
			_wallContacts[i] = 0;
		}
	}

	void resize(const int & nwalls) {
		_wallForces.resize(nwalls);
		_wallContacts.resize(nwalls);
	}

	void resize(const int & nwalls, const int & nwallsFinite) {
		_wallForces.resize(nwalls);
		_wallContacts.resize(nwalls);
	}

	void operator+=(const WallState2d & w) {
		for (size_t i = 0; i < _wallForces.size(); i++) {
			_wallForces[i] += w._wallForces[i];
			_wallContacts[i] += w._wallContacts[i];
		}	
	}

	vector<size_t> 		_wallContacts; 			// number of contacts on infinite walls
	vector<Vector2d> 	_wallForces; 			// forces on infinite walls
};


struct CData {
	
	void operator+=(const CData & c) {
		_cpairs.insert( _cpairs.end(),	c._cpairs.begin(),	c._cpairs.end());
		_forces.insert( _forces.end(),	c._forces.begin(),	c._forces.end());
		_normals.insert(_normals.end(),	c._normals.begin(),	c._normals.end());
		_clocs.insert(_clocs.end(),		c._clocs.begin(),		c._clocs.end());
	}

	// maybe add more variables like branch vectors in the future
	vector<Vector2i> _cpairs;				// contact pairs
	vector<Vector2d> _forces;				// contact forces
	vector<Vector2d> _normals;				// contact normals
	vector<Vector2d> _clocs;				// location of contact points
};

#endif /* WORLDSTATES_H_ */
