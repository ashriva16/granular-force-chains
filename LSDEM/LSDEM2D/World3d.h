/*
 * World3d.h
 *
 *  Created on: September 9, 2014
 *      Author: Reid
 */

#ifndef WORLD3D_H_
#define WORLD3D_H_

#include "definitions.h"
#include "WorldStates.h"
#include "Grain3d.h"
#include "Wall3d.h"
#include "WallFrustum3d.h"

// comment!
class World3d {
	
public:
	World3d() {_dt = 0; _gdamping = 0; _ngrains = 0; _maxGrainIdInContact = 0; _nwalls = 0;}
	
	// regular walls only
	World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const double & dt, const double & gdamping, const double & verletRadius):
		_grains(grains), _walls(walls), _dt(dt), _gdamping(gdamping) {
		_ngrains = grains.size();
		_maxGrainIdInContact = grains.size();
		_nwalls = walls.size();
		_nwallsFrustum = 0;
		_globalGrainState.resize(_ngrains);	_globalGrainState.reset(); 
		_globalWallState.resize(_nwalls);	_globalWallState.reset();
		constructNeighborList(verletRadius);
	}

	//Wall frustrums
	 World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const vector<WallFrustum3d> & wallsFrustum,  const double & dt, const double & gdamping, const double & verletRadius):
                _grains(grains), _walls(walls), _wallsFrustum(wallsFrustum), _dt(dt), _gdamping(gdamping) {
    	_ngrains = grains.size();
        _maxGrainIdInContact = grains.size();
		_nwalls = walls.size();
		_nwallsFrustum = wallsFrustum.size();	
        _globalGrainState.resize(_ngrains); _globalGrainState.reset();
        _globalWallState.resize(_nwalls); _globalWallState.reset();
		constructNeighborList(verletRadius);
  	}

	// updates _globalGrainState and _globalWallState
	void computeWorldState() {
		// zero out the global state and define temp variables
		_globalGrainState.reset();
		_globalWallState.reset();
		Vector3d force;
		Vector3d momenti;
		Vector3d momentj;
		Vector3d cmvec; 
		Vector6d stress;
		Vector3d fn, fs;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

		#pragma omp parallel default(none) shared(numprocessors,rank,cout) firstprivate(force, momenti, momentj, cmvec, stress, fn, fs) // reduction(+:threadWallState,threadGrainState) //num_threads(1)
		{
			GrainState3d threadGrainState;
			WallState3d threadWallState;
//			for (size_t repeat = 0; repeat < 1000; repeat++) {
			threadGrainState.resize(_ngrains);
			threadGrainState.reset();
			threadWallState.resize(_nwalls);
			threadWallState.reset();
			size_t ncontacts = 0;
			// go through the grains and compute grain state
			
			#pragma omp for schedule(dynamic, 5) 
			for (size_t i = rank; i < _maxGrainIdInContact; i+=numprocessors) {
				// grain-grain contacts
				for (size_t k = 0; k < _neighbours[i].size(); k++) {
					size_t j = _neighbours[i][k];
					if (_grains[i].bCircleCheck(_grains[j])) { 
						if (_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momentj, momenti)) {
							cmvec = _grains[i].getPosition() - _grains[j].getPosition();
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainForces[j] -= force;
							threadGrainState._grainMoments[i] += momentj;
							threadGrainState._grainMoments[j] += momenti;
							threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
							threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
							threadGrainState._stressVoigt(2) += force(2)*cmvec(2);
							threadGrainState._stressVoigt(3) += 0.5*(force(1)*cmvec(2) + force(2)*cmvec(1));
							threadGrainState._stressVoigt(4) += 0.5*(force(2)*cmvec(0) + force(0)*cmvec(2));
							threadGrainState._stressVoigt(5) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
						}
					}
				}
				// flat walls
				for (size_t j = 0; j < _nwalls; j++) {
					if (_walls[j].bCircleContact(_grains[i])) {
						if (_walls[j].findWallForceMomentFriction(_grains[i], force, momenti, ncontacts, stress, _dt) ) {
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainMoments[i] += momenti;
							threadWallState._wallForces[j] -= force;
							threadWallState._wallContacts[j] += ncontacts;
							threadGrainState._stressVoigt += stress;
						}
					}
				} // end loop over wall j
				// wall frustums
				for (size_t j = 0; j < _nwallsFrustum; j++) {
					if (_wallsFrustum[j].bCircleContact(_grains[i])) {
						if (_wallsFrustum[j].findWallForceMoment(_grains[i], force, momenti, ncontacts, stress, _dt) ) {
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainMoments[i] += momenti;
							threadGrainState._stressVoigt += stress;  
						}
					}
				}
			} // end loop over grains
			
			#pragma omp critical 
			{
				_globalWallState += threadWallState; // FIXME
				_globalGrainState += threadGrainState; // FIXME
			} // FIXME
		
		} // closes openmp parallel section
//		cout << "time to compute world state: rank: " << rank << ", time " << omp_get_wtime() - s2 << endl;
		// MPI calls  sendbuf       recvbuff                                   count       				type               op       comm
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainForces[0].data(),  				_ngrains*3,  			MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainMoments[0].data(), 				_ngrains*3,  			MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._stressVoigt.data(),     				6,           			MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallContacts.data(),     				_nwalls,     			MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallForces[0].data(),    				3*_nwalls,   			MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		
		_globalGrainState._stressVoigt /= findVolume();

	} // end computeWorldState method
	
	CData computeCstate() const {
		CData cDataRank;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank==0) {
//			cout << "1" << endl;
		}
		#pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
		{
			CData cDataThread;
			#pragma omp for schedule(dynamic, 5)
			for (size_t i = rank; i < _ngrains; i+=numprocessors) {
				for (size_t j = i+1; j < _ngrains; j++) {
					if (_grains[i].bCircleCheck(_grains[j])) { 
						CData cDataContact = _grains[i].findContactData(_grains[j]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // close grain subloop
			} // end loop over grains
			#pragma omp critical
			{
				cDataRank += cDataThread;
			}
		} // closes openmp parallel section
		return cDataRank;
	} // end computeCState method
	
	// Same but for contact with walls
	CData computeWallCstate() const {
		CData cDataRank;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		#pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
		{
			CData cDataThread;
			#pragma omp for schedule(dynamic, 5)
			for (size_t i = rank; i < _ngrains; i+=numprocessors) {
				for (size_t j = 0; j < _nwalls; j++) {
					if (_walls[j].bCircleContact(_grains[i])) {
						CData cDataContact = _walls[j].findContactData(_grains[i]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // end loop over wall j
			} // end loop over grains
			#pragma omp critical
			{
				cDataRank += cDataThread;
			}
		} // closes openmp parallel section
		return cDataRank;
	} // end computeWallCState method

	void applyBodyForce(Vector3d bodyForce) {
		for (size_t i = 0; i < _ngrains; i++) {
			_globalGrainState._grainForces[i] += bodyForce;
		}
	}
	
	void applyAcceleration(Vector3d acceleration) {
		for (size_t i = 0; i < _ngrains; i++) {
			_globalGrainState._grainForces[i] += _grains[i].getMass()*acceleration;
		}
	}
	 
	// take a timestep for each grain based off of the world's _globalGrainState
	void grainTimestep() {
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
//		double s2 = omp_get_wtime();
//		for (size_t repeat = 0; repeat < 1000; repeat++) {
		// update center of mass, rotation, etc
		#pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
		for (size_t i = 0; i < _ngrains; i++) {
			_grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gdamping, _dt);
//			_grains[i].updatePoints();
		}
//		 update points for grains on this rank
		#pragma omp parallel for default(none) shared(rank,numprocessors) schedule(static,1) //num_threads(12)
		for (size_t i = rank; i < _ngrains; i+=numprocessors) {
			_grains[i].updatePoints();
		}
	}
	
	// take a time step for part of the grains, while the rest remain at rest
	void grainTimestepPartial(const vector<size_t> activeGrainIndices) {
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

		// update center of mass
		#pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
		for (size_t i = 0; i < activeGrainIndices.size(); i++) {
			size_t idx = activeGrainIndices[i];
			_grains[idx].takeTimestep(_globalGrainState._grainForces[idx], _globalGrainState._grainMoments[idx], _gdamping, _dt);
		}
//		 update points for grains on this rank
		#pragma omp parallel for default(none) shared(rank,numprocessors) schedule(static,1) //num_threads(12)
		for (size_t i = rank; i < _ngrains; i+=numprocessors) {
			_grains[i].updatePoints();
		}
	}

	void constructNeighborList(const double & verletRadius) {
		_neighbours.clear();
		for (size_t i = 0; i < _ngrains; i++) {
			vector<size_t> neighborsI;
			for (size_t j = i+1; j < _ngrains; j++) {
				Vector3d dist = _grains[j].getPosition() - _grains[i].getPosition();
				if (dist.norm() < verletRadius) {
				 	neighborsI.push_back(j);
				}
			}
			_neighbours.push_back(neighborsI);
		}
	}

	void updateNeighborList(const double & verletRadius) {
		for (size_t i = 0; i < _ngrains; i++) {
			vector<size_t> neighborsI;
			for (size_t j = i+1; j < _ngrains; j++) {
				Vector3d dist = _grains[j].getPosition() - _grains[i].getPosition();
				if (dist.norm() < verletRadius) {
				 	neighborsI.push_back(j);
				}
			}
			_neighbours[i] = neighborsI;
		}
	}

	void dragTopWall() {
		Vector3d force;
		Vector3d momenti;
		size_t ncontacts=0;
		Vector6d stress;
		Vector3d threadVelocity(0,0,0);
		Vector3d totalVelocity(0,0,0);
		size_t threadGcount = 0;
		size_t totalGcount = 0;
		#pragma omp parallel default(none) firstprivate(force,momenti,ncontacts,stress,threadVelocity,totalVelocity,threadGcount,totalGcount)
		{
			#pragma omp for schedule(dynamic,5)
			for (size_t i = 0; i < _ngrains; i+=1) {
				if (_walls[1].bCircleContact(_grains[i])) {
					if (_walls[1].findWallForceMomentFriction(_grains[i], force, momenti, ncontacts, stress, _dt) ) {
						threadVelocity += _grains[i].getVelocity();
						threadGcount++;
					}
				}
			}
		}
		#pragma omp critical
		{
			totalVelocity += threadVelocity;
			totalGcount += threadGcount;
		}
		double fac = 1.0;
		Vector3d projectedVelocity = totalVelocity - _walls[1].getNormal()*totalVelocity.dot(_walls[1].getNormal());
		_walls[1].changeVelocity( projectedVelocity*fac/(double)totalGcount );
	}
	
	void wallTimestep(const double & force, const size_t & wallId) {
		_walls[wallId].takeTimestep(force - _globalWallState._wallForces[wallId].norm(), _globalWallState._wallContacts[wallId]);
	}

	void topWallMoveWithGrains(const double & amount, const double & height) {
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
		_walls[1].moveWall(Vector3d(0,0,amount)); 
		#pragma omp parallel for schedule(static,1) // num_threads(2)
		for (size_t i = rank; i < _ngrains; i+=numprocessors) {
			_grains[i].moveGrain(Vector3d(0,0,amount*_grains[i].getPosition()(2)/height ) );
		}
	}
	
	double findVolume() {
		double volume = 1;
		// 6 walls
		if (_walls.size() == 6) {
			volume = (_walls[1].getPosition()(2) - _walls[0].getPosition()(2))*
						(_walls[3].getPosition()(1) - _walls[2].getPosition()(1))*
						(_walls[5].getPosition()(0) - _walls[4].getPosition()(0));
		}
		else if (_nwallsFrustum == 1) {
			volume = _wallsFrustum[0].getVolume();
		}
		else if (_nwallsFrustum == 2) {
			volume = _wallsFrustum[1].getVolume();
		}
		// else if (_walls.size() == 1) {
		// 	double minX = DBL_MAX;
		// 	double minY = DBL_MAX;
		// 	double minZ = DBL_MAX;
		// 	double maxX = -DBL_MAX;
		// 	double maxY = -DBL_MAX;
		// 	double maxZ = -DBL_MAX;
		// 	for (size_t i = 0; i < _grains.size(); i++) {
		// 		maxX = max(maxX,_grains[i].getPosition()(0));
		// 		maxY = max(maxX,_grains[i].getPosition()(1));
		// 		maxZ = max(maxX,_grains[i].getPosition()(2));
		// 		minX = min(minX,_grains[i].getPosition()(0));
		// 		minY = min(minY,_grains[i].getPosition()(1));
		// 		minZ = min(minZ,_grains[i].getPosition()(2));
		// 	}
		// 	volume = (maxX-minX)*(maxY-minY)*(maxZ-minZ);
		// }
		return volume;
	}
	
	void changeDt(const double & dt) {
		_dt = dt;
	}
	void changeGdamping(const double & gdamping) {
		_gdamping = gdamping;
	}
	// changes kn of grains, walls, membrane
	void changeKnAll(const double & kn) {
		for (size_t i = 0; i < _grains.size(); i++) {
			_grains[i].changeKn(kn);
		}
		for (size_t i = 0; i < _walls.size(); i++) {
			_walls[i].changeKn(kn);
		}
		for (size_t i = 0; i < _wallsFrustum.size(); i++) {
			_wallsFrustum[i].changeKn(kn);
		}
	}
	// changes ks of grains, walls
	void changeKsAll(const double & ks) {
		for (size_t i = 0; i < _grains.size(); i++) {
			_grains[i].changeKs(ks);
		}
		for (size_t i = 0; i < _walls.size(); i++) {
			_walls[i].changeKs(ks);
		}
	}
	void changeMaxGrainIdInContact(const size_t & maxGrainId) {
		_maxGrainIdInContact = maxGrainId;
	}
	// const get methods
	const double getWallPenetration(const size_t & wallId) const {
		double penetration = 0; // cumulative penetration amount
		for (size_t i = 0; i < _grains.size(); i++) {
			penetration += _walls[wallId].getMaxPenetration(_grains[i]);
		}
		return penetration;
	}
	const vector<Grain3d> & getGrains() const {
		return _grains;
	}
	const vector<Wall3d> & getWalls() const {
		return _walls;
	}
	const GrainState3d & getGrainState() const {
		return _globalGrainState;
	}
	const WallState3d & getWallState() const {
		return _globalWallState;
	}
		
//	 non const get methods (members can be changed such as wall positions)
	vector<Wall3d> & getWallsNonConst() {
		return _walls;
	}
	vector<WallFrustum3d> & getWallsFrustumNonConst() {
		return _wallsFrustum;
	}
	vector<Grain3d> & getGrainsNonConst() {
		return _grains;
	}

	void outputContacts() {
		// Initialize
		size_t nContactsGrainGrain = 0;
		for (size_t i = 0; i < _ngrains; i++){
			for (size_t j = i+1; j < _ngrains; j++){
				// Normal contacts
				if (_grains[i].bCircleCheck(_grains[j])){
				   if(_grains[i].checkPenetration(_grains[j])) {
				   		nContactsGrainGrain += 1;
				    }
				}			
			} // end loop over grain j
	    }
		printf("Number of grain-grain contacts = %lu\n",nContactsGrainGrain);
	}

	double getKineticEnergy() {
        double ke = 0;
        for (size_t i = 0; i < _ngrains; i++) {
            ke += _grains[i].computeKineticEnergy();
        }
        return ke;
	}

	void outputKineticEnergy() {
        double ke = 0;
        for (size_t i = 0; i < _ngrains; i++) {
            ke += _grains[i].computeKineticEnergy();
        }
        printf("Kinetic energy = %.4f\n", ke);
	}

	void outputStress() {
		cout << "Voigt stress:" << endl;
		cout <<  _globalGrainState._stressVoigt << endl;
	}

private:
	vector<Grain3d>				_grains;
	vector<Wall3d>				_walls;
	vector<WallFrustum3d> 		_wallsFrustum;
	GrainState3d				_globalGrainState;	// grain state of entire assembly
	WallState3d					_globalWallState;		// wall state of entire assembly
	double 						_dt;			 // time increment
	double 						_gdamping; 		// global damping
	size_t 						_ngrains;		// number of grains
	size_t						_maxGrainIdInContact; // max id of grains for which contact check is performed (useful when a lot of grains are in ballistic mode) 
	size_t 						_nwalls;		// number of walls
	size_t 						_nwallsFrustum;
	vector<vector<size_t>>  	_neighbours; 		// neighbors within Verlet radius
};


#endif
