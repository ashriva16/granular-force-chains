/*
 * World2d.h
 *
 *  Created on: October 25, 2016
 */

#ifndef WORLD2D_H_
#define WORLD2D_H_

#include "definitions.h"
#include "WorldStates.h"
#include "Grain2d.h"
#include "Wall2d.h"

class World2d {

public:
	World2d() {_dt = 0; _gDamping = 0.; _nwalls = 0;}

	// flat infinite walls + grains
	World2d(const vector<Grain2d> & grains, const vector<Wall2d> & walls, const double & dt, const double & gDamping, const double & verletRadius):
		_grains(grains), _walls(walls), _dt(dt), _gDamping(gDamping) {
		_ngrains = grains.size();
		_nwalls = walls.size();
		_maxGrainIdInContact = grains.size();
		_globalGrainState.resize(_ngrains);	_globalGrainState.reset();
		_globalWallState.resize(_nwalls);	_globalWallState.reset();
		constructNeighborList(verletRadius);
	}

	// computes the snapshot of the world (all grain forces/moments, all wall forces, total stress) by updating _globalGrainState and _globalWallState
	void computeWorldState() {
        
		// reset the global states
		_globalGrainState.reset();
		_globalWallState.reset();

		// define temp variables
		Vector2d force;					// force on a grain
		force.fill(0);					// initialization
		double momenti = 0.;			// moment on grain i
		double momentj = 0.;			// momnet on grain j
		Vector2d cmvec; 				// branch vector
		Vector3d stress;				// stress in the assembly
		Vector2d fn, fs;				// normal and shear forces
        
        // define parallel processing stuff (MPI is initialized at the mainfile)
		int numprocessors, rank;
	    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 		

	// TODO: Check if introducing the num_threads explicitly is faster or not
	#pragma omp parallel default(none) shared(numprocessors, rank, cout) firstprivate(force, momenti, momentj, cmvec, stress, fn, fs)// num_threads(12)
	{		
		// Initialize local states in each thread to later add all together
		GrainState2d threadGrainState;
		WallState2d threadWallState;
		threadGrainState.resize(_ngrains);
		threadGrainState.reset();
		threadWallState.resize(_nwalls);
		threadWallState.reset();
		size_t ncontacts = 0;

		// Each thread grubs 40 chunks of iterations dynamically
		#pragma omp for schedule(dynamic, 40)
		for (size_t i = rank; i < _maxGrainIdInContact; i+=numprocessors) {
			for (size_t k = 0; k < _neighbours[i].size(); k++) {
				size_t j = _neighbours[i][k];
				if (_grains[i].bcircleGrainIntersection(_grains[j])){   
				  if(_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momenti, momentj)) {
						cmvec = _grains[i].getPosition() - _grains[j].getPosition();
						threadGrainState._grainForces[i] += force;
						threadGrainState._grainForces[j] -= force;
						threadGrainState._grainMoments[i] += momenti;
						threadGrainState._grainMoments[j] += momentj;
						threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
						threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
						threadGrainState._stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
				    }
				}
			} // end loop over grain j
			// flat infinite wall-grain contact detection
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
	    } // end loop over grains
        // Assemble the global state
		#pragma omp critical
		{
			_globalWallState += threadWallState;
			_globalGrainState += threadGrainState;
		}   		    
	} // close opemp parallel section

	// MPI calls  sendbuf       recvbuff                                   		count       		type               op       comm
	MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainForces[0].data(),  		_ngrains*2,  		MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainMoments.data(),    		_ngrains,  			MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._stressVoigt.data(),     		3,           		MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallContacts.data(),     		_nwalls,     		MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallForces[0].data(),     	2*_nwalls,   		MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);

	// Finish computing the stress in the assembly by dividing over the area
	// We assume here that 4 walls surround all the particles - else better to compute the area in a postprocessing stage
	if (_walls.size() == 4) {
			_globalGrainState._stressVoigt /= (_walls[1].getPosition()(0) - _walls[0].getPosition()(0))*(_walls[2].getPosition()(1) - _walls[1].getPosition()(1));
		}
	} // end computeWorldState


	// Gets world state at the given snapshot for output purposes
	CData computeCstate() const {

		CData cDataRank;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		#pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
		{
			CData cDataThread;
			#pragma omp for schedule(dynamic, 5)
			for (size_t i = rank; i < _ngrains; i+=numprocessors) {
				for (size_t j = i+1; j < _ngrains; j++) {
					if (_grains[i].bcircleGrainIntersection(_grains[j])) { 
						CData cDataContact = _grains[i].findContactData(_grains[j]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // close grain subloop
				// Add here particle-wall contact data if implemented in the wall class in the future  
			} // end loop over grains
			#pragma omp critical
			{
				cDataRank += cDataThread;
			}
		} // closes openmp parallel section
		return cDataRank;
	} // end computeCstate method

	void applyBodyForce(Vector2d bodyForce) {
		for (size_t i = 0; i < _ngrains; i++) {
			_globalGrainState._grainForces[i] += bodyForce;
		}
	}
	
	// Need to prescribe negative acceleration for gravity
	void applyAcceleration(Vector2d acceleration) {
		for (size_t i = 0; i < _maxGrainIdInContact; i++) {
			_globalGrainState._grainForces[i] += _grains[i].getMass()*acceleration;
		}
	}

	// take a timestep for each grain based off of the world's _globalGrainState
	void grainTimestep() {

		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

		#pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
		for (size_t i = 0; i < _ngrains; i++) {
			_grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gDamping, _dt);
		}
	}

	//move an entire row of particles (from ID startID to endID) to have a yvalue of finalYPos
	//and velocity g*(finalYPos - origYPos)

	void moveRow(const size_t & startID, const size_t & endID, const double & dropHeight, const double & earthG, const size_t & step){
		for (size_t iD = startID; iD < endID; iD++){
			const Vector2d newPos(_grains[iD].getPosition()[0], dropHeight);
			Vector2d velocity(0 , -1.0*earthG*_dt*step);
			_grains[iD].changePos(newPos); 
	//		_grains[iD].changeVelocity(velocity);	
		}
	
	}

	void constructNeighborList(const double & verletRadius) {
		_neighbours.clear();
		for (size_t i = 0; i < _ngrains; i++) {
			vector<size_t> neighborsI;
			for (size_t j = i+1; j < _ngrains; j++) {
				Vector2d dist = _grains[j].getPosition() - _grains[i].getPosition();
				if (dist.norm() < verletRadius) {
				 	neighborsI.push_back(j);
				}
			}
			_neighbours.push_back(neighborsI);
		}
	}

	void updateNeighborList(const double & verletRadius) {
		for (size_t i = 0; i < _maxGrainIdInContact; i++) {
			vector<size_t> neighborsI;
			for (size_t j = i+1; j < _maxGrainIdInContact; j++) {
				Vector2d dist = _grains[j].getPosition() - _grains[i].getPosition();
				if (dist.norm() < verletRadius) {
				 	neighborsI.push_back(j);
				}
			}
			_neighbours[i] = neighborsI;
		}
	}

	// move grain
	void moveGrain(const size_t & grainid, const Vector2d & amount) {
		_grains[grainid].moveGrain(amount);
	}

	// move infinite wall
	void moveWall(const size_t & wallid, const Vector2d & amount) {
		_walls[wallid].moveWall(amount);
	}

	//timestep wall
	void wallTimestep(const double & force, const size_t & wallId) {
                _walls[wallId].takeTimestep(force - _globalWallState._wallForces[wallId].norm(), _globalWallState._wallContacts[wallId]);
        }


	// change friction for infinite walls
	void changeMuWall(const size_t & wallid, const double & newmu) {
		_walls[wallid].changeMu(newmu);
	}

	// change friction for particles
	void changeMuGrains(const double & newmu) {
		for (size_t i = 0; i < _ngrains; ++i) {
			_grains[i].changeMu(newmu);
		}
	}
	void changeMaxGrainIdInContact(const size_t & maxGrainId) {
		_maxGrainIdInContact = maxGrainId;
	}
	// take a timestep for more advanced walls in the future
	// void takeTimestepWall () {
	// }

	void changeDt(const double & dt) {
		_dt = dt;
	}
	void changeGdamping(const double & gDamping) {
		_gDamping = gDamping;
	}

    // const get methods
	const vector<Grain2d> & getGrains() const {
		return _grains;
	}
	const vector<Wall2d> & getWalls() const {
		return _walls;
	}
	const GrainState2d & getGrainState() const {
		return _globalGrainState;
	}
	const WallState2d & getWallState() const {
		return _globalWallState;
	}

	//	 non const get methods (members can be changed such as wall positions)
	// IMPORTANT FOR WALL MOVEMENT in walls with friction SEE World3d.h


	const vector<Vector2d> & getGrainForces() const {
		return _globalGrainState._grainForces;
	}
	const vector<double> & getGrainMoments() const {
		return _globalGrainState._grainMoments;
	}
	// const vector<double> & getGrainForceNorms() const {
	// 	return _globalGrainState._grainForceNorms;
	// }
	// const vector<int> & getGrainContacts() const{
	// 	return _globalGrainState._grainContacts;
	// }

	void outputStress() {
		cout << "sigxx = " << _globalGrainState._stressVoigt(0) << endl;
		cout << "sigyy = " << _globalGrainState._stressVoigt(1) << endl;
		cout << "sigxy = " << _globalGrainState._stressVoigt(2) << endl;
	}

private:
    vector<Grain2d>   	 _grains;				// vector of grain objects
	vector<Wall2d> 		 _walls;			 	// vector of infinite wall objects
    double 				 _dt;					// time increment
    double 			 	 _gDamping;				// global damping
	GrainState2d		 _globalGrainState;		// grain state of entire assembly
	WallState2d			 _globalWallState;		// wall state of entire assembly
	size_t 				 _ngrains;				// number of grains
	size_t				 _maxGrainIdInContact; // max id of grains for which contact check is performed (useful when a lot of grains are in ballistic mode) 
	size_t 				 _nwalls;				// number of infinite walls
	vector<vector< size_t > >  	_neighbours; 		// neighbors within Verlet radius
};

#endif /* World2d_H_ */
