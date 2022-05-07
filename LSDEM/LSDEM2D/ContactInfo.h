/*
 * ContactInfo.h
 *
 *  Created on: April 4th, 2016
 *      Author: Liuchi
 */

#ifndef CONTACTINFO_H_
#define CONTACTINFO_H_

// TODO: add other fields here
 #include "definitions.h"

struct ContactInfo {
	ContactInfo() {}
	ContactInfo( const vector<Vector2i> & cidtemp, const vector<Vector8d> & ctemp, const vector<Vector6d> & wallctemp, const vector<int> & wallcid):
	_grainPBFN(ctemp), _contactPair(cidtemp), _wallId(wallcid), _wallPFN(wallctemp){}
	// _grainPBFN stands for ContactPosition, BranchVector, ContactForce and ContactNormal respectively

	vector<Vector8d> _grainPBFN;
	vector<Vector2i> _contactPair;
    vector<Vector6d> _wallPFN;
    vector<int>      _wallId;
    
    	
	
	void resize() {

		_grainPBFN.resize(0);
		_contactPair.resize(0);
		_wallPFN.resize(0);
		_wallId.resize(0);
		   
	}
	
	void clear() {

		_grainPBFN.clear();
		_contactPair.clear();
		_wallPFN.clear();
		_wallId.clear();
		
	}
	
};



/*
for (unsigned int i = 0; i < trelax+trot, t++) {
	drumWorld.computeWorldState();
	if (t == trelax-1) {
		cinfo = drumWorld.computeContactInfo();
		for (unsigned int Infoindex  = 0; Infoindex < cinfo._grainPBFN.size(); Infoindex++) {
			fprintf(staticContacts, "%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %d %d\n", cinfo._grainPBFN[Infoindex](0),
				                                                                               cinfo._grainPBFN[Infoindex](1),
				                                                                               cinfo._grainPBFN[Infoindex](2),
				                                                                               cinfo._grainPBFN[Infoindex](3),
				                                                                               cinfo._grainPBFN[Infoindex](4),
				                                                                               cinfo._grainPBFN[Infoindex](5),
				                                                                               cinfo._grainPBFN[Infoindex](6),
				                                                                               cinfo._grainPBFN[Infoindex](7),
				                                                                               cinfo._contactPair[Infoindex](0),
				                                                                               cinfo._contactPair[Infoindex](1));

		}
		for (unsigned int windex = 0; windex < cinfo._wallPBFN.size(); windex++) {
			fprintf(staticWallContacts, "%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %d %d\n", cinfo._wallPBFN[windex](0),
				                                                                                   cinfo._wallPBFN[windex](1),
				                                                                                   cinfo._wallPBFN[windex](2),
				                                                                                   cinfo._wallPBFN[windex](3),
				                                                                                   cinfo._wallPBFN[windex](4),
				                                                                                   cinfo._wallPBFN[windex](5),
				                                                                                   cinfo._wallPBFN[windex](6),
				                                                                                   cinfo._wallPBFN[windex](7),
				                                                                                   cinfo._wallContactPair[windex](0),
				                                                                                   cinfo._wallContactPair[windex](1));
		}
		cinfo.clear();
		fclose(staticContacts);
		fclose(staticWallContacts);

		if (rank == 0) {
			for(unsigned int i = 0; i < drumWorld.getGrains().size(); i++) {
				fprintf(staticpositions, "%4.6f %4.6f %4.6f %d\n", drumWorld.getGrains()[i].getPosition()(0),
					                                               drumWorld.getGrains()[i].getPosition()(0),
					                                               drumWorld.getGrains()[i].getTheta(),
					                                               drumWorld.getGrainContacts()[i]);
				fprintf(staticvelocities, "%4.6f %4.6f %4.6f\n", drumWorld.getGrains()[i].getVelocity()(0),
					                                             drumWorld.getGrains()[i].getVelocity()(1),
					                                             drumWorld.getGrains()[i].getOmega());
			}
			fclose(staticpositions);
			fclose(staticvelocities);
		}
	} // end of static relaxation output

	if (t >= trelax && t%outfile == 0) {
		cinfo = drumWorld.computeContactInfo();
		sprintf(tempfame,"/home/lli/Projects/LSDEMDrum/%s/%d/Contacts_%s_%d_%d.dat", sr, round, sr, tstep, rank);
		Contacts = fopen(tempfname, "w");
		sprintf(tempfname, "/home/lli/Projects/LSDEMDrum/%s/%d/wallContacts_%s_%d_%d.dat", sr, round, sr, tstep, rank);
		wallContacts = fopen(tempfname, "w");
		for (unsigned int tindex = 0; tindex < cinfo._grainPBFN.size(); tindex++) {
			fprintf(Contacts, "%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %d %d\n", cinfo._grainPBFN[tindex](0),
				                                                                         cinfo._grainPBFN[tindex](1),
				                                                                         cinfo._grainPBFN[tindex](2),
				                                                                         cinfo._grainPBFN[tindex](3),
				                                                                         cinfo._grainPBFN[tindex](4),
				                                                                         cinfo._grainPBFN[tindex](5),
				                                                                         cinfo._grainPBFN[tindex](6),
				                                                                         cinfo._grainPBFN[tindex](7),
				                                                                         cinfo._contactPair[tindex](0),
				                                                                         cinfo._contactPair[tindex](1));
		}
		for (unsigned int wtindex = 0; wtindex < cinfo._wallPBFN.size(); wtindex++) {
			fprintf(wallContacts, "%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %4.6f %d %d\n", cinfo._wallPBFN[wtindex](0),
				                                                                             cinfo._wallPBFN[wtindex](1),
				                                                                             cinfo._wallPBFN[wtindex](2),
				                                                                             cinfo._wallPBFN[wtindex](3),
				                                                                             cinfo._wallPBFN[wtindex](4),
				                                                                             cinfo._wallPBFN[wtindex](5),
				                                                                             cinfo._wallPBFN[wtindex](6),
				                                                                             cinfo._wallPBFN[wtindex](7),
				                                                                             cinfo._wallContactPair[wtindex](0),
				                                                                             cinfo._wallContactPair[wtindex](1));
		}
        cinfo.clear();
        fclose(Contacts);
        fclose(wallContacts);
		if (rank == 0) {
			if (t%output == 0) {
				printf( "timestep %d of %d (%4.2f%% complete)\n", t, tmax, 100*double(t)/double(tmax));
				fflush(stdout);
                for (size_t index = 0; index < drumWorld.getGrains().size(); index++) {
                    ke += drumWorld.getGrains()[index].computeKineticEnergy(rotatingspeed);
                }
                cout << "Total kinetic energy is " << ke << endl;
			}
            sprintf(tempfname, "/home/lli/Projects/LSDEMDrum/%s/%d/Positions_%s_%d.dat", sr, round, sr, tstep);
            positions = fopen(tempfname,"w");
            sprintf(tempfname, "/home/lli/Projects/LSDEMDrum/%s/%d/Velocities_%s_%d.dat", sr, round, sr, tstep);
            velocities = fopen(tempfname,"w");
            sprintf(tempfname, "/home/lli/Projects/LSDEMDrum/%s/%d/WallPos_%s_%d.dat", sr, round, sr, tstep);
            wallPos = fopen(tempfname,"w");
            sprintf(tempfname, "/home/lli/Projects/LSDEMDrum/%s/%d/WallVel_%s_%d.dat", sr, round, sr, tstep);
            wallVel = fopen(tempfname,"w");                 
            sprintf(tempfname,"/home/lli/Projects/LSDEMDrum/%s/%d/GrainForces_%s_%d.dat", sr, round, sr, tstep);
            grainForces = fopen(tempfname,"w");
            for (size_t i = 0; i < drumWorld.getGrains().size(); i++) {
                fprintf(positions, "%4.6f %4.6f %4.6f %d\n", drumWorld.getGrains()[i].getPosition()(0), 
                                                             drumWorld.getGrains()[i].getPosition()(1), 
                                                             drumWorld.getGrains()[i].getTheta(),
                                                             drumWorld.getGrainContacts()[i]);
                fprintf(velocities, "%4.6f %4.6f %4.6f\n", drumWorld.getGrains()[i].getVelocity()(0), 
                                                           drumWorld.getGrains()[i].getVelocity()(1), 
                                                           drumWorld.getGrains()[i].getOmega());
                fprintf(grainForces, "%4.6f %4.6f %4.6f\n", drumWorld.getGrainForces()[i](0),
                                                            drumWorld.getGrainForces()[i](1),
                                                            drumWorld.getGrainForces()[i](2));
            }

            for (size_t index = 0; index < drumWorld.getWallGrains().size(); index++) {
                fprintf(wallPos, "%4.6f %4.6f %4.6f\n", drumWorld.getWallGrains()[index].getPosition()(0),
                                                        drumWorld.getWallGrains()[index].getPosition()(1),
                                                        drumWorld.getWallGrains()[index].getTheta());
                fprintf(wallVel, "%4.6f %4.6f %4.6f\n", drumWorld.getWallGrains()[index].getVelocity()(0),
                                                        drumWorld.getWallGrains()[index].getVelocity()(1),
                                                        drumWorld.getWallGrains()[index].getOmega());
            }
            fclose(positions);
            fclose(velocities);
            fclose(wallPos);
            fclose(wallVel);
            fclose(grainForces);
		}

	} // end of periodic output
	if (t >= trelax) {
		drumWorld.takeTimestepWall();
	}
	drumWorld.takeTimestepGrain();
}
*/
#endif /* CONTACTINFO_H_ */