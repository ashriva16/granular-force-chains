/*
 * 		maintest.cpp
 *
 *  	Created on: October 27, 2015
 */
 
#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
#include "readInputFile.h"
#include "World2d.h"
#include "Wall2d.h"

int main(int argc, char * argv[]) {

   bool debug = 0;  
//   if (debug) {
  //  omp_set_dynamic(0);     // Explicitly disable dynamic teams
//}

    // MPI initialization 
    int numprocessors, rank; 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    char name[50]; int count;
    MPI_Get_processor_name(name, &count);

    size_t numInpFiles 	= 	9 ; // number of input files  
    size_t numIntParam	= 	5 ; // number of integer input parameters
    size_t numDecParam	= 	6 ; // number of decimal/double parameters
    vector<double>  parValues ;
    vector<size_t>  intParValues(numIntParam)  ; 
    string * inputFiles;
    inputFiles		= 	new string [numInpFiles] ;
    getParameters(argv[1], numInpFiles, numIntParam, numDecParam, 
				inputFiles, intParValues, parValues ) ; 
	  
    // Get morphology, init. position and init. velocity input files
    string file_morph = inputFiles[0];
    string file_pos = inputFiles[1];
    string file_vel = inputFiles[2];
    string morph_dir = inputFiles[3];
    
    // Generate grains
    vector<Grain2d> grains = generateGrainsFromFiles(file_morph, morph_dir, file_pos, file_vel);
    size_t ngrains = grains.size();

    // Change particle properties 
    double rho = 2.65;//3.125e-7;                  // g/pixel^3
    double mu = 0.5;
    double kn = 1e11;
    double ks = 1e11; 
    double cresN = 0.4; //coefficient of restitution 
    double cresS = 0.5;
    double maxVel = 200; //maximum random starting velocity in a given direction
    bool pluviate = intParValues[3]; //in pluviation mode or not?  
    size_t maxGrainIdInContact; 
  	
    size_t numRow     = intParValues[4]; //number of grains per row
    double dropHeight = parValues[4];    //what height to bring particles to before dropping them  
    size_t distBetweenRows = 50; //distance between rows  
    double verletRadius = 100;

    // Generate a square of flat walls (box)
    //      ------(2)-----
    //      |            |
    //      |            |
    //     (3)           (1)
    //      |            |
    //      |            |
    //      ------(0)-----

    // Infinite walls    
    double knWall = 1000000.;
    double muWall = 0.55;         
    double lX = parValues[5];      
    double lY = parValues[3];       

    vector<Wall2d> walls(4);

    // wall            normal               position           knWall   muWall  ID
    walls[0] = Wall2d(Vector2d(0.,1.),  Vector2d(0.,0.),       kn,  mu, 0);
    walls[1] = Wall2d(Vector2d(-1.,0.), Vector2d(lX,0.),       kn,  mu, 1);
    walls[2] = Wall2d(Vector2d(0.,-1.), Vector2d(lX,lY),       kn,  mu, 2);
    walls[3] = Wall2d(Vector2d(1.,0.),  Vector2d(0.,lY),       kn,  mu, 3);
    
    // set up file IO
    string posRotFile = inputFiles[4];
    string velocityFile = inputFiles[5];
    string cInfoDir = inputFiles[6];
    string shearHistFileInName = inputFiles[7];
    string shearHistFileOutName = inputFiles[8];
    

    FILE * positions        = fopen(posRotFile.c_str(),"w");
    FILE * velocities       = fopen(velocityFile.c_str(),"w");
    
    if (shearHistFileInName.substr(shearHistFileInName.length() - 4).compare("None") != 0 ){ //read shear histories if file name isn't 'None'
        readShearHistFile(shearHistFileInName, grains);
    }
  
    // Global parameters
   
    double A                           = parValues[0]; //0.0003; //amplitude of shake 
    double T                           = parValues[1]; //50000; //Time period of shake (in steps)       
    double forcez                      = parValues[2];
    double gDamping                     = 0.00;             // global damping
    double dt                           = 0.000001;           // time step 
    double dropLocation                 = grains[0].getPosition()[1];//height grains are dropped from
    static const unsigned int nTot     = intParValues[0];            // total time steps  
    const size_t nout                  = intParValues[1];  // get output data every nTout time steps 
    const size_t pout                  = 1000; //print data
    const size_t stopStep              = intParValues[2]; //time step to stop shaking     
    cout << "stopStep" << stopStep << endl;
    // Create world
    World2d world(grains, walls, dt, gDamping, verletRadius);

    // Clear up memory from grains since world now has a copy
    grains.clear();

    // Start clock
    double duration;
    double start = omp_get_wtime();

    double earthG = 980;                                // pixel/s^2
    Vector2d gravVec(0*earthG, -1.0*earthG);
    double omega = 2*M_PI/(T);
    double shift = 0;
    const size_t nUpdateNN = 1000; //how often to update the nearest neighbor list
    size_t startID, endID;
    // Start time integration
    for (size_t step = 0; step < nTot; ++step)
    {


	size_t distTraveled = (size_t) (0.5*earthG*pow(dt*step,2));
	if ((distTraveled % distBetweenRows == 0) && pluviate){ //update timestepped grains
		startID = (distTraveled/distBetweenRows)*numRow; endID = startID + numRow;
		if (endID <= ngrains){  
			world.moveRow(startID,endID,((distTraveled/distBetweenRows) + 1)*dropHeight,earthG,step); //move next row of grains
			world.changeMaxGrainIdInContact(endID);}
	}


	if (step % nUpdateNN == 0) //update nearest neighbor list  
		world.updateNeighborList(verletRadius);


        // Compute world state
        world.computeWorldState();

	//make el shaky shake 
	if (A != 0) shift = A*sin(omega*step);
	if (step < stopStep) {
		world.moveWall(0,Vector2d(0,shift)); //shake bottom wall
		world.wallTimestep(forcez, 2); //compress top wall
	}

        // Get output
        if (step % nout == 0)
        { 
            //cout << "step = " << step/nTout << endl;
            // Print positions and velocities to file
            for (size_t i = 0; i < ngrains; i++) { 
                fprintf(positions,  "%4.8f %4.8f %4.8f\n", world.getGrains()[i].getPosition()(0), world.getGrains()[i].getPosition()(1), world.getGrains()[i].getTheta());
                fprintf(velocities, "%4.8f %4.8f %4.8f\n", world.getGrains()[i].getVelocity()(0), world.getGrains()[i].getVelocity()(1), world.getGrains()[i].getOmega());
            }
        }

	if (step % pout == 0){
            duration = omp_get_wtime() - start;
            printf( "Timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", 
		int(step)+1, int(nTot), 100*double(step+1)/double(nTot), duration/60., duration/60.*nTot/(step+1) - duration/60.);
            fflush(stdout);
		
	    cout << "top wall pos" << world.getWalls()[2].getPosition() << endl;  
	    cout << world.getGrains()[10].getPosition() << endl;
	    world.outputStress();	
	cout << "shift: " << shift << endl;
cout << "dropLoc" << distTraveled << "and" << distBetweenRows << endl << endl;
        } // end output step

        // Apply gravity
        world.applyAcceleration(gravVec);

        // Take a time step
        world.grainTimestep();

    } // end time integration

    fclose(positions);
    fclose(velocities);
     
    // Print simulation time and save shear histories
    if (rank ==0) {
        duration = omp_get_wtime() - start;
        printf("Time taken: %dh, %dm, %2.2fs\n", int(floor(duration/3600.)), -int(floor(duration/3600.))+int(floor(fmod(duration/60., 60.))), -floor(fmod(duration/60., 60.))*60. + fmod(duration, 3600.) );
        
        FILE * shearhist   = fopen(shearHistFileOutName.c_str() , "w" ) ; 

        // shear history write out: 
        for (size_t g = 0 ; g < ngrains ; g++)
        {
            Grain2d grain = world.getGrains()[g] ;  
            size_t np = grain.getPointList().size() ; 

            vector<double> nodeSh = grain.getNodeShears() ; 
            vector<size_t> nodeContact = grain.getNodeContact() ; 

            double TOL = 1e-6 ; 

            for (size_t p = 0 ; p < np ; p++)
            {
                if (abs(nodeSh[p]) > TOL)
                {
                    fprintf(shearhist, "%ld %ld %ld %.16f\n", g , p ,  
                        nodeContact[p] , nodeSh[p]) ;
                }
            }
        }
        fclose(shearhist);
    
    }
   
    std::stringstream fname;
    fname << cInfoDir + "cinfo_" << 0 << ".dat";
    cout << "fname: " << fname.str() << endl;
    FILE * cinfo = fopen(fname.str().c_str(), "w");
       
    CData cState = world.computeCstate();
    for (size_t i = 0; i < cState._clocs.size(); i++) {
        fprintf(cinfo, "%d %d ",cState._cpairs[i](0), cState._cpairs[i](1) ); // grains
        fprintf(cinfo, "%.2f %.2f ",cState._forces[i](0), cState._forces[i](1)); // forces
        fprintf(cinfo, "%.3f %.3f ",cState._normals[i](0), cState._normals[i](1)); // normals
        fprintf(cinfo, "%.2f %.2f\n",cState._clocs[i](0), cState._clocs[i](1)); // locations
    }  

    fclose(cinfo);

 
    // Close MPI
    MPI_Finalize();

    cout << "Program Completed" << endl;

    return 0;
}
