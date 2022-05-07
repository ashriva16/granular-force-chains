/*
 * readInputFile.h
 *
 *  Created on: Jun 3, 2014
 *      Modified by: Liuchi
 */

#ifndef READINPUTFILE_H_
#define READINPUTFILE_H_

#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"

// creates a vector of grain objects from input files
vector<Grain2d> generateGrainsFromFiles(string morphologyMaterialIDs,
										string morphologyDir,
	                                   	string Positions,
	                                   	string Velocities) {

	string line;
	string line_property;
	string line_velocity;
	string line_position;
	string partial;
	string propertyfile;
	istringstream iss;
	ifstream file_information(morphologyMaterialIDs.c_str());		// construct an ifstream and open the grain morphology file
	ifstream file_position(Positions.c_str());						// construct another ifstream and open the grain position file
	ifstream file_velocity(Velocities.c_str());						// construct another ifstream and open the grain velocity file
	
	getline(file_information, line);								// read a line (first) of the morphology file (number of particles)
    iss.str(line);													// copy line string to into iss (we basically bind iss to the line we just read)
    getline(iss, partial, ' ');										// extracts characters until delimiter ' ' is found; the latter is extracted and discarded 
	size_t numberOfGrains = atoi(line.c_str());						// converts the string to an integer
	iss.clear();													// clear the error state of the stream (e.g. end-of-file -> no error)
    char tempfilename[100];
    // Initialize the vector of grain objects
	vector<Grain2d> grainList(numberOfGrains);
	// temp stuff
	Vector2d point;
	Vector2d position;
	double theta;
	Vector2d velocity;
	double omega;

	// Go through each grain 
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {

        // Read morphology index for each particle - each index has each own property .dat file
        getline(file_information, line);
        iss.str(line);
        getline(iss, partial, ' ');
		int morphologyID = atoi(partial.c_str());
        iss.clear();
        sprintf(tempfilename, "grainproperty%d.dat", morphologyID);
        propertyfile = morphologyDir + tempfilename;
        ifstream file_gp(propertyfile.c_str());

        // mass
        getline(file_gp, line_property);
        double mass = atof(line_property.c_str());
	
        // moment of inertia
        getline(file_gp, line_property);
		double momentOfInertia = atof(line_property.c_str());

		// cmLset (center of mass)
		getline(file_gp, line_property);
		Vector2d cmLset;
		iss.str(line_property);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		iss.clear();

		// number of points on the grain surface (INTEGER)
		getline(file_gp, line_property);
		int npoints = atoi(line_property.c_str());

		// the point coordinates
		getline(file_gp, line_property);
		vector<Vector2d> pointList(npoints);
		iss.str(line_property);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// bounding box radius
		getline(file_gp, line_property);
		double bboxRadius = atof(line_property.c_str());

		// level set dimensions (INTEGERS)
		getline(file_gp, line_property);
		iss.str(line_property);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file_gp, line_property);
		vector<float> lsetvec(xdim*ydim);
		iss.str(line_property);
		for (int i = 0; i < xdim*ydim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// kn
		getline(file_gp, line_property);
		double kn = atof(line_property.c_str());

        // ks
		getline(file_gp, line_property);
		double ks = atof(line_property.c_str());

        // mu
		getline(file_gp, line_property);
		double mu = atof(line_property.c_str());

	//cresN 
		getline(file_gp, line);
		double cresN = atof(line.c_str());
	//cresS 
		getline(file_gp, line);
		double cresS = atof(line.c_str());
	
        // Clear string for property file ready for next grain - do we need this?
        propertyfile.clear();
        line_property.clear();

	    // Read position and theta from position file
	    getline(file_position, line_position);
		iss.str(line_position);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		theta = atof(partial.c_str());
		iss.clear();

		// Read velocity and omega from velocity file
		getline(file_velocity, line_velocity);
		iss.str(line_velocity);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		omega = atof(partial.c_str());
		iss.clear();
	   	
		// Create level set object  - should we we use new?
		Levelset2d lset(lsetvec, xdim, ydim);

		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, velocity, momentOfInertia, theta, omega, cmLset, pointList, bboxRadius, lset, grainidx, morphologyID, kn, ks, mu, cresN, cresS);
	}

	return grainList;
} // end generateGrainsFromFiles


// creates a vector of grain objects from a single input file (e.g. from level set imaging)
vector<Grain2d> generateGrainsFromFile(string filename) {

	string line;
	string partial;
	istringstream iss;
	ifstream file(filename.c_str());								

	getline(file, line);
	size_t numberOfGrains = atoi(line.c_str());

    // Initialize the vector of grain objects
	vector<Grain2d> grainList(numberOfGrains);

	// temp stuff
	Vector2d point;
	Vector2d position;
	Vector2d velocity;

	// morphologyID is of no importance here since each grain is special; just pass any ID to the constructor
	size_t morphologyID = 0;

	// Go through each grain 
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {

        // mass
        getline(file, line);
        double mass = atof(line.c_str());

		// position
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		iss.clear();

		// velocity
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		iss.clear();

        // moment of inertia
        getline(file, line);
		double momentOfInertia = atof(line.c_str());

		// theta
		getline(file, line);
		double theta = atof(line.c_str());

		// omega
		getline(file, line);
		double omega = atof(line.c_str());

		// cmLset (center of mass)
		getline(file, line);
		Vector2d cmLset;
		iss.str(line);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		iss.clear();

		// number of points on the grain surface (INTEGER)
		getline(file, line);
		int npoints = atoi(line.c_str());

		// the point coordinates
		getline(file, line);
		vector<Vector2d> pointList(npoints);
		iss.str(line);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// bounding box radius
		getline(file, line);
		double bboxRadius = atof(line.c_str());

		// level set dimensions (INTEGERS)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file, line);
		vector<float> lsetvec(xdim*ydim);
		iss.str(line);
		for (int i = 0; i < xdim*ydim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

	// kn
		getline(file, line);
		double kn = atof(line.c_str());

        // ks
		getline(file, line);
		double ks = atof(line.c_str());

        // mu
		getline(file, line);
		double mu = atof(line.c_str());

	//cresN 
		getline(file, line);
		double cresN = atof(line.c_str());
	//cresS 
		getline(file, line);
		double cresS = atof(line.c_str());
		
		// Create level set object  - should we we use new?
		Levelset2d lset(lsetvec, xdim, ydim);

		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, velocity, momentOfInertia, theta, omega, cmLset, pointList, bboxRadius, lset, grainidx, morphologyID, kn, ks, mu, cresS, cresN);
	}

	return grainList;
} // end generateGrainsFromFile


vector<Vector2d> readPositionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string  line;
	string partial;
	istringstream iss;
	vector<Vector2d> positions;
	positions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		iss.clear();
	}
	return positions;
} // end readPositionFile


vector<double> readRotationFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> rotations;
	rotations.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		rotations[grainidx] = atof(partial.c_str());
		iss.clear();
	}
	return rotations;
} // end readRotationFile


void readShearHistFile(string filename, vector<Grain2d> & grains) {
	ifstream file(filename.c_str());

	if (file.good() == false)
	{
		cout << "here" << endl;
		std::cerr << "file does not exist" << endl ;
		exit(0) ; 
	}

	string   line;
	string 	partial;
	istringstream iss;
	size_t grainidx;
	size_t nodeidx;
	size_t contactingidx;
	double shear;
	while (true) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		grainidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		nodeidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		contactingidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		shear = atof(partial.c_str());
		grains[grainidx].changeShearHist(nodeidx, contactingidx, shear);
		iss.clear();
		if (file.eof()) {
			break;
		}
	}
}

// For imposing strain boundary conditions through wall displacements
vector<double> readWallPositionFile(string filename, size_t nwalls) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> wallPositions;
	wallPositions.resize(nwalls);
	for (size_t wallidx = 0; wallidx < nwalls; wallidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		wallPositions[wallidx] = atof(partial.c_str());
		iss.clear();
	}
	return wallPositions;
} // end readWallPositionFile


vector<Matrix<double, 5, 1> > readContactParameters(string filename, size_t contacttype) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Matrix<double, 5, 1> > contactParameters;
	contactParameters.resize(contacttype);
	for (size_t grainidx = 0; grainidx < contacttype; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		contactParameters[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](3) = atof(partial.c_str());
		getline(iss, partial, ' ');
		contactParameters[grainidx](4) = atof(partial.c_str());
		iss.clear();
	}
	return contactParameters;	
} // end readContactParameters


vector<double> readMaximumValues(string filename) {
    // return the minimum, maximum positions of all grains: x_min, x_max, y_min, y_max, max_radius;
    ifstream file(filename.c_str());
    string line;
    vector<double> extremevalues;
    for (size_t index = 0; index < 5; index++) {
        getline(file, line);
        extremevalues.push_back(atof(line.c_str()));
    }
    return extremevalues;
} // end readMaximumValues


vector<Vector3d> readCorrectionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector3d> corrections;
	corrections.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		corrections[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		corrections[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		corrections[grainidx](2) = atof(partial.c_str());
		iss.clear();

	}
	return corrections;
} // end readCorrectionFile

//get parameters 
void getParameters(const char* filename, const size_t & numInpFiles, const size_t &  numIntParam, size_t & numDecParam,  
		string * inputFiles, vector<size_t> & intParVal , vector<double> & parValues)
{
	ifstream infile(filename);
	if (!infile.is_open()) 
	{
		std::cout<<"Could not open parameter file" << filename << endl ;
		exit(0) ; 
	}
	cout << filename << endl;

	for (size_t i=0; i<numInpFiles ; i++) 
	{
		string ch;
		getline(infile,ch);
		size_t sl = ch.size();
		size_t pre = ch.find("=") + 2;
		inputFiles[i] =  ch.substr(pre,sl-pre);			
	}

	string dumy1 , dumy2 ; 
	for (size_t p=0 ; p < numIntParam ; p++) 
		infile >> dumy1 >> dumy2 >> intParVal[p] ; 

	cout << intParVal[1] << endl; 
	cout << intParVal[0] << endl;
	//numDecParam += intParVal[1] - intParVal[0] -1 ;  
	parValues.resize(numDecParam) ; 
	cout << numDecParam << endl;
	for (size_t p=0; p < numDecParam ;  p++)
		infile >> dumy1 >> dumy2 >> parValues[p] ; 
}


#endif // READINPUTFILE_H_ 
