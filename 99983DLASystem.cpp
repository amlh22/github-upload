//
//  DLASystem.cpp
//

#include "DLASystem.h"

// colors
namespace colours {
	GLfloat blue[] = { 0.1, 0.3, 0.9, 1.0 };   // blue
	GLfloat red[] = { 1.0, 0.2, 0.1, 0.2 };   // red
	GLfloat green[] = { 0.3, 0.6, 0.3, 1.0 };     // green
	GLfloat paleGrey[] = { 0.7, 0.7, 0.7, 1.0 };     // green
	GLfloat darkGrey[] = { 0.2, 0.2, 0.2, 1.0 };     // green
}


// this function gets called every step,
//   if there is an active particle then it gets moved,
//   if not then add a particle
void DLASystem::Update() {
	if (lastParticleIsActive == 1)
		moveLastParticle();
	else if (numParticles < endNum) {
		addParticleOnAddCircle();
		setParticleActive();
	}
	if (lastParticleIsActive == 0 || slowNotFast == 1)
		glutPostRedisplay(); //Tell GLUT that the display has changed
}


void DLASystem::clearParticles() {
	// delete particles and the particle list
	for (int i = 0; i<numParticles; i++) {
		delete particleList[i];
	}
	particleList.clear();
	numParticles = 0;
}

// remove any existing particles and setup initial condition
void DLASystem::Reset() {
	// stop running
	running = 0;

	clearParticles();

	lastParticleIsActive = 0;

	// set the grid to zero
	for (int i = 0; i<gridSize; i++) {
		for (int j = 0; j<gridSize; j++) {
			grid[i][j] = 0;
		}
	}

	// setup initial condition and parameters
	addCircle = 5;
	killCircle = 2.0*addCircle;
	clusterRadius = 0.0;
	// add a single particle at the origin, disable when enable square
	double pos[] = { 0.0, 0.0 };
	addParticle(pos);

	// set the view
	int InitialViewSize = 40;
	setViewSize(InitialViewSize);

	// Enable this section to start with seed particles creating a square

	/*for (double i =-InitialViewSize/4; i <= InitialViewSize/4; i++) {
	double pos[] = { i, -InitialViewSize/4 };
	addParticle(pos);

	double posy[] = { -InitialViewSize/4,i };
	addParticle(posy);

	double posx[] = { i,InitialViewSize/4 };
	addParticle(posx);

	double posy2[] = { InitialViewSize/4,i };
	addParticle(posy2);
	}*/

	// write initial values for Nc and R to the file
	writeToFile(numParticles, clusterRadius);
}

// check the given position (x, y) to determine if it lies inside the boundaries of our grid.
// returns true if position is valid, false otherwise.
bool DLASystem::isPositionValid(int x, int y) {
	if (x >= 0 && x < gridSize && y >= 0 && y < gridSize)
		return true;
	else
		return false;
}

// set the value of a grid cell for a particular position
// note the position has the initial particle at (0,0)
// but this corresponds to the middle of the grid array ie grid[ halfGrid ][ halfGrid ]
void DLASystem::setGrid(double pos[], int val) {
	int halfGrid = gridSize / 2;
	int x = (int)(pos[0] + halfGrid);
	int y = (int)(pos[1] + halfGrid);

	// make sure (x, y) lies within the boundaries of our grid to prevent corrupting the memory and crashing
	if (!isPositionValid(x, y)) {
		cout << "Cannot set cell. Invalid coordinates: (" << x << ", " << y << ")" << endl;
		return;
	}

	grid[x][y] = val;
}

// read the grid cell for a given position
int DLASystem::readGrid(double pos[]) {
	int halfGrid = gridSize / 2;
	int x = (int)(pos[0] + halfGrid);
	int y = (int)(pos[1] + halfGrid);

	// make sure (x, y) lies within the boundaries of our grid to prevent corrupting the memory and crashing
	if (!isPositionValid(x, y)) {
		cout << "Cannot set cell. Invalid coordinates: (" << x << ", " << y << ")" << endl;
		return 0;  // this should not happen.
				   // return 0 to indicate the cell is not occupied. (In reality, this cell does not exist).
	}

	return grid[x][y];
}

// check if the cluster is big enough and we should stop:
// to be safe, we need the killCircle to be at least 2 less than the edge of the grid
int DLASystem::checkStop() {
	if (killCircle + 2 >= gridSize / 2) {
		pauseRunning();
		cout << "STOP" << endl;
		glutPostRedisplay(); // update display
		return 1;
	}
	else return 0;
}

// return the current number of particles in the cluster
int DLASystem::getNumParticles() {
	return numParticles;
}

// get sticking probability
double DLASystem::getStickProbability() {
	return stickProbability;
}

// set sticking probability
void DLASystem::setStickProbability(double stickProbability) {
	double p = stickProbability;
	// make sure the probability is always in the range [0, 1]
	if (p < 0)
		p = 0;
	else if (p > 1)
		p = 1;

	this->stickProbability = p;
}

// add a particle to the system at a specific position
void DLASystem::addParticle(double pos[]) {
	// create a new particle
	Particle * p = new Particle(pos);
	// push_back means "add this to the end of the list"
	particleList.push_back(p);
	numParticles++;

	// pos coordinates should be -gridSize/2 < x < gridSize/2
	setGrid(pos, 1);
}

// add a particle to the system at a random position on the addCircle
// if we hit an occupied site then we do nothing except print a message
// (this should never happen)
void DLASystem::addParticleOnAddCircle() {
	double pos[2];
	double theta = rgen.random01() * 2 * M_PI;
	pos[0] = floor(addCircle * cos(theta));
	pos[1] = floor(addCircle * sin(theta));
	if (readGrid(pos) == 0)
		addParticle(pos);
	else
		cout << "FAIL " << pos[0] << " " << pos[1] << endl;
}

// send back the position of a neighbour of a given grid cell
// NOTE: there is no check that the neighbour is inside the grid,
// this has to be done separately...
void DLASystem::setPosNeighbour(double setpos[], double pos[], int val) {
	if (attachDiag == 0 && val > 3) {
		// Trying to find a diagonal neighbour while diagonal sticking is disabled.
		// This is an error and should not happen. Print error message and return.
		cout << "Invalid neighbour value (" << val << "). attachDiag = " << attachDiag << "." << endl;
		return;
	}

	switch (val) {
	case 0:  // Right
		setpos[0] = pos[0] + 1.0;
		setpos[1] = pos[1];
		break;
	case 1:  // Left
		setpos[0] = pos[0] - 1.0;
		setpos[1] = pos[1];
		break;
	case 2:  // Up
		setpos[0] = pos[0];
		setpos[1] = pos[1] + 1.0;
		break;
	case 3:  // Down
		setpos[0] = pos[0];
		setpos[1] = pos[1] - 1.0;
		break;
	case 4:  // Up and Right
		setpos[0] = pos[0] + 1.0;
		setpos[1] = pos[1] + 1.0;
		break;
	case 5:  // Up and Left
		setpos[0] = pos[0] - 1.0;
		setpos[1] = pos[1] + 1.0;
		break;
	case 6:  // Down and Right
		setpos[0] = pos[0] + 1.0;
		setpos[1] = pos[1] - 1.0;
		break;
	case 7:  // Down and Left
		setpos[0] = pos[0] - 1.0;
		setpos[1] = pos[1] - 1.0;
		break;
	}
}

// if the view is smaller than the kill circle then increase the view area (zoom out)
void DLASystem::updateViewSize() {
	double mult = 1.2;
	if (viewSize < 2.0*killCircle) {
		setViewSize(viewSize * mult);
	}
}

// set the view to be the size of the add circle (ie zoom in on the cluster)
void DLASystem::viewAddCircle() {
	setViewSize(2.0*addCircle);  // factor of 2 is to go from radius to diameter
}

// return the current radius of the particle cluster
double DLASystem::getClusterRadius() {
	return clusterRadius;
}

// when we add a particle to the cluster, we should update the cluster radius
// and the sizes of the addCircle and the killCircle
void DLASystem::updateClusterRadius(double pos[]) {

	double rr = distanceFromOrigin(pos);
	if (rr > clusterRadius) {
		clusterRadius = rr;
		// this is how big addCircle is supposed to be:
		//   either 20% more than cluster radius, or at least 2 bigger.
		double check = clusterRadius * addRatio;
		if (check < clusterRadius + 2)
			check = clusterRadius + 2;
		// if it is smaller then update everything...
		if (addCircle < check) {
			addCircle = check;
			killCircle = killRatio * addCircle;
			updateViewSize();
		}
		checkStop();
	}
}

// make a random move of the last particle in the particleList
void DLASystem::moveLastParticle() {
	int rr = rgen.randomInt(4);  // pick a random number in the range 0-3, which direction do we hop?
	double newpos[2];

	Particle *lastP = particleList[numParticles - 1];

	setPosNeighbour(newpos, lastP->pos, rr);

	if (distanceFromOrigin(newpos) > killCircle) {
		//cout << "#deleting particle" << endl;
		setGrid(lastP->pos, 0);
		particleList.pop_back();  // remove particle from particleList
		numParticles--;
		setParticleInactive();
	}
	// check if destination is empty
	else if (readGrid(newpos) == 0) {
		setGrid(lastP->pos, 0);  // set the old grid site to empty
								 // update the position
		particleList[numParticles - 1]->pos[0] = newpos[0];
		particleList[numParticles - 1]->pos[1] = newpos[1];
		setGrid(lastP->pos, 1);  // set the new grid site to be occupied

								 // check if we stick
		if (checkStick()) {
			//cout << "stick" << endl;
			setParticleInactive();  // make the particle inactive (stuck)
			updateClusterRadius(lastP->pos);  // update the cluster radius, addCircle, etc.

			// print Nc and R to the output file for Nc = 1, 10, 100, 200, 300, ..., 1000.
			if (numParticles == 1 || numParticles <= 10 || numParticles % 100 == 0) {
				writeToFile(numParticles, clusterRadius);
			}

			if (numParticles >= endNum) {
				// close the file if the simulation has ended
				if (logfile.is_open())
					logfile.close();
			}
		}
	}
	else {
		// if we get to here then we are trying to move to an occupied site
		// (this should never happen as long as the sticking probability is 1.0)
		cout << "particle bounced off at position (" << lastP->pos[0] << ", " << lastP->pos[1] << ")" << endl;
		
		if (!bounceEnabled) {
			// we should not get here if bouncing is disabled.
			cout << "(this should not happen with sticking probability = 1)";
		}
	}
}

// check if the last particle should stick (to a neighbour)
int DLASystem::checkStick() {
	Particle *lastP = particleList[numParticles - 1];
	int result = 0;

	// loop over neighbours

	// if diagonal sticking is disabled, we have 4 neighbours to check: Right, Left, Up and Down.
	// if diagonal sticking is enabled, we have 4 additional diagonal neighbours to check (8 in total).
	int numNeighbours = (attachDiag == 0) ? 4 : 8;

	for (int i = 0; i < numNeighbours; i++) {
		double checkpos[2];
		setPosNeighbour(checkpos, lastP->pos, i);
		// if the neighbour is occupied...
		if (readGrid(checkpos) == 1) {
			if (!shouldBounceOff()) {
				// stick only if particle does not bounce off.
				result = 1;
			}
			else {
				// if particle touches the cluster but decides to bounce off, treat this as if it didn't
				// touch the cluster at all. in this case, particle will continue moving until it sticks
				// or gets too far away from the cluster and killed. there is nothing extra that needs
				// to be done here.
			}
		}
	}


	return result;
}

// decide if a particle that touches the cluster should bounce off
bool DLASystem::shouldBounceOff() {
	if (!bounceEnabled)  // if bouncing is not enabled, it should always stick
		return false;
	if (stickProbability == 1.0)
		return false;  // if p = 1, it should always stick

	// this is the upper cutoff that is used in probability calculation. its value shouldn't affect
	// the calculation as long as it is sufficiently large. however, setting a too low value can
	// affect the precision of the calculation negatively.
	int m = 1000;

	// this is the random number. value returned by the randomInt() function does not include its
	// upper bound, so we add +1 to make it inclusive. i.e. we want a number in the range [0, m].
	int x = rgen.randomInt(m + 1);

	double ratio = (double)x / (double)m;
	bool shouldStick = (ratio < stickProbability);
	bool shouldBounceOff = !shouldStick;
	return shouldBounceOff;
}

// print the number of particles (Nc) and size (R) of the cluster to a file.
void DLASystem::writeToFile(int Nc, double R) {
	if (!logfile.is_open()) {  // if the file is not already opened, open it
		logfile.open(OUTPUT_FILE_NAME, ios::app);  // ios:trunc means truncate (i.e. delete the file
													 // if it already exists and create a new one)
		if (!logfile.is_open()) {
			// if something goes wrong and we can't open the file, print an error message
			

			char errorMessage[128];  // char buffer of arbitrary size to hold the error message
			strerror_s(errorMessage, 128, errno);  // errno is a global int that holds the last error number.
												   // here, we convert the error number to a human-readable string
												   // to figure out what kind of error this is.
			cout << "Error opening file (" << errno << ") : " << errorMessage << endl;
			if (errno == 13) {
				// this is a common "Permission denied" error. print some extra info to help evaluator fix the problem.
				cout << "Please make sure you have write permissions and folder in which ";
				cout << "executable is running is not set to read-only." << endl;
			}
			return;
		}

		logfile << "Nc\tR" << endl;  // print a simple header which consists of Nc and R,e
									 // separated by the tab character
	}

	// print Nc and R, separated by the tab character
	logfile << numParticles << "\t" << clusterRadius << endl;
}

// constructor
DLASystem::DLASystem(Window *set_win) {
	cout << "creating system, gridSize " << gridSize << endl;
	win = set_win;
	numParticles = 0;
	endNum = 4000;

	// make sure diagonal sticking is disabled by default
	attachDiag = 0;

	// value that controls the "stickiness" of particles. if "bounceEnabled" is false,
	// this value will be ignored.
	// this is used in "Exercise 3".
	stickProbability = 1;

	// allocate memory for the grid, remember to free the memory in destructor
	grid = new int*[gridSize];
	for (int i = 0; i<gridSize; i++) {
		grid[i] = new int[gridSize];
	}
	slowNotFast = 1;
	// reset initial parameters
	Reset();

	addRatio = 1.5;   // how much bigger the addCircle should be, compared to cluster radius
	killRatio = 2.0;   // how much bigger is the killCircle, compared to the addCircle

					   // this opens a logfile, if we want to...
					   //logfile.open("opfile.txt");
}

// destructor
DLASystem::~DLASystem() {
	// strictly we should not print inside the destructor but never mind...
	cout << "deleting system" << endl;
	// delete the particles
	clearParticles();
	// delete the grid
	for (int i = 0; i<gridSize; i++)
		delete[] grid[i];
	delete[] grid;

	if (logfile.is_open())
		logfile.close();
}



// this draws the system
void DLASystem::DrawSquares() {

	// draw the particles
	double halfSize = 0.5;
	for (int p = 0; p<numParticles; p++) {
		double *vec = particleList[p]->pos;
		glPushMatrix();
		if (p == numParticles - 1 && lastParticleIsActive == 1)
			glColor4fv(colours::red);
		else if (p == 0)
			glColor4fv(colours::green);
		else
			glColor4fv(colours::blue);
		glRectd(drawScale*(vec[0] - halfSize),
			drawScale*(vec[1] - halfSize),
			drawScale*(vec[0] + halfSize),
			drawScale*(vec[1] + halfSize));
		glPopMatrix();
	}

	// print some information (at top left)
	// this ostringstream is a way to create a string with numbers and words (similar to cout << ... )
	ostringstream str;
	str << "num " << numParticles << " size " << clusterRadius << " | ";
	
	if (bounceEnabled) {
		str << "sticking probability = " << getStickProbability();
	}
	else {
		str << "bouncing disabled";
	}

	// print the string
	win->displayString(str, -0.9, 0.9, colours::red);

	// if we are paused then print this (at bottom left)
	if (running == 0) {
		ostringstream pauseStr;
		pauseStr << "paused";
		win->displayString(pauseStr, -0.9, -0.9, colours::red);
	}

}
