/*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#include <sbpl/discrete_space_information/environment_navxythetavsteer.h>
#include <sbpl/planners/planner.h>
#include <sbpl/utils/2Dgridsearch.h>
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
#include <../../../Software/ACADOtoolkit/examples/ocp/rocket_with_templates_ocp.hpp>
#include <string.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

#if TIME_DEBUG
static clock_t time3_addallout = 0;
static clock_t time_gethash = 0;
static clock_t time_createhash = 0;
static clock_t time_getsuccs = 0;
#endif

static long int checks = 0;

//function prototypes

/* Constructor and destructor */

EnvironmentNAVXYTHETAVSTEER::EnvironmentNAVXYTHETAVSTEER(){
	EnvNAVXYTHETAVSTEERCfg.obsthresh = NAVXYTHETAVSTEER_DEFAULTOBSTHRESHOLD;
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh = EnvNAVXYTHETAVSTEERCfg.obsthresh; 
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh = -1; 

	grid2Dsearchfromstart = NULL;
	grid2Dsearchfromgoal = NULL;
	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;
	
	iteration = 0;

	EnvNAVXYTHETAVSTEER.bInitialized = false;

	EnvNAVXYTHETAVSTEERCfg.NumThetaDirs = NAVXYTHETAVSTEER_DEFAULTTHETADIRS;

	//no memory allocated in cfg yet
	EnvNAVXYTHETAVSTEERCfg.Grid2D = NULL;
	EnvNAVXYTHETAVSTEERCfg.ActionsV = NULL;
	EnvNAVXYTHETAVSTEERCfg.PredActionsV = NULL;
}

EnvironmentNAVXYTHETAVSTEER::~EnvironmentNAVXYTHETAVSTEER(){
	SBPL_PRINTF("destroying XYTHETAVSTEER\n");
	if (grid2Dsearchfromstart != NULL) delete grid2Dsearchfromstart;
	grid2Dsearchfromstart = NULL;

	if (grid2Dsearchfromgoal != NULL) delete grid2Dsearchfromgoal;
	grid2Dsearchfromgoal = NULL;

	if (EnvNAVXYTHETAVSTEERCfg.Grid2D != NULL) {
		for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++)
			delete[] EnvNAVXYTHETAVSTEERCfg.Grid2D[x];
		delete[] EnvNAVXYTHETAVSTEERCfg.Grid2D;
		EnvNAVXYTHETAVSTEERCfg.Grid2D = NULL;
	}

	//delete actions
	if (EnvNAVXYTHETAVSTEERCfg.ActionsV != NULL) {
		for (int ind = 0; ind < (int)EnvNAVXYTHETAVSTEERCfg.ActionsV->size(); ind++)
			delete EnvNAVXYTHETAVSTEERCfg.ActionsV->at(ind);
		delete EnvNAVXYTHETAVSTEERCfg.ActionsV;
		EnvNAVXYTHETAVSTEERCfg.ActionsV = NULL;
	}
	if (EnvNAVXYTHETAVSTEERCfg.PredActionsV != NULL) {
		delete[] EnvNAVXYTHETAVSTEERCfg.PredActionsV;
		EnvNAVXYTHETAVSTEERCfg.PredActionsV = NULL;
	}
	
	//delete the states themselves first
	for (int i = 0; i < (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size(); i++) {
		delete EnvNAVXYTHETAVSTEER.StateID2DataTable.at(i);
		EnvNAVXYTHETAVSTEER.StateID2DataTable.at(i) = NULL;
	}
	EnvNAVXYTHETAVSTEER.StateID2DataTable.clear();

	//delete hashtable
	if (EnvNAVXYTHETAVSTEER.Data2StateIDHashTable != NULL) {
		delete[] EnvNAVXYTHETAVSTEER.Data2StateIDHashTable;
		EnvNAVXYTHETAVSTEER.Data2StateIDHashTable = NULL;
	}
}

//-------------------problem specific and local functions---------------------

static unsigned int inthash(unsigned int key)
{
	key += (key << 12);
	key ^= (key >> 22);
	key += (key << 4);
	key ^= (key >> 9);
	key += (key << 10);
	key ^= (key >> 2);
	key += (key << 7);
	key ^= (key >> 12);
	return key;
}

//examples of hash functions: map state coordinates onto a hash value
//#define GETHASHBIN(X, Y) (Y*WIDTH_Y+X) 
//here we have state coord: <X1, X2, X3, X4>

/*
* PAY ATTENTION TO UNSIGNED!!! I WAIT TO CHANGE IT, BECAUSE
* I HAVE TO DECIDE IF DISCRETIZE VELOCITY OR NOT!!!
*/
unsigned int EnvironmentNAVXYTHETAVSTEER::GETHASHBIN(unsigned int x, unsigned int y, unsigned int theta, unsigned int v, unsigned int steer)
{
	return inthash((inthash(x) + (inthash(y) << 1) + (inthash(theta) << 2) + (inthash(v) << 3) + (inthash(steer) << 4))) &
		(EnvNAVXYTHETAVSTEER.HashTableSize - 1);
}

void EnvironmentNAVXYTHETAVSTEER::PrintHashTableHist()
{
	int s0 = 0, s1 = 0, s50 = 0, s100 = 0, s200 = 0, s300 = 0, slarge = 0;

	for (int j = 0; j < (int)EnvNAVXYTHETAVSTEER.HashTableSize; j++) {
		if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() == 0)
			s0++;
		else if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() < 50)
			s1++;
		else if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() < 100)
			s50++;
		else if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() < 200)
			s100++;
		else if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() < 300)
			s200++;
		else if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[j].size() < 400)
			s300++;
		else
			slarge++;
	}
	SBPL_PRINTF("hash table histogram: 0:%d, <50:%d, <100:%d, <200:%d, <300:%d, <400:%d >400:%d\n", s0, s1, s50, s100,
				s200, s300, slarge);
}

void EnvironmentNAVXYTHETAVSTEER::ReadConfiguration(FILE* fCfg)
{
	//read in the configuration of environment and initialize  EnvCfg structure
	char sTemp[1024], sTemp1[1024];
	int dTemp;
	int x, y;

	//discretization(cells)
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (discretization)\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "discretization(cells):");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format (discretization)\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (discretization)\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EnvWidth_c = atoi(sTemp);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (discretization)\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EnvHeight_c = atoi(sTemp);

	// Scan for optional NumThetaDirs parameter. Check for following obsthresh.
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "NumThetaDirs:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	EnvNAVXYTHETAVSTEERCfg.NumThetaDirs = atoi(sTemp);

	//numv: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (obsthresh)\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "NumV:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.numV = atoi(sTemp);
	
	//velocities: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (obsthresh)\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "Velocities:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}
	
	/* SUBSTITUTE VECTOR WITH POINTER */
	for(int i=0;i<EnvNAVXYTHETAVSTEERCfg.numV;i++){
		if (fscanf(fCfg, "%s", sTemp) != 1) {
			SBPL_ERROR("ERROR: ran out of env file early\n");
			throw new SBPL_Exception();
		}
		
		EnvNAVXYTHETAVSTEERCfg.velocities.push_back(atof(sTemp));
	}
	
	// Scan for optional NumSteers parameter.
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "NumSteers:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	EnvNAVXYTHETAVSTEERCfg.numSteers = atoi(sTemp);
	
	//obsthresh: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (obsthresh)\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "obsthresh:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}

	// obsthresh
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.obsthresh = atoi(sTemp);
	SBPL_PRINTF("obsthresh = %d\n", EnvNAVXYTHETAVSTEERCfg.obsthresh);

	//cost_inscribed_thresh: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "cost_inscribed_thresh:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh = atoi(sTemp);
	SBPL_PRINTF("cost_inscribed_thresh = %d\n", EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh);

	//cost_possibly_circumscribed_thresh: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "cost_possibly_circumscribed_thresh:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh = atoi(sTemp);
	SBPL_PRINTF("cost_possibly_circumscribed_thresh = %d\n", EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh);

	//cellsize
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "cellsize(meters):");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.cellsize_m = atof(sTemp);

	//start(meters,rads,m/s,rads):
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.StartX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.StartY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.StartTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.StartV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.velocities);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.StartSteer = ContSteer2Disc(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.numSteers);

	if (EnvNAVXYTHETAVSTEERCfg.StartX_c < 0 || EnvNAVXYTHETAVSTEERCfg.StartX_c >= EnvNAVXYTHETAVSTEERCfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.StartY_c < 0 || EnvNAVXYTHETAVSTEERCfg.StartY_c >= EnvNAVXYTHETAVSTEERCfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.StartTheta < 0 || EnvNAVXYTHETAVSTEERCfg.StartTheta >= EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.StartV < 0 || EnvNAVXYTHETAVSTEERCfg.StartV >= EnvNAVXYTHETAVSTEERCfg.numV) {
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.StartSteer < 0 || EnvNAVXYTHETAVSTEERCfg.StartSteer >= EnvNAVXYTHETAVSTEERCfg.numSteers) {
		SBPL_ERROR("ERROR: illegal start coordinates for steer\n");
		throw new SBPL_Exception();
	}

	//end(meters,rads,m/s):
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EndX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EndY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EndTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EndV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.velocities);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVSTEERCfg.EndSteer = ContSteer2Disc(atof(sTemp), EnvNAVXYTHETAVSTEERCfg.numSteers);

	if (EnvNAVXYTHETAVSTEERCfg.EndX_c < 0 || EnvNAVXYTHETAVSTEERCfg.EndX_c >= EnvNAVXYTHETAVSTEERCfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.EndY_c < 0 || EnvNAVXYTHETAVSTEERCfg.EndY_c >= EnvNAVXYTHETAVSTEERCfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.EndTheta < 0 || EnvNAVXYTHETAVSTEERCfg.EndTheta >= EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.EndV < 0 || EnvNAVXYTHETAVSTEERCfg.EndV >= EnvNAVXYTHETAVSTEERCfg.numV) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.EndSteer < 0 || EnvNAVXYTHETAVSTEERCfg.EndSteer >= EnvNAVXYTHETAVSTEERCfg.numSteers) {
		SBPL_ERROR("ERROR: illegal start coordinates for steer\n");
		throw new SBPL_Exception();
	}

	//allocate the 2D environment
	EnvNAVXYTHETAVSTEERCfg.Grid2D = new unsigned char*[EnvNAVXYTHETAVSTEERCfg.EnvWidth_c];
	for (x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVSTEERCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETAVSTEERCfg.EnvHeight_c];
	}

	//environment:
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	for (y = EnvNAVXYTHETAVSTEERCfg.EnvHeight_c-1; y >= 0; y--)
		for (x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
			if (fscanf(fCfg, "%d", &dTemp) != 1) {
				SBPL_ERROR("ERROR: incorrect format of config file\n");
				throw new SBPL_Exception();
			}
			EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] = dTemp;
		}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c, EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.StartSteer)){
		SBPL_ERROR("ERROR: invalid start state\n");
		throw new SBPL_Exception();
	}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c, EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.EndV, 0)){
		SBPL_ERROR("ERROR: invalid end state\n");
		throw new SBPL_Exception();
	}
}

EnvNAVXYTHETAVSTEERHashEntry_t* EnvironmentNAVXYTHETAVSTEER::GetHashEntry(unsigned int x, unsigned int y, unsigned int theta, unsigned int v, unsigned int steer)
{
	//clock_t currenttime = clock();

	int binid = GETHASHBIN(x, y, theta, v, steer);

#if DEBUG
	if ((int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid].size() > 500)
	{
		SBPL_PRINTF("WARNING: Hash table has a bin %d (X1=%d X2=%d X3=%d X4=%d) of size %d\n",
					binid, X1, X2, X3, X4, (int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid].size());

		PrintHashTableHist();
	}
#endif

	//iterate over the states in the bin and select the perfect match
	for (int ind = 0; ind < (int)EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid].size(); ind++) {
		if (EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind]->x == x &&
			EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind]->y == y &&
			EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind]->theta == theta &&
			EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind]->v == v &&
			EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind]->steer == steer)
		{
			//time_gethash += clock()-currenttime;
			return EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[binid][ind];
		}
	}

	//time_gethash += clock()-currenttime;

	return NULL;
}

EnvNAVXYTHETAVSTEERHashEntry_t* EnvironmentNAVXYTHETAVSTEER::CreateNewHashEntry(unsigned int x, unsigned int y, unsigned int theta,
													unsigned int v, unsigned int steer)
{
	int i;

	//clock_t currenttime = clock();

	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = new EnvNAVXYTHETAVSTEERHashEntry_t;

	HashEntry->x = x;
	HashEntry->y = y;
	HashEntry->theta = theta;
	HashEntry->v = v;
	HashEntry->steer = steer;
	HashEntry->iteration = 0;

	HashEntry->stateID = EnvNAVXYTHETAVSTEER.StateID2DataTable.size();

	//insert into the tables
	EnvNAVXYTHETAVSTEER.StateID2DataTable.push_back(HashEntry);

	//get the hash table bin
	i = GETHASHBIN(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->steer);

	//insert the entry into the bin
	EnvNAVXYTHETAVSTEER.Data2StateIDHashTable[i].push_back(HashEntry);

	/*
	* REVIEW BELOW OPERATIONS. THEY WAS ALREADY IMPLEMENTED
	*/
	//insert into and initialize the mappings
	int* entry = new int[NUMOFINDICES_STATEID2IND];
	StateID2IndexMapping.push_back(entry);
	for (i = 0; i < NUMOFINDICES_STATEID2IND; i++) {
		StateID2IndexMapping[HashEntry->stateID][i] = -1;
	}

	if (HashEntry->stateID != (int)StateID2IndexMapping.size() - 1) {
		SBPL_ERROR("ERROR in Env... function: last state has incorrect stateID\n");
		throw new SBPL_Exception();
	}

	//time_createhash += clock()-currenttime;

	return HashEntry;
}

void EnvironmentNAVXYTHETAVSTEER::CreateStartandGoalStates()
{
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry;

	//create start state
	unsigned int x = 0;
	unsigned int y = 0;
	unsigned int theta = 0;
	unsigned int v = 0;
	unsigned int steer = 0;
	HashEntry = CreateNewHashEntry(x, y, theta, v, steer);
	EnvNAVXYTHETAVSTEER.startstateid = HashEntry->stateID;

	//create goal state
	x = y = theta = v = steer = 1;
	HashEntry = CreateNewHashEntry(x, y, theta, v, steer);
	EnvNAVXYTHETAVSTEER.goalstateid = HashEntry->stateID;
}

void EnvironmentNAVXYTHETAVSTEER::ComputeReplanningDataforAction(EnvNAVXYTHETAVSTEERAction_t* action){
	int j;

	//iterate over all the cells involved in the action
	sbpl_xy_theta_v_steer_cell_t startcell5d, endcell5d;
	for (int i = 0; i < (int)action->intersectingcellsV.size(); i++) {
		//compute the translated affected search Pose - what state has an
		//outgoing action whose intersecting cell is at 0,0
		startcell5d.steer = action->startsteer;
		startcell5d.v = action->startv;
		startcell5d.theta = action->starttheta;
		startcell5d.x = -action->intersectingcellsV.at(i).x;
		startcell5d.y = -action->intersectingcellsV.at(i).y;

		//compute the translated affected search Pose - what state has an
		//incoming action whose intersecting cell is at 0,0
		endcell5d.steer = action->endsteer;
		endcell5d.v = action->endv;
		endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		endcell5d.x = startcell5d.x + action->dX;
		endcell5d.y = startcell5d.y + action->dY;

		//store the cells if not already there
		for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
			if (affectedsuccstatesV.at(j) == endcell5d) break;
		}
		if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell5d);

		for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
			if (affectedpredstatesV.at(j) == startcell5d) break;
		}
		if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell5d);
	}//over intersecting cells

	//---intersecting cell = origin
	//compute the translated affected search Pose - what state has an outgoing action whose intersecting cell is at 0,0
	startcell5d.steer = action->startsteer;
	startcell5d.v = action->startv;
	startcell5d.theta = action->starttheta;
	startcell5d.x = -0;
	startcell5d.y = -0;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell5d.steer = action->endsteer;
	endcell5d.v = action->endv;
	endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	endcell5d.x = startcell5d.x + action->dX;
	endcell5d.y = startcell5d.y + action->dY;

	//store the cells if not already there
	for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
		if (affectedsuccstatesV.at(j) == endcell5d) break;
	}
	if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell5d);

	for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
		if (affectedpredstatesV.at(j) == startcell5d) break;
	}
	if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell5d);

	//---intersecting cell = outcome state
	//compute the translated affected search Pose - what state has an outgoing action whose intersecting cell is at 0,0
	startcell5d.steer = action->startsteer;
	startcell5d.v = action->startv;
	startcell5d.theta = action->starttheta;
	startcell5d.x = -action->dX;
	startcell5d.y = -action->dY;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell5d.steer = action->endsteer;
	endcell5d.v = action->endv;
	endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	endcell5d.x = startcell5d.x + action->dX;
	endcell5d.y = startcell5d.y + action->dY;

	for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
		if (affectedsuccstatesV.at(j) == endcell5d) break;
	}
	if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell5d);

	for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
		if (affectedpredstatesV.at(j) == startcell5d) break;
	}
	if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell5d);
}

void EnvironmentNAVXYTHETAVSTEER::ComputeReplanningData()
{
	//iterate over all actions
	//Steer angles
	for(int sind = 0;sind < EnvNAVXYTHETAVSTEERCfg.numSteers;sind++){
		//velocities
		for(int vind = 0; vind < EnvNAVXYTHETAVSTEERCfg.numV; vind++){
			//orientations
			for (int tind = 0; tind < EnvNAVXYTHETAVSTEERCfg.NumThetaDirs; tind++) {
				int idx = sind*EnvNAVXYTHETAVSTEERCfg.numV*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+vind*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+tind;
				//actions
				for (int aind = 0; aind < (int)(EnvNAVXYTHETAVSTEERCfg.ActionsV->at(idx)->size()); aind++) {
					//compute replanning data for this action 
					ComputeReplanningDataforAction(&EnvNAVXYTHETAVSTEERCfg.ActionsV->at(idx)->at(aind));
				}
			}
		}
	}
}

void EnvironmentNAVXYTHETAVSTEER::PrecomputeActionswithCompleteMotionPrimitive(vector<SBPL_xythetavsteer_mprimitive>* motionprimitiveV){
	SBPL_PRINTF("Pre-computing action data using motion primitives for every pair velocity/angle...\n");
	EnvNAVXYTHETAVSTEERCfg.ActionsV = new vector<vector<EnvNAVXYTHETAVSTEERAction_t> *>();
	EnvNAVXYTHETAVSTEERCfg.PredActionsV = new vector<EnvNAVXYTHETAVSTEERActionIndex_t> [EnvNAVXYTHETAVSTEERCfg.numSteers*EnvNAVXYTHETAVSTEERCfg.numV*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs];
	
	vector<sbpl_2Dcell_t> footprint;

	/*
	* if (motionprimitiveV->size() % EnvNAVXYTHETALATCfg.NumThetaDirs != 0) {
	* SBPL_ERROR("ERROR: motionprimitives should be uniform across actions\n");
	* throw new SBPL_Exception();
	* }
	*/

	//EnvNAVXYTHETAVCfg.actionwidth = ((int)motionprimitiveV->size()) / EnvNAVXYTHETALATCfg.NumThetaDirs;

	//Allocate vectors
	for(int sind = 0; sind < EnvNAVXYTHETAVSTEERCfg.numSteers; sind++)
		for(int vind = 0; vind < EnvNAVXYTHETAVSTEERCfg.numV; vind++)
			for (int tind = 0; tind < EnvNAVXYTHETAVSTEERCfg.NumThetaDirs; tind++){
				EnvNAVXYTHETAVSTEERCfg.ActionsV->push_back(new vector<EnvNAVXYTHETAVSTEERAction_t>());
			}
	
	//iterate over source angles|speeds|angles
	int maxnumofactions = 0;
	for(int sind = 0; sind < EnvNAVXYTHETAVSTEERCfg.numSteers; sind++){
		for(int vind = 0; vind < EnvNAVXYTHETAVSTEERCfg.numV; vind++){
			for (int tind = 0; tind < EnvNAVXYTHETAVSTEERCfg.NumThetaDirs; tind++) {
				SBPL_PRINTF("pre-computing for pair (angle, speed, angle) (%d,%d,%d) out of (%d,%d,%d)\n", sind, vind, tind, EnvNAVXYTHETAVSTEERCfg.numSteers, EnvNAVXYTHETAVSTEERCfg.numV, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);

				//compute current index
				int vector_index = sind*EnvNAVXYTHETAVSTEERCfg.numV*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+vind*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+tind;
				
				//EnvNAVXYTHETAVCfg.ActionsV->at(vector_index) = new vector<EnvNAVXYTHETAVAction_t>();

				//compute sourcepose
				sbpl_xy_theta_v_steer_pt_t sourcepose;
				sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
				sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
				sourcepose.theta = DiscTheta2ContNotUnif(tind, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
				sourcepose.v = DiscV2Cont(vind, EnvNAVXYTHETAVSTEERCfg.velocities);
				sourcepose.steer = DiscSteer2Cont(sind, EnvNAVXYTHETAVSTEERCfg.numSteers);

				//iterate over motion primitives
				int numofactions = 0;
				int aind = -1;
				for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
					//find a motion primitive for this angle
					if (motionprimitiveV->at(mind).start_theta_disc != tind || motionprimitiveV->at(mind).start_v_disc != vind || motionprimitiveV->at(mind).start_steer_disc != sind) continue;

					aind++;
					numofactions++;
					EnvNAVXYTHETAVSTEERAction_t element_to_add;
					
					//action index
					element_to_add.aind = aind;

					//start angle
					element_to_add.starttheta = tind;
					
					//start velocity
					element_to_add.startv = vind;
					
					//start steer angle
					element_to_add.startsteer = sind;

					//compute dislocation
					element_to_add.endtheta = motionprimitiveV->at(mind).endcell.theta;
					element_to_add.endv = motionprimitiveV->at(mind).endcell.v;
					element_to_add.endsteer = motionprimitiveV->at(mind).endcell.steer;
					element_to_add.dX = motionprimitiveV->at(mind).endcell.x;
					element_to_add.dY = motionprimitiveV->at(mind).endcell.y;

					//compute and store interm points as well as intersecting cells
					element_to_add.intersectingcellsV.clear();
					element_to_add.intermptV.clear();
					element_to_add.interm3DcellsV.clear();

					sbpl_xy_theta_v_steer_cell_t previnterm3Dcell;
					previnterm3Dcell.x = 0;
					previnterm3Dcell.y = 0;

					// Compute all the intersected cells for this action (intermptV and interm3DcellsV)
					for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
						sbpl_xy_theta_v_steer_pt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
						element_to_add.intermptV.push_back(intermpt);

						// also compute the intermediate discrete cells if not there already
						sbpl_xy_theta_v_steer_pt_t pose;
						pose.x = intermpt.x + sourcepose.x;
						pose.y = intermpt.y + sourcepose.y;
						pose.theta = intermpt.theta;
						pose.v = intermpt.v;
						pose.steer = intermpt.steer;

						sbpl_xy_theta_v_steer_cell_t intermediate2dCell;
						intermediate2dCell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
						intermediate2dCell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETAVSTEERCfg.cellsize_m);

						// add unique cells to the list
						if (element_to_add.interm3DcellsV.size() == 0 || intermediate2dCell.x
							!= previnterm3Dcell.x || intermediate2dCell.y != previnterm3Dcell.y) {
							element_to_add.interm3DcellsV.push_back(intermediate2dCell);
						}

						previnterm3Dcell = intermediate2dCell;
					}

					//compute linear and angular time
					double linear_distance = 0;
					//double medium_velocity = element_to_add.intermptV[0].v;
					//double xmax=element_to_add.intermptV[0].x, xmin=element_to_add.intermptV[0].x, ymax=element_to_add.intermptV[0].y, ymin=element_to_add.intermptV[0].y;
					for (unsigned int i = 1; i<element_to_add.intermptV.size(); i++) {
						double x0 = element_to_add.intermptV[i - 1].x;
						double y0 = element_to_add.intermptV[i - 1].y;
						double x1 = element_to_add.intermptV[i].x;
						double y1 = element_to_add.intermptV[i].y;
						double dx = x1 - x0;
						double dy = y1 - y0;
						linear_distance += sqrt(dx * dx + dy * dy);
						//medium_velocity += fabs(element_to_add.intermptV[i].v);
						
						/*if(element_to_add.intermptV[i].x > xmax)
							xmax = element_to_add.intermptV[i].x;
						
						if(element_to_add.intermptV[i].x < xmin)
							xmin = element_to_add.intermptV[i].x;
						
						if(element_to_add.intermptV[i].y > ymax)
							ymax = element_to_add.intermptV[i].y;
						
						if(element_to_add.intermptV[i].y < ymin)
							ymin = element_to_add.intermptV[i].y;*/
					}
					
					//medium_velocity /= (int)element_to_add.intermptV.size();
					//double linear_time = linear_distance / medium_velocity;
					
					/*
					* Have some sense in ackermann vehicle?
					double angular_distance =
							fabs(computeMinUnsignedAngleDiff(DiscTheta2Cont(EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).endtheta,
																			EnvNAVXYTHETAVCfg.NumThetaDirs),
															DiscTheta2Cont(EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).starttheta,
																			EnvNAVXYTHETAVCfg.NumThetaDirs)));
					double angular_time = angular_distance / ((PI_CONST / 4.0) /
										EnvNAVXYTHETALATCfg.timetoturn45degsinplace_secs);
					*/
					
					//make the cost the max of the two times
					element_to_add.cost = NAVXYTHETAVSTEER_COSTMULT_MTOMM * motionprimitiveV->at(mind).additionalactioncostmult;
					//element_to_add.cost = NAVXYTHETAV_COSTMULT_MTOMM*(linear_distance/medium_velocity);
					//double linear_distance_action = element_to_add.dX*element_to_add.dX + element_to_add.dY*element_to_add.dY;
					//element_to_add.cost = (int)((linear_distance/linear_distance_action) * NAVXYTHETAV_COSTMULT_MTOMM * motionprimitiveV->at(mind).additionalactioncostmult);
					//element_to_add.cost = linear_time * NAVXYTHETAV_COSTMULT_MTOMM * motionprimitiveV->at(mind).additionalactioncostmult;
					//element_to_add.cost = ((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin))*1000;
					
					//if(element_to_add.cost <= 0)
						//printf("%d\n",element_to_add.cost);
					
							//(int)(ceil(NAVXYTHETALAT_COSTMULT_MTOMM * max(linear_time, angular_time)));
					//use any additional cost multiplier
					//EnvNAVXYTHETALATCfg.ActionsV[tind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;

					//now compute the intersecting cells for this motion (including ignoring the source footprint)
					get_2d_motion_cells(EnvNAVXYTHETAVSTEERCfg.FootprintPolygon, motionprimitiveV->at(mind).intermptV, &element_to_add.intersectingcellsV, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	/*
		#if DEBUG
					SBPL_FPRINTF(fDeb,
								"action action_index=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) endv=%3d"
								"cost=%4d (mprimID %3d: %3d %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
								vector_index,
								aind,
								EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).dX,
								EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).dY,
								EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).endtheta,
								EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).intermptV[0].theta * 180 / PI_CONST,
								EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).intermptV[EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).intermptV.size() - 1].theta * 180 / PI_CONST, EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).cost,
								motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
								motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta, motionprimitiveV->at(mind).endcell.v,
								(int)EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).interm3DcellsV.size(),
								(int)EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind).intersectingcellsV.size());
		#endif
	*/
					//add to the list of backward actions
					int possteer = (element_to_add.endsteer >= 0)?element_to_add.endsteer:element_to_add.endsteer+EnvNAVXYTHETAVSTEERCfg.numSteers;
					int postheta = (element_to_add.endtheta >= 0)?element_to_add.endtheta:element_to_add.endtheta+EnvNAVXYTHETAVSTEERCfg.NumThetaDirs;
					int targetindex = possteer*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs*EnvNAVXYTHETAVSTEERCfg.numV+element_to_add.endv*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+postheta;
					EnvNAVXYTHETAVSTEERActionIndex_t idx;
					idx.vector_index = vector_index;
					idx.aind = aind;
					EnvNAVXYTHETAVSTEERCfg.ActionsV->at(vector_index)->push_back(element_to_add);
					EnvNAVXYTHETAVSTEERCfg.PredActionsV[targetindex].push_back(idx);
				}

				if (maxnumofactions < numofactions) maxnumofactions = numofactions;
			}
		}
	}

	//at this point we don't allow nonuniform number of actions -> WHY???
	/*
	if (motionprimitiveV->size() != (size_t)(EnvNAVXYTHETALATCfg.NumThetaDirs * maxnumofactions)) {
		SBPL_ERROR("ERROR: nonuniform number of actions is not supported "
				"(maxnumofactions=%d while motprims=%d thetas=%d\n",
				maxnumofactions, (unsigned int)motionprimitiveV->size(), EnvNAVXYTHETALATCfg.NumThetaDirs);
		throw new SBPL_Exception();
	}
	*/
	
	//now compute replanning data
	ComputeReplanningData();

	SBPL_PRINTF("done pre-computing action data based on motion primitives\n");
}

void EnvironmentNAVXYTHETAVSTEER::InitializeEnvConfig(vector<SBPL_xythetavsteer_mprimitive>* motionprimitiveV)
{
	//aditional to configuration file initialization of EnvCfg if necessary
	
	/*
	//dXY dirs
	EnvNAVXYTHETALATCfg.dXY[0][0] = -1;
	EnvNAVXYTHETALATCfg.dXY[0][1] = -1;
	EnvNAVXYTHETALATCfg.dXY[1][0] = -1;
	EnvNAVXYTHETALATCfg.dXY[1][1] = 0;
	EnvNAVXYTHETALATCfg.dXY[2][0] = -1;
	EnvNAVXYTHETALATCfg.dXY[2][1] = 1;
	EnvNAVXYTHETALATCfg.dXY[3][0] = 0;
	EnvNAVXYTHETALATCfg.dXY[3][1] = -1;
	EnvNAVXYTHETALATCfg.dXY[4][0] = 0;
	EnvNAVXYTHETALATCfg.dXY[4][1] = 1;
	EnvNAVXYTHETALATCfg.dXY[5][0] = 1;
	EnvNAVXYTHETALATCfg.dXY[5][1] = -1;
	EnvNAVXYTHETALATCfg.dXY[6][0] = 1;
	EnvNAVXYTHETALATCfg.dXY[6][1] = 0;
	EnvNAVXYTHETALATCfg.dXY[7][0] = 1;
	EnvNAVXYTHETALATCfg.dXY[7][1] = 1;
	*/

	sbpl_xy_theta_pt_t temppose;
	temppose.x = 0.0;
	temppose.y = 0.0;
	temppose.theta = 0.0;
	vector<sbpl_2Dcell_t> footprint;
	get_2d_footprint_cells(EnvNAVXYTHETAVSTEERCfg.FootprintPolygon, &footprint, temppose, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	SBPL_PRINTF("number of cells in footprint of the robot = %d\n", (unsigned int)footprint.size());

	for (vector<sbpl_2Dcell_t>::iterator it = footprint.begin(); it != footprint.end(); ++it) {
		SBPL_PRINTF("Footprint cell at (%d, %d)\n", it->x, it->y);
	}

#if DEBUG
	SBPL_FPRINTF(fDeb, "footprint cells (size=%d):\n", (int)footprint.size());
	for(int i = 0; i < (int) footprint.size(); i++)
	{
		SBPL_FPRINTF(fDeb, "%d %d (cont: %.3f %.3f)\n", footprint.at(i).x, footprint.at(i).y,
					DISCXY2CONT(footprint.at(i).x, EnvNAVXYTHETAVSTEERCfg.cellsize_m),
					DISCXY2CONT(footprint.at(i).y, EnvNAVXYTHETAVSTEERCfg.cellsize_m));
	}
#endif
	
	PrecomputeActionswithCompleteMotionPrimitive(motionprimitiveV);
}

bool EnvironmentNAVXYTHETAVSTEER::InitGeneral(vector<SBPL_xythetavsteer_mprimitive>* motionprimitiveV){
	//Initialize other parameters of the environment
	InitializeEnvConfig(motionprimitiveV);

	//initialize Environment
	InitializeEnvironment();

	//pre-compute heuristics
	ComputeHeuristicValues();

	return true;
}

void EnvironmentNAVXYTHETAVSTEER::InitializeEnvironment()
{
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry;

	int maxsize = EnvNAVXYTHETAVSTEERCfg.EnvWidth_c * EnvNAVXYTHETAVSTEERCfg.EnvHeight_c * EnvNAVXYTHETAVSTEERCfg.NumThetaDirs * EnvNAVXYTHETAVSTEERCfg.numV * EnvNAVXYTHETAVSTEERCfg.numSteers;

	SBPL_PRINTF("environment stores states in hashtable\n");

	//initialize the map from Data to StateID
	//Maximum hash table size
	EnvNAVXYTHETAVSTEER.HashTableSize = EnvNAVXYTHETAVSTEERCfg.NumThetaDirs * EnvNAVXYTHETAVSTEERCfg.EnvWidth_c * EnvNAVXYTHETAVSTEERCfg.EnvHeight_c * 8 * 8;/*EnvNAVXYTHETAVCfg.numV;*/ //should be power of two - REVIEW -> USE DYNAMIC PARAMETERS
	EnvNAVXYTHETAVSTEER.Data2StateIDHashTable = new vector<EnvNAVXYTHETAVSTEERHashEntry_t*> [EnvNAVXYTHETAVSTEER.HashTableSize];

	//initialize the map from StateID to Data
	EnvNAVXYTHETAVSTEER.StateID2DataTable.clear();

	//create start state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c,
										EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.StartSteer)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c,
												EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.StartSteer);
	}
	EnvNAVXYTHETAVSTEER.startstateid = HashEntry->stateID;

	//create goal state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c,
										EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.EndSteer)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c,
												EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.EndSteer);
	}
	EnvNAVXYTHETAVSTEER.goalstateid = HashEntry->stateID;

	//initialized
	EnvNAVXYTHETAVSTEER.bInitialized = true;
}

bool EnvironmentNAVXYTHETAVSTEER::ReadMotionPrimitives(FILE* fMotPrims){
	char sTemp[1024], sExpected[1024];
	float fTemp;
	int dTemp;
	int totalNumofActions = 0;

	SBPL_PRINTF("Reading in motion primitives...");

	//read in the resolution
	strcpy(sExpected, "resolution_m:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
	if (fabs(fTemp - EnvNAVXYTHETAVSTEERCfg.cellsize_m) > ERR_EPS) {
		SBPL_ERROR("ERROR: invalid resolution %f (instead of %f) in the dynamics file\n", fTemp,
				EnvNAVXYTHETAVSTEERCfg.cellsize_m);
		return false;
	}

	//read in the angular resolution
	strcpy(sExpected, "numberofangles:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%d", &dTemp) == 0) return false;
	if (dTemp != EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		return false;
	}
	
	//read in the velocity resolution
	strcpy(sExpected, "numberofvelocities:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%d", &dTemp) == 0) return false;
	if (dTemp != EnvNAVXYTHETAVSTEERCfg.numV) {
		SBPL_ERROR("ERROR: invalid velocity resolution %d velocities (instead of %d velocities) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVSTEERCfg.numV);
		return false;
	}
	
	//read in the velocity values
	strcpy(sExpected, "velocities:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	for(int i=0;i<EnvNAVXYTHETAVSTEERCfg.numV;i++){
		if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
		if (fTemp != EnvNAVXYTHETAVSTEERCfg.velocities[i]) {
			SBPL_ERROR("ERROR: invalid velocity value %f velocity (instead of %f velocity) in the motion primitives file\n",
					fTemp, EnvNAVXYTHETAVSTEERCfg.velocities[i]);
			return false;
		}
	}
	
	//read in the velocity resolution
	strcpy(sExpected, "numberofsteeringangles:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%d", &dTemp) == 0) return false;
	if (dTemp != EnvNAVXYTHETAVSTEERCfg.numSteers) {
		SBPL_ERROR("ERROR: invalid steering angles resolution %d (instead of %d steering angles) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVSTEERCfg.numSteers);
		return false;
	}

	//read in the total number of actions
	strcpy(sExpected, "totalnumberofprimitives:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%d", &totalNumofActions) == 0) {
		return false;
	}

	for (int i = 0; i < totalNumofActions; i++) {
		SBPL_xythetavsteer_mprimitive motprim;

		if (ReadSingleMotionPrimitive(&motprim, fMotPrims) == false) return false;

		EnvNAVXYTHETAVSTEERCfg.mprimV.push_back(motprim);
	}
	SBPL_PRINTF("done reading of motion primitives");

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::ReadSingleMotionPrimitive(SBPL_xythetavsteer_mprimitive* pMotPrim, FILE* fIn){
	char sTemp[1024];
	int dTemp;
	double fTemp;
	char sExpected[1024];
	int numofIntermPoses;

	//read in actionID
	strcpy(sExpected, "primID:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%d", &pMotPrim->motprimID) != 1) return false;

	//read in start cell
	strcpy(sExpected, "start_pose_disc:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading start x\n");
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading start y\n");
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading start angle\n");
		return false;
	}
	pMotPrim->start_theta_disc = dTemp;
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading start v\n");
		return false;
	}
	pMotPrim->start_v_disc = dTemp;
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading start v\n");
		return false;
	}
	pMotPrim->start_steer_disc = dTemp;

	//read in end cell
	strcpy(sExpected, "end_pose_disc:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}

	if (ReadSingleCell(&pMotPrim->endcell, fIn) == false) {
		SBPL_ERROR("ERROR: failed to read in endsearchpose\n");
		return false;
	}

	//read in action cost
	strcpy(sExpected, "additionalactioncostmult:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%lf", &fTemp) != 1) return false;
	pMotPrim->additionalactioncostmult = fTemp;

	//read in intermediate poses
	strcpy(sExpected, "intermediateposes:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%d", &numofIntermPoses) != 1) return false;
	for (int i = 0; i < numofIntermPoses; i++) {
		sbpl_xy_theta_v_steer_pt_t intermpose;
		if (ReadSinglePose(&intermpose, fIn) == false) {
			SBPL_ERROR("ERROR: failed to read in intermediate poses\n");
			return false;
		}
		pMotPrim->intermptV.push_back(intermpose);
	}

	//check that the last pose corresponds correctly to the last pose
	sbpl_xy_theta_v_steer_pt_t sourcepose;
	sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	sourcepose.theta = DiscTheta2ContNotUnif(pMotPrim->start_theta_disc, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	sourcepose.v = DiscV2Cont(pMotPrim->start_v_disc, EnvNAVXYTHETAVSTEERCfg.velocities);
	sourcepose.steer = DiscSteer2Cont(pMotPrim->start_steer_disc, EnvNAVXYTHETAVSTEERCfg.numSteers);
	double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
	double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
	double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;
	double mp_endv_ms = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v;
	double mp_endsteer_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].steer;
	int endx_disc = CONTXY2DISC(mp_endx_m, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int endy_disc = CONTXY2DISC(mp_endy_m, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int endtheta_disc = ContTheta2DiscNotUnif(mp_endtheta_rad, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	int endv_disc = ContV2Disc(mp_endv_ms, EnvNAVXYTHETAVSTEERCfg.velocities);
	int endsteer_disc = ContSteer2Disc(mp_endsteer_rad, EnvNAVXYTHETAVSTEERCfg.numSteers);
	if (endx_disc != pMotPrim->endcell.x || endy_disc != pMotPrim->endcell.y || endtheta_disc != pMotPrim->endcell.theta || endv_disc != pMotPrim->endcell.v || endsteer_disc != pMotPrim->endcell.steer) {
		SBPL_ERROR( "ERROR: incorrect primitive %d with startangle=%d, startv=%d and startsteer=%d "
				"last interm point %f %f %f %f %f does not match end pose %d %d %d %d %d\n",
				pMotPrim->motprimID, pMotPrim->start_theta_disc, pMotPrim->start_v_disc, pMotPrim->start_steer_disc,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].steer,
				pMotPrim->endcell.x, pMotPrim->endcell.y,
				pMotPrim->endcell.theta, pMotPrim->endcell.v, pMotPrim->endcell.steer);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::ReadSingleCell(sbpl_xy_theta_v_steer_cell_t* cell, FILE* fIn){
	char sTemp[60];

	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->x = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->y = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->theta = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->v = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->steer = atoi(sTemp);
	
	//normalize the angle
	cell->theta = NORMALIZEDISCTHETA(cell->theta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::ReadSinglePose(sbpl_xy_theta_v_steer_pt_t* pose, FILE* fIn){
	char sTemp[60];

	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->x = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->y = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->theta = round_two_decimal(atof(sTemp));
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->v = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->steer = atof(sTemp);

	pose->theta = normalizeAngle(pose->theta);
	pose->steer = pose->steer;

	return true;
}

/*void EnvironmentNAVXYTHETAV::AddAllOutcomes(unsigned int SourceX1, unsigned int SourceX2, unsigned int SourceX3,
									unsigned int SourceX4, CMDPACTION* action, int cost)
{
	EnvXXXHashEntry_t* OutHashEntry;
	int i;
	float CumProb = 0.0;

	//iterate over outcomes
	for (i = 0; i < 2; i++) {
		unsigned int newX1 = SourceX1 + i;
		unsigned int newX2 = SourceX2 + i;
		unsigned int newX3 = SourceX3 + i;
		unsigned int newX4 = SourceX4 + i;

		//add the outcome
		if ((OutHashEntry = GetHashEntry(newX1, newX2, newX3, newX4)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX1, newX2, newX3, newX4);
		}
		float Prob = 0.5; //probability of the outcome
		action->AddOutcome(OutHashEntry->stateID, cost, Prob);
		CumProb += Prob;

	} //while

	if (CumProb != 1.0) {
		SBPL_ERROR("ERROR in EnvXXX... function: prob. of all action outcomes=%f\n", CumProb);
		throw new SBPL_Exception();
	}
}*/

//------------------------------------------------------------------------------

//------------------------------Heuristic computation--------------------------

void EnvironmentNAVXYTHETAVSTEER::EnsureHeuristicsUpdated(bool bGoalHeuristics){
	if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
		grid2Dsearchfromstart->search(EnvNAVXYTHETAVSTEERCfg.Grid2D, EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c,
									EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeStartHeuristics = false;
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVSTEERCfg.EndX_c,
																				EnvNAVXYTHETAVSTEERCfg.EndY_c)/EnvNAVXYTHETAVSTEERCfg.velocities.at(EnvNAVXYTHETAVSTEERCfg.velocities.size()-1)));
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.EndX_c,
																				EnvNAVXYTHETAVCfg.EndY_c)));*/
	}

	if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
		grid2Dsearchfromgoal->search(EnvNAVXYTHETAVSTEERCfg.Grid2D, EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c,
									EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeGoalHeuristics = false;
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVSTEERCfg.StartX_c,
																				EnvNAVXYTHETAVSTEERCfg.StartY_c)/EnvNAVXYTHETAVSTEERCfg.velocities.at(EnvNAVXYTHETAVSTEERCfg.velocities.size()-1)));
		
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.StartX_c,
																				EnvNAVXYTHETAVCfg.StartY_c)));*/
	}
}

void EnvironmentNAVXYTHETAVSTEER::ComputeHeuristicValues()
{
	//whatever necessary pre-computation of heuristic values is done here
	SBPL_PRINTF("Precomputing heuristics\n");
	
	//allocated 2D grid searches
	grid2Dsearchfromstart = new SBPL2DGridSearch(EnvNAVXYTHETAVSTEERCfg.EnvWidth_c, EnvNAVXYTHETAVSTEERCfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	grid2Dsearchfromgoal = new SBPL2DGridSearch(EnvNAVXYTHETAVSTEERCfg.EnvWidth_c, EnvNAVXYTHETAVSTEERCfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVSTEERCfg.cellsize_m);

	//set OPEN type to sliding buckets
	grid2Dsearchfromstart->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	grid2Dsearchfromgoal->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	
	SBPL_PRINTF("done\n");
}

void EnvironmentNAVXYTHETAVSTEER::PrintHeuristicValues(){
	#ifndef ROS
		const char* heur = "heur.txt";
	#endif
	FILE* fHeur = SBPL_FOPEN(heur, "w");
	if (fHeur == NULL) {
		SBPL_ERROR("ERROR: could not open debug file to write heuristic\n");
		throw new SBPL_Exception();
	}
	SBPL2DGridSearch* grid2Dsearch = NULL;
	
	if(grid2Dsearchfromgoal != NULL){
		grid2Dsearch = grid2Dsearchfromgoal;
			SBPL_FPRINTF(fHeur, "goal heuristics:\n");
	}
	else if(grid2Dsearchfromstart != NULL){
		grid2Dsearch = grid2Dsearchfromstart;
		SBPL_FPRINTF(fHeur, "start heuristics:\n");
	}
	else{
		SBPL_FCLOSE(fHeur);
		return;
	}

	for (int y = 0; y < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c; y++) {
		for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
			if (grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y) < INFINITECOST)
				SBPL_FPRINTF(fHeur, "%5d ", grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y));
			else
				SBPL_FPRINTF(fHeur, "XXXXX ");
		}
		SBPL_FPRINTF(fHeur, "\n");
	}
	
	SBPL_FCLOSE(fHeur);
}

//-----------interface with outside functions-----------------------------------

bool EnvironmentNAVXYTHETAVSTEER::InitializeEnv(const char* sEnvFile, const vector<sbpl_2Dpt_t>& perimeterptsV, const char* sMotPrimFile)
{
	EnvNAVXYTHETAVSTEERCfg.FootprintPolygon = perimeterptsV;

	FILE* fCfg = fopen(sEnvFile, "r");
	if (fCfg == NULL) {
		SBPL_ERROR("ERROR: unable to open %s\n", sEnvFile);
		throw new SBPL_Exception();
	}
	ReadConfiguration(fCfg);
	fclose(fCfg);

	if (sMotPrimFile != NULL) {
		FILE* fMotPrim = fopen(sMotPrimFile, "r");
		if (fMotPrim == NULL) {
			SBPL_ERROR("ERROR: unable to open %s\n", sMotPrimFile);
			throw new SBPL_Exception();
		}
		if (ReadMotionPrimitives(fMotPrim) == false) {
			SBPL_ERROR("ERROR: failed to read in motion primitive file\n");
			throw new SBPL_Exception();
		}
		InitGeneral(&EnvNAVXYTHETAVSTEERCfg.mprimV);
		fclose(fMotPrim);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

	SBPL_PRINTF("size of env: %d by %d\n", EnvNAVXYTHETAVSTEERCfg.EnvWidth_c, EnvNAVXYTHETAVSTEERCfg.EnvHeight_c);

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::InitializeEnv(const char* sEnvFile){
	return false;
}

bool EnvironmentNAVXYTHETAVSTEER::SetEnvParameter(const char* parameter, int value){
	if (EnvNAVXYTHETAVSTEER.bInitialized == true) {
		SBPL_ERROR("ERROR: all parameters must be set before initialization of the environment\n");
		return false;
	}

	SBPL_PRINTF("setting parameter %s to %d\n", parameter, value);

	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh = (unsigned char)value;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh = value;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVSTEERCfg.obsthresh = (unsigned char)value;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		return false;
	}

	return true;
}

int EnvironmentNAVXYTHETAVSTEER::GetEnvParameter(const char* parameter)
{
	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		return (int)EnvNAVXYTHETAVSTEERCfg.obsthresh;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		throw new SBPL_Exception();
	}
}

bool EnvironmentNAVXYTHETAVSTEER::InitializeMDPCfg(MDPConfig *MDPCfg)
{
	//initialize MDPCfg with the start and goal ids
	MDPCfg->goalstateid = EnvNAVXYTHETAVSTEER.goalstateid;
	MDPCfg->startstateid = EnvNAVXYTHETAVSTEER.startstateid;

	return true;
}

int EnvironmentNAVXYTHETAVSTEER::GetFromToHeuristic(int FromStateID, int ToStateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if(FromStateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size() || ToStateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size()){
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVSTEER... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	//get X, Y for the state
	EnvNAVXYTHETAVSTEERHashEntry_t* FromHashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[FromStateID];
	EnvNAVXYTHETAVSTEERHashEntry_t* ToHashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[ToStateID];

	//TODO - check if one of the gridsearches already computed and then use it.

	return (int)(NAVXYTHETAVSTEER_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y)/EnvNAVXYTHETAVSTEERCfg.velocities.at(EnvNAVXYTHETAVSTEERCfg.velocities.size()-1));
	//return (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y));
}

int EnvironmentNAVXYTHETAVSTEER::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[stateID];
	//computes distances from start state that is grid2D, so it is EndX_c EndY_c
	int h2D = grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y); 
	int hEuclid = (int)(NAVXYTHETAVSTEER_COSTMULT_MTOMM * EuclideanDistance_m(HashEntry->x, HashEntry->y,
																		EnvNAVXYTHETAVSTEERCfg.EndX_c,
																		EnvNAVXYTHETAVSTEERCfg.EndY_c));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETAVSTEERCfg.velocities.at(EnvNAVXYTHETAVSTEERCfg.velocities.size()-1));
	//return (int)(((double)__max(h2D, hEuclid)));
}

int EnvironmentNAVXYTHETAVSTEER::GetStartHeuristic(int stateID)
{
	
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[stateID];
	int h2D = grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y);
	int hEuclid = (int)(NAVXYTHETAVSTEER_COSTMULT_MTOMM * EuclideanDistance_m(EnvNAVXYTHETAVSTEERCfg.StartX_c,
																		EnvNAVXYTHETAVSTEERCfg.StartY_c, HashEntry->x,
																		HashEntry->y));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETAVSTEERCfg.velocities.at(EnvNAVXYTHETAVSTEERCfg.velocities.size()-1));
	//return (int)(((double)__max(h2D, hEuclid)));
}

void EnvironmentNAVXYTHETAVSTEER::SetAllActionsandAllOutcomes(CMDPSTATE* state){
	int cost;
	
#if DEBUG
	if(state->StateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size())
	{
		SBPL_ERROR("ERROR in Env... function: stateID illegal\n");
		throw new SBPL_Exception();
	}

	if((int)state->Actions.size() != 0)
	{
		SBPL_ERROR("ERROR in Env_setAllActionsandAllOutcomes: actions already exist for the state\n");
		throw new SBPL_Exception();
	}
#endif

	//goal state should be absorbing
	if (state->StateID == EnvNAVXYTHETAVSTEER.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[state->StateID];
	
	int vector_index = HashEntry->steer*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs*EnvNAVXYTHETAVSTEERCfg.numV + HashEntry->v*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+HashEntry->theta;

	//iterate through actions
	for (int aind = 0; aind < EnvNAVXYTHETAVSTEERCfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVSTEERAction_t* nav5daction = &EnvNAVXYTHETAVSTEERCfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav5daction->dX;
		int newY = HashEntry->y + nav5daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav5daction->endtheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		int newV = nav5daction->endv;
		int newSteer = nav5daction->endsteer;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->steer, nav5daction);
		if (cost >= INFINITECOST) continue;

		//add the action
		CMDPACTION* action = state->AddAction(aind);

#if TIME_DEBUG
		clock_t currenttime = clock();
#endif

		EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV, newSteer)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV, newSteer);
		}
		action->AddOutcome(OutHashEntry->stateID, cost, 1.0);

#if TIME_DEBUG
		time3_addallout += clock()-currenttime;
#endif
	}
}

void EnvironmentNAVXYTHETAVSTEER::SetAllPreds(CMDPSTATE* state)
{
	//implement this if the planner needs access to predecessors

	SBPL_ERROR("ERROR in EnvNAVXYTHETAVSTEER... function: SetAllPreds is undefined\n");
	throw new SBPL_Exception();
}

void EnvironmentNAVXYTHETAVSTEER::GetSuccs(int SourceStateID, vector<int>* SuccIDV, vector<int>* CostV)
{
	GetSuccs(SourceStateID, SuccIDV, CostV, NULL);
}

void EnvironmentNAVXYTHETAVSTEER::GetSuccs(int sourceStateID, vector<int>* succIDV, vector<int>* costV, vector<EnvNAVXYTHETAVSTEERAction_t*>* actionindV){
	int aind;

#if TIME_DEBUG
	clock_t currenttime = clock();
#endif

	//clear the successor array
	succIDV->clear();
	costV->clear();
	//succIDV->reserve(EnvNAVXYTHETAVCfg.actionwidth);
	//costV->reserve(EnvNAVXYTHETALATCfg.actionwidth);
	if (actionindV != NULL) {
		actionindV->clear();
		//actionindV->reserve(EnvNAVXYTHETALATCfg.actionwidth);
	}

	//goal state should be absorbing
	if (sourceStateID == EnvNAVXYTHETAVSTEER.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[sourceStateID];
	
	int vector_index = (unsigned int)HashEntry->steer*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs*EnvNAVXYTHETAVSTEERCfg.numV + (unsigned int)HashEntry->v*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//iterate through actions
	for (aind = 0; aind < EnvNAVXYTHETAVSTEERCfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVSTEERAction_t* nav5daction = &EnvNAVXYTHETAVSTEERCfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav5daction->dX;
		int newY = HashEntry->y + nav5daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav5daction->endtheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		int newV = nav5daction->endv;
		int newSteer = nav5daction->endsteer;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		int cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->steer, nav5daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV, newSteer)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV, newSteer);
		}

		succIDV->push_back(OutHashEntry->stateID);
		costV->push_back(cost);
		if (actionindV != NULL) actionindV->push_back(nav5daction);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentNAVXYTHETAVSTEER::GetPreds(int TargetStateID, vector<int>* PredIDV, vector<int>* CostV)
{
	//TODO- to support tolerance, need:
	// a) generate preds for goal state based on all possible goal state variable settings,
	// b) change goal check condition in gethashentry 
	// c) change getpredsofchangedcells and getsuccsofchangedcells functions

	int aind;

#if TIME_DEBUG
	clock_t currenttime = clock();
#endif

	//get X, Y for the state
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[TargetStateID];
	int targetindex = (unsigned int)HashEntry->steer*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs*EnvNAVXYTHETAVSTEERCfg.numV + (unsigned int)HashEntry->v*EnvNAVXYTHETAVSTEERCfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//clear the preds array
	PredIDV->clear();
	CostV->clear();
	PredIDV->reserve(EnvNAVXYTHETAVSTEERCfg.PredActionsV[targetindex].size());
	CostV->reserve(EnvNAVXYTHETAVSTEERCfg.PredActionsV[targetindex].size());

	//iterate through actions
	vector<EnvNAVXYTHETAVSTEERActionIndex_t>* actionsV = &EnvNAVXYTHETAVSTEERCfg.PredActionsV[targetindex];
	for (aind = 0; aind < (int)actionsV->size(); aind++) {

		EnvNAVXYTHETAVSTEERAction_t* nav5daction = &EnvNAVXYTHETAVSTEERCfg.ActionsV->at(actionsV->at(aind).vector_index)->at(actionsV->at(aind).aind);

		int predX = HashEntry->x - nav5daction->dX;
		int predY = HashEntry->y - nav5daction->dY;
		int predTheta = nav5daction->starttheta;
		int predV = nav5daction->startv;
		int predSteer = nav5daction->startsteer;

		//skip the invalid cells
		if (!IsValidCell(predX, predY)) continue;

		//get cost
		int cost = GetActionCost(predX, predY, predTheta, predV, predSteer, nav5daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(predX, predY, predTheta, predV, predSteer)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(predX, predY, predTheta, predV, predSteer);
		}

		PredIDV->push_back(OutHashEntry->stateID);
		CostV->push_back(cost);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

int EnvironmentNAVXYTHETAVSTEER::SizeofCreatedEnv()
{
	return (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size();
}

void EnvironmentNAVXYTHETAVSTEER::PrintState(int stateID, bool bVerbose, FILE* fOut /*=NULL*/)
{
#if DEBUG
	if(stateID >= (int)EnvNAVXYTHETAVSTEER.StateID2DataTable.size())
	{
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVSTEER... function: stateID illegal (2)\n");
		throw new SBPL_Exception();
	}
#endif

	if (fOut == NULL) fOut = stdout;

	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[stateID];

	if (stateID == EnvNAVXYTHETAVSTEER.goalstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a goal state\n");
	}
	
	if (stateID == EnvNAVXYTHETAVSTEER.startstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a start state\n");
	}

	SBPL_FPRINTF(fOut, "x=%d y=%d theta=%d v=%d steer=%d\n", HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->steer);
}

void EnvironmentNAVXYTHETAVSTEER::PrintEnv_Config(FILE* fOut)
{
	if(fOut != NULL){
		fprintf(fOut, "discretization(cells): %d %d\n", EnvNAVXYTHETAVSTEERCfg.EnvWidth_c, EnvNAVXYTHETAVSTEERCfg.EnvHeight_c);
		fprintf(fOut, "NumThetaDirs: %d\n", EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		fprintf(fOut, "NumV: %d\n", EnvNAVXYTHETAVSTEERCfg.numV);
		fprintf(fOut, "Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVSTEERCfg.numV; i++)
			fprintf(fOut, " %f", EnvNAVXYTHETAVSTEERCfg.velocities.at(i));
		fprintf(fOut, "NumSteers: %d\n", EnvNAVXYTHETAVSTEERCfg.numSteers);
		fprintf(fOut, "\n");
		fprintf(fOut, "obsthresh: %d\n", EnvNAVXYTHETAVSTEERCfg.obsthresh);
		fprintf(fOut, "cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "cellsize(meters): %f\n", EnvNAVXYTHETAVSTEERCfg.cellsize_m);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "start(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.velocities), DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.StartSteer, EnvNAVXYTHETAVSTEERCfg.numSteers));
		fprintf(fOut, "end(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.velocities), DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.numSteers));
		fprintf(fOut, "environment:\n");
		
		for (int y = EnvNAVXYTHETAVSTEERCfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
				fprintf(fOut, "%d ", EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y]);
			}
			
			fprintf(fOut, "\n");
		}
	}
	else{
		printf("discretization(cells): %d %d\n", EnvNAVXYTHETAVSTEERCfg.EnvWidth_c, EnvNAVXYTHETAVSTEERCfg.EnvHeight_c);
		printf("NumThetaDirs: %d\n", EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
		printf("NumV: %d\n", EnvNAVXYTHETAVSTEERCfg.numV);
		printf("Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVSTEERCfg.numV; i++)
			printf(" %f", EnvNAVXYTHETAVSTEERCfg.velocities.at(i));
		fprintf(fOut, "NumSteers: %d\n", EnvNAVXYTHETAVSTEERCfg.numSteers);
		printf("\n");
		printf("obsthresh: %d\n", EnvNAVXYTHETAVSTEERCfg.obsthresh);
		printf("cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh);
		printf("cellsize(meters): %f\n", EnvNAVXYTHETAVSTEERCfg.cellsize_m);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh);
		printf("start(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.velocities), DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.StartSteer, EnvNAVXYTHETAVSTEERCfg.numSteers));
		printf("end(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.velocities), DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.numSteers));
		printf("environment:\n");
		
		for (int y = EnvNAVXYTHETAVSTEERCfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
				printf("%d ", EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y]);
			}
			
			printf("\n");
		}
	}
}

bool EnvironmentNAVXYTHETAVSTEER::InitializeEnv(int width, int height, int numthetadirs, int numv, int numsteers, vector<double> velocities,
							const unsigned char* mapdata,
							double startx, double starty, double starttheta, double startv, double startsteer,
							double goalx, double goaly, double goaltheta, double goalv, double goalsteer,
							const vector<sbpl_2Dpt_t>& perimeterptsV, double cellsize_m,
							unsigned char obsthresh, unsigned char cost_inscribed_thresh,
							unsigned char cost_possibly_circumscribed_thresh, const char* sMotPrimFile){
	SBPL_PRINTF("env: initialize with width=%d height=%d start=%.3f %.3f %.3f %.3f %.3f "
				"goalx=%.3f %.3f %.3f %.3f %.3f cellsize=%.3f numthetadirs=%d numv=%d numsteers=%d, obsthresh=%d cost_inscribed_thresh=%d cost_possibly_circumscribed_thresh=%d\n",
				width, height, startx, starty, starttheta, startv, startsteer, goalx, goaly, goaltheta, goalv, goalsteer, cellsize_m, numthetadirs, numv, numsteers,
				obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh);

	SBPL_PRINTF("NOTE: goaltol parameters currently unused\n");

	SBPL_PRINTF("perimeter has size=%d\n", (unsigned int)perimeterptsV.size());

	for (int i = 0; i < (int)perimeterptsV.size(); i++) {
		SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, perimeterptsV.at(i).x, perimeterptsV.at(i).y);
	}

	EnvNAVXYTHETAVSTEERCfg.obsthresh = obsthresh;
	EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh = cost_inscribed_thresh;
	EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh = cost_possibly_circumscribed_thresh;
	EnvNAVXYTHETAVSTEERCfg.NumThetaDirs = numthetadirs;
	EnvNAVXYTHETAVSTEERCfg.numV = numv;
	EnvNAVXYTHETAVSTEERCfg.numSteers = numsteers;
	
	EnvNAVXYTHETAVSTEERCfg.velocities.clear();
	for(int i=0;i<numv;i++){
		EnvNAVXYTHETAVSTEERCfg.velocities.push_back(velocities.at(i));
	}
	
	//TODO - need to set the tolerance as well

	SetConfiguration(width, height, mapdata, CONTXY2DISC(startx, cellsize_m), CONTXY2DISC(starty, cellsize_m),
					ContTheta2DiscNotUnif(starttheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs), ContV2Disc(startv, EnvNAVXYTHETAVSTEERCfg.velocities), ContSteer2Disc(startsteer, EnvNAVXYTHETAVSTEERCfg.numSteers),
					CONTXY2DISC(goalx, cellsize_m), CONTXY2DISC(goaly, cellsize_m), ContTheta2DiscNotUnif(goaltheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs),
					ContV2Disc(goalv, EnvNAVXYTHETAVSTEERCfg.velocities), ContSteer2Disc(goalsteer, EnvNAVXYTHETAVSTEERCfg.numSteers), cellsize_m, perimeterptsV);

	if (sMotPrimFile != NULL) {
		FILE* fMotPrim = fopen(sMotPrimFile, "r");
		if (fMotPrim == NULL) {
			SBPL_ERROR("ERROR: unable to open %s\n", sMotPrimFile);
			throw new SBPL_Exception();
		}

		if (ReadMotionPrimitives(fMotPrim) == false) {
			SBPL_ERROR("ERROR: failed to read in motion primitive file\n");
			throw new SBPL_Exception();
		}
		fclose(fMotPrim);
	}

	if (EnvNAVXYTHETAVSTEERCfg.mprimV.size() != 0) {
		InitGeneral(&EnvNAVXYTHETAVSTEERCfg.mprimV);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::UpdateCost(int x, int y, unsigned char newcost){
	EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] = newcost;

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

bool EnvironmentNAVXYTHETAVSTEER::SetMap(const unsigned char* mapdata){
	for (int xind = 0; xind < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; xind++) {
		for (int yind = 0; yind < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c; yind++) {
			EnvNAVXYTHETAVSTEERCfg.Grid2D[xind][EnvNAVXYTHETAVSTEERCfg.EnvHeight_c-1-yind] = mapdata[xind + yind * EnvNAVXYTHETAVSTEERCfg.EnvWidth_c];
		}
	}

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

void EnvironmentNAVXYTHETAVSTEER::GetPredsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV, vector<int> *preds_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_steer_cell_t affectedcell;
	EnvNAVXYTHETAVSTEERHashEntry_t* affectedHashEntry;

	//increment iteration for processing savings
	iteration++;

	for (int i = 0; i < (int)changedcellsV->size(); i++) {
		cell = changedcellsV->at(i);

		//now iterate over all states that could potentially be affected
		for (int sind = 0; sind < (int)affectedpredstatesV.size(); sind++) {
			affectedcell = affectedpredstatesV.at(sind);

			//translate to correct for the offset
			affectedcell.x = affectedcell.x + cell.x;
			affectedcell.y = affectedcell.y + cell.y;

			//insert only if it was actually generated
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v, affectedcell.steer);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

void EnvironmentNAVXYTHETAVSTEER::GetSuccsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV,
										vector<int> *succs_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_steer_cell_t affectedcell;
	EnvNAVXYTHETAVSTEERHashEntry_t* affectedHashEntry;

	SBPL_ERROR("ERROR: getsuccs is not supported currently\n");
	throw new SBPL_Exception();

	//increment iteration for processing savings
	iteration++;

	//TODO - check
	for (int i = 0; i < (int)changedcellsV->size(); i++) {
		cell = changedcellsV->at(i);

		//now iterate over all states that could potentially be affected
		for (int sind = 0; sind < (int)affectedsuccstatesV.size(); sind++) {
			affectedcell = affectedsuccstatesV.at(sind);

			//translate to correct for the offset
			affectedcell.x = affectedcell.x + cell.x;
			affectedcell.y = affectedcell.y + cell.y;

			//insert only if it was actually generated
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v, affectedcell.steer);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

bool EnvironmentNAVXYTHETAVSTEER::IsObstacle(int x, int y){
	return (EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] >= EnvNAVXYTHETAVSTEERCfg.obsthresh);
}

bool EnvironmentNAVXYTHETAVSTEER::IsValidConfiguration(int x, int y, int theta, int v, int steer){
	vector<sbpl_2Dcell_t> footprint;
	sbpl_xy_theta_v_steer_pt_t pose;

	//compute continuous pose
	pose.x = DISCXY2CONT(x, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	pose.y = DISCXY2CONT(y, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	pose.theta = DiscTheta2ContNotUnif(theta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	pose.v = DiscV2Cont(v, EnvNAVXYTHETAVSTEERCfg.velocities);
	pose.steer = DiscSteer2Cont(steer, EnvNAVXYTHETAVSTEERCfg.numSteers);

	//compute footprint cells
	get_2d_footprint_cells(EnvNAVXYTHETAVSTEERCfg.FootprintPolygon, &footprint, pose, EnvNAVXYTHETAVSTEERCfg.cellsize_m);

	//iterate over all footprint cells
	for (int find = 0; find < (int)footprint.size(); find++) {
		int x = footprint.at(find).x;
		int y = footprint.at(find).y;

		if (x < 0 || x >= EnvNAVXYTHETAVSTEERCfg.EnvWidth_c || y < 0 || y >= EnvNAVXYTHETAVSTEERCfg.EnvHeight_c ||
			EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] >= EnvNAVXYTHETAVSTEERCfg.obsthresh)
		{
			return false;
		}
	}
	
	if(v < 0 || v >= EnvNAVXYTHETAVSTEERCfg.numV){
		return false;
	}
	
	if(steer < 0 || steer >= EnvNAVXYTHETAVSTEERCfg.numSteers){
		return false;
	}

	return true;
}

void EnvironmentNAVXYTHETAVSTEER::GetEnvParams(int *size_x, int *size_y, int * numthetadirs, int * numv, int * numsteers, vector<double> * velocities,
							double* startx, double* starty, double* starttheta, double * startv, double * startsteer,
							double* goalx, double* goaly, double* goaltheta, double * goalv, double * goalsteer, double* cellsize_m,
							unsigned char* cost_inscribed_thresh, unsigned char* cost_possibly_circumscribed_thresh,
							unsigned char* obsthresh, vector<SBPL_xythetavsteer_mprimitive>* motionprimitiveV){
	
	*size_x = EnvNAVXYTHETAVSTEERCfg.EnvWidth_c;
	*size_y = EnvNAVXYTHETAVSTEERCfg.EnvHeight_c;

	*numthetadirs = EnvNAVXYTHETAVSTEERCfg.NumThetaDirs;
	*numv = EnvNAVXYTHETAVSTEERCfg.numV;
	for(int i;i<EnvNAVXYTHETAVSTEERCfg.velocities.size();i++){
		velocities->push_back(EnvNAVXYTHETAVSTEERCfg.velocities.at(i));
	}
	*numsteers = EnvNAVXYTHETAVSTEERCfg.numSteers;
	
	*startx = DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	*starty = DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.StartY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	*starttheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.StartTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	*startv = DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.StartV, EnvNAVXYTHETAVSTEERCfg.velocities);
	*startsteer = DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.StartSteer, EnvNAVXYTHETAVSTEERCfg.numSteers);
	*goalx = DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	*goaly = DISCXY2CONT(EnvNAVXYTHETAVSTEERCfg.EndY_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	*goaltheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVSTEERCfg.EndTheta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	*goalv = DiscV2Cont(EnvNAVXYTHETAVSTEERCfg.EndV, EnvNAVXYTHETAVSTEERCfg.velocities);
	*goalsteer = DiscSteer2Cont(EnvNAVXYTHETAVSTEERCfg.EndSteer, EnvNAVXYTHETAVSTEERCfg.numSteers);

	*cellsize_m = EnvNAVXYTHETAVSTEERCfg.cellsize_m;
	
	*cost_inscribed_thresh = EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh;
	*cost_possibly_circumscribed_thresh = EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh;
	*obsthresh = EnvNAVXYTHETAVSTEERCfg.obsthresh;

	*motionprimitiveV = EnvNAVXYTHETAVSTEERCfg.mprimV;
}

const EnvNAVXYTHETAVSTEERConfig_t* EnvironmentNAVXYTHETAVSTEER::GetEnvNavConfig(){
	return &EnvNAVXYTHETAVSTEERCfg;
}

void EnvironmentNAVXYTHETAVSTEER::PrintTimeStat(FILE* fOut){
	/*TO BE IMPLEMENTED*/
}

unsigned char EnvironmentNAVXYTHETAVSTEER::GetMapCost(int x, int y){
	return EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y];
}

bool EnvironmentNAVXYTHETAVSTEER::IsWithinMapCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c);
}

bool EnvironmentNAVXYTHETAVSTEER::PoseContToDisc(double px, double py, double pth, double pv, double pst, int &ix, int &iy, int &ith, int &iv, int &ist) const{
	ix = CONTXY2DISC(px, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	iy = CONTXY2DISC(py, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	ith = ContTheta2DiscNotUnif(pth, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	iv = ContV2Disc(pv, EnvNAVXYTHETAVSTEERCfg.velocities);
	ist = ContSteer2Disc(pst, EnvNAVXYTHETAVSTEERCfg.numV);
	
	return (pth >= -2 * PI_CONST) && (pth <= 2 * PI_CONST) && (ix >= 0 && ix < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVSTEERCfg.numV) && (ist >= 0) && (ist < EnvNAVXYTHETAVSTEERCfg.numSteers);
}

bool EnvironmentNAVXYTHETAVSTEER::PoseDiscToCont(int ix, int iy, int ith, int iv, int ist, double &px, double &py, double &pth, double &pv, double &pst) const{
	px = DISCXY2CONT(ix, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	py = DISCXY2CONT(iy, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	pth = normalizeAngle(DiscTheta2ContNotUnif(ith, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs));
	pv = DiscV2Cont(iv, EnvNAVXYTHETAVSTEERCfg.velocities);
	pst = DiscSteer2Cont(ist, EnvNAVXYTHETAVSTEERCfg.numSteers);
	
	return (ith >= 0) && (ith < EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) && (ix >= 0 && ix < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVSTEERCfg.numV) && (ist >= 0) && (ist < EnvNAVXYTHETAVSTEERCfg.numSteers);
}

int EnvironmentNAVXYTHETAVSTEER::SetStart(double x, double y, double theta, double v, double steer){
	int startx = CONTXY2DISC(x, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int starty = CONTXY2DISC(y, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int starttheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	int startv = ContV2Disc(v, EnvNAVXYTHETAVSTEERCfg.velocities);
	int startsteer = ContSteer2Disc(steer, EnvNAVXYTHETAVSTEERCfg.numSteers);

	if (!IsWithinMapCell(startx, starty)) {
		SBPL_ERROR("ERROR: trying to set a start cell %d %d that is outside of map\n", startx, starty);
		return -1;
	}

	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d)\n", x, y, theta, v, steer, startx, starty, starttheta, startv, startsteer);

	if (!IsValidConfiguration(startx, starty, starttheta, startv, startsteer)) {
		SBPL_PRINTF("WARNING: start configuration %d %d %d %d %d is invalid\n", startx, starty, starttheta, startv, startsteer);
	}

	EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(startx, starty, starttheta, startv, startsteer)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(startx, starty, starttheta, startv, startsteer);
	}

	if (EnvNAVXYTHETAVSTEER.startstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true;
		//because termination condition can be not all states TODO - make it dependent on term. condition
		bNeedtoRecomputeGoalHeuristics = true; 
	}

	//set start
	EnvNAVXYTHETAVSTEER.startstateid = OutHashEntry->stateID;
	EnvNAVXYTHETAVSTEERCfg.StartX_c = startx;
	EnvNAVXYTHETAVSTEERCfg.StartY_c = starty;
	EnvNAVXYTHETAVSTEERCfg.StartTheta = starttheta;
	EnvNAVXYTHETAVSTEERCfg.StartV = startv;
	EnvNAVXYTHETAVSTEERCfg.StartSteer = startsteer;

	return EnvNAVXYTHETAVSTEER.startstateid;
}

int EnvironmentNAVXYTHETAVSTEER::SetGoal(double x, double y, double theta, double v, double steer){
	int goalx = CONTXY2DISC(x, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int goaly = CONTXY2DISC(y, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
	int goaltheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs);
	int goalv = ContV2Disc(v, EnvNAVXYTHETAVSTEERCfg.velocities);
	int goalsteer = ContSteer2Disc(steer, EnvNAVXYTHETAVSTEERCfg.numSteers);

	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d)\n", x, y, theta, v, steer, goalx, goaly, goaltheta, goalv, goalsteer);

	if (!IsWithinMapCell(goalx, goaly)) {
		SBPL_ERROR("ERROR: trying to set a goal cell %d %d that is outside of map\n", goalx, goaly);
		return -1;
	}

	if (!IsValidConfiguration(goalx, goaly, goaltheta, goalv, goalsteer)) {
		SBPL_PRINTF("WARNING: goal configuration is invalid\n");
	}

	EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(goalx, goaly, goaltheta, goalv, goalsteer)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(goalx, goaly, goaltheta, goalv, goalsteer);
	}

	//need to recompute start heuristics?
	if (EnvNAVXYTHETAVSTEER.goalstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true; //because termination condition may not plan all the way to the new goal
		bNeedtoRecomputeGoalHeuristics = true; //because goal heuristics change
	}

	EnvNAVXYTHETAVSTEER.goalstateid = OutHashEntry->stateID;

	EnvNAVXYTHETAVSTEERCfg.EndX_c = goalx;
	EnvNAVXYTHETAVSTEERCfg.EndY_c = goaly;
	EnvNAVXYTHETAVSTEERCfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVSTEERCfg.EndV = goalv;
	EnvNAVXYTHETAVSTEERCfg.EndSteer = goalsteer;

	return EnvNAVXYTHETAVSTEER.goalstateid;
}

void EnvironmentNAVXYTHETAVSTEER::GetCoordFromState(int stateID, int& x, int& y, int& theta, int& v, int& steer) const{
	EnvNAVXYTHETAVSTEERHashEntry_t* HashEntry = EnvNAVXYTHETAVSTEER.StateID2DataTable[stateID];
	x = HashEntry->x;
	y = HashEntry->y;
	theta = HashEntry->theta;
	v = HashEntry->v;
	steer = HashEntry->steer;
}

int EnvironmentNAVXYTHETAVSTEER::GetStateFromCoord(int x, int y, int theta, int v, int steer){
	EnvNAVXYTHETAVSTEERHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(x, y, theta, v, steer)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(x, y, theta, v, steer);
	}
	return OutHashEntry->stateID;
}

void EnvironmentNAVXYTHETAVSTEER::ConvertStateIDPathintoXYThetaVPath(vector<int>* stateIDPath, vector<sbpl_xy_theta_v_steer_pt_t>* xythetavPath){
	vector<EnvNAVXYTHETAVSTEERAction_t*> actionV;
	vector<int> CostV;
	vector<int> SuccIDV;
	int targetx_c, targety_c, targettheta_c, targetv_c, targetsteer_c;
	int sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c;

	SBPL_PRINTF("checks=%ld\n", checks);

	xythetavPath->clear();

#if DEBUG
	SBPL_FPRINTF(fDeb, "converting stateid path into coordinates:\n");
#endif

	for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
		int sourceID = stateIDPath->at(pind);
		int targetID = stateIDPath->at(pind + 1);

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c);
#endif

		//get successors and pick the target via the cheapest action
		SuccIDV.clear();
		CostV.clear();
		actionV.clear();
		GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

		int bestcost = INFINITECOST;
		int bestsind = -1;

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c);
		GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c, targetsteer_c);
		SBPL_FPRINTF(fDeb, "looking for %d %d %d %d %d -> %d %d %d %d %d (numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c,
					targetx_c, targety_c, targettheta_c, targetv_c, targetsteer_c, (int)SuccIDV.size());
#endif

		for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
#if DEBUG
			int x_c, y_c, theta_c, v_c steer_c;
			GetCoordFromState(SuccIDV[sind], x_c, y_c, theta_c, v_c, steer_c);
			SBPL_FPRINTF(fDeb, "succ: %d %d %d %d %d\n", x_c, y_c, theta_c, v_c, steer_c);
#endif

			if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
				bestcost = CostV[sind];
				bestsind = sind;
			}
		}
		if (bestsind == -1) {
			SBPL_ERROR("ERROR: successor not found for transition:\n");
			GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c);
			GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c, targetsteer_c);
			SBPL_PRINTF("%d %d %d %d %d -> %d %d %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c, targetx_c, targety_c,
						targettheta_c, targetv_c, targetsteer_c);
			throw new SBPL_Exception();
		}

		//now push in the actual path
		//int sourcex_c, sourcey_c, sourcetheta_c;
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourcesteer_c);
		double x, y;
		x = DISCXY2CONT(sourcex_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
		y = DISCXY2CONT(sourcey_c, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
		
		for (int ipind = 0; ipind < ((int)actionV[bestsind]->intermptV.size()) - 1; ipind++) {
			//translate appropriately
			sbpl_xy_theta_v_steer_pt_t intermpt = actionV[bestsind]->intermptV[ipind];
			intermpt.x += x;
			intermpt.y += y;

#if DEBUG
			int nx = CONTXY2DISC(intermpt.x, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
			int ny = CONTXY2DISC(intermpt.y, EnvNAVXYTHETAVSTEERCfg.cellsize_m);
			SBPL_FPRINTF(fDeb, "%.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d cost=%d) ",
				intermpt.x, intermpt.y, intermpt.theta, intermpt.v, intermpt.steer
				nx, ny, ContTheta2DiscNotUnif(intermpt.theta, EnvNAVXYTHETAVSTEERCfg.NumThetaDirs),
				ContV2Disc(intermpt.v, EnvNAVXYTHETAVSTEERCfg.velocities), ContSteer2Disc(intermpt.steerm EnvNAVXYTHETAVSTEERCfg.numSteers), EnvNAVXYTHETAVSTEERCfg.Grid2D[nx][ny]);
			if(ipind == 0) SBPL_FPRINTF(fDeb, "first (heur=%d)\n", GetStartHeuristic(sourceID));
			else SBPL_FPRINTF(fDeb, "\n");
#endif

			//store
			xythetavPath->push_back(intermpt);
		}
	}
}

/* Useful functions */
double EnvironmentNAVXYTHETAVSTEER::EuclideanDistance_m(int x1, int y1, int x2, int y2)
{
	int sqdist = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return EnvNAVXYTHETAVSTEERCfg.cellsize_m * sqrt((double)sqdist);
}

bool EnvironmentNAVXYTHETAVSTEER::IsValidCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c &&
			EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] < EnvNAVXYTHETAVSTEERCfg.obsthresh);
}

int EnvironmentNAVXYTHETAVSTEER::GetActionCost(int sourceX, int sourceY, int sourceTheta, int sourceV, int sourceSteer, EnvNAVXYTHETAVSTEERAction_t* action){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_steer_cell_t interm5Dcell;
	int i;

	//TODO - go over bounding box (minpt and maxpt) to test validity and skip
	//testing boundaries below, also order intersect cells so that the four
	//farthest pts go first

	if (!IsValidCell(sourceX, sourceY)) return INFINITECOST;
	if (!IsValidCell(sourceX + action->dX, sourceY + action->dY)) return INFINITECOST;

	if (EnvNAVXYTHETAVSTEERCfg.Grid2D[sourceX + action->dX][sourceY + action->dY] >=
		EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh) 
	{
		return INFINITECOST;
	}

	//need to iterate over discretized center cells and compute cost based on them
	unsigned char maxcellcost = 0;
	for (i = 0; i < (int)action->interm3DcellsV.size(); i++) {
		interm5Dcell = action->interm3DcellsV.at(i);
		interm5Dcell.x = interm5Dcell.x + sourceX;
		interm5Dcell.y = interm5Dcell.y + sourceY;

		if (interm5Dcell.x < 0 || interm5Dcell.x >= EnvNAVXYTHETAVSTEERCfg.EnvWidth_c || interm5Dcell.y < 0
			|| interm5Dcell.y >= EnvNAVXYTHETAVSTEERCfg.EnvHeight_c) return INFINITECOST;

		maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVSTEERCfg.Grid2D[interm5Dcell.x][interm5Dcell.y]);

		//check that the robot is NOT in the cell at which there is no valid orientation
		if (maxcellcost >= EnvNAVXYTHETAVSTEERCfg.cost_inscribed_thresh) return INFINITECOST;
	}

	//check collisions that for the particular footprint orientation along the action
	if (EnvNAVXYTHETAVSTEERCfg.FootprintPolygon.size() > 1 && (int)maxcellcost >=
		EnvNAVXYTHETAVSTEERCfg.cost_possibly_circumscribed_thresh)
	{
		checks++;

		for (i = 0; i < (int)action->intersectingcellsV.size(); i++) {
			//get the cell in the map
			cell = action->intersectingcellsV.at(i);
			cell.x = cell.x + sourceX;
			cell.y = cell.y + sourceY;

			//check validity
			if (!IsValidCell(cell.x, cell.y)) return INFINITECOST;

			//if(EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y] > currentmaxcost)
			////cost computation changed: cost = max(cost of centers of the
			//robot along action)
			//	currentmaxcost = EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y];
			//	//intersecting cells are only used for collision checking
		}
	}

	//to ensure consistency of h2D:
	maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVSTEERCfg.Grid2D[sourceX][sourceY]);
	int currentmaxcost =
			(int)__max(maxcellcost, EnvNAVXYTHETAVSTEERCfg.Grid2D[sourceX + action->dX][sourceY + action->dY]);
			
	return action->cost * (currentmaxcost + 1); //use cell cost as multiplicative factor
}

void EnvironmentNAVXYTHETAVSTEER::SetConfiguration(int width, int height, const unsigned char * mapdata, int startx, int starty,
								int starttheta, int startv, int startsteer, int goalx, int goaly, int goaltheta, int goalv, int goalsteer,
							double cellsize_m, const vector<sbpl_2Dpt_t>& perimeterptsV){
	EnvNAVXYTHETAVSTEERCfg.EnvWidth_c = width;
	EnvNAVXYTHETAVSTEERCfg.EnvHeight_c = height;
	EnvNAVXYTHETAVSTEERCfg.StartX_c = startx;
	EnvNAVXYTHETAVSTEERCfg.StartY_c = starty;
	EnvNAVXYTHETAVSTEERCfg.StartTheta = starttheta;
	EnvNAVXYTHETAVSTEERCfg.StartSteer = startsteer;

	if(!IsWithinMapCell(EnvNAVXYTHETAVSTEERCfg.StartX_c, EnvNAVXYTHETAVSTEERCfg.StartY_c)){
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.StartTheta < 0 || EnvNAVXYTHETAVSTEERCfg.StartTheta >= EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVSTEERCfg.StartV < 0 || EnvNAVXYTHETAVSTEERCfg.StartV >= EnvNAVXYTHETAVSTEERCfg.numV){
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVSTEERCfg.StartSteer < 0 || EnvNAVXYTHETAVSTEERCfg.StartSteer >= EnvNAVXYTHETAVSTEERCfg.numSteers){
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVSTEERCfg.EndX_c = goalx;
	EnvNAVXYTHETAVSTEERCfg.EndY_c = goaly;
	EnvNAVXYTHETAVSTEERCfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVSTEERCfg.EndV = goalv;
	EnvNAVXYTHETAVSTEERCfg.EndSteer = goalsteer;

	if(!IsWithinMapCell(EnvNAVXYTHETAVSTEERCfg.EndX_c, EnvNAVXYTHETAVSTEERCfg.EndY_c)){
		SBPL_ERROR("ERROR: illegal goal coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVSTEERCfg.EndTheta < 0 || EnvNAVXYTHETAVSTEERCfg.EndTheta >= EnvNAVXYTHETAVSTEERCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVSTEERCfg.EndV < 0 || EnvNAVXYTHETAVSTEERCfg.EndV >= EnvNAVXYTHETAVSTEERCfg.numV){
		SBPL_ERROR("ERROR: illegal goal coordinates for v\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVSTEERCfg.EndSteer < 0 || EnvNAVXYTHETAVSTEERCfg.EndSteer >= EnvNAVXYTHETAVSTEERCfg.numSteers){
		SBPL_ERROR("ERROR: illegal goal coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVSTEERCfg.FootprintPolygon = perimeterptsV;

	EnvNAVXYTHETAVSTEERCfg.cellsize_m = cellsize_m;
	
	//allocate the 2D environment
	EnvNAVXYTHETAVSTEERCfg.Grid2D = new unsigned char*[EnvNAVXYTHETAVSTEERCfg.EnvWidth_c];
	for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVSTEERCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETAVSTEERCfg.EnvHeight_c];
	}

	//environment:
	if (0 == mapdata) {
		for (int y = 0; y < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVSTEERCfg.Grid2D[x][y] = 0;
			}
		}
	}
	else {
		for (int y = 0; y < EnvNAVXYTHETAVSTEERCfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVSTEERCfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVSTEERCfg.Grid2D[x][EnvNAVXYTHETAVSTEERCfg.EnvHeight_c-1-y] = mapdata[x + y * width];
			}
		}
	}
}

//------------------------------------------------------------------------------