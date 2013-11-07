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

#include <sbpl/discrete_space_information/environment_navxythetav.h>
#include <sbpl/planners/planner.h>
#include <sbpl/utils/2Dgridsearch.h>
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
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

EnvironmentNAVXYTHETAV::EnvironmentNAVXYTHETAV(){
	EnvNAVXYTHETAVCfg.obsthresh = NAVXYTHETAV_DEFAULTOBSTHRESHOLD;
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVCfg.cost_inscribed_thresh = EnvNAVXYTHETAVCfg.obsthresh; 
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh = -1; 

	grid2Dsearchfromstart = NULL;
	grid2Dsearchfromgoal = NULL;
	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;
	
	iteration = 0;

	EnvNAVXYTHETAV.bInitialized = false;

	EnvNAVXYTHETAVCfg.NumThetaDirs = NAVXYTHETAV_DEFAULTTHETADIRS;

	//no memory allocated in cfg yet
	EnvNAVXYTHETAVCfg.Grid2D = NULL;
	EnvNAVXYTHETAVCfg.ActionsV = NULL;
	EnvNAVXYTHETAVCfg.PredActionsV = NULL;
}

EnvironmentNAVXYTHETAV::~EnvironmentNAVXYTHETAV(){
	SBPL_PRINTF("destroying XYTHETAV\n");
	if (grid2Dsearchfromstart != NULL) delete grid2Dsearchfromstart;
	grid2Dsearchfromstart = NULL;

	if (grid2Dsearchfromgoal != NULL) delete grid2Dsearchfromgoal;
	grid2Dsearchfromgoal = NULL;

	if (EnvNAVXYTHETAVCfg.Grid2D != NULL) {
		for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++)
			delete[] EnvNAVXYTHETAVCfg.Grid2D[x];
		delete[] EnvNAVXYTHETAVCfg.Grid2D;
		EnvNAVXYTHETAVCfg.Grid2D = NULL;
	}

	//delete actions
	if (EnvNAVXYTHETAVCfg.ActionsV != NULL) {
		for (int ind = 0; ind < (int)EnvNAVXYTHETAVCfg.ActionsV->size(); ind++)
			delete EnvNAVXYTHETAVCfg.ActionsV->at(ind);
		delete EnvNAVXYTHETAVCfg.ActionsV;
		EnvNAVXYTHETAVCfg.ActionsV = NULL;
	}
	if (EnvNAVXYTHETAVCfg.PredActionsV != NULL) {
		delete[] EnvNAVXYTHETAVCfg.PredActionsV;
		EnvNAVXYTHETAVCfg.PredActionsV = NULL;
	}
	
	//delete the states themselves first
	for (int i = 0; i < (int)EnvNAVXYTHETAV.StateID2DataTable.size(); i++) {
		delete EnvNAVXYTHETAV.StateID2DataTable.at(i);
		EnvNAVXYTHETAV.StateID2DataTable.at(i) = NULL;
	}
	EnvNAVXYTHETAV.StateID2DataTable.clear();

	//delete hashtable
	if (EnvNAVXYTHETAV.Data2StateIDHashTable != NULL) {
		delete[] EnvNAVXYTHETAV.Data2StateIDHashTable;
		EnvNAVXYTHETAV.Data2StateIDHashTable = NULL;
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
unsigned int EnvironmentNAVXYTHETAV::GETHASHBIN(unsigned int x, unsigned int y, unsigned int theta, unsigned int v)
{
	return inthash((inthash(x) + (inthash(y) << 1) + (inthash(theta) << 2) + (inthash(v) << 3))) &
		(EnvNAVXYTHETAV.HashTableSize - 1);
}

void EnvironmentNAVXYTHETAV::PrintHashTableHist()
{
	int s0 = 0, s1 = 0, s50 = 0, s100 = 0, s200 = 0, s300 = 0, slarge = 0;

	for (int j = 0; j < (int)EnvNAVXYTHETAV.HashTableSize; j++) {
		if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() == 0)
			s0++;
		else if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() < 50)
			s1++;
		else if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() < 100)
			s50++;
		else if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() < 200)
			s100++;
		else if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() < 300)
			s200++;
		else if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[j].size() < 400)
			s300++;
		else
			slarge++;
	}
	SBPL_PRINTF("hash table histogram: 0:%d, <50:%d, <100:%d, <200:%d, <300:%d, <400:%d >400:%d\n", s0, s1, s50, s100,
				s200, s300, slarge);
}

void EnvironmentNAVXYTHETAV::ReadConfiguration(FILE* fCfg)
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
	EnvNAVXYTHETAVCfg.EnvWidth_c = atoi(sTemp);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (discretization)\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.EnvHeight_c = atoi(sTemp);

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
	
	EnvNAVXYTHETAVCfg.NumThetaDirs = atoi(sTemp);

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
	EnvNAVXYTHETAVCfg.numV = atoi(sTemp);
	
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
	for(int i=0;i<EnvNAVXYTHETAVCfg.numV;i++){
		if (fscanf(fCfg, "%s", sTemp) != 1) {
			SBPL_ERROR("ERROR: ran out of env file early\n");
			throw new SBPL_Exception();
		}
		
		EnvNAVXYTHETAVCfg.velocities.push_back(atof(sTemp));
	}
	
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
	EnvNAVXYTHETAVCfg.obsthresh = atoi(sTemp);
	SBPL_PRINTF("obsthresh = %d\n", EnvNAVXYTHETAVCfg.obsthresh);

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
	EnvNAVXYTHETAVCfg.cost_inscribed_thresh = atoi(sTemp);
	SBPL_PRINTF("cost_inscribed_thresh = %d\n", EnvNAVXYTHETAVCfg.cost_inscribed_thresh);

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
	EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh = atoi(sTemp);
	SBPL_PRINTF("cost_possibly_circumscribed_thresh = %d\n", EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh);

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
	EnvNAVXYTHETAVCfg.cellsize_m = atof(sTemp);

	//start(meters,rads,m/s):
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.StartX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.StartY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.StartTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVCfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.StartV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVCfg.velocities);

	if (EnvNAVXYTHETAVCfg.StartX_c < 0 || EnvNAVXYTHETAVCfg.StartX_c >= EnvNAVXYTHETAVCfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.StartY_c < 0 || EnvNAVXYTHETAVCfg.StartY_c >= EnvNAVXYTHETAVCfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.StartTheta < 0 || EnvNAVXYTHETAVCfg.StartTheta >= EnvNAVXYTHETAVCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.StartV < 0 || EnvNAVXYTHETAVCfg.StartV >= EnvNAVXYTHETAVCfg.numV) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
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
	EnvNAVXYTHETAVCfg.EndX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.EndY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVCfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.EndTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVCfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.EndV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVCfg.velocities);

	if (EnvNAVXYTHETAVCfg.EndX_c < 0 || EnvNAVXYTHETAVCfg.EndX_c >= EnvNAVXYTHETAVCfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.EndY_c < 0 || EnvNAVXYTHETAVCfg.EndY_c >= EnvNAVXYTHETAVCfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.EndTheta < 0 || EnvNAVXYTHETAVCfg.EndTheta >= EnvNAVXYTHETAVCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.EndV < 0 || EnvNAVXYTHETAVCfg.EndV >= EnvNAVXYTHETAVCfg.numV) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}

	//allocate the 2D environment
	EnvNAVXYTHETAVCfg.Grid2D = new unsigned char*[EnvNAVXYTHETAVCfg.EnvWidth_c];
	for (x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETAVCfg.EnvHeight_c];
	}

	//environment:
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	for (y = EnvNAVXYTHETAVCfg.EnvHeight_c-1; y >= 0; y--)
		for (x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
			if (fscanf(fCfg, "%d", &dTemp) != 1) {
				SBPL_ERROR("ERROR: incorrect format of config file\n");
				throw new SBPL_Exception();
			}
			EnvNAVXYTHETAVCfg.Grid2D[x][y] = dTemp;
		}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c, EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.StartV)){
		SBPL_ERROR("ERROR: invalid start state\n");
		throw new SBPL_Exception();
	}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c, EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.EndV)){
		SBPL_ERROR("ERROR: invalid end state\n");
		throw new SBPL_Exception();
	}
}

EnvNAVXYTHETAVHashEntry_t* EnvironmentNAVXYTHETAV::GetHashEntry(unsigned int x, unsigned int y, unsigned int theta, unsigned int v)
{
	//clock_t currenttime = clock();

	int binid = GETHASHBIN(x, y, theta, v);

#if DEBUG
	if ((int)EnvNAVXYTHETAV.Data2StateIDHashTable[binid].size() > 500)
	{
		SBPL_PRINTF("WARNING: Hash table has a bin %d (X1=%d X2=%d X3=%d X4=%d) of size %d\n",
					binid, X1, X2, X3, X4, (int)EnvNAVXYTHETAV.Data2StateIDHashTable[binid].size());

		PrintHashTableHist();
	}
#endif

	//iterate over the states in the bin and select the perfect match
	for (int ind = 0; ind < (int)EnvNAVXYTHETAV.Data2StateIDHashTable[binid].size(); ind++) {
		if (EnvNAVXYTHETAV.Data2StateIDHashTable[binid][ind]->x == x &&
			EnvNAVXYTHETAV.Data2StateIDHashTable[binid][ind]->y == y &&
			EnvNAVXYTHETAV.Data2StateIDHashTable[binid][ind]->theta == theta &&
			EnvNAVXYTHETAV.Data2StateIDHashTable[binid][ind]->v == v)
		{
			//time_gethash += clock()-currenttime;
			return EnvNAVXYTHETAV.Data2StateIDHashTable[binid][ind];
		}
	}

	//time_gethash += clock()-currenttime;

	return NULL;
}

EnvNAVXYTHETAVHashEntry_t* EnvironmentNAVXYTHETAV::CreateNewHashEntry(unsigned int x, unsigned int y, unsigned int theta,
													unsigned int v)
{
	int i;

	//clock_t currenttime = clock();

	EnvNAVXYTHETAVHashEntry_t* HashEntry = new EnvNAVXYTHETAVHashEntry_t;

	HashEntry->x = x;
	HashEntry->y = y;
	HashEntry->theta = theta;
	HashEntry->v = v;
	HashEntry->iteration = 0;

	HashEntry->stateID = EnvNAVXYTHETAV.StateID2DataTable.size();

	//insert into the tables
	EnvNAVXYTHETAV.StateID2DataTable.push_back(HashEntry);

	//get the hash table bin
	i = GETHASHBIN(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v);

	//insert the entry into the bin
	EnvNAVXYTHETAV.Data2StateIDHashTable[i].push_back(HashEntry);

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

void EnvironmentNAVXYTHETAV::CreateStartandGoalStates()
{
	EnvNAVXYTHETAVHashEntry_t* HashEntry;

	//create start state
	unsigned int x = 0;
	unsigned int y = 0;
	unsigned int theta = 0;
	unsigned int v = 0;
	HashEntry = CreateNewHashEntry(x, y, theta, v);
	EnvNAVXYTHETAV.startstateid = HashEntry->stateID;

	//create goal state
	x = y = theta = v = 1;
	HashEntry = CreateNewHashEntry(x, y, theta, v);
	EnvNAVXYTHETAV.goalstateid = HashEntry->stateID;
}

void EnvironmentNAVXYTHETAV::ComputeReplanningDataforAction(EnvNAVXYTHETAVAction_t* action){
	int j;

	//iterate over all the cells involved in the action
	sbpl_xy_theta_v_cell_t startcell4d, endcell4d;
	for (int i = 0; i < (int)action->intersectingcellsV.size(); i++) {
		//compute the translated affected search Pose - what state has an
		//outgoing action whose intersecting cell is at 0,0
		startcell4d.v = action->startv;
		startcell4d.theta = action->starttheta;
		startcell4d.x = -action->intersectingcellsV.at(i).x;
		startcell4d.y = -action->intersectingcellsV.at(i).y;

		//compute the translated affected search Pose - what state has an
		//incoming action whose intersecting cell is at 0,0
		endcell4d.v = action->endv;
		endcell4d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
		endcell4d.x = startcell4d.x + action->dX;
		endcell4d.y = startcell4d.y + action->dY;

		//store the cells if not already there
		for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
			if (affectedsuccstatesV.at(j) == endcell4d) break;
		}
		if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell4d);

		for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
			if (affectedpredstatesV.at(j) == startcell4d) break;
		}
		if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell4d);
	}//over intersecting cells

	//---intersecting cell = origin
	//compute the translated affected search Pose - what state has an outgoing action whose intersecting cell is at 0,0
	startcell4d.v = action->startv;
	startcell4d.theta = action->starttheta;
	startcell4d.x = -0;
	startcell4d.y = -0;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell4d.v = action->endv;
	endcell4d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	endcell4d.x = startcell4d.x + action->dX;
	endcell4d.y = startcell4d.y + action->dY;

	//store the cells if not already there
	for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
		if (affectedsuccstatesV.at(j) == endcell4d) break;
	}
	if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell4d);

	for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
		if (affectedpredstatesV.at(j) == startcell4d) break;
	}
	if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell4d);

	//---intersecting cell = outcome state
	//compute the translated affected search Pose - what state has an outgoing action whose intersecting cell is at 0,0
	startcell4d.v = action->startv;
	startcell4d.theta = action->starttheta;
	startcell4d.x = -action->dX;
	startcell4d.y = -action->dY;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell4d.v = action->endv;
	endcell4d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	endcell4d.x = startcell4d.x + action->dX;
	endcell4d.y = startcell4d.y + action->dY;

	for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
		if (affectedsuccstatesV.at(j) == endcell4d) break;
	}
	if (j == (int)affectedsuccstatesV.size()) affectedsuccstatesV.push_back(endcell4d);

	for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
		if (affectedpredstatesV.at(j) == startcell4d) break;
	}
	if (j == (int)affectedpredstatesV.size()) affectedpredstatesV.push_back(startcell4d);
}

void EnvironmentNAVXYTHETAV::ComputeReplanningData()
{
	//iterate over all actions
	//velocities
	for(int vind = 0; vind < EnvNAVXYTHETAVCfg.numV; vind++){
		//orientations
		for (int tind = 0; tind < EnvNAVXYTHETAVCfg.NumThetaDirs; tind++) {
			//actions
			for (int aind = 0; aind < (int)(EnvNAVXYTHETAVCfg.ActionsV->at(vind*EnvNAVXYTHETAVCfg.NumThetaDirs+tind)->size()); aind++) {
				//compute replanning data for this action 
				ComputeReplanningDataforAction(&EnvNAVXYTHETAVCfg.ActionsV->at(vind*EnvNAVXYTHETAVCfg.NumThetaDirs+tind)->at(aind));
			}
		}
	}
}

void EnvironmentNAVXYTHETAV::PrecomputeActionswithCompleteMotionPrimitive(vector<SBPL_xythetav_mprimitive>* motionprimitiveV){
	SBPL_PRINTF("Pre-computing action data using motion primitives for every pair velocity/angle...\n");
	EnvNAVXYTHETAVCfg.ActionsV = new vector<vector<EnvNAVXYTHETAVAction_t> *>();
	EnvNAVXYTHETAVCfg.PredActionsV = new vector<EnvNAVXYTHETAVActionIndex_t> [EnvNAVXYTHETAVCfg.numV*EnvNAVXYTHETAVCfg.NumThetaDirs];
	
	vector<sbpl_2Dcell_t> footprint;

	/*
	* if (motionprimitiveV->size() % EnvNAVXYTHETALATCfg.NumThetaDirs != 0) {
	* SBPL_ERROR("ERROR: motionprimitives should be uniform across actions\n");
	* throw new SBPL_Exception();
	* }
	*/

	//EnvNAVXYTHETAVCfg.actionwidth = ((int)motionprimitiveV->size()) / EnvNAVXYTHETALATCfg.NumThetaDirs;

	//Allocate vectors
	for(int vind = 0; vind < EnvNAVXYTHETAVCfg.numV; vind++)
		for (int tind = 0; tind < EnvNAVXYTHETAVCfg.NumThetaDirs; tind++){
			EnvNAVXYTHETAVCfg.ActionsV->push_back(new vector<EnvNAVXYTHETAVAction_t>());
		}
	
	//iterate over source angles
	int maxnumofactions = 0;
	for(int vind = 0; vind < EnvNAVXYTHETAVCfg.numV; vind++){
		for (int tind = 0; tind < EnvNAVXYTHETAVCfg.NumThetaDirs; tind++) {
			SBPL_PRINTF("pre-computing for pair (speed, angle) (%d,%d) out of (%d,%d)\n", vind, tind, EnvNAVXYTHETAVCfg.numV, EnvNAVXYTHETAVCfg.NumThetaDirs);

			//compute current index
			int vector_index = vind*EnvNAVXYTHETAVCfg.NumThetaDirs+tind;
			
			//EnvNAVXYTHETAVCfg.ActionsV->at(vector_index) = new vector<EnvNAVXYTHETAVAction_t>();

			//compute sourcepose
			sbpl_xy_theta_v_pt_t sourcepose;
			sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVCfg.cellsize_m);
			sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVCfg.cellsize_m);
			sourcepose.theta = DiscTheta2ContNotUnif(tind, EnvNAVXYTHETAVCfg.NumThetaDirs);
			sourcepose.v = DiscV2Cont(vind, EnvNAVXYTHETAVCfg.velocities);

			//iterate over motion primitives
			int numofactions = 0;
			int aind = -1;
			for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
				//find a motion primitive for this angle
				if (motionprimitiveV->at(mind).start_theta_disc != tind || motionprimitiveV->at(mind).start_v_disc != vind) continue;

				aind++;
				numofactions++;
				EnvNAVXYTHETAVAction_t element_to_add;
				
				//action index
				element_to_add.aind = aind;

				//start angle
				element_to_add.starttheta = tind;
				
				//start velocity
				element_to_add.startv = vind;

				//compute dislocation
				element_to_add.endtheta = motionprimitiveV->at(mind).endcell.theta;
				element_to_add.endv = motionprimitiveV->at(mind).endcell.v;
				element_to_add.dX = motionprimitiveV->at(mind).endcell.x;
				element_to_add.dY = motionprimitiveV->at(mind).endcell.y;

				//compute and store interm points as well as intersecting cells
				element_to_add.intersectingcellsV.clear();
				element_to_add.intermptV.clear();
				element_to_add.interm3DcellsV.clear();

				sbpl_xy_theta_v_cell_t previnterm3Dcell;
				previnterm3Dcell.x = 0;
				previnterm3Dcell.y = 0;

				// Compute all the intersected cells for this action (intermptV and interm3DcellsV)
				for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
					sbpl_xy_theta_v_pt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
					element_to_add.intermptV.push_back(intermpt);

					// also compute the intermediate discrete cells if not there already
					sbpl_xy_theta_v_pt_t pose;
					pose.x = intermpt.x + sourcepose.x;
					pose.y = intermpt.y + sourcepose.y;
					pose.theta = intermpt.theta;
					pose.v = intermpt.v;

					sbpl_xy_theta_v_cell_t intermediate2dCell;
					intermediate2dCell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETAVCfg.cellsize_m);
					intermediate2dCell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETAVCfg.cellsize_m);

					// add unique cells to the list
					if (element_to_add.interm3DcellsV.size() == 0 || intermediate2dCell.x
						!= previnterm3Dcell.x || intermediate2dCell.y != previnterm3Dcell.y) {
						element_to_add.interm3DcellsV.push_back(intermediate2dCell);
					}

					previnterm3Dcell = intermediate2dCell;
				}

				//compute linear and angular time
				double linear_distance = 0;
				double medium_velocity = element_to_add.intermptV[0].v;
				//double xmax=element_to_add.intermptV[0].x, xmin=element_to_add.intermptV[0].x, ymax=element_to_add.intermptV[0].y, ymin=element_to_add.intermptV[0].y;
				for (unsigned int i = 1; i<element_to_add.intermptV.size(); i++) {
					double x0 = element_to_add.intermptV[i - 1].x;
					double y0 = element_to_add.intermptV[i - 1].y;
					double x1 = element_to_add.intermptV[i].x;
					double y1 = element_to_add.intermptV[i].y;
					double dx = x1 - x0;
					double dy = y1 - y0;
					linear_distance += sqrt(dx * dx + dy * dy);
					medium_velocity += fabs(element_to_add.intermptV[i].v);
					
					/*if(element_to_add.intermptV[i].x > xmax)
						xmax = element_to_add.intermptV[i].x;
					
					if(element_to_add.intermptV[i].x < xmin)
						xmin = element_to_add.intermptV[i].x;
					
					if(element_to_add.intermptV[i].y > ymax)
						ymax = element_to_add.intermptV[i].y;
					
					if(element_to_add.intermptV[i].y < ymin)
						ymin = element_to_add.intermptV[i].y;*/
				}
				
				medium_velocity /= (int)element_to_add.intermptV.size();
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
				element_to_add.cost = NAVXYTHETAV_COSTMULT_MTOMM * motionprimitiveV->at(mind).additionalactioncostmult;
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
				get_2d_motion_cells(EnvNAVXYTHETAVCfg.FootprintPolygon, motionprimitiveV->at(mind).intermptV, &element_to_add.intersectingcellsV, EnvNAVXYTHETAVCfg.cellsize_m);
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
				int postheta = (element_to_add.endtheta >= 0)?element_to_add.endtheta:element_to_add.endtheta+EnvNAVXYTHETAVCfg.NumThetaDirs;
				int targetindex = element_to_add.endv*EnvNAVXYTHETAVCfg.NumThetaDirs+postheta;
				EnvNAVXYTHETAVActionIndex_t idx;
				idx.vector_index = vector_index;
				idx.aind = aind;
				EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->push_back(element_to_add);
				EnvNAVXYTHETAVCfg.PredActionsV[targetindex].push_back(idx);
			}

			if (maxnumofactions < numofactions) maxnumofactions = numofactions;
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

void EnvironmentNAVXYTHETAV::InitializeEnvConfig(vector<SBPL_xythetav_mprimitive>* motionprimitiveV)
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
	get_2d_footprint_cells(EnvNAVXYTHETAVCfg.FootprintPolygon, &footprint, temppose, EnvNAVXYTHETAVCfg.cellsize_m);
	SBPL_PRINTF("number of cells in footprint of the robot = %d\n", (unsigned int)footprint.size());

	for (vector<sbpl_2Dcell_t>::iterator it = footprint.begin(); it != footprint.end(); ++it) {
		SBPL_PRINTF("Footprint cell at (%d, %d)\n", it->x, it->y);
	}

#if DEBUG
	SBPL_FPRINTF(fDeb, "footprint cells (size=%d):\n", (int)footprint.size());
	for(int i = 0; i < (int) footprint.size(); i++)
	{
		SBPL_FPRINTF(fDeb, "%d %d (cont: %.3f %.3f)\n", footprint.at(i).x, footprint.at(i).y,
					DISCXY2CONT(footprint.at(i).x, EnvNAVXYTHETAVCfg.cellsize_m),
					DISCXY2CONT(footprint.at(i).y, EnvNAVXYTHETAVCfg.cellsize_m));
	}
#endif
	
	PrecomputeActionswithCompleteMotionPrimitive(motionprimitiveV);
}

bool EnvironmentNAVXYTHETAV::InitGeneral(vector<SBPL_xythetav_mprimitive>* motionprimitiveV){
	//Initialize other parameters of the environment
	InitializeEnvConfig(motionprimitiveV);

	//initialize Environment
	InitializeEnvironment();

	//pre-compute heuristics
	ComputeHeuristicValues();

	return true;
}

void EnvironmentNAVXYTHETAV::InitializeEnvironment()
{
	EnvNAVXYTHETAVHashEntry_t* HashEntry;

	int maxsize = EnvNAVXYTHETAVCfg.EnvWidth_c * EnvNAVXYTHETAVCfg.EnvHeight_c * EnvNAVXYTHETAVCfg.NumThetaDirs * EnvNAVXYTHETAVCfg.numV;

	SBPL_PRINTF("environment stores states in hashtable\n");

	//initialize the map from Data to StateID
	//Maximum hash table size
	EnvNAVXYTHETAV.HashTableSize = EnvNAVXYTHETAVCfg.NumThetaDirs * EnvNAVXYTHETAVCfg.EnvWidth_c * EnvNAVXYTHETAVCfg.EnvHeight_c * 8;/*EnvNAVXYTHETAVCfg.numV;*/ //should be power of two - REVIEW -> USE DYNAMIC PARAMETERS
	EnvNAVXYTHETAV.Data2StateIDHashTable = new vector<EnvNAVXYTHETAVHashEntry_t*> [EnvNAVXYTHETAV.HashTableSize];

	//initialize the map from StateID to Data
	EnvNAVXYTHETAV.StateID2DataTable.clear();

	//create start state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c,
										EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.StartV)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c,
												EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.StartV);
	}
	EnvNAVXYTHETAV.startstateid = HashEntry->stateID;

	//create goal state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c,
										EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.EndV)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c,
												EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.EndV);
	}
	EnvNAVXYTHETAV.goalstateid = HashEntry->stateID;

	//initialized
	EnvNAVXYTHETAV.bInitialized = true;
}

bool EnvironmentNAVXYTHETAV::ReadMotionPrimitives(FILE* fMotPrims){
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
	if (fabs(fTemp - EnvNAVXYTHETAVCfg.cellsize_m) > ERR_EPS) {
		SBPL_ERROR("ERROR: invalid resolution %f (instead of %f) in the dynamics file\n", fTemp,
				EnvNAVXYTHETAVCfg.cellsize_m);
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
	if (dTemp != EnvNAVXYTHETAVCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVCfg.NumThetaDirs);
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
	if (dTemp != EnvNAVXYTHETAVCfg.numV) {
		SBPL_ERROR("ERROR: invalid velocity resolution %d velocities (instead of %d velocities) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVCfg.numV);
		return false;
	}
	
	//read in the velocity values
	strcpy(sExpected, "velocities:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	for(int i=0;i<EnvNAVXYTHETAVCfg.numV;i++){
		if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
		if (fTemp != EnvNAVXYTHETAVCfg.velocities[i]) {
			SBPL_ERROR("ERROR: invalid velocity value %f velocity (instead of %f velocity) in the motion primitives file\n",
					fTemp, EnvNAVXYTHETAVCfg.velocities[i]);
			return false;
		}
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
		SBPL_xythetav_mprimitive motprim;

		if (ReadSingleMotionPrimitive(&motprim, fMotPrims) == false) return false;

		EnvNAVXYTHETAVCfg.mprimV.push_back(motprim);
	}
	SBPL_PRINTF("done reading of motion primitives");

	return true;
}

bool EnvironmentNAVXYTHETAV::ReadSingleMotionPrimitive(SBPL_xythetav_mprimitive* pMotPrim, FILE* fIn){
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
		sbpl_xy_theta_v_pt_t intermpose;
		if (ReadSinglePose(&intermpose, fIn) == false) {
			SBPL_ERROR("ERROR: failed to read in intermediate poses\n");
			return false;
		}
		pMotPrim->intermptV.push_back(intermpose);
	}

	//check that the last pose corresponds correctly to the last pose
	sbpl_xy_theta_v_pt_t sourcepose;
	sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVCfg.cellsize_m);
	sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVCfg.cellsize_m);
	sourcepose.theta = DiscTheta2ContNotUnif(pMotPrim->start_theta_disc, EnvNAVXYTHETAVCfg.NumThetaDirs);
	sourcepose.v = DiscV2Cont(pMotPrim->start_v_disc, EnvNAVXYTHETAVCfg.velocities);
	double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
	double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
	double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;
	double mp_endv_ms = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v;
	int endx_disc = CONTXY2DISC(mp_endx_m, EnvNAVXYTHETAVCfg.cellsize_m);
	int endy_disc = CONTXY2DISC(mp_endy_m, EnvNAVXYTHETAVCfg.cellsize_m);
	int endtheta_disc = ContTheta2DiscNotUnif(mp_endtheta_rad, EnvNAVXYTHETAVCfg.NumThetaDirs);
	int endv_disc = ContV2Disc(mp_endv_ms, EnvNAVXYTHETAVCfg.velocities);
	if (endx_disc != pMotPrim->endcell.x || endy_disc != pMotPrim->endcell.y || endtheta_disc != pMotPrim->endcell.theta || endv_disc != pMotPrim->endcell.v) {
		SBPL_ERROR( "ERROR: incorrect primitive %d with startangle=%d and startv=%d "
				"last interm point %f %f %f %f does not match end pose %d %d %d %d\n",
				pMotPrim->motprimID, pMotPrim->start_theta_disc, pMotPrim->start_v_disc,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v,
				pMotPrim->endcell.x, pMotPrim->endcell.y,
				pMotPrim->endcell.theta, pMotPrim->endcell.v);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAV::ReadSingleCell(sbpl_xy_theta_v_cell_t* cell, FILE* fIn){
	char sTemp[60];

	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->x = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->y = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->theta = atoi(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	cell->v = atoi(sTemp);
	
	//normalize the angle
	cell->theta = NORMALIZEDISCTHETA(cell->theta, EnvNAVXYTHETAVCfg.NumThetaDirs);

	return true;
}

bool EnvironmentNAVXYTHETAV::ReadSinglePose(sbpl_xy_theta_v_pt_t* pose, FILE* fIn){
	char sTemp[60];

	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->x = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->y = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->theta = atof(sTemp);
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	pose->v = atof(sTemp);

	pose->theta = normalizeAngle(pose->theta);

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

void EnvironmentNAVXYTHETAV::EnsureHeuristicsUpdated(bool bGoalHeuristics){
	if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
		grid2Dsearchfromstart->search(EnvNAVXYTHETAVCfg.Grid2D, EnvNAVXYTHETAVCfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c,
									EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeStartHeuristics = false;
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.EndX_c,
																				EnvNAVXYTHETAVCfg.EndY_c)/NAVXYTHETAV_DEFAULTMEDIUMVELOCITY));
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.EndX_c,
																				EnvNAVXYTHETAVCfg.EndY_c)));*/
	}

	if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
		grid2Dsearchfromgoal->search(EnvNAVXYTHETAVCfg.Grid2D, EnvNAVXYTHETAVCfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c,
									EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeGoalHeuristics = false;
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.StartX_c,
																				EnvNAVXYTHETAVCfg.StartY_c)/NAVXYTHETAV_DEFAULTMEDIUMVELOCITY));
		
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.StartX_c,
																				EnvNAVXYTHETAVCfg.StartY_c)));*/
	}
}

void EnvironmentNAVXYTHETAV::ComputeHeuristicValues()
{
	//whatever necessary pre-computation of heuristic values is done here
	SBPL_PRINTF("Precomputing heuristics\n");
	
	//allocated 2D grid searches
	grid2Dsearchfromstart = new SBPL2DGridSearch(EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVCfg.cellsize_m);
	grid2Dsearchfromgoal = new SBPL2DGridSearch(EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVCfg.cellsize_m);

	//set OPEN type to sliding buckets
	grid2Dsearchfromstart->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	grid2Dsearchfromgoal->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	
	SBPL_PRINTF("done\n");
}

void EnvironmentNAVXYTHETAV::PrintHeuristicValues(){
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

	for (int y = 0; y < EnvNAVXYTHETAVCfg.EnvHeight_c; y++) {
		for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
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

bool EnvironmentNAVXYTHETAV::InitializeEnv(const char* sEnvFile, const vector<sbpl_2Dpt_t>& perimeterptsV, const char* sMotPrimFile)
{
	EnvNAVXYTHETAVCfg.FootprintPolygon = perimeterptsV;

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
		InitGeneral(&EnvNAVXYTHETAVCfg.mprimV);
		fclose(fMotPrim);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

	SBPL_PRINTF("size of env: %d by %d\n", EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c);

	return true;
}

bool EnvironmentNAVXYTHETAV::InitializeEnv(const char* sEnvFile){
	return false;
}

bool EnvironmentNAVXYTHETAV::SetEnvParameter(const char* parameter, int value){
	if (EnvNAVXYTHETAV.bInitialized == true) {
		SBPL_ERROR("ERROR: all parameters must be set before initialization of the environment\n");
		return false;
	}

	SBPL_PRINTF("setting parameter %s to %d\n", parameter, value);

	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVCfg.cost_inscribed_thresh = (unsigned char)value;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh = value;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVCfg.obsthresh = (unsigned char)value;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		return false;
	}

	return true;
}

int EnvironmentNAVXYTHETAV::GetEnvParameter(const char* parameter)
{
	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVCfg.cost_inscribed_thresh;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		return (int)EnvNAVXYTHETAVCfg.obsthresh;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		throw new SBPL_Exception();
	}
}

bool EnvironmentNAVXYTHETAV::InitializeMDPCfg(MDPConfig *MDPCfg)
{
	//initialize MDPCfg with the start and goal ids
	MDPCfg->goalstateid = EnvNAVXYTHETAV.goalstateid;
	MDPCfg->startstateid = EnvNAVXYTHETAV.startstateid;

	return true;
}

int EnvironmentNAVXYTHETAV::GetFromToHeuristic(int FromStateID, int ToStateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if(FromStateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size() || ToStateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size()){
		SBPL_ERROR("ERROR in EnvNAVXYTHETALAT... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	//get X, Y for the state
	EnvNAVXYTHETAVHashEntry_t* FromHashEntry = EnvNAVXYTHETAV.StateID2DataTable[FromStateID];
	EnvNAVXYTHETAVHashEntry_t* ToHashEntry = EnvNAVXYTHETAV.StateID2DataTable[ToStateID];

	//TODO - check if one of the gridsearches already computed and then use it.

	return (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y)/NAVXYTHETAV_DEFAULTMEDIUMVELOCITY);
	//return (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y));
}

int EnvironmentNAVXYTHETAV::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[stateID];
	//computes distances from start state that is grid2D, so it is EndX_c EndY_c
	int h2D = grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y); 
	int hEuclid = (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(HashEntry->x, HashEntry->y,
																		EnvNAVXYTHETAVCfg.EndX_c,
																		EnvNAVXYTHETAVCfg.EndY_c));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / NAVXYTHETAV_DEFAULTMEDIUMVELOCITY);
	//return (int)(((double)__max(h2D, hEuclid)));
}

int EnvironmentNAVXYTHETAV::GetStartHeuristic(int stateID)
{
	
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[stateID];
	int h2D = grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y);
	int hEuclid = (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(EnvNAVXYTHETAVCfg.StartX_c,
																		EnvNAVXYTHETAVCfg.StartY_c, HashEntry->x,
																		HashEntry->y));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / NAVXYTHETAV_DEFAULTMEDIUMVELOCITY);
	//return (int)(((double)__max(h2D, hEuclid)));
}

void EnvironmentNAVXYTHETAV::SetAllActionsandAllOutcomes(CMDPSTATE* state){
	int cost;
	
#if DEBUG
	if(state->StateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size())
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
	if (state->StateID == EnvNAVXYTHETAV.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[state->StateID];
	
	int vector_index = HashEntry->v*EnvNAVXYTHETAVCfg.NumThetaDirs+HashEntry->theta;

	//iterate through actions
	for (int aind = 0; aind < EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVAction_t* nav4daction = &EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav4daction->dX;
		int newY = HashEntry->y + nav4daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav4daction->endtheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
		int newV = nav4daction->endv;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, nav4daction);
		if (cost >= INFINITECOST) continue;

		//add the action
		CMDPACTION* action = state->AddAction(aind);

#if TIME_DEBUG
		clock_t currenttime = clock();
#endif

		EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV);
		}
		action->AddOutcome(OutHashEntry->stateID, cost, 1.0);

#if TIME_DEBUG
		time3_addallout += clock()-currenttime;
#endif
	}
}

void EnvironmentNAVXYTHETAV::SetAllPreds(CMDPSTATE* state)
{
	//implement this if the planner needs access to predecessors

	SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: SetAllPreds is undefined\n");
	throw new SBPL_Exception();
}

void EnvironmentNAVXYTHETAV::GetSuccs(int SourceStateID, vector<int>* SuccIDV, vector<int>* CostV)
{
	GetSuccs(SourceStateID, SuccIDV, CostV, NULL);
}

void EnvironmentNAVXYTHETAV::GetSuccs(int sourceStateID, vector<int>* succIDV, vector<int>* costV, vector<EnvNAVXYTHETAVAction_t*>* actionindV){
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
	if (sourceStateID == EnvNAVXYTHETAV.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[sourceStateID];
	
	int vector_index = (unsigned int)HashEntry->v*EnvNAVXYTHETAVCfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//iterate through actions
	for (aind = 0; aind < EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVAction_t* nav4daction = &EnvNAVXYTHETAVCfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav4daction->dX;
		int newY = HashEntry->y + nav4daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav4daction->endtheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
		int newV = nav4daction->endv;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		int cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, nav4daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV);
		}

		succIDV->push_back(OutHashEntry->stateID);
		costV->push_back(cost);
		if (actionindV != NULL) actionindV->push_back(nav4daction);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentNAVXYTHETAV::GetPreds(int TargetStateID, vector<int>* PredIDV, vector<int>* CostV)
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
	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[TargetStateID];
	int targetindex = (unsigned int)HashEntry->v*EnvNAVXYTHETAVCfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//clear the preds array
	PredIDV->clear();
	CostV->clear();
	PredIDV->reserve(EnvNAVXYTHETAVCfg.PredActionsV[targetindex].size());
	CostV->reserve(EnvNAVXYTHETAVCfg.PredActionsV[targetindex].size());

	//iterate through actions
	vector<EnvNAVXYTHETAVActionIndex_t>* actionsV = &EnvNAVXYTHETAVCfg.PredActionsV[targetindex];
	for (aind = 0; aind < (int)actionsV->size(); aind++) {

		EnvNAVXYTHETAVAction_t* nav4daction = &EnvNAVXYTHETAVCfg.ActionsV->at(actionsV->at(aind).vector_index)->at(actionsV->at(aind).aind);

		int predX = HashEntry->x - nav4daction->dX;
		int predY = HashEntry->y - nav4daction->dY;
		int predTheta = nav4daction->starttheta;
		int predV = nav4daction->startv;

		//skip the invalid cells
		if (!IsValidCell(predX, predY)) continue;

		//get cost
		int cost = GetActionCost(predX, predY, predTheta, predV, nav4daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(predX, predY, predTheta, predV)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(predX, predY, predTheta, predV);
		}

		PredIDV->push_back(OutHashEntry->stateID);
		CostV->push_back(cost);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

int EnvironmentNAVXYTHETAV::SizeofCreatedEnv()
{
	return (int)EnvNAVXYTHETAV.StateID2DataTable.size();
}

void EnvironmentNAVXYTHETAV::PrintState(int stateID, bool bVerbose, FILE* fOut /*=NULL*/)
{
#if DEBUG
	if(stateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size())
	{
		SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal (2)\n");
		throw new SBPL_Exception();
	}
#endif

	if (fOut == NULL) fOut = stdout;

	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[stateID];

	if (stateID == EnvNAVXYTHETAV.goalstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a goal state\n");
	}
	
	if (stateID == EnvNAVXYTHETAV.startstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a start state\n");
	}

	SBPL_FPRINTF(fOut, "x=%d y=%d theta=%d v=%d\n", HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v);
}

void EnvironmentNAVXYTHETAV::PrintEnv_Config(FILE* fOut)
{
	if(fOut != NULL){
		fprintf(fOut, "discretization(cells): %d %d\n", EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c);
		fprintf(fOut, "NumThetaDirs: %d\n", EnvNAVXYTHETAVCfg.NumThetaDirs);
		fprintf(fOut, "NumV: %d\n", EnvNAVXYTHETAVCfg.numV);
		fprintf(fOut, "Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVCfg.numV; i++)
			fprintf(fOut, " %f", EnvNAVXYTHETAVCfg.velocities.at(i));
		fprintf(fOut, "\n");
		fprintf(fOut, "obsthresh: %d\n", EnvNAVXYTHETAVCfg.obsthresh);
		fprintf(fOut, "cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_inscribed_thresh);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "cellsize(meters): %f\n", EnvNAVXYTHETAVCfg.cellsize_m);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "start(meters,rads,m/s): %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVCfg.StartY_c, EnvNAVXYTHETAVCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVCfg.StartV, EnvNAVXYTHETAVCfg.velocities));
		fprintf(fOut, "end(meters,rads,m/s): %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVCfg.EndY_c, EnvNAVXYTHETAVCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVCfg.EndV, EnvNAVXYTHETAVCfg.velocities));
		fprintf(fOut, "environment:\n");
		
		for (int y = EnvNAVXYTHETAVCfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
				fprintf(fOut, "%d ", EnvNAVXYTHETAVCfg.Grid2D[x][y]);
			}
			
			fprintf(fOut, "\n");
		}
	}
	else{
		printf("discretization(cells): %d %d\n", EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c);
		printf("NumThetaDirs: %d\n", EnvNAVXYTHETAVCfg.NumThetaDirs);
		printf("NumV: %d\n", EnvNAVXYTHETAVCfg.numV);
		printf("Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVCfg.numV; i++)
			printf(" %f", EnvNAVXYTHETAVCfg.velocities.at(i));
		printf("\n");
		printf("obsthresh: %d\n", EnvNAVXYTHETAVCfg.obsthresh);
		printf("cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_inscribed_thresh);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh);
		printf("cellsize(meters): %f\n", EnvNAVXYTHETAVCfg.cellsize_m);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh);
		printf("start(meters,rads,m/s): %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVCfg.StartY_c, EnvNAVXYTHETAVCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVCfg.StartV, EnvNAVXYTHETAVCfg.velocities));
		printf("end(meters,rads,m/s): %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVCfg.EndY_c, EnvNAVXYTHETAVCfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVCfg.EndV, EnvNAVXYTHETAVCfg.velocities));
		printf("environment:\n");
		
		for (int y = EnvNAVXYTHETAVCfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
				printf("%d ", EnvNAVXYTHETAVCfg.Grid2D[x][y]);
			}
			
			printf("\n");
		}
	}
}

bool EnvironmentNAVXYTHETAV::InitializeEnv(int width, int height, int numthetadirs, int numv, vector<double> velocities,
							const unsigned char* mapdata,
							double startx, double starty, double starttheta, double startv,
							double goalx, double goaly, double goaltheta, double goalv,
							const vector<sbpl_2Dpt_t>& perimeterptsV, double cellsize_m,
							unsigned char obsthresh, unsigned char cost_inscribed_thresh,
							unsigned char cost_possibly_circumscribed_thresh, const char* sMotPrimFile){
	SBPL_PRINTF("env: initialize with width=%d height=%d start=%.3f %.3f %.3f %.3f "
				"goalx=%.3f %.3f %.3f %.3f cellsize=%.3f numthetadirs=%d numv=%d, obsthresh=%d cost_inscribed_thresh=%d cost_possibly_circumscribed_thresh=%d\n",
				width, height, startx, starty, starttheta, startv, goalx, goaly, goaltheta, goalv, cellsize_m, numthetadirs, numv,
				obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh);

	SBPL_PRINTF("NOTE: goaltol parameters currently unused\n");

	SBPL_PRINTF("perimeter has size=%d\n", (unsigned int)perimeterptsV.size());

	for (int i = 0; i < (int)perimeterptsV.size(); i++) {
		SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, perimeterptsV.at(i).x, perimeterptsV.at(i).y);
	}

	EnvNAVXYTHETAVCfg.obsthresh = obsthresh;
	EnvNAVXYTHETAVCfg.cost_inscribed_thresh = cost_inscribed_thresh;
	EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh = cost_possibly_circumscribed_thresh;
	EnvNAVXYTHETAVCfg.NumThetaDirs = numthetadirs;
	EnvNAVXYTHETAVCfg.numV = numv;
	
	for(int i=0;i<numv;i++){
		EnvNAVXYTHETAVCfg.velocities.at(i) = velocities.at(i);
	}

	//TODO - need to set the tolerance as well

	SetConfiguration(width, height, mapdata, CONTXY2DISC(startx, cellsize_m), CONTXY2DISC(starty, cellsize_m),
					ContTheta2DiscNotUnif(starttheta, EnvNAVXYTHETAVCfg.NumThetaDirs), ContV2Disc(startv, EnvNAVXYTHETAVCfg.velocities),
					CONTXY2DISC(goalx, cellsize_m), CONTXY2DISC(goaly, cellsize_m), ContTheta2DiscNotUnif(goaltheta, EnvNAVXYTHETAVCfg.NumThetaDirs),
					ContV2Disc(goalv, EnvNAVXYTHETAVCfg.velocities), cellsize_m, perimeterptsV);

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

	if (EnvNAVXYTHETAVCfg.mprimV.size() != 0) {
		InitGeneral(&EnvNAVXYTHETAVCfg.mprimV);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAV::UpdateCost(int x, int y, unsigned char newcost){
	EnvNAVXYTHETAVCfg.Grid2D[x][y] = newcost;

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

bool EnvironmentNAVXYTHETAV::SetMap(const unsigned char* mapdata){
	for (int xind = 0; xind < EnvNAVXYTHETAVCfg.EnvWidth_c; xind++) {
		for (int yind = 0; yind < EnvNAVXYTHETAVCfg.EnvHeight_c; yind++) {
			EnvNAVXYTHETAVCfg.Grid2D[xind][EnvNAVXYTHETAVCfg.EnvHeight_c-1-yind] = mapdata[xind + yind * EnvNAVXYTHETAVCfg.EnvWidth_c];
		}
	}

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

void EnvironmentNAVXYTHETAV::GetPredsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV, vector<int> *preds_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_cell_t affectedcell;
	EnvNAVXYTHETAVHashEntry_t* affectedHashEntry;

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
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

void EnvironmentNAVXYTHETAV::GetSuccsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV,
										vector<int> *succs_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_cell_t affectedcell;
	EnvNAVXYTHETAVHashEntry_t* affectedHashEntry;

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
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

bool EnvironmentNAVXYTHETAV::IsObstacle(int x, int y){
	return (EnvNAVXYTHETAVCfg.Grid2D[x][y] >= EnvNAVXYTHETAVCfg.obsthresh);
}

bool EnvironmentNAVXYTHETAV::IsValidConfiguration(int x, int y, int theta, int v){
	vector<sbpl_2Dcell_t> footprint;
	sbpl_xy_theta_v_pt_t pose;

	//compute continuous pose
	pose.x = DISCXY2CONT(x, EnvNAVXYTHETAVCfg.cellsize_m);
	pose.y = DISCXY2CONT(y, EnvNAVXYTHETAVCfg.cellsize_m);
	pose.theta = DiscTheta2ContNotUnif(theta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	pose.v = DiscV2Cont(v, EnvNAVXYTHETAVCfg.velocities);

	//compute footprint cells
	get_2d_footprint_cells(EnvNAVXYTHETAVCfg.FootprintPolygon, &footprint, pose, EnvNAVXYTHETAVCfg.cellsize_m);

	//iterate over all footprint cells
	for (int find = 0; find < (int)footprint.size(); find++) {
		int x = footprint.at(find).x;
		int y = footprint.at(find).y;

		if (x < 0 || x >= EnvNAVXYTHETAVCfg.EnvWidth_c || y < 0 || y >= EnvNAVXYTHETAVCfg.EnvHeight_c ||
			EnvNAVXYTHETAVCfg.Grid2D[x][y] >= EnvNAVXYTHETAVCfg.obsthresh)
		{
			return false;
		}
	}
	
	if(v < 0 || v >= EnvNAVXYTHETAVCfg.numV){
		return false;
	}

	return true;
}

void EnvironmentNAVXYTHETAV::GetEnvParams(int *size_x, int *size_y, int * numthetadirs, int * numv, vector<double> * velocities,
							double* startx, double* starty, double* starttheta, double * startv, 
							double* goalx, double* goaly, double* goaltheta, double * goalv, double* cellsize_m,
							unsigned char* cost_inscribed_thresh, unsigned char* cost_possibly_circumscribed_thresh,
							unsigned char* obsthresh, vector<SBPL_xythetav_mprimitive>* motionprimitiveV){
	
	*size_x = EnvNAVXYTHETAVCfg.EnvWidth_c;
	*size_y = EnvNAVXYTHETAVCfg.EnvHeight_c;

	*numthetadirs = EnvNAVXYTHETAVCfg.NumThetaDirs;
	*numv = EnvNAVXYTHETAVCfg.numV;
	for(int i;i<EnvNAVXYTHETAVCfg.velocities.size();i++){
		velocities->push_back(EnvNAVXYTHETAVCfg.velocities.at(i));
	}
	
	*startx = DISCXY2CONT(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.cellsize_m);
	*starty = DISCXY2CONT(EnvNAVXYTHETAVCfg.StartY_c, EnvNAVXYTHETAVCfg.cellsize_m);
	*starttheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.StartTheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	*startv = DiscV2Cont(EnvNAVXYTHETAVCfg.StartV, EnvNAVXYTHETAVCfg.velocities);
	*goalx = DISCXY2CONT(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.cellsize_m);
	*goaly = DISCXY2CONT(EnvNAVXYTHETAVCfg.EndY_c, EnvNAVXYTHETAVCfg.cellsize_m);
	*goaltheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVCfg.EndTheta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	*goalv = DiscV2Cont(EnvNAVXYTHETAVCfg.EndV, EnvNAVXYTHETAVCfg.velocities);

	*cellsize_m = EnvNAVXYTHETAVCfg.cellsize_m;
	
	*cost_inscribed_thresh = EnvNAVXYTHETAVCfg.cost_inscribed_thresh;
	*cost_possibly_circumscribed_thresh = EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh;
	*obsthresh = EnvNAVXYTHETAVCfg.obsthresh;

	*motionprimitiveV = EnvNAVXYTHETAVCfg.mprimV;
}

const EnvNAVXYTHETAVConfig_t* EnvironmentNAVXYTHETAV::GetEnvNavConfig(){
	return &EnvNAVXYTHETAVCfg;
}

void EnvironmentNAVXYTHETAV::PrintTimeStat(FILE* fOut){
	/*TO BE IMPLEMENTED*/
}

unsigned char EnvironmentNAVXYTHETAV::GetMapCost(int x, int y){
	return EnvNAVXYTHETAVCfg.Grid2D[x][y];
}

bool EnvironmentNAVXYTHETAV::IsWithinMapCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVCfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVCfg.EnvHeight_c);
}

bool EnvironmentNAVXYTHETAV::PoseContToDisc(double px, double py, double pth, double pv, int &ix, int &iy, int &ith, int &iv) const{
	ix = CONTXY2DISC(px, EnvNAVXYTHETAVCfg.cellsize_m);
	iy = CONTXY2DISC(py, EnvNAVXYTHETAVCfg.cellsize_m);
	ith = ContTheta2DiscNotUnif(pth, EnvNAVXYTHETAVCfg.NumThetaDirs);
	iv = ContV2Disc(pv, EnvNAVXYTHETAVCfg.velocities);
	
	return (pth >= -2 * PI_CONST) && (pth <= 2 * PI_CONST) && (ix >= 0 && ix < EnvNAVXYTHETAVCfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVCfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVCfg.numV);
}

bool EnvironmentNAVXYTHETAV::PoseDiscToCont(int ix, int iy, int ith, int iv, double &px, double &py, double &pth, double &pv) const{
	px = DISCXY2CONT(ix, EnvNAVXYTHETAVCfg.cellsize_m);
	py = DISCXY2CONT(iy, EnvNAVXYTHETAVCfg.cellsize_m);
	pth = normalizeAngle(DiscTheta2ContNotUnif(ith, EnvNAVXYTHETAVCfg.NumThetaDirs));
	pv = DiscV2Cont(iv, EnvNAVXYTHETAVCfg.velocities);
	
	return (ith >= 0) && (ith < EnvNAVXYTHETAVCfg.NumThetaDirs) && (ix >= 0 && ix < EnvNAVXYTHETAVCfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVCfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVCfg.numV);
}

int EnvironmentNAVXYTHETAV::SetStart(double x, double y, double theta, double v){
	int startx = CONTXY2DISC(x, EnvNAVXYTHETAVCfg.cellsize_m);
	int starty = CONTXY2DISC(y, EnvNAVXYTHETAVCfg.cellsize_m);
	int starttheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	int startv = ContV2Disc(v, EnvNAVXYTHETAVCfg.velocities);

	if (!IsWithinMapCell(startx, starty)) {
		SBPL_ERROR("ERROR: trying to set a start cell %d %d that is outside of map\n", startx, starty);
		return -1;
	}

	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f (%d %d %d %d)\n", x, y, theta, v, startx, starty, starttheta, startv);

	if (!IsValidConfiguration(startx, starty, starttheta, startv)) {
		SBPL_PRINTF("WARNING: start configuration %d %d %d %d is invalid\n", startx, starty, starttheta, startv);
	}

	EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(startx, starty, starttheta, startv)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(startx, starty, starttheta, startv);
	}

	if (EnvNAVXYTHETAV.startstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true;
		//because termination condition can be not all states TODO - make it dependent on term. condition
		bNeedtoRecomputeGoalHeuristics = true; 
	}

	//set start
	EnvNAVXYTHETAV.startstateid = OutHashEntry->stateID;
	EnvNAVXYTHETAVCfg.StartX_c = startx;
	EnvNAVXYTHETAVCfg.StartY_c = starty;
	EnvNAVXYTHETAVCfg.StartTheta = starttheta;
	EnvNAVXYTHETAVCfg.StartV = startv;

	return EnvNAVXYTHETAV.startstateid;
}

int EnvironmentNAVXYTHETAV::SetGoal(double x, double y, double theta, double v){
	int goalx = CONTXY2DISC(x, EnvNAVXYTHETAVCfg.cellsize_m);
	int goaly = CONTXY2DISC(y, EnvNAVXYTHETAVCfg.cellsize_m);
	int goaltheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVCfg.NumThetaDirs);
	int goalv = ContV2Disc(v, EnvNAVXYTHETAVCfg.velocities);

	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f (%d %d %d %d)\n", x, y, theta, v, goalx, goaly, goaltheta, goalv);

	if (!IsWithinMapCell(goalx, goaly)) {
		SBPL_ERROR("ERROR: trying to set a goal cell %d %d that is outside of map\n", goalx, goaly);
		return -1;
	}

	if (!IsValidConfiguration(goalx, goaly, goaltheta, goalv)) {
		SBPL_PRINTF("WARNING: goal configuration is invalid\n");
	}

	EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(goalx, goaly, goaltheta, goalv)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(goalx, goaly, goaltheta, goalv);
	}

	//need to recompute start heuristics?
	if (EnvNAVXYTHETAV.goalstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true; //because termination condition may not plan all the way to the new goal
		bNeedtoRecomputeGoalHeuristics = true; //because goal heuristics change
	}

	EnvNAVXYTHETAV.goalstateid = OutHashEntry->stateID;

	EnvNAVXYTHETAVCfg.EndX_c = goalx;
	EnvNAVXYTHETAVCfg.EndY_c = goaly;
	EnvNAVXYTHETAVCfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVCfg.EndV = goalv;

	return EnvNAVXYTHETAV.goalstateid;
}

void EnvironmentNAVXYTHETAV::GetCoordFromState(int stateID, int& x, int& y, int& theta, int& v) const{
	EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[stateID];
	x = HashEntry->x;
	y = HashEntry->y;
	theta = HashEntry->theta;
	v = HashEntry->v;
}

int EnvironmentNAVXYTHETAV::GetStateFromCoord(int x, int y, int theta, int v){
	EnvNAVXYTHETAVHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(x, y, theta, v)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(x, y, theta, v);
	}
	return OutHashEntry->stateID;
}

void EnvironmentNAVXYTHETAV::ConvertStateIDPathintoXYThetaVPath(vector<int>* stateIDPath, vector<sbpl_xy_theta_v_pt_t>* xythetavPath){
	vector<EnvNAVXYTHETAVAction_t*> actionV;
	vector<int> CostV;
	vector<int> SuccIDV;
	int targetx_c, targety_c, targettheta_c, targetv_c;
	int sourcex_c, sourcey_c, sourcetheta_c, sourcev_c;

	SBPL_PRINTF("checks=%ld\n", checks);

	xythetavPath->clear();

#if DEBUG
	SBPL_FPRINTF(fDeb, "converting stateid path into coordinates:\n");
#endif

	for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
		int sourceID = stateIDPath->at(pind);
		int targetID = stateIDPath->at(pind + 1);

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c);
#endif

		//get successors and pick the target via the cheapest action
		SuccIDV.clear();
		CostV.clear();
		actionV.clear();
		GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

		int bestcost = INFINITECOST;
		int bestsind = -1;

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c);
		GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c);
		SBPL_FPRINTF(fDeb, "looking for %d %d %d %d -> %d %d %d %d (numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c,
					targetx_c, targety_c, targettheta_c, targetv_c, (int)SuccIDV.size());
#endif

		for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
#if DEBUG
			int x_c, y_c, theta_c, v_c;
			GetCoordFromState(SuccIDV[sind], x_c, y_c, theta_c, v_c);
			SBPL_FPRINTF(fDeb, "succ: %d %d %d %d\n", x_c, y_c, theta_c, v_c);
#endif

			if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
				bestcost = CostV[sind];
				bestsind = sind;
			}
		}
		if (bestsind == -1) {
			SBPL_ERROR("ERROR: successor not found for transition:\n");
			GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c);
			GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c);
			SBPL_PRINTF("%d %d %d %d -> %d %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, targetx_c, targety_c,
						targettheta_c, targetv_c);
			throw new SBPL_Exception();
		}

		//now push in the actual path
		//int sourcex_c, sourcey_c, sourcetheta_c;
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c);
		double x, y;
		x = DISCXY2CONT(sourcex_c, EnvNAVXYTHETAVCfg.cellsize_m);
		y = DISCXY2CONT(sourcey_c, EnvNAVXYTHETAVCfg.cellsize_m);
		
		for (int ipind = 0; ipind < ((int)actionV[bestsind]->intermptV.size()) - 1; ipind++) {
			//translate appropriately
			sbpl_xy_theta_v_pt_t intermpt = actionV[bestsind]->intermptV[ipind];
			intermpt.x += x;
			intermpt.y += y;

#if DEBUG
			int nx = CONTXY2DISC(intermpt.x, EnvNAVXYTHETAVCfg.cellsize_m);
			int ny = CONTXY2DISC(intermpt.y, EnvNAVXYTHETAVCfg.cellsize_m);
			SBPL_FPRINTF(fDeb, "%.3f %.3f %.3f %.3f (%d %d %d %d cost=%d) ",
				intermpt.x, intermpt.y, intermpt.theta, intermpt.v,
				nx, ny, ContTheta2DiscNotUnif(intermpt.theta, EnvNAVXYTHETALATCfg.NumThetaDirs),
				ConvV2Disc(intermpt.v, EnvNAVXYTHETAVCfg.velocities), EnvNAVXYTHETALATCfg.Grid2D[nx][ny]);
			if(ipind == 0) SBPL_FPRINTF(fDeb, "first (heur=%d)\n", GetStartHeuristic(sourceID));
			else SBPL_FPRINTF(fDeb, "\n");
#endif

			//store
			xythetavPath->push_back(intermpt);
		}
	}
}

/* Useful functions */
double EnvironmentNAVXYTHETAV::EuclideanDistance_m(int x1, int y1, int x2, int y2)
{
	int sqdist = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return EnvNAVXYTHETAVCfg.cellsize_m * sqrt((double)sqdist);
}

bool EnvironmentNAVXYTHETAV::IsValidCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVCfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVCfg.EnvHeight_c &&
			EnvNAVXYTHETAVCfg.Grid2D[x][y] < EnvNAVXYTHETAVCfg.obsthresh);
}

int EnvironmentNAVXYTHETAV::GetActionCost(int sourceX, int sourceY, int sourceTheta, int sourceV, EnvNAVXYTHETAVAction_t* action){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_cell_t interm4Dcell;
	int i;

	//TODO - go over bounding box (minpt and maxpt) to test validity and skip
	//testing boundaries below, also order intersect cells so that the four
	//farthest pts go first

	if (!IsValidCell(sourceX, sourceY)) return INFINITECOST;
	if (!IsValidCell(sourceX + action->dX, sourceY + action->dY)) return INFINITECOST;

	if (EnvNAVXYTHETAVCfg.Grid2D[sourceX + action->dX][sourceY + action->dY] >=
		EnvNAVXYTHETAVCfg.cost_inscribed_thresh) 
	{
		return INFINITECOST;
	}

	//need to iterate over discretized center cells and compute cost based on them
	unsigned char maxcellcost = 0;
	for (i = 0; i < (int)action->interm3DcellsV.size(); i++) {
		interm4Dcell = action->interm3DcellsV.at(i);
		interm4Dcell.x = interm4Dcell.x + sourceX;
		interm4Dcell.y = interm4Dcell.y + sourceY;

		if (interm4Dcell.x < 0 || interm4Dcell.x >= EnvNAVXYTHETAVCfg.EnvWidth_c || interm4Dcell.y < 0
			|| interm4Dcell.y >= EnvNAVXYTHETAVCfg.EnvHeight_c) return INFINITECOST;

		maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVCfg.Grid2D[interm4Dcell.x][interm4Dcell.y]);

		//check that the robot is NOT in the cell at which there is no valid orientation
		if (maxcellcost >= EnvNAVXYTHETAVCfg.cost_inscribed_thresh) return INFINITECOST;
	}

	//check collisions that for the particular footprint orientation along the action
	if (EnvNAVXYTHETAVCfg.FootprintPolygon.size() > 1 && (int)maxcellcost >=
		EnvNAVXYTHETAVCfg.cost_possibly_circumscribed_thresh)
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
	maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVCfg.Grid2D[sourceX][sourceY]);
	int currentmaxcost =
			(int)__max(maxcellcost, EnvNAVXYTHETAVCfg.Grid2D[sourceX + action->dX][sourceY + action->dY]);
			
	return action->cost * (currentmaxcost + 1); //use cell cost as multiplicative factor
}

void EnvironmentNAVXYTHETAV::SetConfiguration(int width, int height, const unsigned char * mapdata, int startx, int starty,
								int starttheta, int startv, int goalx, int goaly, int goaltheta, int goalv,
							double cellsize_m, const vector<sbpl_2Dpt_t>& perimeterptsV){
	EnvNAVXYTHETAVCfg.EnvWidth_c = width;
	EnvNAVXYTHETAVCfg.EnvHeight_c = height;
	EnvNAVXYTHETAVCfg.StartX_c = startx;
	EnvNAVXYTHETAVCfg.StartY_c = starty;
	EnvNAVXYTHETAVCfg.StartTheta = starttheta;

	if(!IsWithinMapCell(EnvNAVXYTHETAVCfg.StartX_c, EnvNAVXYTHETAVCfg.StartY_c)){
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.StartTheta < 0 || EnvNAVXYTHETAVCfg.StartTheta >= EnvNAVXYTHETAVCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVCfg.StartV < 0 || EnvNAVXYTHETAVCfg.StartV >= EnvNAVXYTHETAVCfg.numV){
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVCfg.EndX_c = goalx;
	EnvNAVXYTHETAVCfg.EndY_c = goaly;
	EnvNAVXYTHETAVCfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVCfg.EndV = goalv;

	if(!IsWithinMapCell(EnvNAVXYTHETAVCfg.EndX_c, EnvNAVXYTHETAVCfg.EndY_c)){
		SBPL_ERROR("ERROR: illegal goal coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVCfg.EndTheta < 0 || EnvNAVXYTHETAVCfg.EndTheta >= EnvNAVXYTHETAVCfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVCfg.EndV < 0 || EnvNAVXYTHETAVCfg.EndV >= EnvNAVXYTHETAVCfg.numV){
		SBPL_ERROR("ERROR: illegal goal coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVCfg.FootprintPolygon = perimeterptsV;

	EnvNAVXYTHETAVCfg.cellsize_m = cellsize_m;
	
	//allocate the 2D environment
	EnvNAVXYTHETAVCfg.Grid2D = new unsigned char*[EnvNAVXYTHETAVCfg.EnvWidth_c];
	for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVCfg.Grid2D[x] = new unsigned char[EnvNAVXYTHETAVCfg.EnvHeight_c];
	}

	//environment:
	if (0 == mapdata) {
		for (int y = 0; y < EnvNAVXYTHETAVCfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVCfg.Grid2D[x][y] = 0;
			}
		}
	}
	else {
		for (int y = 0; y < EnvNAVXYTHETAVCfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVCfg.Grid2D[x][EnvNAVXYTHETAVCfg.EnvHeight_c-1-y] = mapdata[x + y * width];
			}
		}
	}
}

//------------------------------------------------------------------------------