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

#include <sbpl/discrete_space_information/environment_navxythetavomega.h>
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

EnvironmentNAVXYTHETAVOMEGA::EnvironmentNAVXYTHETAVOMEGA(){
	EnvNAVXYTHETAVOMEGACfg.obsthresh = NAVXYTHETAVOMEGA_DEFAULTOBSTHRESHOLD;
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh = EnvNAVXYTHETAVOMEGACfg.obsthresh; 
	//the value that pretty much makes it disabled
	EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh = -1; 

	grid2Dsearchfromstart = NULL;
	grid2Dsearchfromgoal = NULL;
	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;
	
	iteration = 0;

	EnvNAVXYTHETAVOMEGA.bInitialized = false;

	EnvNAVXYTHETAVOMEGACfg.NumThetaDirs = NAVXYTHETAVOMEGA_DEFAULTTHETADIRS;

	//no memory allocated in cfg yet
	EnvNAVXYTHETAVOMEGACfg.Grid2D = NULL;
	EnvNAVXYTHETAVOMEGACfg.ActionsV = NULL;
	EnvNAVXYTHETAVOMEGACfg.PredActionsV = NULL;
}

EnvironmentNAVXYTHETAVOMEGA::~EnvironmentNAVXYTHETAVOMEGA(){
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("destroying XYTHETAVOMEGA\n");
#endif
	if (grid2Dsearchfromstart != NULL) delete grid2Dsearchfromstart;
	grid2Dsearchfromstart = NULL;

	if (grid2Dsearchfromgoal != NULL) delete grid2Dsearchfromgoal;
	grid2Dsearchfromgoal = NULL;

	if (EnvNAVXYTHETAVOMEGACfg.Grid2D != NULL) {
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++)
			delete[] EnvNAVXYTHETAVOMEGACfg.Grid2D[x];
		delete[] EnvNAVXYTHETAVOMEGACfg.Grid2D;
		EnvNAVXYTHETAVOMEGACfg.Grid2D = NULL;
	}

	//delete actions
	if (EnvNAVXYTHETAVOMEGACfg.ActionsV != NULL) {
		for (int ind = 0; ind < (int)EnvNAVXYTHETAVOMEGACfg.ActionsV->size(); ind++)
			delete EnvNAVXYTHETAVOMEGACfg.ActionsV->at(ind);
		delete EnvNAVXYTHETAVOMEGACfg.ActionsV;
		EnvNAVXYTHETAVOMEGACfg.ActionsV = NULL;
	}
	if (EnvNAVXYTHETAVOMEGACfg.PredActionsV != NULL) {
		delete[] EnvNAVXYTHETAVOMEGACfg.PredActionsV;
		EnvNAVXYTHETAVOMEGACfg.PredActionsV = NULL;
	}
	
	//delete the states themselves first
	for (int i = 0; i < (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size(); i++) {
		delete EnvNAVXYTHETAVOMEGA.StateID2DataTable.at(i);
		EnvNAVXYTHETAVOMEGA.StateID2DataTable.at(i) = NULL;
	}
	EnvNAVXYTHETAVOMEGA.StateID2DataTable.clear();

	//delete hashtable
	if (EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable != NULL) {
		delete[] EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable;
		EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable = NULL;
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
unsigned int EnvironmentNAVXYTHETAVOMEGA::GETHASHBIN(unsigned int x, unsigned int y, unsigned int theta, unsigned int v, unsigned int omega)
{
	return inthash((inthash(x) + (inthash(y) << 1) + (inthash(theta) << 2) + (inthash(v) << 3) + (inthash(omega) << 4))) &
		(EnvNAVXYTHETAVOMEGA.HashTableSize - 1);
}

void EnvironmentNAVXYTHETAVOMEGA::PrintHashTableHist()
{
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	int s0 = 0, s1 = 0, s50 = 0, s100 = 0, s200 = 0, s300 = 0, slarge = 0;

	for (int j = 0; j < (int)EnvNAVXYTHETAVOMEGA.HashTableSize; j++) {
		if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() == 0)
			s0++;
		else if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() < 50)
			s1++;
		else if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() < 100)
			s50++;
		else if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() < 200)
			s100++;
		else if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() < 300)
			s200++;
		else if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[j].size() < 400)
			s300++;
		else
			slarge++;
	}
	SBPL_PRINTF("hash table histogram: 0:%d, <50:%d, <100:%d, <200:%d, <300:%d, <400:%d >400:%d\n", s0, s1, s50, s100,
				s200, s300, slarge);
#endif
}

void EnvironmentNAVXYTHETAVOMEGA::ReadConfiguration(FILE* fCfg)
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
	EnvNAVXYTHETAVOMEGACfg.EnvWidth_c = atoi(sTemp);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (discretization)\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EnvHeight_c = atoi(sTemp);

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
	
	EnvNAVXYTHETAVOMEGACfg.NumThetaDirs = atoi(sTemp);

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
	EnvNAVXYTHETAVOMEGACfg.numV = atoi(sTemp);
	
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
	for(int i=0;i<EnvNAVXYTHETAVOMEGACfg.numV;i++){
		if (fscanf(fCfg, "%s", sTemp) != 1) {
			SBPL_ERROR("ERROR: ran out of env file early\n");
			throw new SBPL_Exception();
		}
		
		EnvNAVXYTHETAVOMEGACfg.velocities.push_back(atof(sTemp));
	}
	
	// Scan for optional NumOmegas parameter.
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "NumOmegas:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	EnvNAVXYTHETAVOMEGACfg.numOmega = atoi(sTemp);
	
	//angular velocities: 
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early (obsthresh)\n");
		throw new SBPL_Exception();
	}
	strcpy(sTemp1, "Omegas:");
	if (strcmp(sTemp1, sTemp) != 0) {
		SBPL_ERROR("ERROR: configuration file has incorrect format\n");
		SBPL_PRINTF("Expected %s got %s\n", sTemp1, sTemp);
		SBPL_PRINTF("see existing examples of env files for the right format of heading\n");
		throw new SBPL_Exception();
	}
	
	/* SUBSTITUTE VECTOR WITH POINTER */
	for(int i=0;i<EnvNAVXYTHETAVOMEGACfg.numOmega;i++){
		if (fscanf(fCfg, "%s", sTemp) != 1) {
			SBPL_ERROR("ERROR: ran out of env file early\n");
			throw new SBPL_Exception();
		}
		
		EnvNAVXYTHETAVOMEGACfg.omegas.push_back(atof(sTemp));
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
	EnvNAVXYTHETAVOMEGACfg.obsthresh = atoi(sTemp);
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("obsthresh = %d\n", EnvNAVXYTHETAVOMEGACfg.obsthresh);
#endif

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
	EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh = atoi(sTemp);
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("cost_inscribed_thresh = %d\n", EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh);
#endif

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
	EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh = atoi(sTemp);
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("cost_possibly_circumscribed_thresh = %d\n", EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh);
#endif

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
	EnvNAVXYTHETAVOMEGACfg.cellsize_m = atof(sTemp);

	//start(meters,rads,m/s,rad/s):
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.StartX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.StartY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.StartTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.StartV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.velocities);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.StartOmega = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.omegas);

	if (EnvNAVXYTHETAVOMEGACfg.StartX_c < 0 || EnvNAVXYTHETAVOMEGACfg.StartX_c >= EnvNAVXYTHETAVOMEGACfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.StartY_c < 0 || EnvNAVXYTHETAVOMEGACfg.StartY_c >= EnvNAVXYTHETAVOMEGACfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.StartTheta < 0 || EnvNAVXYTHETAVOMEGACfg.StartTheta >= EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.StartV < 0 || EnvNAVXYTHETAVOMEGACfg.StartV >= EnvNAVXYTHETAVOMEGACfg.numV) {
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.StartOmega < 0 || EnvNAVXYTHETAVOMEGACfg.StartOmega >= EnvNAVXYTHETAVOMEGACfg.numOmega) {
		SBPL_ERROR("ERROR: illegal start coordinates for omega\n");
		throw new SBPL_Exception();
	}

	//end(meters,rads,m/s,rad/s):
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EndX_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EndY_c = CONTXY2DISC(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EndTheta = ContTheta2DiscNotUnif(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EndV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.velocities);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVOMEGACfg.EndOmega = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVOMEGACfg.omegas);

	if (EnvNAVXYTHETAVOMEGACfg.EndX_c < 0 || EnvNAVXYTHETAVOMEGACfg.EndX_c >= EnvNAVXYTHETAVOMEGACfg.EnvWidth_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.EndY_c < 0 || EnvNAVXYTHETAVOMEGACfg.EndY_c >= EnvNAVXYTHETAVOMEGACfg.EnvHeight_c) {
		SBPL_ERROR("ERROR: illegal end coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.EndTheta < 0 || EnvNAVXYTHETAVOMEGACfg.EndTheta >= EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.EndV < 0 || EnvNAVXYTHETAVOMEGACfg.EndV >= EnvNAVXYTHETAVOMEGACfg.numV) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.EndOmega < 0 || EnvNAVXYTHETAVOMEGACfg.EndOmega >= EnvNAVXYTHETAVOMEGACfg.numOmega) {
		SBPL_ERROR("ERROR: illegal start coordinates for omega\n");
		throw new SBPL_Exception();
	}

	//allocate the 2D environment
	EnvNAVXYTHETAVOMEGACfg.Grid2D = new double*[EnvNAVXYTHETAVOMEGACfg.EnvWidth_c];
	for (x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVOMEGACfg.Grid2D[x] = new double[EnvNAVXYTHETAVOMEGACfg.EnvHeight_c];
	}

	//environment:
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	
	double fTemp;
	for (y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1; y >= 0; y--)
		for (x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
			if (fscanf(fCfg, "%lf", &fTemp) != 1) {
				SBPL_ERROR("ERROR: incorrect format of config file\n");
				throw new SBPL_Exception();
			}
			EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] = fTemp;
		}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c, EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.StartOmega)){
		SBPL_ERROR("ERROR: invalid start state\n");
		throw new SBPL_Exception();
	}
	
	if(!IsValidConfiguration(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c, EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.EndOmega)){
		SBPL_ERROR("ERROR: invalid end state\n");
		throw new SBPL_Exception();
	}
}

EnvNAVXYTHETAVOMEGAHashEntry_t* EnvironmentNAVXYTHETAVOMEGA::GetHashEntry(unsigned int x, unsigned int y, unsigned int theta, unsigned int v, unsigned int omega)
{
	//clock_t currenttime = clock();

	int binid = GETHASHBIN(x, y, theta, v, omega);

#if DEBUG
	if ((int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid].size() > 500)
	{
		SBPL_PRINTF("WARNING: Hash table has a bin %d (X1=%d X2=%d X3=%d X4=%d) of size %d\n",
					binid, X1, X2, X3, X4, (int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid].size());

		PrintHashTableHist();
	}
#endif

	//iterate over the states in the bin and select the perfect match
	for (int ind = 0; ind < (int)EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid].size(); ind++) {
		if (EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind]->x == x &&
			EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind]->y == y &&
			EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind]->theta == theta &&
			EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind]->v == v &&
			EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind]->omega == omega)
		{
			//time_gethash += clock()-currenttime;
			return EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[binid][ind];
		}
	}

	//time_gethash += clock()-currenttime;

	return NULL;
}

EnvNAVXYTHETAVOMEGAHashEntry_t* EnvironmentNAVXYTHETAVOMEGA::CreateNewHashEntry(unsigned int x, unsigned int y, unsigned int theta,
													unsigned int v, unsigned int omega)
{
	int i;

	//clock_t currenttime = clock();

	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = new EnvNAVXYTHETAVOMEGAHashEntry_t;

	HashEntry->x = x;
	HashEntry->y = y;
	HashEntry->theta = theta;
	HashEntry->v = v;
	HashEntry->omega = omega;
	HashEntry->iteration = 0;

	HashEntry->stateID = EnvNAVXYTHETAVOMEGA.StateID2DataTable.size();

	//insert into the tables
	EnvNAVXYTHETAVOMEGA.StateID2DataTable.push_back(HashEntry);

	//get the hash table bin
	i = GETHASHBIN(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->omega);

	//insert the entry into the bin
	EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable[i].push_back(HashEntry);

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

void EnvironmentNAVXYTHETAVOMEGA::CreateStartandGoalStates()
{
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry;

	//create start state
	unsigned int x = 0;
	unsigned int y = 0;
	unsigned int theta = 0;
	unsigned int v = 0;
	unsigned int omega = 0;
	HashEntry = CreateNewHashEntry(x, y, theta, v, omega);
	EnvNAVXYTHETAVOMEGA.startstateid = HashEntry->stateID;

	//create goal state
	x = y = theta = v = omega = 1;
	HashEntry = CreateNewHashEntry(x, y, theta, v, omega);
	EnvNAVXYTHETAVOMEGA.goalstateid = HashEntry->stateID;
}

void EnvironmentNAVXYTHETAVOMEGA::ComputeReplanningDataforAction(EnvNAVXYTHETAVOMEGAAction_t* action){
	int j;

	//iterate over all the cells involved in the action
	sbpl_xy_theta_v_omega_cell_t startcell5d, endcell5d;
	for (int i = 0; i < (int)action->intersectingcellsV.size(); i++) {
		//compute the translated affected search Pose - what state has an
		//outgoing action whose intersecting cell is at 0,0
		startcell5d.omega = action->startomega;
		startcell5d.v = action->startv;
		startcell5d.theta = action->starttheta;
		startcell5d.x = -action->intersectingcellsV.at(i).x;
		startcell5d.y = -action->intersectingcellsV.at(i).y;

		//compute the translated affected search Pose - what state has an
		//incoming action whose intersecting cell is at 0,0
		endcell5d.omega = action->endomega;
		endcell5d.v = action->endv;
		endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
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
	startcell5d.omega = action->startomega;
	startcell5d.v = action->startv;
	startcell5d.theta = action->starttheta;
	startcell5d.x = -0;
	startcell5d.y = -0;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell5d.omega = action->endomega;
	endcell5d.v = action->endv;
	endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
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
	startcell5d.omega = action->startomega;
	startcell5d.v = action->startv;
	startcell5d.theta = action->starttheta;
	startcell5d.x = -action->dX;
	startcell5d.y = -action->dY;

	//compute the translated affected search Pose - what state has an incoming action whose intersecting cell is at 0,0
	endcell5d.omega = action->endomega;
	endcell5d.v = action->endv;
	endcell5d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
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

void EnvironmentNAVXYTHETAVOMEGA::ComputeReplanningData()
{
	//iterate over all actions
	//Omega speed
	for(int sind = 0;sind < EnvNAVXYTHETAVOMEGACfg.numOmega;sind++){
		//velocities
		for(int vind = 0; vind < EnvNAVXYTHETAVOMEGACfg.numV; vind++){
			//orientations
			for (int tind = 0; tind < EnvNAVXYTHETAVOMEGACfg.NumThetaDirs; tind++) {
				int idx = sind*EnvNAVXYTHETAVOMEGACfg.numV*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+vind*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+tind;
				//actions
				for (int aind = 0; aind < (int)(EnvNAVXYTHETAVOMEGACfg.ActionsV->at(idx)->size()); aind++) {
					//compute replanning data for this action 
					ComputeReplanningDataforAction(&EnvNAVXYTHETAVOMEGACfg.ActionsV->at(idx)->at(aind));
				}
			}
		}
	}
}

void EnvironmentNAVXYTHETAVOMEGA::PrecomputeActionswithCompleteMotionPrimitive(vector<SBPL_xythetavomega_mprimitive>* motionprimitiveV){
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("Pre-computing action data using motion primitives for every pair velocity/angle...\n");
#endif
	EnvNAVXYTHETAVOMEGACfg.ActionsV = new vector<vector<EnvNAVXYTHETAVOMEGAAction_t> *>();
	EnvNAVXYTHETAVOMEGACfg.PredActionsV = new vector<EnvNAVXYTHETAVOMEGAActionIndex_t> [EnvNAVXYTHETAVOMEGACfg.numOmega*EnvNAVXYTHETAVOMEGACfg.numV*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs];
	
	vector<sbpl_2Dcell_t> footprint;

	/*
	* if (motionprimitiveV->size() % EnvNAVXYTHETALATCfg.NumThetaDirs != 0) {
	* SBPL_ERROR("ERROR: motionprimitives should be uniform across actions\n");
	* throw new SBPL_Exception();
	* }
	*/

	//EnvNAVXYTHETAVCfg.actionwidth = ((int)motionprimitiveV->size()) / EnvNAVXYTHETALATCfg.NumThetaDirs;

	//Allocate vectors
	for(int sind = 0; sind < EnvNAVXYTHETAVOMEGACfg.numOmega; sind++)
		for(int vind = 0; vind < EnvNAVXYTHETAVOMEGACfg.numV; vind++)
			for (int tind = 0; tind < EnvNAVXYTHETAVOMEGACfg.NumThetaDirs; tind++){
				EnvNAVXYTHETAVOMEGACfg.ActionsV->push_back(new vector<EnvNAVXYTHETAVOMEGAAction_t>());
			}
	
	//iterate over source angles|speeds|angles
	int maxnumofactions = 0;
	for(int sind = 0; sind < EnvNAVXYTHETAVOMEGACfg.numOmega; sind++){
		for(int vind = 0; vind < EnvNAVXYTHETAVOMEGACfg.numV; vind++){
			for (int tind = 0; tind < EnvNAVXYTHETAVOMEGACfg.NumThetaDirs; tind++) {
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
				SBPL_PRINTF("pre-computing for pair (angle, speed, angle) (%d,%d,%d) out of (%d,%d,%d)\n", sind, vind, tind, EnvNAVXYTHETAVOMEGACfg.numOmega, EnvNAVXYTHETAVOMEGACfg.numV, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
#endif

				//compute current index
				int vector_index = sind*EnvNAVXYTHETAVOMEGACfg.numV*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+vind*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+tind;
				
				//EnvNAVXYTHETAVCfg.ActionsV->at(vector_index) = new vector<EnvNAVXYTHETAVAction_t>();

				//compute sourcepose
				sbpl_xy_theta_v_omega_pt_t sourcepose;
				sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
				sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
				sourcepose.theta = DiscTheta2ContNotUnif(tind, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
				sourcepose.v = DiscV2Cont(vind, EnvNAVXYTHETAVOMEGACfg.velocities);
				sourcepose.omega = DiscV2Cont(sind, EnvNAVXYTHETAVOMEGACfg.omegas);

				//iterate over motion primitives
				int numofactions = 0;
				int aind = -1;
				for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
					//find a motion primitive for this angle
					if (motionprimitiveV->at(mind).start_theta_disc != tind || motionprimitiveV->at(mind).start_v_disc != vind || motionprimitiveV->at(mind).start_omega_disc != sind) continue;

					aind++;
					numofactions++;
					EnvNAVXYTHETAVOMEGAAction_t element_to_add;
					
					//action index
					element_to_add.aind = aind;

					//start angle
					element_to_add.starttheta = tind;
					
					//start velocity
					element_to_add.startv = vind;
					
					//start omega speed
					element_to_add.startomega = sind;

					//compute dislocation
					element_to_add.endtheta = motionprimitiveV->at(mind).endcell.theta;
					element_to_add.endv = motionprimitiveV->at(mind).endcell.v;
					element_to_add.endomega = motionprimitiveV->at(mind).endcell.omega;
					element_to_add.dX = motionprimitiveV->at(mind).endcell.x;
					element_to_add.dY = motionprimitiveV->at(mind).endcell.y;

					//compute and store interm points as well as intersecting cells
					element_to_add.intersectingcellsV.clear();
					element_to_add.intermptV.clear();
					element_to_add.interm3DcellsV.clear();

					sbpl_xy_theta_v_omega_cell_t previnterm3Dcell;
					previnterm3Dcell.x = 0;
					previnterm3Dcell.y = 0;

					// Compute all the intersected cells for this action (intermptV and interm3DcellsV)
					for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
						sbpl_xy_theta_v_omega_pt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
						element_to_add.intermptV.push_back(intermpt);

						// also compute the intermediate discrete cells if not there already
						sbpl_xy_theta_v_omega_pt_t pose;
						pose.x = intermpt.x + sourcepose.x;
						pose.y = intermpt.y + sourcepose.y;
						pose.theta = intermpt.theta;
						pose.v = intermpt.v;
						pose.omega = intermpt.omega;

						sbpl_xy_theta_v_omega_cell_t intermediate2dCell;
						intermediate2dCell.x = CONTXY2DISC(pose.x, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
						intermediate2dCell.y = CONTXY2DISC(pose.y, EnvNAVXYTHETAVOMEGACfg.cellsize_m);

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
					element_to_add.cost = NAVXYTHETAVOMEGA_COSTMULT_MTOMM * motionprimitiveV->at(mind).additionalactioncostmult;
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
					get_2d_motion_cells(EnvNAVXYTHETAVOMEGACfg.FootprintPolygon, motionprimitiveV->at(mind).intermptV, &element_to_add.intersectingcellsV, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
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
					int postheta = (element_to_add.endtheta >= 0)?element_to_add.endtheta:element_to_add.endtheta+EnvNAVXYTHETAVOMEGACfg.NumThetaDirs;
					int targetindex = element_to_add.endomega*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs*EnvNAVXYTHETAVOMEGACfg.numV+element_to_add.endv*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+postheta;
					EnvNAVXYTHETAVOMEGAActionIndex_t idx;
					idx.vector_index = vector_index;
					idx.aind = aind;
					EnvNAVXYTHETAVOMEGACfg.ActionsV->at(vector_index)->push_back(element_to_add);
					EnvNAVXYTHETAVOMEGACfg.PredActionsV[targetindex].push_back(idx);
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

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("done pre-computing action data based on motion primitives\n");
#endif
}

void EnvironmentNAVXYTHETAVOMEGA::InitializeEnvConfig(vector<SBPL_xythetavomega_mprimitive>* motionprimitiveV)
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
	get_2d_footprint_cells(EnvNAVXYTHETAVOMEGACfg.FootprintPolygon, &footprint, temppose, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("number of cells in footprint of the robot = %d\n", (unsigned int)footprint.size());

	for (vector<sbpl_2Dcell_t>::iterator it = footprint.begin(); it != footprint.end(); ++it) {
		SBPL_PRINTF("Footprint cell at (%d, %d)\n", it->x, it->y);
	}
#endif

#if DEBUG
	SBPL_FPRINTF(fDeb, "footprint cells (size=%d):\n", (int)footprint.size());
	for(int i = 0; i < (int) footprint.size(); i++)
	{
		SBPL_FPRINTF(fDeb, "%d %d (cont: %.3f %.3f)\n", footprint.at(i).x, footprint.at(i).y,
					DISCXY2CONT(footprint.at(i).x, EnvNAVXYTHETAVOMEGACfg.cellsize_m),
					DISCXY2CONT(footprint.at(i).y, EnvNAVXYTHETAVOMEGACfg.cellsize_m));
	}
#endif
	
	PrecomputeActionswithCompleteMotionPrimitive(motionprimitiveV);
}

bool EnvironmentNAVXYTHETAVOMEGA::InitGeneral(vector<SBPL_xythetavomega_mprimitive>* motionprimitiveV){
	//Initialize other parameters of the environment
	InitializeEnvConfig(motionprimitiveV);

	//initialize Environment
	InitializeEnvironment();

	//pre-compute heuristics
	ComputeHeuristicValues();

	return true;
}

void EnvironmentNAVXYTHETAVOMEGA::InitializeEnvironment()
{
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry;

	int maxsize = EnvNAVXYTHETAVOMEGACfg.EnvWidth_c * EnvNAVXYTHETAVOMEGACfg.EnvHeight_c * EnvNAVXYTHETAVOMEGACfg.NumThetaDirs * EnvNAVXYTHETAVOMEGACfg.numV * EnvNAVXYTHETAVOMEGACfg.numOmega;

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("environment stores states in hashtable\n");
#endif

	//initialize the map from Data to StateID
	//Maximum hash table size
	EnvNAVXYTHETAVOMEGA.HashTableSize = EnvNAVXYTHETAVOMEGACfg.NumThetaDirs * EnvNAVXYTHETAVOMEGACfg.EnvWidth_c * EnvNAVXYTHETAVOMEGACfg.EnvHeight_c * 8 * 8;/*EnvNAVXYTHETAVCfg.numV;*/ //should be power of two - REVIEW -> USE DYNAMIC PARAMETERS
	EnvNAVXYTHETAVOMEGA.Data2StateIDHashTable = new vector<EnvNAVXYTHETAVOMEGAHashEntry_t*> [EnvNAVXYTHETAVOMEGA.HashTableSize];

	//initialize the map from StateID to Data
	EnvNAVXYTHETAVOMEGA.StateID2DataTable.clear();

	//create start state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c,
										EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.StartOmega)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c,
												EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.StartOmega);
	}
	EnvNAVXYTHETAVOMEGA.startstateid = HashEntry->stateID;

	//create goal state
	if ((HashEntry = GetHashEntry(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c,
										EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.EndOmega)) == NULL) {
		//have to create a new entry
		HashEntry = CreateNewHashEntry(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c,
												EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.EndOmega);
	}
	EnvNAVXYTHETAVOMEGA.goalstateid = HashEntry->stateID;

	//initialized
	EnvNAVXYTHETAVOMEGA.bInitialized = true;
}

bool EnvironmentNAVXYTHETAVOMEGA::ReadMotionPrimitives(FILE* fMotPrims){
	char sTemp[1024], sExpected[1024];
	float fTemp;
	int dTemp;
	int totalNumofActions = 0;

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("Reading in motion primitives...");
#endif

	//read in the resolution
	strcpy(sExpected, "resolution_m:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
	if (fabs(fTemp - EnvNAVXYTHETAVOMEGACfg.cellsize_m) > ERR_EPS) {
		SBPL_ERROR("ERROR: invalid resolution %f (instead of %f) in the dynamics file\n", fTemp,
				EnvNAVXYTHETAVOMEGACfg.cellsize_m);
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
	if (dTemp != EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
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
	if (dTemp != EnvNAVXYTHETAVOMEGACfg.numV) {
		SBPL_ERROR("ERROR: invalid velocity resolution %d velocities (instead of %d velocities) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVOMEGACfg.numV);
		return false;
	}
	
	//read in the velocity values
	strcpy(sExpected, "velocities:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	for(int i=0;i<EnvNAVXYTHETAVOMEGACfg.numV;i++){
		if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
		if (fTemp != EnvNAVXYTHETAVOMEGACfg.velocities[i]) {
			SBPL_ERROR("ERROR: invalid velocity value %f velocity (instead of %f velocity) in the motion primitives file\n",
					fTemp, EnvNAVXYTHETAVOMEGACfg.velocities[i]);
			return false;
		}
	}
	
	//read in the angular velocity resolution
	strcpy(sExpected, "numberofomegas:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fMotPrims, "%d", &dTemp) == 0) return false;
	if (dTemp != EnvNAVXYTHETAVOMEGACfg.numOmega) {
		SBPL_ERROR("ERROR: invalid angular velocity resolution %d (instead of %d angular velocity) in the motion primitives file\n",
				dTemp, EnvNAVXYTHETAVOMEGACfg.numOmega);
		return false;
	}
	
	//read in the angular velocity values
	strcpy(sExpected, "omegas:");
	if (fscanf(fMotPrims, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	for(int i=0;i<EnvNAVXYTHETAVOMEGACfg.numOmega;i++){
		if (fscanf(fMotPrims, "%f", &fTemp) == 0) return false;
		if (fTemp != EnvNAVXYTHETAVOMEGACfg.omegas[i]) {
			SBPL_ERROR("ERROR: invalid velocity value %f velocity (instead of %f velocity) in the motion primitives file\n",
					fTemp, EnvNAVXYTHETAVOMEGACfg.omegas[i]);
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
		SBPL_xythetavomega_mprimitive motprim;

		if (ReadSingleMotionPrimitive(&motprim, fMotPrims) == false) return false;

		EnvNAVXYTHETAVOMEGACfg.mprimV.push_back(motprim);
	}
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("done reading of motion primitives");
#endif

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::ReadSingleMotionPrimitive(SBPL_xythetavomega_mprimitive* pMotPrim, FILE* fIn){
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
		SBPL_ERROR("ERROR reading start omega\n");
		return false;
	}
	pMotPrim->start_omega_disc = dTemp;

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
		sbpl_xy_theta_v_omega_pt_t intermpose;
		if (ReadSinglePose(&intermpose, fIn) == false) {
			SBPL_ERROR("ERROR: failed to read in intermediate poses\n");
			return false;
		}
		pMotPrim->intermptV.push_back(intermpose);
	}

	//check that the last pose corresponds correctly to the last pose
	sbpl_xy_theta_v_omega_pt_t sourcepose;
	sourcepose.x = DISCXY2CONT(0, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	sourcepose.y = DISCXY2CONT(0, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	sourcepose.theta = DiscTheta2ContNotUnif(pMotPrim->start_theta_disc, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	sourcepose.v = DiscV2Cont(pMotPrim->start_v_disc, EnvNAVXYTHETAVOMEGACfg.velocities);
	sourcepose.omega = DiscV2Cont(pMotPrim->start_omega_disc, EnvNAVXYTHETAVOMEGACfg.omegas);
	double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
	double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
	double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;
	double mp_endv_ms = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v;
	double mp_endomega_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].omega;
	int endx_disc = CONTXY2DISC(mp_endx_m, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int endy_disc = CONTXY2DISC(mp_endy_m, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int endtheta_disc = ContTheta2DiscNotUnif(mp_endtheta_rad, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	int endv_disc = ContV2Disc(mp_endv_ms, EnvNAVXYTHETAVOMEGACfg.velocities);
	int endomega_disc = ContV2Disc(mp_endomega_rad, EnvNAVXYTHETAVOMEGACfg.omegas);
	if (endx_disc != pMotPrim->endcell.x || endy_disc != pMotPrim->endcell.y || endtheta_disc != pMotPrim->endcell.theta || endv_disc != pMotPrim->endcell.v || endomega_disc != pMotPrim->endcell.omega) {
		SBPL_ERROR( "ERROR: incorrect primitive %d with startangle=%d, startv=%d and startomega=%d "
				"last interm point %f %f %f %f %f does not match end pose %d %d %d %d %d\n",
				pMotPrim->motprimID, pMotPrim->start_theta_disc, pMotPrim->start_v_disc, pMotPrim->start_omega_disc,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v,
				pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].omega,
				pMotPrim->endcell.x, pMotPrim->endcell.y,
				pMotPrim->endcell.theta, pMotPrim->endcell.v, pMotPrim->endcell.omega);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::ReadSingleCell(sbpl_xy_theta_v_omega_cell_t* cell, FILE* fIn){
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
	cell->omega = atoi(sTemp);
	
	//normalize the angle
	cell->theta = NORMALIZEDISCTHETA(cell->theta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::ReadSinglePose(sbpl_xy_theta_v_omega_pt_t* pose, FILE* fIn){
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
	pose->omega = atof(sTemp);

	pose->theta = normalizeAngle(pose->theta);
	pose->omega = pose->omega;

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

void EnvironmentNAVXYTHETAVOMEGA::EnsureHeuristicsUpdated(bool bGoalHeuristics){
	if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
		//Create temporary map
		unsigned char ** map;
	
		map = new unsigned char*[EnvNAVXYTHETAVOMEGACfg.EnvWidth_c];
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
			map[x] = new unsigned char[EnvNAVXYTHETAVOMEGACfg.EnvHeight_c];
		}
		
		for (int y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1; y >= 0; y--)
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				map[x][y] = (int)(EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y]);
			}
		
		grid2Dsearchfromstart->search(map, EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c,
									EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeStartHeuristics = false;
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVOMEGACfg.EndX_c,
																				EnvNAVXYTHETAVOMEGACfg.EndY_c)/EnvNAVXYTHETAVOMEGACfg.velocities.at(EnvNAVXYTHETAVOMEGACfg.velocities.size()-1)));
#endif
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.EndX_c,
																				EnvNAVXYTHETAVCfg.EndY_c)));*/
		
		//Destroy temporary map
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++)
			delete[] map[x];
		delete[] map;
	}

	if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
		//Create temporary map
		unsigned char ** map;
	
		map = new unsigned char*[EnvNAVXYTHETAVOMEGACfg.EnvWidth_c];
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
			map[x] = new unsigned char[EnvNAVXYTHETAVOMEGACfg.EnvHeight_c];
		}
		
		for (int y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1; y >= 0; y--)
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				map[x][y] = (int)(EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y]);
			}
			
		grid2Dsearchfromgoal->search(map, EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh,
									EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c,
									EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c,
									SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
		bNeedtoRecomputeGoalHeuristics = false;
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
		SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVOMEGACfg.StartX_c,
																				EnvNAVXYTHETAVOMEGACfg.StartY_c)/EnvNAVXYTHETAVOMEGACfg.velocities.at(EnvNAVXYTHETAVOMEGACfg.velocities.size()-1)));
#endif
		
		/*SBPL_PRINTF("2dsolcost_infullunits=%d\n",
					(int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETAVCfg.StartX_c,
																				EnvNAVXYTHETAVCfg.StartY_c)));*/
		
		//Destroy temporary map
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++)
			delete[] map[x];
		delete[] map;
	}
}

void EnvironmentNAVXYTHETAVOMEGA::ComputeHeuristicValues()
{
	//whatever necessary pre-computation of heuristic values is done here
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("Precomputing heuristics\n");
#endif
	
	//allocated 2D grid searches
	grid2Dsearchfromstart = new SBPL2DGridSearch(EnvNAVXYTHETAVOMEGACfg.EnvWidth_c, EnvNAVXYTHETAVOMEGACfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	grid2Dsearchfromgoal = new SBPL2DGridSearch(EnvNAVXYTHETAVOMEGACfg.EnvWidth_c, EnvNAVXYTHETAVOMEGACfg.EnvHeight_c,
												(float)EnvNAVXYTHETAVOMEGACfg.cellsize_m);

	//set OPEN type to sliding buckets
	grid2Dsearchfromstart->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	grid2Dsearchfromgoal->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
	
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("done\n");
#endif
}

void EnvironmentNAVXYTHETAVOMEGA::PrintHeuristicValues(){
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
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

	for (int y = 0; y < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c; y++) {
		for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
			if (grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y) < INFINITECOST)
				SBPL_FPRINTF(fHeur, "%5d ", grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y));
			else
				SBPL_FPRINTF(fHeur, "XXXXX ");
		}
		SBPL_FPRINTF(fHeur, "\n");
	}
	
	SBPL_FCLOSE(fHeur);
#endif
}

//-----------interface with outside functions-----------------------------------

bool EnvironmentNAVXYTHETAVOMEGA::InitializeEnv(const char* sEnvFile, const vector<sbpl_2Dpt_t>& perimeterptsV, const char* sMotPrimFile)
{
	EnvNAVXYTHETAVOMEGACfg.FootprintPolygon = perimeterptsV;

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
		InitGeneral(&EnvNAVXYTHETAVOMEGACfg.mprimV);
		fclose(fMotPrim);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("size of env: %d by %d\n", EnvNAVXYTHETAVOMEGACfg.EnvWidth_c, EnvNAVXYTHETAVOMEGACfg.EnvHeight_c);
#endif

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::InitializeEnv(const char* sEnvFile){
	return false;
}

bool EnvironmentNAVXYTHETAVOMEGA::SetEnvParameter(const char* parameter, int value){
	if (EnvNAVXYTHETAVOMEGA.bInitialized == true) {
		SBPL_ERROR("ERROR: all parameters must be set before initialization of the environment\n");
		return false;
	}

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("setting parameter %s to %d\n", parameter, value);
#endif

	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh = (unsigned char)value;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh = value;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		if (value < 0 || value > 255) {
			SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
			return false;
		}
		EnvNAVXYTHETAVOMEGACfg.obsthresh = (unsigned char)value;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		return false;
	}

	return true;
}

int EnvironmentNAVXYTHETAVOMEGA::GetEnvParameter(const char* parameter)
{
	if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh;
	}
	else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
		return (int)EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh;
	}
	else if (strcmp(parameter, "cost_obsthresh") == 0) {
		return (int)EnvNAVXYTHETAVOMEGACfg.obsthresh;
	}
	else {
		SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
		throw new SBPL_Exception();
	}
}

bool EnvironmentNAVXYTHETAVOMEGA::InitializeMDPCfg(MDPConfig *MDPCfg)
{
	//initialize MDPCfg with the start and goal ids
	MDPCfg->goalstateid = EnvNAVXYTHETAVOMEGA.goalstateid;
	MDPCfg->startstateid = EnvNAVXYTHETAVOMEGA.startstateid;

	return true;
}

int EnvironmentNAVXYTHETAVOMEGA::GetFromToHeuristic(int FromStateID, int ToStateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if(FromStateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size() || ToStateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size()){
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVOMEGA... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	//get X, Y for the state
	EnvNAVXYTHETAVOMEGAHashEntry_t* FromHashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[FromStateID];
	EnvNAVXYTHETAVOMEGAHashEntry_t* ToHashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[ToStateID];

	//TODO - check if one of the gridsearches already computed and then use it.

	return (int)(NAVXYTHETAVOMEGA_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y)/EnvNAVXYTHETAVOMEGACfg.velocities.at(EnvNAVXYTHETAVOMEGACfg.velocities.size()-1));
	//return (int)(NAVXYTHETAV_COSTMULT_MTOMM * EuclideanDistance_m(FromHashEntry->x, FromHashEntry->y, ToHashEntry->x, ToHashEntry->y));
}

int EnvironmentNAVXYTHETAVOMEGA::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVOMEGA... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[stateID];
	//computes distances from start state that is grid2D, so it is EndX_c EndY_c
	int h2D = grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y); 
	int hEuclid = (int)(NAVXYTHETAVOMEGA_COSTMULT_MTOMM * EuclideanDistance_m(HashEntry->x, HashEntry->y,
																		EnvNAVXYTHETAVOMEGACfg.EndX_c,
																		EnvNAVXYTHETAVOMEGACfg.EndY_c));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETAVOMEGACfg.velocities.at(EnvNAVXYTHETAVOMEGACfg.velocities.size()-1));
	//return (int)(((double)__max(h2D, hEuclid)));
}

int EnvironmentNAVXYTHETAVOMEGA::GetStartHeuristic(int stateID)
{
	
#if USE_HEUR==0
	return 0;
#endif

#if DEBUG
	if (stateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size()) {
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVOMEGA... function: stateID illegal\n");
		throw new SBPL_Exception();
	}
#endif

	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[stateID];
	int h2D = grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(HashEntry->x, HashEntry->y);
	int hEuclid = (int)(NAVXYTHETAVOMEGA_COSTMULT_MTOMM * EuclideanDistance_m(EnvNAVXYTHETAVOMEGACfg.StartX_c,
																		EnvNAVXYTHETAVOMEGACfg.StartY_c, HashEntry->x,
																		HashEntry->y));

	//define this function if it is used in the planner (heuristic backward search would use it)
	return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETAVOMEGACfg.velocities.at(EnvNAVXYTHETAVOMEGACfg.velocities.size()-1));
	//return (int)(((double)__max(h2D, hEuclid)));
}

void EnvironmentNAVXYTHETAVOMEGA::SetAllActionsandAllOutcomes(CMDPSTATE* state){
	int cost;
	
#if DEBUG
	if(state->StateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size())
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
	if (state->StateID == EnvNAVXYTHETAVOMEGA.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[state->StateID];
	
	int vector_index = HashEntry->omega*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs*EnvNAVXYTHETAVOMEGACfg.numV + HashEntry->v*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+HashEntry->theta;

	//iterate through actions
	for (int aind = 0; aind < EnvNAVXYTHETAVOMEGACfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVOMEGAAction_t* nav5daction = &EnvNAVXYTHETAVOMEGACfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav5daction->dX;
		int newY = HashEntry->y + nav5daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav5daction->endtheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
		int newV = nav5daction->endv;
		int newOmega = nav5daction->endomega;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->omega, nav5daction);
		if (cost >= INFINITECOST) continue;

		//add the action
		CMDPACTION* action = state->AddAction(aind);

#if TIME_DEBUG
		clock_t currenttime = clock();
#endif

		EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV, newOmega)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV, newOmega);
		}
		action->AddOutcome(OutHashEntry->stateID, cost, 1.0);

#if TIME_DEBUG
		time3_addallout += clock()-currenttime;
#endif
	}
}

void EnvironmentNAVXYTHETAVOMEGA::SetAllPreds(CMDPSTATE* state)
{
	//implement this if the planner needs access to predecessors

	SBPL_ERROR("ERROR in EnvNAVXYTHETAVOMEGA... function: SetAllPreds is undefined\n");
	throw new SBPL_Exception();
}

void EnvironmentNAVXYTHETAVOMEGA::GetSuccs(int SourceStateID, vector<int>* SuccIDV, vector<int>* CostV)
{
	GetSuccs(SourceStateID, SuccIDV, CostV, NULL);
}

void EnvironmentNAVXYTHETAVOMEGA::GetSuccs(int sourceStateID, vector<int>* succIDV, vector<int>* costV, vector<EnvNAVXYTHETAVOMEGAAction_t*>* actionindV){
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
	if (sourceStateID == EnvNAVXYTHETAVOMEGA.goalstateid) return;

	//get X, Y for the state
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[sourceStateID];
	
	int vector_index = (unsigned int)HashEntry->omega*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs*EnvNAVXYTHETAVOMEGACfg.numV + (unsigned int)HashEntry->v*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//iterate through actions
	for (aind = 0; aind < EnvNAVXYTHETAVOMEGACfg.ActionsV->at(vector_index)->size(); aind++) {
		EnvNAVXYTHETAVOMEGAAction_t* nav5daction = &EnvNAVXYTHETAVOMEGACfg.ActionsV->at(vector_index)->at(aind);
		int newX = HashEntry->x + nav5daction->dX;
		int newY = HashEntry->y + nav5daction->dY;
		int newTheta = NORMALIZEDISCTHETA(nav5daction->endtheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
		int newV = nav5daction->endv;
		int newOmega = nav5daction->endomega;

		//skip the invalid cells
		if (!IsValidCell(newX, newY)) continue;

		//get cost
		int cost = GetActionCost(HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->omega, nav5daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(newX, newY, newTheta, newV, newOmega)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(newX, newY, newTheta, newV, newOmega);
		}

		succIDV->push_back(OutHashEntry->stateID);
		costV->push_back(cost);
		if (actionindV != NULL) actionindV->push_back(nav5daction);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentNAVXYTHETAVOMEGA::GetPreds(int TargetStateID, vector<int>* PredIDV, vector<int>* CostV)
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
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[TargetStateID];
	int targetindex = (unsigned int)HashEntry->omega*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs*EnvNAVXYTHETAVOMEGACfg.numV + (unsigned int)HashEntry->v*EnvNAVXYTHETAVOMEGACfg.NumThetaDirs+(unsigned int)HashEntry->theta;

	//clear the preds array
	PredIDV->clear();
	CostV->clear();
	PredIDV->reserve(EnvNAVXYTHETAVOMEGACfg.PredActionsV[targetindex].size());
	CostV->reserve(EnvNAVXYTHETAVOMEGACfg.PredActionsV[targetindex].size());

	//iterate through actions
	vector<EnvNAVXYTHETAVOMEGAActionIndex_t>* actionsV = &EnvNAVXYTHETAVOMEGACfg.PredActionsV[targetindex];
	for (aind = 0; aind < (int)actionsV->size(); aind++) {

		EnvNAVXYTHETAVOMEGAAction_t* nav5daction = &EnvNAVXYTHETAVOMEGACfg.ActionsV->at(actionsV->at(aind).vector_index)->at(actionsV->at(aind).aind);

		int predX = HashEntry->x - nav5daction->dX;
		int predY = HashEntry->y - nav5daction->dY;
		int predTheta = nav5daction->starttheta;
		int predV = nav5daction->startv;
		int predOmega = nav5daction->startomega;

		//skip the invalid cells
		if (!IsValidCell(predX, predY)) continue;

		//get cost
		int cost = GetActionCost(predX, predY, predTheta, predV, predOmega, nav5daction);
		if (cost >= INFINITECOST) continue;

		EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
		if ((OutHashEntry = GetHashEntry(predX, predY, predTheta, predV, predOmega)) == NULL) {
			//have to create a new entry
			OutHashEntry = CreateNewHashEntry(predX, predY, predTheta, predV, predOmega);
		}

		PredIDV->push_back(OutHashEntry->stateID);
		CostV->push_back(cost);
	}

#if TIME_DEBUG
	time_getsuccs += clock()-currenttime;
#endif
}

int EnvironmentNAVXYTHETAVOMEGA::SizeofCreatedEnv()
{
	return (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size();
}

void EnvironmentNAVXYTHETAVOMEGA::PrintState(int stateID, bool bVerbose, FILE* fOut /*=NULL*/)
{
#if DEBUG
	if(stateID >= (int)EnvNAVXYTHETAVOMEGA.StateID2DataTable.size())
	{
		SBPL_ERROR("ERROR in EnvNAVXYTHETAVOMEGA... function: stateID illegal (2)\n");
		throw new SBPL_Exception();
	}
#endif
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	if (fOut == NULL) fOut = stdout;

	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[stateID];

	
	if (stateID == EnvNAVXYTHETAVOMEGA.goalstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a goal state\n");
	}
	
	if (stateID == EnvNAVXYTHETAVOMEGA.startstateid && bVerbose) {
		SBPL_FPRINTF(fOut, "the state is a start state\n");
	}

	SBPL_FPRINTF(fOut, "x=%d y=%d theta=%d v=%d omega=%d\n", HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v, HashEntry->omega);
#endif
}

void EnvironmentNAVXYTHETAVOMEGA::PrintEnv_Config(FILE* fOut)
{
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	if(fOut != NULL){
		fprintf(fOut, "discretization(cells): %d %d\n", EnvNAVXYTHETAVOMEGACfg.EnvWidth_c, EnvNAVXYTHETAVOMEGACfg.EnvHeight_c);
		fprintf(fOut, "NumThetaDirs: %d\n", EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
		fprintf(fOut, "NumV: %d\n", EnvNAVXYTHETAVOMEGACfg.numV);
		fprintf(fOut, "Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVOMEGACfg.numV; i++)
			fprintf(fOut, " %f", EnvNAVXYTHETAVOMEGACfg.velocities.at(i));
		fprintf(fOut, "NumOmegas: %d\n", EnvNAVXYTHETAVOMEGACfg.numOmega);
		fprintf(fOut, "Omegas:");
		for(int i = 0; i<EnvNAVXYTHETAVOMEGACfg.numOmega; i++)
			fprintf(fOut, " %f", EnvNAVXYTHETAVOMEGACfg.omegas.at(i));
		fprintf(fOut, "\n");
		fprintf(fOut, "obsthresh: %d\n", EnvNAVXYTHETAVOMEGACfg.obsthresh);
		fprintf(fOut, "cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "cellsize(meters): %f\n", EnvNAVXYTHETAVOMEGACfg.cellsize_m);
		fprintf(fOut, "cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh);
		fprintf(fOut, "start(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.velocities), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartOmega, EnvNAVXYTHETAVOMEGACfg.omegas));
		fprintf(fOut, "end(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.velocities), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndOmega, EnvNAVXYTHETAVOMEGACfg.omegas));
		fprintf(fOut, "environment:\n");
		
		for (int y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				fprintf(fOut, "%lf ", EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y]);
			}
			
			fprintf(fOut, "\n");
		}
	}
	else{
		printf("discretization(cells): %d %d\n", EnvNAVXYTHETAVOMEGACfg.EnvWidth_c, EnvNAVXYTHETAVOMEGACfg.EnvHeight_c);
		printf("NumThetaDirs: %d\n", EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
		printf("NumV: %d\n", EnvNAVXYTHETAVOMEGACfg.numV);
		printf("Velocities:");
		for(int i = 0; i<EnvNAVXYTHETAVOMEGACfg.numV; i++)
			printf(" %f", EnvNAVXYTHETAVOMEGACfg.velocities.at(i));
		printf("NumOmegas: %d\n", EnvNAVXYTHETAVOMEGACfg.numOmega);
		printf("Omegas:");
		for(int i = 0; i<EnvNAVXYTHETAVOMEGACfg.numOmega; i++)
			printf(" %f", EnvNAVXYTHETAVOMEGACfg.omegas.at(i));
		printf("\n");
		printf("obsthresh: %d\n", EnvNAVXYTHETAVOMEGACfg.obsthresh);
		printf("cost_inscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh);
		printf("cellsize(meters): %f\n", EnvNAVXYTHETAVOMEGACfg.cellsize_m);
		printf("cost_possibly_circumscribed_thresh: %d\n", EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh);
		printf("start(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.velocities), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartOmega, EnvNAVXYTHETAVOMEGACfg.omegas));
		printf("end(meters,rads,m/s): %f %f %f %f %f\n", DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m), DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.velocities), DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndOmega, EnvNAVXYTHETAVOMEGACfg.omegas));
		printf("environment:\n");
		
		for (int y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1; y >= 0; y--){
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				printf("%lf ", EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y]);
			}
			
			printf("\n");
		}
	}
#endif
}

bool EnvironmentNAVXYTHETAVOMEGA::InitializeEnv(int width, int height, int numthetadirs, int numv, int numomega, vector<double> velocities, vector<double> omegas,
							const double* mapdata,
							double startx, double starty, double starttheta, double startv, double startomega,
							double goalx, double goaly, double goaltheta, double goalv, double goalomega,
							const vector<sbpl_2Dpt_t>& perimeterptsV, double cellsize_m,
							unsigned char obsthresh, unsigned char cost_inscribed_thresh,
							unsigned char cost_possibly_circumscribed_thresh, const char* sMotPrimFile){
#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("env: initialize with width=%d height=%d start=%.3f %.3f %.3f %.3f %.3f "
				"goalx=%.3f %.3f %.3f %.3f %.3f cellsize=%.3f numthetadirs=%d numv=%d numomega=%d, obsthresh=%d cost_inscribed_thresh=%d cost_possibly_circumscribed_thresh=%d\n",
				width, height, startx, starty, starttheta, startv, startomega, goalx, goaly, goaltheta, goalv, goalomega, cellsize_m, numthetadirs, numv, numomega,
				obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh);

	SBPL_PRINTF("NOTE: goaltol parameters currently unused\n");

	SBPL_PRINTF("perimeter has size=%d\n", (unsigned int)perimeterptsV.size());

	for (int i = 0; i < (int)perimeterptsV.size(); i++) {
		SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, perimeterptsV.at(i).x, perimeterptsV.at(i).y);
	}
#endif

	EnvNAVXYTHETAVOMEGACfg.obsthresh = obsthresh;
	EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh = cost_inscribed_thresh;
	EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh = cost_possibly_circumscribed_thresh;
	EnvNAVXYTHETAVOMEGACfg.NumThetaDirs = numthetadirs;
	EnvNAVXYTHETAVOMEGACfg.numV = numv;
	EnvNAVXYTHETAVOMEGACfg.numOmega = numomega;
	
	EnvNAVXYTHETAVOMEGACfg.velocities.clear();
	for(int i=0;i<numv;i++){
		EnvNAVXYTHETAVOMEGACfg.velocities.push_back(velocities.at(i));
	}
	
	EnvNAVXYTHETAVOMEGACfg.omegas.clear();
	for(int i=0;i<numomega;i++){
		EnvNAVXYTHETAVOMEGACfg.omegas.push_back(omegas.at(i));
	}
	
	//TODO - need to set the tolerance as well

	SetConfiguration(width, height, mapdata, CONTXY2DISC(startx, cellsize_m), CONTXY2DISC(starty, cellsize_m),
					ContTheta2DiscNotUnif(starttheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs), ContV2Disc(startv, EnvNAVXYTHETAVOMEGACfg.velocities), ContV2Disc(startomega, EnvNAVXYTHETAVOMEGACfg.omegas),
					CONTXY2DISC(goalx, cellsize_m), CONTXY2DISC(goaly, cellsize_m), ContTheta2DiscNotUnif(goaltheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs),
					ContV2Disc(goalv, EnvNAVXYTHETAVOMEGACfg.velocities), ContV2Disc(goalomega, EnvNAVXYTHETAVOMEGACfg.omegas), cellsize_m, perimeterptsV);

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

	if (EnvNAVXYTHETAVOMEGACfg.mprimV.size() != 0) {
		InitGeneral(&EnvNAVXYTHETAVOMEGACfg.mprimV);
	}
	else{
		//InitGeneral( NULL);
		return false;
	}

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::UpdateCost(int x, int y, double newcost){
	EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] = newcost;

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

bool EnvironmentNAVXYTHETAVOMEGA::SetMap(const double * mapdata){
	for (int xind = 0; xind < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; xind++) {
		for (int yind = 0; yind < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c; yind++) {
			EnvNAVXYTHETAVOMEGACfg.Grid2D[xind][EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1-yind] = mapdata[xind + yind * EnvNAVXYTHETAVOMEGACfg.EnvWidth_c];
		}
	}

	bNeedtoRecomputeStartHeuristics = true;
	bNeedtoRecomputeGoalHeuristics = true;

	return true;
}

void EnvironmentNAVXYTHETAVOMEGA::GetPredsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV, vector<int> *preds_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_omega_cell_t affectedcell;
	EnvNAVXYTHETAVOMEGAHashEntry_t* affectedHashEntry;

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
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v, affectedcell.omega);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

void EnvironmentNAVXYTHETAVOMEGA::GetSuccsOfChangedEdges(vector<sbpl_2Dcell_t> const * changedcellsV,
										vector<int> *succs_of_changededgesIDV){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_omega_cell_t affectedcell;
	EnvNAVXYTHETAVOMEGAHashEntry_t* affectedHashEntry;

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
			affectedHashEntry = GetHashEntry(affectedcell.x, affectedcell.y, affectedcell.theta, affectedcell.v, affectedcell.omega);
			if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
				succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
				affectedHashEntry->iteration = iteration; //mark as already inserted
			}
		}
	}
}

bool EnvironmentNAVXYTHETAVOMEGA::IsObstacle(int x, int y){
	return (EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] >= EnvNAVXYTHETAVOMEGACfg.obsthresh);
}

bool EnvironmentNAVXYTHETAVOMEGA::IsValidConfiguration(int x, int y, int theta, int v, int omega){
	vector<sbpl_2Dcell_t> footprint;
	sbpl_xy_theta_v_omega_pt_t pose;

	//compute continuous pose
	pose.x = DISCXY2CONT(x, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	pose.y = DISCXY2CONT(y, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	pose.theta = DiscTheta2ContNotUnif(theta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	pose.v = DiscV2Cont(v, EnvNAVXYTHETAVOMEGACfg.velocities);
	pose.omega = DiscV2Cont(omega, EnvNAVXYTHETAVOMEGACfg.omegas);

	//compute footprint cells
	get_2d_footprint_cells(EnvNAVXYTHETAVOMEGACfg.FootprintPolygon, &footprint, pose, EnvNAVXYTHETAVOMEGACfg.cellsize_m);

	//iterate over all footprint cells
	for (int find = 0; find < (int)footprint.size(); find++) {
		int x = footprint.at(find).x;
		int y = footprint.at(find).y;

		if (x < 0 || x >= EnvNAVXYTHETAVOMEGACfg.EnvWidth_c || y < 0 || y >= EnvNAVXYTHETAVOMEGACfg.EnvHeight_c ||
			EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] >= EnvNAVXYTHETAVOMEGACfg.obsthresh)
		{
			return false;
		}
	}
	
	if(v < 0 || v >= EnvNAVXYTHETAVOMEGACfg.numV){
		return false;
	}
	
	if(omega < 0 || omega >= EnvNAVXYTHETAVOMEGACfg.numOmega){
		return false;
	}

	return true;
}

void EnvironmentNAVXYTHETAVOMEGA::GetEnvParams(int *size_x, int *size_y, int * numthetadirs, int * numv, int * numomega, vector<double> * velocities, vector<double> * omegas,
							double* startx, double* starty, double* starttheta, double * startv, double * startomega,
							double* goalx, double* goaly, double* goaltheta, double * goalv, double * goalomega, double* cellsize_m,
							unsigned char* cost_inscribed_thresh, unsigned char* cost_possibly_circumscribed_thresh,
							unsigned char* obsthresh, vector<SBPL_xythetavomega_mprimitive>* motionprimitiveV){
	
	*size_x = EnvNAVXYTHETAVOMEGACfg.EnvWidth_c;
	*size_y = EnvNAVXYTHETAVOMEGACfg.EnvHeight_c;

	*numthetadirs = EnvNAVXYTHETAVOMEGACfg.NumThetaDirs;
	*numv = EnvNAVXYTHETAVOMEGACfg.numV;
	for(int i;i<EnvNAVXYTHETAVOMEGACfg.velocities.size();i++){
		velocities->push_back(EnvNAVXYTHETAVOMEGACfg.velocities.at(i));
	}
	*numomega = EnvNAVXYTHETAVOMEGACfg.numOmega;
	for(int i;i<EnvNAVXYTHETAVOMEGACfg.omegas.size();i++){
		omegas->push_back(EnvNAVXYTHETAVOMEGACfg.omegas.at(i));
	}
	
	*startx = DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	*starty = DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.StartY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	*starttheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.StartTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	*startv = DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartV, EnvNAVXYTHETAVOMEGACfg.velocities);
	*startomega = DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.StartOmega, EnvNAVXYTHETAVOMEGACfg.omegas);
	*goalx = DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	*goaly = DISCXY2CONT(EnvNAVXYTHETAVOMEGACfg.EndY_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	*goaltheta = DiscTheta2ContNotUnif(EnvNAVXYTHETAVOMEGACfg.EndTheta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	*goalv = DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndV, EnvNAVXYTHETAVOMEGACfg.velocities);
	*goalomega = DiscV2Cont(EnvNAVXYTHETAVOMEGACfg.EndOmega, EnvNAVXYTHETAVOMEGACfg.omegas);

	*cellsize_m = EnvNAVXYTHETAVOMEGACfg.cellsize_m;
	
	*cost_inscribed_thresh = EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh;
	*cost_possibly_circumscribed_thresh = EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh;
	*obsthresh = EnvNAVXYTHETAVOMEGACfg.obsthresh;

	*motionprimitiveV = EnvNAVXYTHETAVOMEGACfg.mprimV;
}

const EnvNAVXYTHETAVOMEGAConfig_t* EnvironmentNAVXYTHETAVOMEGA::GetEnvNavConfig(){
	return &EnvNAVXYTHETAVOMEGACfg;
}

void EnvironmentNAVXYTHETAVOMEGA::PrintTimeStat(FILE* fOut){
	/*TO BE IMPLEMENTED*/
}

double EnvironmentNAVXYTHETAVOMEGA::GetMapCost(int x, int y){
	return EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y];
}

bool EnvironmentNAVXYTHETAVOMEGA::IsWithinMapCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c);
}

bool EnvironmentNAVXYTHETAVOMEGA::PoseContToDisc(double px, double py, double pth, double pv, double pom, int &ix, int &iy, int &ith, int &iv, int &iom) const{
	ix = CONTXY2DISC(px, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	iy = CONTXY2DISC(py, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	ith = ContTheta2DiscNotUnif(pth, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	iv = ContV2Disc(pv, EnvNAVXYTHETAVOMEGACfg.velocities);
	iom = ContV2Disc(pom, EnvNAVXYTHETAVOMEGACfg.omegas);
	
	return (pth >= -2 * PI_CONST) && (pth <= 2 * PI_CONST) && (ix >= 0 && ix < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVOMEGACfg.numV) && (iom >= 0) && (iom < EnvNAVXYTHETAVOMEGACfg.numOmega);
}

bool EnvironmentNAVXYTHETAVOMEGA::PoseDiscToCont(int ix, int iy, int ith, int iv, int iom, double &px, double &py, double &pth, double &pv, double &pom) const{
	px = DISCXY2CONT(ix, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	py = DISCXY2CONT(iy, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	pth = normalizeAngle(DiscTheta2ContNotUnif(ith, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs));
	pv = DiscV2Cont(iv, EnvNAVXYTHETAVOMEGACfg.velocities);
	pom = DiscV2Cont(iom, EnvNAVXYTHETAVOMEGACfg.omegas);
	
	return (ith >= 0) && (ith < EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) && (ix >= 0 && ix < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c && iy >= 0 && iy < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c) && (iv >= 0) && (iv < EnvNAVXYTHETAVOMEGACfg.numV) && (iom >= 0) && (iom < EnvNAVXYTHETAVOMEGACfg.numOmega);
}

int EnvironmentNAVXYTHETAVOMEGA::SetStart(double x, double y, double theta, double v, double omega){
	int startx = CONTXY2DISC(x, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int starty = CONTXY2DISC(y, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int starttheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	int startv = ContV2Disc(v, EnvNAVXYTHETAVOMEGACfg.velocities);
	int startomega = ContV2Disc(omega, EnvNAVXYTHETAVOMEGACfg.omegas);

	if (!IsWithinMapCell(startx, starty)) {
		SBPL_ERROR("ERROR: trying to set a start cell %d %d that is outside of map\n", startx, starty);
		return -1;
	}

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d)\n", x, y, theta, v, omega, startx, starty, starttheta, startv, startomega);
#endif

	if (!IsValidConfiguration(startx, starty, starttheta, startv, startomega)) {
		SBPL_PRINTF("WARNING: start configuration %d %d %d %d %d is invalid\n", startx, starty, starttheta, startv, startomega);
	}

	EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(startx, starty, starttheta, startv, startomega)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(startx, starty, starttheta, startv, startomega);
	}

	if (EnvNAVXYTHETAVOMEGA.startstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true;
		//because termination condition can be not all states TODO - make it dependent on term. condition
		bNeedtoRecomputeGoalHeuristics = true; 
	}

	//set start
	EnvNAVXYTHETAVOMEGA.startstateid = OutHashEntry->stateID;
	EnvNAVXYTHETAVOMEGACfg.StartX_c = startx;
	EnvNAVXYTHETAVOMEGACfg.StartY_c = starty;
	EnvNAVXYTHETAVOMEGACfg.StartTheta = starttheta;
	EnvNAVXYTHETAVOMEGACfg.StartV = startv;
	EnvNAVXYTHETAVOMEGACfg.StartOmega = startomega;

	return EnvNAVXYTHETAVOMEGA.startstateid;
}

int EnvironmentNAVXYTHETAVOMEGA::SetGoal(double x, double y, double theta, double v, double omega){
	int goalx = CONTXY2DISC(x, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int goaly = CONTXY2DISC(y, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
	int goaltheta = ContTheta2DiscNotUnif(theta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs);
	int goalv = ContV2Disc(v, EnvNAVXYTHETAVOMEGACfg.velocities);
	int goalomega = ContV2Disc(omega, EnvNAVXYTHETAVOMEGACfg.omegas);

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("env: setting start to %.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d)\n", x, y, theta, v, omega, goalx, goaly, goaltheta, goalv, goalomega);
#endif

	if (!IsWithinMapCell(goalx, goaly)) {
		SBPL_ERROR("ERROR: trying to set a goal cell %d %d that is outside of map\n", goalx, goaly);
		return -1;
	}

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	if (!IsValidConfiguration(goalx, goaly, goaltheta, goalv, goalomega)) {
		SBPL_PRINTF("WARNING: goal configuration is invalid\n");
	}
#endif

	EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(goalx, goaly, goaltheta, goalv, goalomega)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(goalx, goaly, goaltheta, goalv, goalomega);
	}

	//need to recompute start heuristics?
	if (EnvNAVXYTHETAVOMEGA.goalstateid != OutHashEntry->stateID) {
		bNeedtoRecomputeStartHeuristics = true; //because termination condition may not plan all the way to the new goal
		bNeedtoRecomputeGoalHeuristics = true; //because goal heuristics change
	}

	EnvNAVXYTHETAVOMEGA.goalstateid = OutHashEntry->stateID;

	EnvNAVXYTHETAVOMEGACfg.EndX_c = goalx;
	EnvNAVXYTHETAVOMEGACfg.EndY_c = goaly;
	EnvNAVXYTHETAVOMEGACfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVOMEGACfg.EndV = goalv;
	EnvNAVXYTHETAVOMEGACfg.EndOmega = goalomega;

	return EnvNAVXYTHETAVOMEGA.goalstateid;
}

void EnvironmentNAVXYTHETAVOMEGA::GetCoordFromState(int stateID, int& x, int& y, int& theta, int& v, int& omega) const{
	EnvNAVXYTHETAVOMEGAHashEntry_t* HashEntry = EnvNAVXYTHETAVOMEGA.StateID2DataTable[stateID];
	x = HashEntry->x;
	y = HashEntry->y;
	theta = HashEntry->theta;
	v = HashEntry->v;
	omega = HashEntry->omega;
}

int EnvironmentNAVXYTHETAVOMEGA::GetStateFromCoord(int x, int y, int theta, int v, int omega){
	EnvNAVXYTHETAVOMEGAHashEntry_t* OutHashEntry;
	if ((OutHashEntry = GetHashEntry(x, y, theta, v, omega)) == NULL) {
		//have to create a new entry
		OutHashEntry = CreateNewHashEntry(x, y, theta, v, omega);
	}
	return OutHashEntry->stateID;
}

void EnvironmentNAVXYTHETAVOMEGA::ConvertStateIDPathintoXYThetaVOmegaPath(vector<int>* stateIDPath, vector<sbpl_xy_theta_v_omega_pt_t>* xythetavPath){
	vector<EnvNAVXYTHETAVOMEGAAction_t*> actionV;
	vector<int> CostV;
	vector<int> SuccIDV;
	int targetx_c, targety_c, targettheta_c, targetv_c, targetomega_c;
	int sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c;

#if NAVXYTHETAVOMEGA_PERFORMANCE_TEST == 0
	SBPL_PRINTF("checks=%ld\n", checks);
#endif

	xythetavPath->clear();

#if DEBUG
	SBPL_FPRINTF(fDeb, "converting stateid path into coordinates:\n");
#endif

	for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
		int sourceID = stateIDPath->at(pind);
		int targetID = stateIDPath->at(pind + 1);

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c);
#endif

		//get successors and pick the target via the cheapest action
		SuccIDV.clear();
		CostV.clear();
		actionV.clear();
		GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

		int bestcost = INFINITECOST;
		int bestsind = -1;

#if DEBUG
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c);
		GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c, targetomega_c);
		SBPL_FPRINTF(fDeb, "looking for %d %d %d %d %d -> %d %d %d %d %d (numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c,
					targetx_c, targety_c, targettheta_c, targetv_c, targetomega_c, (int)SuccIDV.size());
#endif

		for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
#if DEBUG
			int x_c, y_c, theta_c, v_c omega_c;
			GetCoordFromState(SuccIDV[sind], x_c, y_c, theta_c, v_c, omega_c);
			SBPL_FPRINTF(fDeb, "succ: %d %d %d %d %d\n", x_c, y_c, theta_c, v_c, omega_c);
#endif

			if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
				bestcost = CostV[sind];
				bestsind = sind;
			}
		}
		if (bestsind == -1) {
			SBPL_ERROR("ERROR: successor not found for transition:\n");
			GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c);
			GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c, targetv_c, targetomega_c);
			SBPL_PRINTF("%d %d %d %d %d -> %d %d %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c, targetx_c, targety_c,
						targettheta_c, targetv_c, targetomega_c);
			throw new SBPL_Exception();
		}

		//now push in the actual path
		//int sourcex_c, sourcey_c, sourcetheta_c;
		GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c, sourcev_c, sourceomega_c);
		double x, y;
		x = DISCXY2CONT(sourcex_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
		y = DISCXY2CONT(sourcey_c, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
		
		for (int ipind = 0; ipind < ((int)actionV[bestsind]->intermptV.size()) - 1; ipind++) {
			//translate appropriately
			sbpl_xy_theta_v_omega_pt_t intermpt = actionV[bestsind]->intermptV[ipind];
			intermpt.x += x;
			intermpt.y += y;

#if DEBUG
			int nx = CONTXY2DISC(intermpt.x, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
			int ny = CONTXY2DISC(intermpt.y, EnvNAVXYTHETAVOMEGACfg.cellsize_m);
			SBPL_FPRINTF(fDeb, "%.3f %.3f %.3f %.3f %.3f (%d %d %d %d %d cost=%d) ",
				intermpt.x, intermpt.y, intermpt.theta, intermpt.v, intermpt.omega,
				nx, ny, ContTheta2DiscNotUnif(intermpt.theta, EnvNAVXYTHETAVOMEGACfg.NumThetaDirs),
				ContV2Disc(intermpt.v, EnvNAVXYTHETAVOMEGACfg.velocities), ContV2Disc(intermpt.omega, EnvNAVXYTHETAVOMEGACfg.omegas), EnvNAVXYTHETAVOMEGACfg.Grid2D[nx][ny]);
			if(ipind == 0) SBPL_FPRINTF(fDeb, "first (heur=%d)\n", GetStartHeuristic(sourceID));
			else SBPL_FPRINTF(fDeb, "\n");
#endif

			//store
			xythetavPath->push_back(intermpt);
		}
	}
}

/* Useful functions */
double EnvironmentNAVXYTHETAVOMEGA::EuclideanDistance_m(int x1, int y1, int x2, int y2)
{
	int sqdist = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return EnvNAVXYTHETAVOMEGACfg.cellsize_m * sqrt((double)sqdist);
}

bool EnvironmentNAVXYTHETAVOMEGA::IsValidCell(int x, int y){
	return (x >= 0 && x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c && y >= 0 && y < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c &&
			EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] < EnvNAVXYTHETAVOMEGACfg.obsthresh);
}

int EnvironmentNAVXYTHETAVOMEGA::GetActionCost(int sourceX, int sourceY, int sourceTheta, int sourceV, int sourceOmega, EnvNAVXYTHETAVOMEGAAction_t* action){
	sbpl_2Dcell_t cell;
	sbpl_xy_theta_v_omega_cell_t interm5Dcell;
	int i;

	//TODO - go over bounding box (minpt and maxpt) to test validity and skip
	//testing boundaries below, also order intersect cells so that the four
	//farthest pts go first

	if (!IsValidCell(sourceX, sourceY)) return INFINITECOST;
	if (!IsValidCell(sourceX + action->dX, sourceY + action->dY)) return INFINITECOST;

	if (EnvNAVXYTHETAVOMEGACfg.Grid2D[sourceX + action->dX][sourceY + action->dY] >=
		EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh) 
	{
		return INFINITECOST;
	}

	//need to iterate over discretized center cells and compute cost based on them
	double maxcellcost = 0;
	for (i = 0; i < (int)action->interm3DcellsV.size(); i++) {
		interm5Dcell = action->interm3DcellsV.at(i);
		interm5Dcell.x = interm5Dcell.x + sourceX;
		interm5Dcell.y = interm5Dcell.y + sourceY;

		if (interm5Dcell.x < 0 || interm5Dcell.x >= EnvNAVXYTHETAVOMEGACfg.EnvWidth_c || interm5Dcell.y < 0
			|| interm5Dcell.y >= EnvNAVXYTHETAVOMEGACfg.EnvHeight_c) return INFINITECOST;

		maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVOMEGACfg.Grid2D[interm5Dcell.x][interm5Dcell.y]);

		//check that the robot is NOT in the cell at which there is no valid orientation
		if (maxcellcost >= EnvNAVXYTHETAVOMEGACfg.cost_inscribed_thresh) return INFINITECOST;
	}

	//check collisions that for the particular footprint orientation along the action
	if (EnvNAVXYTHETAVOMEGACfg.FootprintPolygon.size() > 1 && (int)maxcellcost >=
		EnvNAVXYTHETAVOMEGACfg.cost_possibly_circumscribed_thresh)
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
	maxcellcost = __max(maxcellcost, EnvNAVXYTHETAVOMEGACfg.Grid2D[sourceX][sourceY]);
	int currentmaxcost =
			(int)__max(maxcellcost, EnvNAVXYTHETAVOMEGACfg.Grid2D[sourceX + action->dX][sourceY + action->dY]);
			
	return action->cost * (currentmaxcost + 1); //use cell cost as multiplicative factor
}

void EnvironmentNAVXYTHETAVOMEGA::SetConfiguration(int width, int height, const double * mapdata, int startx, int starty,
								int starttheta, int startv, int startomega, int goalx, int goaly, int goaltheta, int goalv, int goalomega,
							double cellsize_m, const vector<sbpl_2Dpt_t>& perimeterptsV){
	EnvNAVXYTHETAVOMEGACfg.EnvWidth_c = width;
	EnvNAVXYTHETAVOMEGACfg.EnvHeight_c = height;
	EnvNAVXYTHETAVOMEGACfg.StartX_c = startx;
	EnvNAVXYTHETAVOMEGACfg.StartY_c = starty;
	EnvNAVXYTHETAVOMEGACfg.StartTheta = starttheta;
	EnvNAVXYTHETAVOMEGACfg.StartV = startv;
	EnvNAVXYTHETAVOMEGACfg.StartOmega = startomega;

	if(!IsWithinMapCell(EnvNAVXYTHETAVOMEGACfg.StartX_c, EnvNAVXYTHETAVOMEGACfg.StartY_c)){
		SBPL_ERROR("ERROR: illegal start coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.StartTheta < 0 || EnvNAVXYTHETAVOMEGACfg.StartTheta >= EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal start coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVOMEGACfg.StartV < 0 || EnvNAVXYTHETAVOMEGACfg.StartV >= EnvNAVXYTHETAVOMEGACfg.numV){
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVOMEGACfg.StartOmega < 0 || EnvNAVXYTHETAVOMEGACfg.StartOmega >= EnvNAVXYTHETAVOMEGACfg.numOmega){
		SBPL_ERROR("ERROR: illegal start coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVOMEGACfg.EndX_c = goalx;
	EnvNAVXYTHETAVOMEGACfg.EndY_c = goaly;
	EnvNAVXYTHETAVOMEGACfg.EndTheta = goaltheta;
	EnvNAVXYTHETAVOMEGACfg.EndV = goalv;
	EnvNAVXYTHETAVOMEGACfg.EndOmega = goalomega;

	if(!IsWithinMapCell(EnvNAVXYTHETAVOMEGACfg.EndX_c, EnvNAVXYTHETAVOMEGACfg.EndY_c)){
		SBPL_ERROR("ERROR: illegal goal coordinates\n");
		throw new SBPL_Exception();
	}
	if (EnvNAVXYTHETAVOMEGACfg.EndTheta < 0 || EnvNAVXYTHETAVOMEGACfg.EndTheta >= EnvNAVXYTHETAVOMEGACfg.NumThetaDirs) {
		SBPL_ERROR("ERROR: illegal goal coordinates for theta\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVOMEGACfg.EndV < 0 || EnvNAVXYTHETAVOMEGACfg.EndV >= EnvNAVXYTHETAVOMEGACfg.numV){
		SBPL_ERROR("ERROR: illegal goal coordinates for v\n");
		throw new SBPL_Exception();
	}
	if(EnvNAVXYTHETAVOMEGACfg.EndOmega < 0 || EnvNAVXYTHETAVOMEGACfg.EndOmega >= EnvNAVXYTHETAVOMEGACfg.numOmega){
		SBPL_ERROR("ERROR: illegal goal coordinates for v\n");
		throw new SBPL_Exception();
	}

	EnvNAVXYTHETAVOMEGACfg.FootprintPolygon = perimeterptsV;

	EnvNAVXYTHETAVOMEGACfg.cellsize_m = cellsize_m;
	
	//allocate the 2D environment
	EnvNAVXYTHETAVOMEGACfg.Grid2D = new double*[EnvNAVXYTHETAVOMEGACfg.EnvWidth_c];
	for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
		EnvNAVXYTHETAVOMEGACfg.Grid2D[x] = new double[EnvNAVXYTHETAVOMEGACfg.EnvHeight_c];
	}

	//environment:
	if (0 == mapdata) {
		for (int y = 0; y < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVOMEGACfg.Grid2D[x][y] = 0;
			}
		}
	}
	else {
		for (int y = 0; y < EnvNAVXYTHETAVOMEGACfg.EnvHeight_c; y++) {
			for (int x = 0; x < EnvNAVXYTHETAVOMEGACfg.EnvWidth_c; x++) {
				EnvNAVXYTHETAVOMEGACfg.Grid2D[x][EnvNAVXYTHETAVOMEGACfg.EnvHeight_c-1-y] = mapdata[x + y * width];
			}
		}
	}
}

//------------------------------------------------------------------------------