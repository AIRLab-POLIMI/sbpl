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
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

//extern clock_t time3_addallout;
//extern clock_t time_gethash;
//extern clock_t time_createhash;

//function prototypes

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

	//start(meters,rads):
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
	EnvNAVXYTHETAVCfg.StartTheta = ContTheta2Disc(atof(sTemp), EnvNAVXYTHETAVCfg.NumThetaDirs);
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

	//end(meters,rads):
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
	EnvNAVXYTHETAVCfg.EndTheta = ContTheta2Disc(atof(sTemp), EnvNAVXYTHETAVCfg.NumThetaDirs);
	if (fscanf(fCfg, "%s", sTemp) != 1) {
		SBPL_ERROR("ERROR: ran out of env file early\n");
		throw new SBPL_Exception();
	}
	EnvNAVXYTHETAVCfg.StartV = ContV2Disc(atof(sTemp), EnvNAVXYTHETAVCfg.velocities);

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
	/* NON MI CONVINCE, VEDERE COSA SUCCEDE */
	for (y = 0; y < EnvNAVXYTHETAVCfg.EnvHeight_c; y++)
		for (x = 0; x < EnvNAVXYTHETAVCfg.EnvWidth_c; x++) {
			if (fscanf(fCfg, "%d", &dTemp) != 1) {
				SBPL_ERROR("ERROR: incorrect format of config file\n");
				throw new SBPL_Exception();
			}
			EnvNAVXYTHETAVCfg.Grid2D[x][y] = dTemp;
		}
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

void EnvironmentNAVXYTHETAV::PrecomputeActionswithCompleteMotionPrimitive(std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV){
	
}

/*
 * SAME AS ABOVE FOR UNSIGNED
 */
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

/*
 * SAME AS ABOVE FOR UNSIGNED
 */
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

void EnvironmentNAVXYTHETAV::InitializeEnvironment()
{

    //initialize the map from Data to StateID
    //Maximum hash table size
    EnvNAVXYTHETAV.HashTableSize = 16 * 1024 * 1024 * 8; //should be power of two - REVIEW -> USE DYNAMIC PARAMETERS
    EnvNAVXYTHETAV.Data2StateIDHashTable = new vector<EnvNAVXYTHETAVHashEntry_t*> [EnvNAVXYTHETAV.HashTableSize];

    //initialize the map from StateID to Coord
    EnvNAVXYTHETAV.StateID2DataTable.clear();

    //create start and goal states
    CreateStartandGoalStates();
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
	strcpy(sExpected, "startpose_disc:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading startangle\n");
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading startangle\n");
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading startangle\n");
		return false;
	}
	pMotPrim->start_theta_disc = dTemp;
	if (fscanf(fIn, "%d", &dTemp) == 0) {
		SBPL_ERROR("ERROR reading startangle\n");
		return false;
	}
	pMotPrim->start_v_disc = dTemp;

	//read in end cell
	strcpy(sExpected, "endpose_disc:");
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
	/*
	strcpy(sExpected, "additionalactioncostmult:");
	if (fscanf(fIn, "%s", sTemp) == 0) return false;
	if (strcmp(sTemp, sExpected) != 0) {
		SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
		return false;
	}
	if (fscanf(fIn, "%d", &dTemp) != 1) return false;
	pMotPrim->additionalactioncostmult = dTemp;
	*/

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
	sourcepose.theta = DiscTheta2Cont(pMotPrim->start_theta_disc, EnvNAVXYTHETAVCfg.NumThetaDirs);
	sourcepose.v = DiscV2Cont(pMotPrim->start_v_disc, EnvNAVXYTHETAVCfg.velocities);
	double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
	double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
	double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;
	double mp_endv_ms = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].v;
	int endx_disc = CONTXY2DISC(mp_endx_m, EnvNAVXYTHETAVCfg.cellsize_m);
	int endy_disc = CONTXY2DISC(mp_endy_m, EnvNAVXYTHETAVCfg.cellsize_m);
	int endtheta_disc = ContTheta2Disc(mp_endtheta_rad, EnvNAVXYTHETAVCfg.NumThetaDirs);
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

void EnvironmentNAVXYTHETAV::ComputeHeuristicValues()
{
    //whatever necessary pre-computation of heuristic values is done here
    SBPL_PRINTF("Precomputing heuristics\n");

    SBPL_PRINTF("done\n");
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
	else
		InitGeneral( NULL);

	SBPL_PRINTF("size of env: %d by %d\n", EnvNAVXYTHETAVCfg.EnvWidth_c, EnvNAVXYTHETAVCfg.EnvHeight_c);

	return true;
}

bool EnvironmentNAVXYTHETAV::InitializeEnv(const char* sEnvFile){
	return false;
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
    if(FromStateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size() ||
       ToStateID >= (int)EnvNAVXYTHETAV.StateID2DataTable.size())
    {
        SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
        throw new SBPL_Exception();
    }
#endif

    //define this function if it is used in the planner

    SBPL_ERROR("ERROR in EnvNAVXYTHETAV.. function: FromToHeuristic is undefined\n");
    throw new SBPL_Exception();

    return 0;
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

    //define this function if it used in the planner (heuristic forward search would use it)

    SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: GetGoalHeuristic is undefined\n");
    throw new SBPL_Exception();
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

    //define this function if it used in the planner (heuristic backward search would use it)

    SBPL_ERROR("ERROR in EnvNAVXYTHETAV.. function: GetStartHeuristic is undefined\n");
    throw new SBPL_Exception();

    return 0;
}

/*
 * THIS FUNCTION MUST BE REVIEWED BECAUSE I HAVE DIFFERENT ACTIONS
 */
void EnvironmentNAVXYTHETAV::SetAllActionsandAllOutcomes(CMDPSTATE* state)
{
/*
#if DEBUG
    if (state->StateID >= (int)EnvNAVXYTHETAV.StateID2CoordTable.size()) {
        SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal\n");
        throw new SBPL_Exception();
    }

    if ((int)state->Actions.size() != 0) {
        SBPL_ERROR("ERROR in Env_setAllActionsandAllOutcomes: actions already exist for the state\n");
        throw new SBPL_Exception();
    }
#endif

    //if it is goal then no successors
    if (state->StateID == EnvNAVXYTHETAV.goalstateid) return;

    //get values for the state
    EnvXXXHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[state->StateID];

    //iterate through the actions for the state
    for (int aind = 0; aind < NAVXYTHETAV_MAXACTIONSWIDTH; aind++) {
        int cost = 1;

        //Add Action
        CMDPACTION* action = state->AddAction(aind);

        //clock_t currenttime = clock();
        //add all the outcomes to the action
		AddAllOutcomes(HashEntry->X1, HashEntry->X2, HashEntry->X3, HashEntry->X4, action, cost);

        //you can break if the number of actual actions is smaller than the maximum possible

        //time3_addallout += clock()-currenttime;
    }*/
}

void EnvironmentNAVXYTHETAV::SetAllPreds(CMDPSTATE* state)
{
    //implement this if the planner needs access to predecessors

    SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: SetAllPreds is undefined\n");
    throw new SBPL_Exception();
}

void EnvironmentNAVXYTHETAV::GetSuccs(int SourceStateID, vector<int>* SuccIDV, vector<int>* CostV)
{
    SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: GetSuccs is undefined\n");
    throw new SBPL_Exception();
}

void EnvironmentNAVXYTHETAV::GetPreds(int TargetStateID, vector<int>* PredIDV, vector<int>* CostV)
{
    SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: GetPreds is undefined\n");
    throw new SBPL_Exception();
}

int EnvironmentNAVXYTHETAV::SizeofCreatedEnv()
{
    return (int)EnvNAVXYTHETAV.StateID2DataTable.size();
}

void EnvironmentNAVXYTHETAV::PrintState(int stateID, bool bVerbose, FILE* fOut /*=NULL*/)
{
#if DEBUG
    if(stateID >= (int)EnvNAVXYTHETAV.StateID2CoordTable.size())
    {
        SBPL_ERROR("ERROR in EnvNAVXYTHETAV... function: stateID illegal (2)\n");
        throw new SBPL_Exception();
    }
#endif

    if (fOut == NULL) fOut = stdout;

    EnvNAVXYTHETAVHashEntry_t* HashEntry = EnvNAVXYTHETAV.StateID2DataTable[stateID];

    if (stateID == EnvNAVXYTHETAV.goalstateid) {
        SBPL_FPRINTF(fOut, "the state is a goal state\n");
    }

    SBPL_FPRINTF(fOut, "X1=%d X2=%d X3=%d X4=%d\n", HashEntry->x, HashEntry->y, HashEntry->theta, HashEntry->v);
}

void EnvironmentNAVXYTHETAV::PrintEnv_Config(FILE* fOut)
{
    //implement this if the planner needs to print out EnvXXX. configuration

    SBPL_ERROR("ERROR in EnvXXX... function: PrintEnv_Config is undefined\n");
    throw new SBPL_Exception();
}

//------------------------------------------------------------------------------
