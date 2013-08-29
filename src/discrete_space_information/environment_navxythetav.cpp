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
}

void EnvironmentNAVXYTHETAV::InitializeEnvConfig()
{
    //aditional to configuration file initialization of EnvCfg if necessary
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

bool EnvironmentNAVXYTHETAV::InitializeEnv(const char* sEnvFile)
{
    FILE* fCfg = fopen(sEnvFile, "r");
    if (fCfg == NULL) {
        SBPL_ERROR("ERROR: unable to open %s\n", sEnvFile);
        throw new SBPL_Exception();
    }
    ReadConfiguration(fCfg);
    fclose(fCfg);

    //Initialize other parameters of the environment
    InitializeEnvConfig();

    //initialize Environment
    InitializeEnvironment();

    //pre-compute heuristics
    ComputeHeuristicValues();

    return true;
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
