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

#ifndef __ENVIRONMENT_NAVXYTHETAV_H_
#define __ENVIRONMENT_NAVXYTHETAV_H_

#include <cstdio>
#include <vector>
#include <sbpl/discrete_space_information/environment.h>

/* Vedere se serve questa define */
#define NAVXYTHETAV_MAXACTIONSWIDTH 9

#define NAVXYTHETAV_DEFAULTOBSTHRESHOLD 254

//class CMDPACTION;

class CMDPSTATE;
class MDPConfig;

/* 
 * Missing the defintion of action. In xytheta is used also an heuristic
 * that execute 2D search in environment. Now, i use euclidean distance
 * after i will change it.
 */

/*
 * The velocity will be discretized or not?
 */

typedef struct ENV_NAVXYTHETAV_CONFIG
{
    //parameters that are read from the configuration file
    int EnvWidth_c;
    int EnvHeight_c;
    int NumThetaDirs;
	int numV;
    int StartX_c;
    int StartY_c;
    int StartTheta;
	int StartV;
    int EndX_c;
    int EndY_c;
    int EndTheta;
	int EndV;
	unsigned char ** Grid2D;
	
	// the value at which and above which cells are obstacles in the maps sent from outside
    // the default is defined above
    unsigned char obsthresh;

    // the value at which and above which until obsthresh (not including it)
    // cells have the nearest obstacle at distance smaller than or equal to
    // the inner circle of the robot. In other words, the robot is definitely
    // colliding with the obstacle, independently of its orientation
    // if no such cost is known, then it should be set to obsthresh (if center
    // of the robot collides with obstacle, then the whole robot collides with
    // it independently of its rotation)
    unsigned char cost_inscribed_thresh;

    // the value at which and above which until cost_inscribed_thresh (not including it) cells
    // **may** have a nearest osbtacle within the distance that is in between
    // the robot inner circle and the robot outer circle
    // any cost below this value means that the robot will NOT collide with any
    // obstacle, independently of its orientation
    // if no such cost is known, then it should be set to 0 or -1 (then no cell
    // cost will be lower than it, and therefore the robot's footprint will
    // always be checked)
    int cost_possibly_circumscribed_thresh; // it has to be integer, because -1 means that it is not provided.
    
    /* Maybe below parameter for ackermann vehicle with this state have no sense */
	//double nominalvel_mpersecs;
    //double timetoturn45degsinplace_secs;
    
    double cellsize_m;

	/* Maybe below parameter is not really useful */
    //int dXY[NAVXYTHETALAT_DXYWIDTH][2];
    
    /* Actions... wait definition before use */
    //array of actions, ActionsV[i][j] - jth action for sourcetheta = i
    //EnvNAVXYTHETALATAction_t** ActionsV; 
    //PredActionsV[i] - vector of pointers to the actions that result in a state with theta = i
    //std::vector<EnvNAVXYTHETALATAction_t*>* PredActionsV; 

    //int actionwidth; //number of motion primitives
    //std::vector<SBPL_xytheta_mprimitive> mprimV;

    //std::vector<sbpl_2Dpt_t> FootprintPolygon;
} EnvNAVXYTHETAVConfig_t;

typedef struct ENVNAVXYTHETAVHASHENTRY
{
    int stateID;
    int x;
    int y;
    int theta;
    int v;
	
	//See if below parameter is useful
	int iteration;
} EnvNAVXYTHETAVHashEntry_t;

typedef struct
{
    int startstateid;
    int goalstateid;

    //hash table of size x_size*y_size (?*possible_v*possible_theta?). Maps from coords to stateId
    int HashTableSize;
    std::vector<EnvNAVXYTHETAVHashEntry_t*>* Data2StateIDHashTable;

    //vector that maps from stateID to coords
    std::vector<EnvNAVXYTHETAVHashEntry_t*> StateID2DataTable;
	
	bool bInitialized;
} EnvironmentNAVXYTHETAV_t;

/** \brief This is 4 dimensions environment for ackermann vehicle.
 *         The 4 dimensions are the position (x,y),
 * 			the heading of the vehicle theta
 * 			and the velocity of the vehicle v
 */
class EnvironmentNAVXYTHETAV : public DiscreteSpaceInformation
{
public:
    /**
     * \brief see comments on the same function in the parent class
     */
    virtual bool InitializeEnv(const char* sEnvFile);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual bool InitializeMDPCfg(MDPConfig *MDPCfg);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetFromToHeuristic(int FromStateID, int ToStateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetGoalHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetStartHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void SetAllActionsandAllOutcomes(CMDPSTATE* state);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void SetAllPreds(CMDPSTATE* state);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void GetSuccs(int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void GetPreds(int TargetStateID, std::vector<int>* PredIDV, std::vector<int>* CostV);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int SizeofCreatedEnv();

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void PrintState(int stateID, bool bVerbose, FILE* fOut = NULL);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void PrintEnv_Config(FILE* fOut);

    ~EnvironmentNAVXYTHETAV() { }

protected:
    //member variables
    EnvNAVXYTHETAVConfig_t EnvNAVXYTHETAVCfg;
    EnvironmentNAVXYTHETAV_t EnvNAVXYTHETAV;

    virtual void ReadConfiguration(FILE* fCfg);

    virtual void InitializeEnvConfig();

    virtual unsigned int GETHASHBIN(int x, int y, int theta, double v);

    virtual void PrintHashTableHist();

    virtual EnvNAVXYTHETAVHashEntry_t* GetHashEntry(int x, int y, int theta, int v);

    virtual EnvNAVXYTHETAVHashEntry_t* CreateNewHashEntry(int x, int y, int theta, int v);

    virtual void CreateStartandGoalStates();

    virtual void InitializeEnvironment();

	/* maybe not useful */
    /*virtual void AddAllOutcomes(int SourceX, int SourceY, int SourceTheta,
                                double SourceV, CMDPACTION* action, int cost);*/

    virtual void ComputeHeuristicValues();
};

#endif

