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
#include <sbpl/utils/utils.h>

#define NAVXYTHETAV_DEFAULTOBSTHRESHOLD 254
#define NAVXYTHETAV_DEFAULTTHETADIRS 16
#define NAVXYTHETAV_DEFAULTMEDIUMVELOCITY 8.0

#define NAVXYTHETAV_COSTMULT_MTOMM 1000

//class CMDPACTION;

class CMDPSTATE;
class MDPConfig;
class SBPL2DGridSearch;

/* 
* Missing the defintion of action. In xytheta is used also an heuristic
* that execute 2D search in environment. Now, i use euclidean distance
* after i will change it.
*/

typedef struct
{
	int motprimID;
	int start_theta_disc;
	int start_v_disc;
	double additionalactioncostmult;
	sbpl_xy_theta_v_cell_t endcell;
	std::vector<sbpl_xy_theta_v_pt_t> intermptV;
} SBPL_xythetav_mprimitive;

typedef struct
{
	unsigned char aind; //index of the action
	int starttheta;
	int startv;
	int dX;
	int dY;
	int endtheta;
	int endv;
	int cost;
	std::vector<sbpl_2Dcell_t> intersectingcellsV;
	std::vector<sbpl_xy_theta_v_pt_t> intermptV;
	std::vector<sbpl_xy_theta_v_cell_t> interm3DcellsV;
} EnvNAVXYTHETAVAction_t;

class sbpl_xythetav_groupaction_t{
private:
	int start_theta;
	int start_v;
	std::vector<EnvNAVXYTHETAVAction_t> actions;
public:
	sbpl_xythetav_groupaction_t(){
	}
};

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
	
	/*
	* the value at which and above which cells are obstacles in the maps sent from outside
	* the default is defined above
	*/
	unsigned char obsthresh;

	/* 
	* the value at which and above which until obsthresh (not including it)
	* cells have the nearest obstacle at distance smaller than or equal to
	* the inner circle of the robot. In other words, the robot is definitely
	* colliding with the obstacle, independently of its orientation
	* if no such cost is known, then it should be set to obsthresh (if center
	* of the robot collides with obstacle, then the whole robot collides with
	* it independently of its rotation)
	*/
	unsigned char cost_inscribed_thresh;

	/*
	* the value at which and above which until cost_inscribed_thresh (not including it) cells
	* **may** have a nearest osbtacle within the distance that is in between
	* the robot inner circle and the robot outer circle
	* any cost below this value means that the robot will NOT collide with any
	* obstacle, independently of its orientation
	* if no such cost is known, then it should be set to 0 or -1 (then no cell
	* cost will be lower than it, and therefore the robot's footprint will
	* always be checked)
	*/
	int cost_possibly_circumscribed_thresh; // it has to be integer, because -1 means that it is not provided.
	
	double cellsize_m;
	std::vector<double> velocities;

	/* Maybe below parameter is not really useful */
	//int dXY[NAVXYTHETALAT_DXYWIDTH][2];
	
	/* FOR ACTIONS MUST BE DEFINED MORE COMPLEX DATA STRUCTURES */
	/* SEE IF USEFUL DATA STRUCTURE CREATED OR IS ENOUGH A VECTOR OF VECTORS */
	/* Actions... wait definition before use */
	//array of actions, ActionsV[i][j] - jth action for sourcetheta+sourcev*NUMTHETA = i
	//EnvNAVXYTHETAVAction_t** ActionsV;
	std::vector<std::vector<EnvNAVXYTHETAVAction_t>* > *ActionsV;
	
	/* SEE IF NEEDED A PREDECESSOR ARRAY */
	//PredActionsV[i] - vector of pointers to the actions that result in a state with theta+v*NUMTHETA = i
	std::vector<EnvNAVXYTHETAVAction_t*>* PredActionsV;

	//int actionwidth; //number of motion primitives
	std::vector<SBPL_xythetav_mprimitive> mprimV;

	std::vector<sbpl_2Dpt_t> FootprintPolygon;
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
	EnvironmentNAVXYTHETAV();
	
	/**
	* \brief initialization of environment from file. See .cfg files for
	*        examples it also takes the perimeter of the robot with respect to some
	*        reference point centered at x=0,y=0, orientation = 0 (along x axis) and speed = 0.
	*        The perimeter is defined in meters as a sequence of vertices of a
	*        polygon defining the perimeter. If vector is of zero size, then robot
	*        is assumed to be point robot (you may want to inflate all obstacles by
	*        its actual radius) Motion primitives file defines the motion primitives
	*        available to the robot
	*/
	virtual bool InitializeEnv(const char* sEnvFile, const std::vector<sbpl_2Dpt_t>& perimeterptsV,
							const char* sMotPrimFile);
	
	/**
	* \brief see comments on the same function in the parent class
	*/
	virtual bool InitializeEnv(const char* sEnvFile);
	
	/**
	* \brief way to set up various parameters - see function body for the list of parameters
	*/
	virtual bool SetEnvParameter(const char* parameter, int value);
	
	/**
	* \brief returns the value of specific parameter - see function body for the list of parameters
	*/
	virtual int GetEnvParameter(const char* parameter);

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
	* \brief returns all successors states, costs of corresponding actions
	*        and pointers to corresponding actions, each of which is a motion
	*        primitive
	*        if actionindV is NULL, then pointers to actions are not returned
	*/
	virtual void GetSuccs(int sourceStateID, std::vector<int>* succIDV, std::vector<int>* costV,
						std::vector<EnvNAVXYTHETAVAction_t*>* actionindV = NULL);

	/**
	* \brief see comments on the same function in the parent class
	*/
	virtual void GetPreds(int TargetStateID, std::vector<int>* PredIDV, std::vector<int>* CostV);
	
	/**
	* \brief see comments on the same function in the parent class
	*/
	virtual void EnsureHeuristicsUpdated(bool bGoalHeuristics);

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
	
	/**
	* \brief initialize environment by parameters.
	*/
	virtual bool InitializeEnv(int width, int height, int numthetadirs, int numv, std::vector<double> velocities,
							const unsigned char* mapdata,
							double startx, double starty, double starttheta, double startv,
							double goalx, double goaly, double goaltheta, double goalv,
							const std::vector<sbpl_2Dpt_t>& perimeterptsV, double cellsize_m,
							unsigned char obsthresh, unsigned char cost_inscribed_thresh,
							unsigned char cost_possibly_circumscribed_thresh, const char* sMotPrimFile);
	
	/**
	* \brief update the traversability of a cell<x,y>
	*/
	virtual bool UpdateCost(int x, int y, unsigned char newcost);
	
	/**
	* \brief re-setting the whole 2D map. 
	* 			The function transforms from linear array mapdata to the 2D matrix used internally: Grid2D[x][y] = mapdata[x+y*width]
	*/
	virtual bool SetMap(const unsigned char* mapdata);
	
	/**
	* \brief this function fill in Predecessor/Successor states of edges whose costs changed
	*        It takes in an array of cells whose traversability changed, and
	*        returns (in vector preds_of_changededgesIDV) the IDs of all
	*        states that have outgoing edges that go through the changed
	*        cells
	*/
	virtual void GetPredsOfChangedEdges(std::vector<sbpl_2Dcell_t> const * changedcellsV,
										std::vector<int> *preds_of_changededgesIDV);
	/**
	* \brief same as GetPredsofChangedEdges, but returns successor states.
	*        Both functions need to be present for incremental search
	*/
	virtual void GetSuccsOfChangedEdges(std::vector<sbpl_2Dcell_t> const * changedcellsV,
										std::vector<int> *succs_of_changededgesIDV);
	
	/**
	* returns true if cell is untraversable
	*/
	virtual bool IsObstacle(int x, int y);
	
	/**
	* \brief returns false if robot intersects obstacles or lies outside of
	*        the map or have an invalid velocity. Note this is pretty expensive
	* 		  operation since it computes the footprint of the robot based on its x,y,theta
	*/
	virtual bool IsValidConfiguration(int x, int y, int theta, int v);
	
	/**
	* \brief returns environment parameters. Useful for creating a copy environment.
	*/
	virtual void GetEnvParams(int *size_x, int *size_y, int * numthetadirs, int * numv, std::vector<double> * velocities,
							double* startx, double* starty, double* starttheta, double * startv, 
							double* goalx, double* goaly, double* goaltheta, double * goalv, double* cellsize_m,
							unsigned char* cost_inscribed_thresh, unsigned char* cost_possibly_circumscribed_thresh,
							unsigned char* obsthresh, std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV);
	
	/**
	* \brief get internal configuration data structure
	*/
	virtual const EnvNAVXYTHETAVConfig_t* GetEnvNavConfig();
	
	/**
	* \brief prints time statistics
	*/
	virtual void PrintTimeStat(FILE* fOut);
	
	/**
	* \brief returns the cost corresponding to the cell <x,y>
	*/
	virtual unsigned char GetMapCost(int x, int y);
	
	/**
	* \brief returns true if cell is within map
	*/
	virtual bool IsWithinMapCell(int x, int y);
	
	/**
	* \brief Transform a pose into discretized form.
	*/
	virtual bool PoseContToDisc(double px, double py, double pth, double pv, int &ix, int &iy, int &ith, int &iv) const;
	
	/** 
	* \brief Transform grid indices into a continuous pose.
	*/
	virtual bool PoseDiscToCont(int ix, int iy, int ith, int iv, double &px, double &py, double &pth, double &pv) const;
	
	/**
	* \brief sets start in meters, radians and m/s
	*/
	virtual int SetStart(double x, double y, double theta, double v);
	
	/**
	* \brief sets goal in meters, radians and m/s
	*/
	virtual int SetGoal(double x, double y, double theta, double v);
	
	/**
	* \brief returns state data of state with ID=stateID
	*/
	virtual void GetCoordFromState(int stateID, int& x, int& y, int& theta, int& v) const;
	
	/**
	* \brief returns stateID for a state with data x,y,theta,v
	*/
	virtual int GetStateFromCoord(int x, int y, int theta, int v);
	
	/** \brief converts a path given by stateIDs into a sequence of coordinates.
	*         Note that since motion primitives are short actions
	*         represented as a sequence of points,
	*         the path returned by this function contains much more points than the
	*         number of points in the input path. The returned coordinates are in
	*         meters,meters,radians,m/s
	*/
	virtual void ConvertStateIDPathintoXYThetaVPath(std::vector<int>* stateIDPath, std::vector<sbpl_xy_theta_v_pt_t>* xythetavPath);

	~EnvironmentNAVXYTHETAV();

protected:
	//member variables
	EnvNAVXYTHETAVConfig_t EnvNAVXYTHETAVCfg;
	EnvironmentNAVXYTHETAV_t EnvNAVXYTHETAV;
	std::vector<sbpl_xy_theta_v_cell_t> affectedsuccstatesV; //arrays of states whose outgoing actions cross cell 0,0
	std::vector<sbpl_xy_theta_v_cell_t> affectedpredstatesV; //arrays of states whose incoming actions cross cell 0,0
	
	//2D search for heuristic computations
	bool bNeedtoRecomputeStartHeuristics; //set whenever grid2Dsearchfromstart needs to be re-executed
	bool bNeedtoRecomputeGoalHeuristics; //set whenever grid2Dsearchfromgoal needs to be re-executed
	SBPL2DGridSearch* grid2Dsearchfromstart; //computes h-values that estimate distances from start x,y to all cells
	SBPL2DGridSearch* grid2Dsearchfromgoal; //computes h-values that estimate distances to goal x,y from all cells
	
	int iteration;

	virtual void ReadConfiguration(FILE* fCfg);

	virtual unsigned int GETHASHBIN(unsigned int x, unsigned int y, unsigned int theta, unsigned int v);

	virtual void PrintHashTableHist();

	virtual EnvNAVXYTHETAVHashEntry_t* GetHashEntry(unsigned int x, unsigned int y, unsigned int theta, unsigned int v);

	virtual EnvNAVXYTHETAVHashEntry_t* CreateNewHashEntry(unsigned int x, unsigned int y, unsigned int theta, unsigned int v);

	virtual void CreateStartandGoalStates();
	
	virtual void ComputeReplanningDataforAction(EnvNAVXYTHETAVAction_t* action);
	virtual void ComputeReplanningData();
	virtual void PrecomputeActionswithCompleteMotionPrimitive(std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV);
	virtual void InitializeEnvConfig(std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV);
	virtual bool InitGeneral(std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV);

	virtual void InitializeEnvironment();
	
	virtual bool ReadMotionPrimitives(FILE* fMotPrims);
	virtual bool ReadSingleMotionPrimitive(SBPL_xythetav_mprimitive* pMotPrim, FILE* fIn);
	virtual bool ReadSingleCell(sbpl_xy_theta_v_cell_t* cell, FILE* fIn);
	virtual bool ReadSinglePose(sbpl_xy_theta_v_pt_t* pose, FILE* fIn);

	/* maybe not useful */
	/*virtual void AddAllOutcomes(int SourceX, int SourceY, int SourceTheta,
								double SourceV, CMDPACTION* action, int cost);*/

	virtual void ComputeHeuristicValues();
	virtual void PrintHeuristicValues();
	
	virtual double EuclideanDistance_m(int x1, int y1, int x2, int y2);
	virtual bool IsValidCell(int x, int y);
	virtual int GetActionCost(int sourceX, int sourceY, int sourceTheta, int sourceV, EnvNAVXYTHETAVAction_t* action);
	virtual void SetConfiguration(int width, int height, const unsigned char * mapdata, int startx, int starty,
								int starttheta, int startv, int goalx, int goaly, int goaltheta, int goalv,
							double cellsize_m, const std::vector<sbpl_2Dpt_t>& perimeterptsV);
};

#endif

