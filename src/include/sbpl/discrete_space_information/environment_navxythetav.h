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

#define NAVXYTHETAV_DEFAULTOBSTHRESHOLD 254 /**< Default obstacle threshold if none is specified */
#define NAVXYTHETAV_DEFAULTTHETADIRS 16 /**< Default orientation directions if none are specified */
#define NAVXYTHETAV_DEFAULTMEDIUMVELOCITY 9.0 /**< Despite the name represents the maximum velocity admissible by a vehicle used */

#define NAVXYTHETAV_COSTMULT_MTOMM 1000 /**< Conversion factor between meters and millimeters */

#define NAVXYTHETAV_PERFORMANCE_TEST GLOBAL_PERFORMANCE_TEST /**< Flag set up if test execution must be done */

//class CMDPACTION;

class CMDPSTATE;
class MDPConfig;
class SBPL2DGridSearch;

/**
 * \brief Data structure to represent an index of an action
 */
typedef struct{
	int vector_index; /**< Subvector index of the action */
	int aind; /**< Index of the action */
} EnvNAVXYTHETAVActionIndex_t;

/**
 * \brief Data structure to represent a lattice primitive
 */
typedef struct
{
	int motprimID; /**< ID of the lattice primitive */
	int start_theta_disc; /**< Start discretized orientation value */
	int start_v_disc; /**< Start discretized velocity value */
	double additionalactioncostmult; /**< Multiplicative cost factor of the primitive */
	sbpl_xy_theta_v_cell_t endcell; /**< Final discretized cell */
	std::vector<sbpl_xy_theta_v_pt_t> intermptV; /**< Intermediate pose */
} SBPL_xythetav_mprimitive;

/**
 * \brief Data structure to represent an action
 */
typedef struct
{
	unsigned char aind; /**< Index of the action */
	int starttheta; /**< Start discretized orientation of the action */
	int startv; /**< Start discretized velocity of the action */
	int dX; /**< Displacement in x direction (in cells) */
	int dY; /**< Displacement in y direction (in cells) */
	int endtheta; /**< Final discretized orientation value */
	int endv; /**< Final discretized velocity value */
	int cost; /**< Cost of the action */
	std::vector<sbpl_2Dcell_t> intersectingcellsV; /**< Cells instersected by the action considering only x and y from the origin */
	std::vector<sbpl_xy_theta_v_pt_t> intermptV; /**< Intermediate poses of the action */
	std::vector<sbpl_xy_theta_v_cell_t> interm3DcellsV; /**< Intermediate states of the action */
} EnvNAVXYTHETAVAction_t;

/**
 * \brief Group of action given a pair of start orientation and start velocity
 */
class sbpl_xythetav_groupaction_t{
private:
	int start_theta;
	int start_v;
	std::vector<EnvNAVXYTHETAVAction_t> actions;
public:
	sbpl_xythetav_groupaction_t(){
	}
};

/**
 * \brief Data structure to store the environment configuration
 */
typedef struct ENV_NAVXYTHETAV_CONFIG
{
	//parameters that are read from the configuration file
	int EnvWidth_c; /**< Width of the environment (in cells) */
	int EnvHeight_c; /**< Height of the environment (in cells) */
	int NumThetaDirs; /**< Number of orientation values */
	int numV; /**< Number of velocities admissed */
	int StartX_c; /**< Start discretized position x */
	int StartY_c; /**< Start discretized position y */
	int StartTheta; /**< Start discretized orientation theta */
	int StartV; /**< Start discretized velocity value */
	int EndX_c; /**< End discretized position x */
	int EndY_c; /**< End discretized position y */
	int EndTheta; /**< End discretized orientation */
	int EndV; /**< End discretized velocity */
	double ** Grid2D; /**< Cost map */
	
	/*
	* the value at which and above which cells are obstacles in the maps sent from outside
	* the default is defined above
	*/
	unsigned char obsthresh; /**< Obstacle threshold */

	/* 
	* the value at which and above which until obsthresh (not including it)
	* cells have the nearest obstacle at distance smaller than or equal to
	* the inner circle of the robot. In other words, the robot is definitely
	* colliding with the obstacle, independently of its orientation
	* if no such cost is known, then it should be set to obsthresh (if center
	* of the robot collides with obstacle, then the whole robot collides with
	* it independently of its rotation)
	*/
	unsigned char cost_inscribed_thresh; /**< Cost inscribed threshold */

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
	int cost_possibly_circumscribed_thresh; /**< Cost possibly circumscribed threshold */ // it has to be integer, because -1 means that it is not provided.
	
	double cellsize_m; /**< Size of a cell in meters */
	std::vector<double> velocities; /**< Vector containing admissible velocities */

	/* Maybe below parameter is not really useful */
	//int dXY[NAVXYTHETALAT_DXYWIDTH][2];
	
	/* FOR ACTIONS MUST BE DEFINED MORE COMPLEX DATA STRUCTURES */
	/* SEE IF USEFUL DATA STRUCTURE CREATED OR IS ENOUGH A VECTOR OF VECTORS */
	/* Actions... wait definition before use */
	//array of actions, ActionsV[i][j] - jth action for sourcetheta+sourcev*NUMTHETA = i
	//EnvNAVXYTHETAVAction_t** ActionsV;
	std::vector<std::vector<EnvNAVXYTHETAVAction_t>* > *ActionsV; /**< All the actions admissible divided by finish index */
	
	/* SEE IF NEEDED A PREDECESSOR ARRAY */
	//PredActionsV[i] - vector of pointers to the actions that result in a state with theta+v*NUMTHETA = i
	std::vector<EnvNAVXYTHETAVActionIndex_t> *PredActionsV; /**< Actions index divided by start state */

	//int actionwidth; //number of motion primitives
	std::vector<SBPL_xythetav_mprimitive> mprimV; /**< All the primitives provided */

	std::vector<sbpl_2Dpt_t> FootprintPolygon; /**< Robot footprint */
} EnvNAVXYTHETAVConfig_t;

/**
 * \brief Data structure to store a hashentry for the environment
 */
typedef struct ENVNAVXYTHETAVHASHENTRY
{
	int stateID; /**< ID of the state */
	int x; /**< x position in cells */
	int y; /**< y position in cells */
	int theta; /**< Discretized orientation value */
	int v; /**< Discretized velocity value */
	
	//See if below parameter is useful
	int iteration;
} EnvNAVXYTHETAVHashEntry_t;

/**
 * \brief Data structure to store important environment informations for planning
 */
typedef struct
{
	int startstateid; /**< ID of the start state */
	int goalstateid; /**< ID of the goal state */

	//hash table of size x_size*y_size (?*possible_v*possible_theta?). Maps from coords to stateId
	int HashTableSize; /**< Size of the hashtable */
	std::vector<EnvNAVXYTHETAVHashEntry_t*>* Data2StateIDHashTable; /**< Hash table that allow to map the state coordinates to an ID */

	//vector that maps from stateID to coords
	std::vector<EnvNAVXYTHETAVHashEntry_t*> StateID2DataTable; /**< Vector that allow to map an ID to a complete state */
	
	bool bInitialized; /**< Flag to know if the environment is already initilized or it must be initilized */
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
	*        reference point centered at x=0,y=0, orientation = 0 (along x axis) and speed = 0 m/s.
	*        The perimeter is defined in meters as a sequence of vertices of a
	*        polygon defining the perimeter. If vector is of zero size, then robot
	*        is assumed to be point robot (you may want to inflate all obstacles by
	*        its actual radius) Motion primitives file defines the motion primitives
	*        available to the robot
	* \param sEnvFile path to the environment configuration file
	* \param perimeterptsV footprint of the robot given as vector of points
	* \param sMotPrimFile path to the motion primitive file
	* \return true if initialization is ok, false otherwise
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
	* \param sourceStateID ID of the state for which the successors are wanted
	* \param[out] succIDV vector in which the IDs of the successor states are stored
	* \param[out] costV vector in which the cost to reach the successors are stored
	* \param[out] actionindV vector in which the actions index to reach te successors are stored; if parameter is NULL they are not stored
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
	* \param width environment width
	* \param heigh environment height
	* \param numthetadirs number of theta used for the discretization
	* \param numv number of velocities used for discretization
	* \param velocities vector containing the velocities used
	* \param mapdata pointer to memory area containing the cost map; the position [0] is the up left corner of the map
	* \param startx start position x of the vehicle in meters
	* \param starty start position y of the vehicle in meters
	* \param starttheta start orientation of the vehicle in rad
	* \param startv start velocity of the vehicle in m/s
	* \param goalx goal position x of the vehicle in meters
	* \param goaly goal position y of the vehicle in meters
	* \param goaltheta goal orientation of the vehicle in rad
	* \param goalv goal velocity of the vehicle in m/s
	* \param perimeterptsV footprint of the robot expressed as vector of points
	* \param cellsize_m size of a cell in meters
	* \param obsthresh obstacle threshold
	* \param cost_inscribed_thresh
	* \param cost_possibly_circumscribed_thresh
	* \param sMotPrimFile path to the motion primitive file
	* \return true if the initialization end well, false otherwise
	*/
	virtual bool InitializeEnv(int width, int height, int numthetadirs, int numv, std::vector<double> velocities,
							const double* mapdata,
							double startx, double starty, double starttheta, double startv,
							double goalx, double goaly, double goaltheta, double goalv,
							const std::vector<sbpl_2Dpt_t>& perimeterptsV, double cellsize_m,
							unsigned char obsthresh, unsigned char cost_inscribed_thresh,
							unsigned char cost_possibly_circumscribed_thresh, const char* sMotPrimFile);
	
	/**
	* \brief update the cost of a cell<x,y>
	* \param x x position of the cell
	* \param y y position of the cell
	* \param newcost updated cost of the cell
	* \return true if the update finish well, false otherwise
	*/
	virtual bool UpdateCost(int x, int y, double newcost);
	
	/**
	* \brief re-setting the whole 2D map. The function transforms from linear array mapdata to 
	* 			the 2D matrix used internally. The position 0 of the array is the up left corner.
	* \param mapdata pointer to memory area where the cost map is stored
	* \return true if the set finish well, false otherwise
	*/
	virtual bool SetMap(const double* mapdata);
	
	/**
	* \brief this function fill in Predecessor/Successor states of edges whose costs changed
	*        It takes in an array of cells whose traversability changed, and
	*        returns (in vector preds_of_changededgesIDV) the IDs of all
	*        states that have outgoing edges that go through the changed
	*        cells
	* \param changedcellsV pointer to the vector containing the cells that have changed the value
	* \param[out] preds_of_changededgesIDV pointer to a vector filled with the ID of the predecessors of the cells that are changed
	*/
	virtual void GetPredsOfChangedEdges(std::vector<sbpl_2Dcell_t> const * changedcellsV,
										std::vector<int> *preds_of_changededgesIDV);
	/**
	* \brief same as GetPredsofChangedEdges, but returns successor states.
	*        Both functions need to be present for incremental search
	* \param changedcellsV pointer to the vector containing the cells that have changed the value
	* \param[out] succs_of_changededgesIDV pointer to a vector filled with the ID of the successors of the cells that are changed
	*/
	virtual void GetSuccsOfChangedEdges(std::vector<sbpl_2Dcell_t> const * changedcellsV,
										std::vector<int> *succs_of_changededgesIDV);
	
	/**
	* \brief this function allow to know if a cell is an obstacle or not
	* \param x x position of the cell
	* \param y y position of the cell
	* \return true if the cell is an obstacle, false otherwise
	*/
	virtual bool IsObstacle(int x, int y);
	
	/**
	* \brief returns false if robot intersects obstacles or lies outside of
	*        the map or have an invalid velocity. Note this is pretty expensive
	* 		  operation since it computes the footprint of the robot based on its x,y,theta
	* \param x x position of the configuration in cells
	* \param y y position of the configuration in cells
	* \param theta discrete orientation of the vehicle in that configuration
	* \param v discrete velocity of the vehicle in that configuration
	* \return true if the configuration is valid, false otherwise
	*/
	virtual bool IsValidConfiguration(int x, int y, int theta, int v);
	
	/**
	* \brief returns environment parameters. Useful for creating a copy environment.
	* \param[out] size_x dimension x of the environment
	* \param[out] size_y dimension y of the environment
	* \param[out] numthetadirs number of orientation discretization values
	* \param[out] numv number of velocities discretization values
	* \param[out] velocities vector of velocities used
	* \param[out] startx start position x of the vehicle in meters
	* \param[out] starty start position y of the vehicle in meters
	* \param[out] starttheta start orientation of the vehicle in rad
	* \param[out] startv start velocity of the vehicle in m/s
	* \param[out] goalx goal position x of the vehicle in meters
	* \param[out] goaly goal position y of the vehicle in meters
	* \param[out] goaltheta goal orientation of the vehicle in rad
	* \param[out] goalv goal velocity of the vehicle in m/s
	* \param[out] cellsize_m size of a cell in meters
	* \param[out] obsthresh obstacle threshold
	* \param[out] cost_inscribed_thresh
	* \param[out] cost_possibly_circumscribed_thresh
	* \param[out] motionprimitiveV vector containing all the motion primitives
	*/
	virtual void GetEnvParams(int *size_x, int *size_y, int * numthetadirs, int * numv, std::vector<double> * velocities,
							double* startx, double* starty, double* starttheta, double * startv, 
							double* goalx, double* goaly, double* goaltheta, double * goalv, double* cellsize_m,
							unsigned char* cost_inscribed_thresh, unsigned char* cost_possibly_circumscribed_thresh,
							unsigned char* obsthresh, std::vector<SBPL_xythetav_mprimitive>* motionprimitiveV);
	
	/**
	* \brief get internal configuration data structure
	* \return the data structure
	*/
	virtual const EnvNAVXYTHETAVConfig_t* GetEnvNavConfig();
	
	/**
	* \brief prints time statistics
	* \param fOut pointer to the file where write the time statistics
	*/
	virtual void PrintTimeStat(FILE* fOut);
	
	/**
	* \brief returns the cost corresponding to the cell <x,y>
	* \param x x cell position
	* \param y y cell position
	* \return the cell cost
	*/
	virtual double GetMapCost(int x, int y);
	
	/**
	* \brief returns true if cell is within map
	* \param x x cell position
	* \param y y cell position
	* \return true if the cell is within the map, false otherwise
	*/
	virtual bool IsWithinMapCell(int x, int y);
	
	/**
	* \brief Transform a pose into discretized form.
	* \param px continuous position x
	* \param py continuous position y
	* \param pth continuous orientation theta
	* \param pv continuous velocity v
	* \param[out] ix discrete position x
	* \param[out] iy discrete position y
	* \param[out] ith discrete orientation theta
	* \param[out] iv discrete velocity v
	* \return true if the conversion is done correctly, false otherwise
	*/
	virtual bool PoseContToDisc(double px, double py, double pth, double pv, int &ix, int &iy, int &ith, int &iv) const;
	
	/** 
	* \brief Transform grid indices into a continuous pose.
	* \param ix discrete position x
	* \param iy discrete position y
	* \param ith discrete orientation theta
	* \param iv discrete velocity v
	* \param[out] px continuous position x
	* \param[out] py continuous position y
	* \param[out] pth continuous orientation theta
	* \param[out] pv continuous velocity v
	* \return true if the conversion is done correctly, false otherwise
	*/
	virtual bool PoseDiscToCont(int ix, int iy, int ith, int iv, double &px, double &py, double &pth, double &pv) const;
	
	/**
	* \brief sets start in meters, radians and m/s
	* \param x x position of the start state in meters
	* \param y y position of the start state in meters
	* \param theta orientation of the start state in rad
	* \param v velocity of the start state in m/s
	* \return the start state ID
	*/
	virtual int SetStart(double x, double y, double theta, double v);
	
	/**
	* \brief sets goal in meters, radians and m/s
	* \param x x position of the goal state in meters
	* \param y y position of the goal state in meters
	* \param theta orientation of the goal state in rad
	* \param v velocity of the goal state in m/s
	* \return the goal state ID
	*/
	virtual int SetGoal(double x, double y, double theta, double v);
	
	/**
	* \brief returns state data of state with ID=stateID
	* \param stateID ID of the state
	* \param[out] x x position of the state with ID stateID in meters
	* \param[out] y y position of the state with ID stateID in meters
	* \param[out] theta orientation of the state with ID stateID in rad
	* \param[out] v velocity of the state with ID stateID in m/s
	*/
	virtual void GetCoordFromState(int stateID, int& x, int& y, int& theta, int& v) const;
	
	/**
	* \brief returns stateID for a state with data x,y,theta,v
	* \param x x position of a state in meters
	* \param y y position of a state in meters
	* \param theta orientation of a state in rad
	* \param v velocity of a state in m/s
	* \return the state ID
	*/
	virtual int GetStateFromCoord(int x, int y, int theta, int v);
	
	/** \brief converts a path given by stateIDs into a sequence of coordinates.
	*         Note that since motion primitives are short actions
	*         represented as a sequence of points,
	*         the path returned by this function contains much more points than the
	*         number of points in the input path. The returned coordinates are in
	*         meters,meters,radians,m/s
	* \param stateIDPath pointer to the path expressed in state IDs
	* \param[out] xythetavPath path of continuous coordinates (x,y,theta,v)
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
	virtual void SetConfiguration(int width, int height, const double * mapdata, int startx, int starty,
								int starttheta, int startv, int goalx, int goaly, int goaltheta, int goalv,
							double cellsize_m, const std::vector<sbpl_2Dpt_t>& perimeterptsV);
};

#endif

