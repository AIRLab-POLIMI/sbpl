/*
* Copyright (c) 2008, Maxim Likhachev
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the University of Pennsylvania nor the names of its
*       contributors may be used to endorse or promote products derived from
*       this software without specific prior written permission.
*
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

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
//#include <sstream>

using namespace std;

#include <sbpl/headers.h>

enum PlannerType
{
	INVALID_PLANNER_TYPE = -1,
	PLANNER_TYPE_ADSTAR,
	PLANNER_TYPE_ARASTAR,
	PLANNER_TYPE_PPCP,
	PLANNER_TYPE_RSTAR,
	PLANNER_TYPE_VI,
	PLANNER_TYPE_ANASTAR,
	PLANNER_TYPE_ANASTARDOUBLE,

	NUM_PLANNER_TYPES
};

std::string PlannerTypeToStr(PlannerType plannerType)
{
	switch (plannerType) {
	case PLANNER_TYPE_ADSTAR:
		return std::string("adstar");
	case PLANNER_TYPE_ARASTAR:
		return std::string("arastar");
	case PLANNER_TYPE_PPCP:
		return std::string("ppcp");
	case PLANNER_TYPE_RSTAR:
		return std::string("rstar");
	case PLANNER_TYPE_VI:
		return std::string("vi");
	case PLANNER_TYPE_ANASTAR:
		return std::string("anastar");
	case PLANNER_TYPE_ANASTARDOUBLE:
		return std::string("anastardouble");
	default:
		return std::string("invalid");
	}
}

PlannerType StrToPlannerType(const char* str)
{
	if (!strcmp(str, "adstar")) {
		return PLANNER_TYPE_ADSTAR;
	}
	else if (!strcmp(str, "arastar")) {
		return PLANNER_TYPE_ARASTAR;
	}
	else if (!strcmp(str, "ppcp")) {
		return PLANNER_TYPE_PPCP;
	}
	else if (!strcmp(str, "rstar")) {
		return PLANNER_TYPE_RSTAR;
	}
	else if (!strcmp(str, "vi")) {
		return PLANNER_TYPE_VI;
	}
	else if (!strcmp(str, "anastar")) {
		return PLANNER_TYPE_ANASTAR;
	}
	else if(!strcmp(str, "anastardouble")) {
		return PLANNER_TYPE_ANASTARDOUBLE;
	}
	else {
		return INVALID_PLANNER_TYPE;
	}
}

enum MainResultType
{
	INVALID_MAIN_RESULT = -1,

	MAIN_RESULT_SUCCESS = 0,
	MAIN_RESULT_FAILURE = 1,
	MAIN_RESULT_INSUFFICIENT_ARGS = 2,
	MAIN_RESULT_INCORRECT_OPTIONS = 3,
	MAIN_RESULT_UNSUPPORTED_ENV = 4,

	NUM_MAIN_RESULTS
};

/*******************************************************************************
* PrintUsage - Prints the proper usage of the sbpl test executable.
*
* @param argv The command-line arguments; used to determine the name of the
*             test executable.
*******************************************************************************/
void PrintUsage(char *argv[])
{
	printf("USAGE: %s [-s] --planner=<planner_t> --search-dir=<search_t> <cfg file> <mot prims>\n",
		argv[0]);
	printf("See '%s -h' for help.\n", argv[0]);
}

/*******************************************************************************
* PrintHelp - Prints a help prompt to the command line when the -h option is
*             used.
*
* @param argv The command line arguments; used to determine the name of the
*             test executable
*******************************************************************************/
void PrintHelp(char** argv)
{
	printf("\n");
	printf("Search-Based Planning Library\n");
	printf("\n");
	printf("    %s -h\n", argv[0]);
	printf("    %s [-s] [--env=<env_t>] [--planner=<planner_t>] [--search-dir=<search_t>] <env cfg> [mot prim]\n",
		argv[0]);
	printf("\n");
	printf("[-s]                      (optional) Find a solution for an example navigation\n");
	printf("                          scenario where the robot only identifies obstacles as\n");
	printf("                          it approaches them.\n");
	printf("--planner=<planner_t>     Select a planner to use for the example.\n");
	printf("                          The default is \"arastar\".\n");
	printf("<planner_t>               One of arastar, adstar, rstar, anastar.\n");
	printf("--search-dir=<search_t>   Select the type of search to run. The default\n");
	printf("                          is \"backwards\".\n");
	printf("<search_t>                One of backward, forward.\n");
	printf("<env cfg>                 Config file representing the environment configuration.\n");
	printf("                          See sbpl/env_examples/ for examples.\n");
	printf("<mot prim>                Motion primitives file for x,y,theta,v lattice\n");
	printf("                          planning.\n");
	printf("\n");
}

/*******************************************************************************
* CheckIsNavigating
* @brief Returns whether the -s option is being used.
*
* @param numOptions The number of options passed through the command line
* @param argv The command-line arguments
* @return whether the -s option was passed in on the cmd line
*******************************************************************************/
bool CheckIsNavigating(int numOptions, char** argv)
{
	for (int i = 1; i < numOptions + 1; i++) {
		if (strcmp(argv[i], "-s") == 0) {
			return true;
		}
	}
	return false;
}

/*******************************************************************************
* CheckSearchDirection -
* @brief Returns the search direction being used
*
* @param numOptions The number of options passed through the command line
* @param argv The command-line arguments
* @return A string representing the search direction; "backward" by default
******************************************************************************/
std::string CheckSearchDirection(int numOptions, char** argv)
{
	int optionLength = strlen("--search-dir=");
	for (int i = 1; i < numOptions + 1; i++) {
		if (strncmp("--search-dir=", argv[i], optionLength) == 0) {
			std::string s(&argv[i][optionLength]);
			return s;
		}
	}
	return std::string("backward");
}

/*******************************************************************************
* CheckPlannerType - Checks for a planner specifier passed in through the
*                    command line. This determines what planner to run in
*                    the example. If none is found, ARA* is assumed.
*
* @param numOptions The number of options passed through the command line
* @param argv The command-line arguments
* @return A string denoting the planner type; "arastar" by default
******************************************************************************/
std::string CheckPlannerType(int numOptions, char** argv)
{
	int optionLength = strlen("--planner=");
	for (int i = 1; i < numOptions + 1; i++) {
		if (strncmp("--planner=", argv[i], optionLength) == 0) {
			std::string s(&argv[i][optionLength]);
			return s;
		}
	}
	return std::string("arastar");
}

/*******************************************************************************
* planxythetav
* @brief An example of planning in four-dimensional space (x,y,theta,v)
*
* @param plannerType The type of planner to be used in this example
* @param envCfgFilename The environment config file. 
* @param motPrimFilename The motion primitives file.
* 
* @return 1 if the planner successfully found a solution; 0 otherwise
*******************************************************************************/
int planxythetavsteer(PlannerType plannerType, char* envCfgFilename, char* motPrimFilename, bool forwardSearch)
{
	int bRet = 0;
	double allocated_time_secs = 1200.0; // in seconds
	double initialEpsilon = 3.0;
	MDPConfig MDPCfg;
	bool bsearchuntilfirstsolution = false;
	bool bforwardsearch = forwardSearch;

	// set the perimeter of the robot (it is given with 0,0,0 robot ref. point for which planning is done)
	vector<sbpl_2Dpt_t> perimeterptsV;
	sbpl_2Dpt_t pt_m;
	double halfwidth = 1.0; //0.3;
	double halflength = 1.0; //0.45;
	pt_m.x = -halflength;
	pt_m.y = -halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = halflength;
	pt_m.y = -halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = halflength;
	pt_m.y = halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = -halflength;
	pt_m.y = halfwidth;
	perimeterptsV.push_back(pt_m);

	// clear the footprint
	//perimeterptsV.clear();

	// Initialize Environment (should be called before initializing anything else)
	EnvironmentNAVXYTHETAVSTEER environment_navxythetavsteer;

	if (!environment_navxythetavsteer.InitializeEnv(envCfgFilename, perimeterptsV, motPrimFilename)) {
		printf("ERROR: InitializeEnv failed\n");
		throw new SBPL_Exception();
	}

	// Initialize MDP Info
	if (!environment_navxythetavsteer.InitializeMDPCfg(&MDPCfg)) {
		printf("ERROR: InitializeMDPCfg failed\n");
		throw new SBPL_Exception();
	}

	// plan a path
	vector<int> solution_stateIDs_V;

	SBPLPlanner* planner = NULL;
	switch (plannerType) {
	case PLANNER_TYPE_ARASTAR:
		printf("Initializing ARAPlanner...\n");
		planner = new ARAPlanner(&environment_navxythetavsteer, bforwardsearch);
		planner->set_initialsolution_eps(initialEpsilon);
		break;
	case PLANNER_TYPE_ADSTAR:
		printf("Initializing ADPlanner...\n");
		planner = new ADPlanner(&environment_navxythetavsteer, bforwardsearch);
		planner->set_initialsolution_eps(initialEpsilon);
		break;
	case PLANNER_TYPE_RSTAR:
		//VIEW WHY
		printf("Invalid configuration: xythetav environment does not support rstar planner...\n");
		//return 0;
		planner = new RSTARPlanner(&environment_navxythetavsteer, bforwardsearch);
		planner->set_initialsolution_eps(initialEpsilon);
		break;
	case PLANNER_TYPE_ANASTAR:
		printf("Initializing anaPlanner...\n");
		planner = new anaPlanner(&environment_navxythetavsteer, bforwardsearch);
		break;
	case PLANNER_TYPE_ANASTARDOUBLE:
		printf("Initializing anaPlanner with key as double type...\n");
		planner = new anaPlannerDouble(&environment_navxythetavsteer, bforwardsearch);
		break;
	default:
		printf("Invalid planner type\n");
		break;
	}

	// set planner properties
	if (planner->set_start(MDPCfg.startstateid) == 0) {
		printf("ERROR: failed to set start state\n");
		throw new SBPL_Exception();
	}
	if (planner->set_goal(MDPCfg.goalstateid) == 0) {
		printf("ERROR: failed to set goal state\n");
		throw new SBPL_Exception();
	}
	
	planner->set_search_mode(bsearchuntilfirstsolution);

	// plan
	printf("start planning...\n");
	bRet = planner->replan(allocated_time_secs, &solution_stateIDs_V);
	printf("done planning\n");
	printf("size of solution=%d\n", (unsigned int)solution_stateIDs_V.size());

	environment_navxythetavsteer.PrintTimeStat(stdout);

	// write solution to sol.txt
	/*for(int contSols=0;contSols<((anaPlannerDouble*)planner)->sols.size();contSols++){
		stringstream nomeSC;
		nomeSC << "solContinuous" << contSols << ".txt";
		FILE* fSolC = fopen(nomeSC.str().c_str(), "w");
		if (fSolC == NULL) {
			printf("ERROR: could not open solution file\n");
			throw new SBPL_Exception();
		}
		// write the continuous solution to file
		vector<sbpl_xy_theta_v_pt_t> xythetavPath;
		environment_navxythetav.ConvertStateIDPathintoXYThetaVPath(&((anaPlannerDouble*)planner)->sols.at(contSols), &xythetavPath);
		printf("solution size=%d\n", (unsigned int)xythetavPath.size());
		for (unsigned int i = 0; i < xythetavPath.size(); i++) {
			fprintf(fSolC, "%.3f %.3f %.3f %.3f\n", xythetavPath.at(i).x, xythetavPath.at(i).y, xythetavPath.at(i).theta, xythetavPath.at(i).v);
		}
		fclose(fSolC);
	}*/
	
	// write solution to sol.txt
	const char* solC = "solContinuous.txt";
	FILE* fSolC = fopen(solC, "w");
	if (fSolC == NULL) {
		printf("ERROR: could not open solution file\n");
		throw new SBPL_Exception();
	}
	// write the continuous solution to file
	vector<sbpl_xy_theta_v_steer_pt_t> xythetavsteerPath;
	environment_navxythetavsteer.ConvertStateIDPathintoXYThetaVPath(&solution_stateIDs_V, &xythetavsteerPath);
	printf("solution size=%d\n", (unsigned int)xythetavsteerPath.size());
	for (unsigned int i = 0; i < xythetavsteerPath.size(); i++) {
		fprintf(fSolC, "%.3f %.3f %.3f %.3f %.3f\n", xythetavsteerPath.at(i).x, xythetavsteerPath.at(i).y, xythetavsteerPath.at(i).theta, xythetavsteerPath.at(i).v, xythetavsteerPath.at(i).steer);
	}
	fclose(fSolC);

	const char* solD = "solDiscrete.txt";
	FILE* fSolD = fopen(solD, "w");
	if (fSolD == NULL) {
		printf("ERROR: could not open solution file\n");
		throw new SBPL_Exception();
	}
	// write the discrete solution to file
	for (size_t i = 0; i < solution_stateIDs_V.size(); i++) {
		int x;
		int y;
		int theta;
		int v;
		int steer;
		environment_navxythetavsteer.GetCoordFromState(solution_stateIDs_V[i], x, y, theta, v, steer);
	
		fprintf(fSolD, "%d %d %d %d %d\t\t%.3f %.3f %.3f %.3f %.3f\n", x, y, theta, v, steer,
		DISCXY2CONT(x, environment_navxythetavsteer.GetEnvNavConfig()->cellsize_m), DISCXY2CONT(y, environment_navxythetavsteer.GetEnvNavConfig()->cellsize_m), DiscTheta2Cont(theta, environment_navxythetavsteer.GetEnvNavConfig()->NumThetaDirs), DiscV2Cont(v, environment_navxythetavsteer.GetEnvNavConfig()->velocities), DiscSteer2Cont(steer, environment_navxythetavsteer.GetEnvNavConfig()->numSteers));
	}
	fclose(fSolD);
	
	environment_navxythetavsteer.PrintTimeStat(stdout);

	if (bRet) {
		// print the solution
		printf("Solution is found\n");
	}
	else {
		printf("Solution does not exist\n");
	}

	fflush(NULL);

	delete planner;

	return bRet;
}

/*******************************************************************************
* planandnavigatexythetav
* @brief An example simulation of how a robot would use (x,y,theta,v) lattice
*        planning.
*
* @param envCfgFilename The environment config file.
*******************************************************************************/
int planandnavigatexythetavsteer(PlannerType plannerType, char* envCfgFilename, char* motPrimFilename, bool forwardSearch)
{
	double allocated_time_secs_foreachplan = 300.0; // in seconds
	double initialEpsilon = 3.0;
	bool bsearchuntilfirstsolution = false;
	bool bforwardsearch = forwardSearch;

	double goaltol_x = 0.001, goaltol_y = 0.001, goaltol_theta = 0.001;

	bool bPrint = false;
	bool bPrintMap = true;

	EnvironmentNAVXYTHETAVSTEER environment_navxythetavsteer;
	EnvironmentNAVXYTHETAVSTEER trueenvironment_navxythetavsteer;

	// set the perimeter of the robot
	// it is given with 0, 0, 0 robot ref. point for which planning is done.
	vector<sbpl_2Dpt_t> perimeterptsV;
	sbpl_2Dpt_t pt_m;
	double halfwidth = 0.01;
	double halflength = 0.01;
	pt_m.x = -halflength;
	pt_m.y = -halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = halflength;
	pt_m.y = -halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = halflength;
	pt_m.y = halfwidth;
	perimeterptsV.push_back(pt_m);
	pt_m.x = -halflength;
	pt_m.y = halfwidth;
	perimeterptsV.push_back(pt_m);

	//	perimeterptsV.clear();

	// initialize true map from the environment file without perimeter or motion primitives
	if (!trueenvironment_navxythetavsteer.InitializeEnv(envCfgFilename, perimeterptsV, motPrimFilename)) {
		printf("ERROR: InitializeEnv failed\n");
		throw new SBPL_Exception();
	}

	// environment parameters
	int size_x = -1, size_y = -1, num_thetas = -1, num_v = -1, num_steers = -1;
	vector<double> velocities;
	
	velocities.push_back(-9.0);
	velocities.push_back(-3.0);
	velocities.push_back(-1.5);
	velocities.push_back(0.0);
	velocities.push_back(1.5);
	velocities.push_back(3.0);
	velocities.push_back(9.0);
	
	double startx = -1, starty = -1, starttheta = -1, startv = -1, startsteer = -1;
	double goalx = -1, goaly = -1, goaltheta = -1, goalv = -1, goalsteer = -1;
	double cellsize_m = 0.0;
	unsigned char obsthresh = 0, cost_inscribed_thresh = 0, cost_possibly_circumscribed_thresh = 0;
	vector<SBPL_xythetavsteer_mprimitive> motionprimitiveV;

	// get environment parameters from the true environment
	trueenvironment_navxythetavsteer.GetEnvParams(&size_x, &size_y, &num_thetas, &num_v, &num_steers, &velocities, &startx, &starty, &starttheta, &startv, &startsteer,
											&goalx, &goaly, &goaltheta, &goalv, &goalsteer, &cellsize_m, &cost_inscribed_thresh,
											&cost_possibly_circumscribed_thresh, &obsthresh, &motionprimitiveV);
	// print the map
	if (bPrintMap) {
		printf("true map:\n");
		for (int y = 0; y < size_y; y++) {
			for (int x = 0; x < size_x; x++) {
				printf("%3d ", trueenvironment_navxythetavsteer.GetMapCost(x, y));
			}
			printf("\n");
		}
		printf("System Pause (return=%d)\n", system("pause"));
	}

	// create an empty map
	unsigned char* map = new unsigned char[size_x * size_y];
	for (int i = 0; i < size_x * size_y; i++) {
		map[i] = 0;
	}

	// check the start and goal obtained from the true environment
	printf("start: %f %f %f %f %f, goal: %f %f %f %f %f\n", startx, starty, starttheta, startv, startsteer, goalx, goaly, goaltheta, goalv, goalsteer);

	// set robot environment parameters (should be done before initialize function is called)
	if (!environment_navxythetavsteer.SetEnvParameter("cost_inscribed_thresh", cost_inscribed_thresh)) {
		printf("ERROR: failed to set parameters\n");
		throw new SBPL_Exception();
	}
	if (!environment_navxythetavsteer.SetEnvParameter("cost_possibly_circumscribed_thresh", cost_possibly_circumscribed_thresh)) {
		printf("ERROR: failed to set parameters\n");
		throw new SBPL_Exception();
	}

	bool envInitialized = environment_navxythetavsteer.InitializeEnv(size_x, size_y, num_thetas, num_v, num_steers, velocities, map,
																	startx, starty, starttheta, startv, startsteer,
																	goalx, goaly, goaltol_theta, goalv, goalsteer,
																	perimeterptsV, cellsize_m, obsthresh, cost_inscribed_thresh,
																	cost_possibly_circumscribed_thresh, motPrimFilename);

	if (!envInitialized) {
		printf("ERROR: InitializeEnv failed\n");
		throw new SBPL_Exception();
	}

	// set start and goal states
	environment_navxythetavsteer.SetStart(startx, starty, starttheta, startv, startsteer);
	environment_navxythetavsteer.SetGoal(goalx, goaly, goaltheta, goalv, goalsteer);

	MDPConfig MDPCfg;

	// initialize MDP info
	if (!environment_navxythetavsteer.InitializeMDPCfg(&MDPCfg)) {
		printf("ERROR: InitializeMDPCfg failed\n");
		throw new SBPL_Exception();
	}

	// create a planner
	vector<int> solution_stateIDs_V;

	SBPLPlanner* planner = NULL;
	switch (plannerType) {
	case PLANNER_TYPE_ARASTAR:
		printf("Initializing ARAPlanner...\n");
		planner = new ARAPlanner(&environment_navxythetavsteer, bforwardsearch);
		break;
	case PLANNER_TYPE_ADSTAR:
		printf("Initializing ADPlanner...\n");
		planner = new ADPlanner(&environment_navxythetavsteer, bforwardsearch);
		break;
	case PLANNER_TYPE_RSTAR:
		printf("Invalid configuration: xytheta environment does not support rstar planner...\n");
		return 0;
	case PLANNER_TYPE_ANASTAR:
		printf("Initializing anaPlanner...\n");
		planner = new anaPlanner(&environment_navxythetavsteer, bforwardsearch);
		break;
	default:
		printf("Invalid planner type\n");
		break;
	}

	// set the start and goal states for the planner and other search variables
	if (planner->set_start(MDPCfg.startstateid) == 0) {
		printf("ERROR: failed to set start state\n");
		throw new SBPL_Exception();
	}
	if (planner->set_goal(MDPCfg.goalstateid) == 0) {
		printf("ERROR: failed to set goal state\n");
		throw new SBPL_Exception();
	}
	planner->set_initialsolution_eps(initialEpsilon);
	planner->set_search_mode(bsearchuntilfirstsolution);

	// compute sensing as a square surrounding the robot with length twice that of the
	// longest motion primitive
	vector<sbpl_2Dcell_t> sensecells;
	double maxMotPrimLengthSquared = 0.0;
	double maxMotPrimLength = 0.0;
	const EnvNAVXYTHETAVSTEERConfig_t* cfg = environment_navxythetavsteer.GetEnvNavConfig();
	for (int i = 0; i < (int)cfg->mprimV.size(); i++) {
		const SBPL_xythetavsteer_mprimitive& mprim = cfg->mprimV.at(i);
		int dx = mprim.endcell.x;
		int dy = mprim.endcell.y;
		if (dx * dx + dy * dy > maxMotPrimLengthSquared) {
			std::cout << "Found a longer motion primitive with dx = " << dx << " and dy = " << dy
				<< " from starttheta = " << (int)mprim.start_theta_disc << ", startv = " << (int)mprim.start_v_disc
				<< " and startsteer = " << (int)mprim.start_steer_disc << std::endl;
			maxMotPrimLengthSquared = dx * dx + dy * dy;
		}
	}
	maxMotPrimLength = sqrt((double)maxMotPrimLengthSquared);
	std::cout << "Maximum motion primitive length: " << maxMotPrimLength << std::endl;

	int sensingRange = (int)ceil(maxMotPrimLength);
	for (int x = -sensingRange; x <= sensingRange; x++) {
		for (int y = -sensingRange; y <= sensingRange; y++) {
			sensecells.push_back(sbpl_2Dcell_t(x, y));
		}
	}

	// create a file to hold the solution vector
	const char* sol = "solNavigating.txt";
	FILE* fSol = fopen(sol, "w");
	if (fSol == NULL) {
		printf("ERROR: could not open solution file\n");
		throw new SBPL_Exception();
	}

	// print the goal pose
	int goalx_c = CONTXY2DISC(goalx, cellsize_m);
	int goaly_c = CONTXY2DISC(goaly, cellsize_m);
	int goaltheta_c = ContTheta2Disc(goaltheta, num_thetas);
	int goalv_c = ContV2Disc(goalv, velocities);
	int goalsteer_c = ContSteer2Disc(goalsteer, num_steers);
	printf("goal_c: %d %d %d %d %d\n", goalx_c, goaly_c, goaltheta_c, goalv_c, goalsteer_c);

	vector<int> preds_of_changededgesIDV;
	vector<sbpl_2Dcell_t> changedcellsV;
	sbpl_2Dcell_t nav2dcell;
	vector<sbpl_xy_theta_v_steer_pt_t> xythetavsteerPath;

	// now comes the main loop
	while (fabs(startx - goalx) > goaltol_x || fabs(starty - goaly) > goaltol_y || fabs(starttheta - goaltheta)
		> goaltol_theta) {
		//simulate sensor data update
		bool bChanges = false;
		preds_of_changededgesIDV.clear();
		changedcellsV.clear();

		// simulate sensing the cells
		for (int i = 0; i < (int)sensecells.size(); i++) {
			int x = CONTXY2DISC(startx, cellsize_m) + sensecells.at(i).x;
			int y = CONTXY2DISC(starty, cellsize_m) + sensecells.at(i).y;

			// ignore if outside the map
			if (x < 0 || x >= size_x || y < 0 || y >= size_y) {
				continue;
			}

			int index = x + y * size_x;
			unsigned char truecost = trueenvironment_navxythetavsteer.GetMapCost(x, y);
			// update the cell if we haven't seen it before
			if (map[index] != truecost) {
				map[index] = truecost;
				environment_navxythetavsteer.UpdateCost(x, y, map[index]);
				printf("setting cost[%d][%d] to %d\n", x, y, map[index]);
				bChanges = true;
				// store the changed cells
				nav2dcell.x = x;
				nav2dcell.y = y;
				changedcellsV.push_back(nav2dcell);
			}
		}

		double TimeStarted = clock();

		// if necessary notify the planner of changes to costmap
		if (bChanges) {
			if (dynamic_cast<ARAPlanner*> (planner) != NULL || dynamic_cast<anaPlannerDouble*> (planner) != NULL) {
				((ARAPlanner*)planner)->costs_changed(); //use by ARA* and ANA* planner (non-incremental)
			}
			else if (dynamic_cast<ADPlanner*> (planner) != NULL) {
				// get the affected states
				environment_navxythetavsteer.GetPredsOfChangedEdges(&changedcellsV, &preds_of_changededgesIDV);
				// let know the incremental planner about them
				//use by AD* planner (incremental)
				((ADPlanner*)planner)->update_preds_of_changededges(&preds_of_changededgesIDV); 
				printf("%d states were affected\n", (int)preds_of_changededgesIDV.size());
			}
		}

		int startx_c = CONTXY2DISC(startx, cellsize_m);
		int starty_c = CONTXY2DISC(starty, cellsize_m);
		int starttheta_c = ContTheta2Disc(starttheta, num_thetas);
		int startv_c = ContV2Disc(startv, velocities);
		int startsteer_c = ContSteer2Disc(startsteer, num_steers);

		// plan a path
		bool bPlanExists = false;

		printf("new planning...\n");
		bPlanExists = (planner->replan(allocated_time_secs_foreachplan, &solution_stateIDs_V) == 1);
		printf("done with the solution of size=%d and sol. eps=%f\n", (unsigned int)solution_stateIDs_V.size(),
			planner->get_solution_eps());
		environment_navxythetavsteer.PrintTimeStat(stdout);

		// write the solution to sol.txt
		fprintf(fSol, "plan time=%.5f eps=%.2f\n", (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC),
				planner->get_solution_eps());
		fflush(fSol);

		xythetavsteerPath.clear();
		environment_navxythetavsteer.ConvertStateIDPathintoXYThetaVPath(&solution_stateIDs_V, &xythetavsteerPath);
		printf("actual path (with intermediate poses) size=%d\n", (unsigned int)xythetavsteerPath.size());
		for (unsigned int i = 0; i < xythetavsteerPath.size(); i++) {
			fprintf(fSol, "%.3f %.3f %.3f %.3f %.3f\n", xythetavsteerPath.at(i).x, xythetavsteerPath.at(i).y, xythetavsteerPath.at(i).theta, xythetavsteerPath.at(i).v, xythetavsteerPath.at(i).steer);
		}
		fprintf(fSol, "*********\n");

		for (int j = 1; j < (int)solution_stateIDs_V.size(); j++) {
			int newx = 0, newy = 0, newtheta = 0, newv = 0, newsteer = 0;
			environment_navxythetavsteer.GetCoordFromState(solution_stateIDs_V[j], newx, newy, newtheta, newv, newsteer);
			fprintf(fSol, "%d %d %d %d %d\n", newx, newy, newtheta, newv, newsteer);
		}
		fflush(fSol);

		// print the map (robot's view of the world and current plan)
		int startindex = startx_c + starty_c * size_x;
		int goalindex = goalx_c + goaly_c * size_x;
		for (int y = 0; bPrintMap && y < size_y; y++) {
			for (int x = 0; x < size_x; x++) {
				int index = x + y * size_x;
				int cost = map[index];
				cost = environment_navxythetavsteer.GetMapCost(x, y);

				// check to see if it is on the path
				bool bOnthePath = false;
				for (int j = 1; j < (int)solution_stateIDs_V.size(); j++) {
					int newx = 0, newy = 0, newtheta = 0, newv = 0, newsteer = 0;
					environment_navxythetavsteer.GetCoordFromState(solution_stateIDs_V[j], newx, newy, newtheta, newv, newsteer);
					if (x == newx && y == newy) bOnthePath = true;
				}

				if (index != startindex && index != goalindex && !bOnthePath) {
					printf("%3d ", cost);
				}
				else if (index == startindex) {
					printf("  X ");
				}
				else if (index == goalindex) {
					printf("  G ");
				}
				else if (bOnthePath) {
					printf("  * ");
				}
				else {
					printf("? ");
				}
			}
			printf("\n");
		}

		// move along the path
		if (bPlanExists && (int)xythetavsteerPath.size() > 1) {
			//get coord of the successor
			int newx, newy, newtheta, newv, newsteer;

			// move until we move into the end of motion primitive
			environment_navxythetavsteer.GetCoordFromState(solution_stateIDs_V[1], newx, newy, newtheta, newv, newsteer);

			printf("moving from %d %d %d %d %d to %d %d %d %d %d\n", startx_c, starty_c, starttheta_c, startv_c, startsteer_c, newx, newy, newtheta, newv, newsteer);

			// this check is weak since true configuration does not know the actual perimeter of the robot
			if (!trueenvironment_navxythetavsteer.IsValidConfiguration(newx, newy, newtheta, newv, newsteer)) {
				printf("ERROR: robot is commanded to move into an invalid configuration "
					"according to true environment\n");
				throw new SBPL_Exception();
			}

			// move
			startx = DISCXY2CONT(newx, cellsize_m);
			starty = DISCXY2CONT(newy, cellsize_m);
			starttheta = DiscTheta2Cont(newtheta, num_thetas);
			startv = DiscV2Cont(newv, velocities);
			startsteer = DiscSteer2Cont(newsteer, num_steers);

			// update the environment
			int newstartstateID = environment_navxythetavsteer.SetStart(startx, starty, starttheta, startv, startsteer);

			// update the planner
			if (planner->set_start(newstartstateID) == 0) {
				printf("ERROR: failed to update robot pose in the planner\n");
				throw new SBPL_Exception();
			}
		}
		else {
			printf("No move is made\n");
		}

		if (bPrint) {
			printf("System Pause (return=%d)\n", system("pause"));
		}
	}

	printf("goal reached!\n");

	fflush(NULL);
	fclose(fSol);

	delete map;
	delete planner;

	return 1;
}

/*******************************************************************************
* main - Parse command line arguments and launch one of the sbpl examples above.
*        Launch sbpl with just the -h option for usage help.
*
* @param argc The number of command-line arguments
* @param argv The command-line arguments
*******************************************************************************/
int main(int argc, char *argv[])
{
	// Print help
	if (argc == 2 && strcmp(argv[1], "-h") == 0) {
		PrintHelp(argv);
		return MAIN_RESULT_SUCCESS;
	}
	else if (argc < 5) {
		PrintUsage(argv);
		return MAIN_RESULT_INSUFFICIENT_ARGS;
	}

	// Find the number of options passed in; cfg and motprim files come after all options
	// options include any string starting with "-"
	int numOptions = 0;
	for (int i = 1; i < argc; i++) {
		if (strncmp("-", argv[i], 1) == 0) {
			numOptions++;
		}
		else {
			break;
		}
	}

	// Check command line arguments to find environment type and whether or not to
	// use one of the navigating examples.
	bool navigating = CheckIsNavigating(numOptions, argv);
	std::string plannerType = CheckPlannerType(numOptions, argv);
	std::string searchDir = CheckSearchDirection(numOptions, argv);

	std::cout << "Environment: xythetavsteer" << "; Planner: " << plannerType << "; Search direction: "
		<< searchDir << std::endl;

	int envArgIdx = numOptions + 1;
	PlannerType planner = StrToPlannerType(plannerType.c_str());
	bool forwardSearch = !strcmp(searchDir.c_str(), "forward");

	char* motPrimFilename = argv[envArgIdx + 1];

	// make sure we've been given valid a valid environment
	if (planner == INVALID_PLANNER_TYPE) {
		PrintUsage(argv);
		return MAIN_RESULT_INCORRECT_OPTIONS;
	}

	// Launch the correct example given the planner and an environment file.
	int plannerRes = 0;
	
	if (navigating) {
			plannerRes = planandnavigatexythetavsteer(planner, argv[envArgIdx], motPrimFilename, forwardSearch);
		}
		else {
			plannerRes = planxythetavsteer(planner, argv[envArgIdx], motPrimFilename, forwardSearch);
		}

	return plannerRes == 1 ? MAIN_RESULT_SUCCESS : MAIN_RESULT_FAILURE;
}
