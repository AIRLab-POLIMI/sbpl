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

#define PERFORMANCETEST_PERFORMANCE_TEST GLOBAL_PERFORMANCE_TEST
#define NUMBER_OF_GOALS 70

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <cstdio>
#include <sstream>

using namespace std;

#include <sbpl/headers.h>

enum EnvironmentType
{
    INVALID_ENV_TYPE = -1,
	ENV_TYPE_XYTHETA,
	ENV_TYPE_XYTHETAV,
	ENV_TYPE_XYTHETAVSTEER,
	ENV_TYPE_XYTHETAV_SAFE,
	ENV_TYPE_XYTHETAVSTEER_SAFE,

    NUM_ENV_TYPES
};

enum PlannerType
{
	INVALID_PLANNER_TYPE = -1,
	PLANNER_TYPE_ADSTAR,
	PLANNER_TYPE_ARASTAR,
	PLANNER_TYPE_ANASTARDOUBLE,

	NUM_PLANNER_TYPES
};

enum MainResultType
{
	INVALID_MAIN_RESULT = -1,

	MAIN_RESULT_SUCCESS = 0,
	MAIN_RESULT_FAILURE = 1,
	MAIN_RESULT_TOO_ARGS = 2,
	MAIN_RESULT_INCORRECT_OPTIONS = 3,
	MAIN_RESULT_UNSUPPORTED_ENV = 4,

	NUM_MAIN_RESULTS
};

typedef struct completepose{
	double x;
	double y;
	double theta;
	double v;
	double steer;
} completepose_t;

/*******************************************************************************
* PrintUsage - Prints the proper usage of the sbpl test executable.
*
* @param argv The command-line arguments; used to determine the name of the
*             test executable.
*******************************************************************************/
void PrintUsage(char *argv[])
{
	/*printf("USAGE: %s [-s] --planner=<planner_t> --search-dir=<search_t> <cfg file> <mot prims>\n",
		argv[0]);
	printf("See '%s -h' for help.\n", argv[0]);*/
	
	printf("USAGE: %s\n", argv[0]);
	printf("See '%s -h' for help.\n", argv[0]);
}

DiscreteSpaceInformation * initEnvironment(DiscreteSpaceInformation * environment, EnvironmentType envType, int width, int height, double * map, int obsthresh, int cost_inscribed_thresh, int cost_possibly_circumscribed_thresh, double nominalvel, double timetoturn, int numtheta, int numV, int numSteers, vector<double> * velocities, MDPConfig * MDPCfg, vector<sbpl_2Dpt_t> * perimeterptsV, completepose_t * startpose, completepose_t * goalpose, char * motionFileName, double cellsize_m){
	switch(envType){
		case ENV_TYPE_XYTHETA:
			printf("Environment: xytheta\n");
			environment = new EnvironmentNAVXYTHETALAT();
			
			((EnvironmentNAVXYTHETALAT *)environment)->SetEnvParameter("cost_inscribed_thresh", cost_inscribed_thresh);
			((EnvironmentNAVXYTHETALAT *)environment)->SetEnvParameter("cost_possibly_circumscribed_thresh", cost_possibly_circumscribed_thresh);
			
			((EnvironmentNAVXYTHETALAT *)environment)->InitializeEnv(width, height, map, startpose->x, startpose->y, startpose->theta, goalpose->x, height-goalpose->y, goalpose->theta, 0, 0, 0, *(perimeterptsV), cellsize_m, nominalvel, timetoturn, obsthresh, motionFileName);
			
			((EnvironmentNAVXYTHETALAT *)environment)->SetStart(startpose->x, startpose->y, startpose->theta);
			((EnvironmentNAVXYTHETALAT *)environment)->SetGoal(goalpose->x, goalpose->y, goalpose->theta);
			
			((EnvironmentNAVXYTHETALAT *)environment)->InitializeMDPCfg(MDPCfg);
			break;
		case ENV_TYPE_XYTHETAV:
		case ENV_TYPE_XYTHETAV_SAFE:
			printf("Environment: xythetav\n");
			environment = new EnvironmentNAVXYTHETAV();
			
			((EnvironmentNAVXYTHETAV *)environment)->InitializeEnv(width, height, numtheta, numV, *(velocities), map, startpose->x, startpose->y, startpose->theta, startpose->v, goalpose->x, goalpose->y, goalpose->theta, goalpose->v, *(perimeterptsV), cellsize_m, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, motionFileName);
			
			((EnvironmentNAVXYTHETAV *)environment)->SetStart(startpose->x, startpose->y, startpose->theta, startpose->v);
			((EnvironmentNAVXYTHETAV *)environment)->SetGoal(goalpose->x, goalpose->y, goalpose->theta, goalpose->v);
			
			((EnvironmentNAVXYTHETAV *)environment)->InitializeMDPCfg(MDPCfg);
			break;
		case ENV_TYPE_XYTHETAVSTEER:
		case ENV_TYPE_XYTHETAVSTEER_SAFE:
			printf("Environment: xythetavsteer\n");
			environment = new EnvironmentNAVXYTHETAVSTEER();
			
			((EnvironmentNAVXYTHETAVSTEER *)environment)->InitializeEnv(width, height, numtheta, numV, numSteers, *(velocities), map, startpose->x, startpose->y, startpose->theta, startpose->v, startpose->steer, goalpose->x, goalpose->y, goalpose->theta, goalpose->v, goalpose->steer, *(perimeterptsV), cellsize_m, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, motionFileName);
			
			((EnvironmentNAVXYTHETAVSTEER *)environment)->SetStart(startpose->x, startpose->y, startpose->theta, startpose->v, startpose->steer);
			((EnvironmentNAVXYTHETAVSTEER *)environment)->SetGoal(goalpose->x, goalpose->y, goalpose->theta, goalpose->v, goalpose->steer);
			
			((EnvironmentNAVXYTHETAVSTEER *)environment)->InitializeMDPCfg(MDPCfg);
			break;
		default:
			printf("Environment: ???\n");
			return 0;
			break;
	}
	
	return environment;
}

SBPLPlanner * initPlanner(SBPLPlanner * planner, PlannerType planType, DiscreteSpaceInformation * environment, bool bforwardsearch, bool bsearchuntilfirstsolution, double initialEpsilon, MDPConfig * MDPCfg){
	switch(planType){
		case PLANNER_TYPE_ARASTAR:
			printf("Planner: ARA*\n");
			planner = new ARAPlanner(environment, bforwardsearch);
			
			planner->set_initialsolution_eps(initialEpsilon);
			planner->set_start(MDPCfg->startstateid);
			planner->set_goal(MDPCfg->goalstateid);
			planner->set_search_mode(bsearchuntilfirstsolution);
			break;
		case PLANNER_TYPE_ANASTARDOUBLE:
			printf("Planner: ANA*\n");
			planner = new anaPlannerDouble(environment, bforwardsearch);
			
			planner->set_start(MDPCfg->startstateid);
			planner->set_goal(MDPCfg->goalstateid);
			planner->set_search_mode(bsearchuntilfirstsolution);
			break;
		case PLANNER_TYPE_ADSTAR:
			printf("Planner: AD*\n");
			planner = new ADPlanner(environment, bforwardsearch);
			
			planner->set_start(MDPCfg->startstateid);
			planner->set_goal(MDPCfg->goalstateid);
			planner->set_search_mode(bsearchuntilfirstsolution);
			planner->set_initialsolution_eps(initialEpsilon);
			break;
		default:
			printf("Planner: ???\n");
			return 0;
			break;
	}
	
	return planner;
}

int saveSolution(DiscreteSpaceInformation * environment, EnvironmentType envType, vector<int> * solution_stateIDs_V){
	int ret = 0;
	int height = 0;
	const char* solC = "solContinuous.txt";
	vector<sbpl_xy_theta_pt_t> xythetaPath;
	vector<sbpl_xy_theta_v_pt_t> xythetavPath;
	vector<sbpl_xy_theta_v_steer_pt_t> xythetavsteerPath;
	FILE* fSolC = fopen(solC, "w");
	if (fSolC == NULL) {
		printf("ERROR: could not open solution file\n");
		throw new SBPL_Exception();
	}
	
	switch(envType){
		case ENV_TYPE_XYTHETA:
			height = ((EnvironmentNAVXYTHETALAT *)environment)->GetEnvNavConfig()->EnvHeight_c;
			((EnvironmentNAVXYTHETALAT *)environment)->ConvertStateIDPathintoXYThetaPath(solution_stateIDs_V, &xythetaPath);
			for (unsigned int i = 0; i < xythetaPath.size(); i++) {
				fprintf(fSolC, "%.3f %.3f %.3f\n", xythetaPath.at(i).x, xythetaPath.at(i).y, xythetaPath.at(i).theta);
			}
			break;
		case ENV_TYPE_XYTHETAV:
		case ENV_TYPE_XYTHETAV_SAFE:
			((EnvironmentNAVXYTHETAV *)environment)->ConvertStateIDPathintoXYThetaVPath(solution_stateIDs_V, &xythetavPath);
			for (unsigned int i = 0; i < xythetavPath.size(); i++) {
				fprintf(fSolC, "%.3f %.3f %.3f %.3f\n", xythetavPath.at(i).x, xythetavPath.at(i).y, xythetavPath.at(i).theta, xythetavPath.at(i).v);
			}
			break;
		case ENV_TYPE_XYTHETAVSTEER:
		case ENV_TYPE_XYTHETAVSTEER_SAFE:
			((EnvironmentNAVXYTHETAVSTEER *)environment)->ConvertStateIDPathintoXYThetaVSteerPath(solution_stateIDs_V, &xythetavsteerPath);
			for (unsigned int i = 0; i < xythetavsteerPath.size(); i++) {
				fprintf(fSolC, "%.3f %.3f %.3f %.3f %.3f\n", xythetavsteerPath.at(i).x, xythetavsteerPath.at(i).y, xythetavsteerPath.at(i).theta, xythetavsteerPath.at(i).v, xythetavsteerPath.at(i).steer);
			}
			break;
		default:
			ret = -1;
			break;
	}
	
	fclose(fSolC);
	
	return ret;
}

int moveFiles(EnvironmentType envType, PlannerType planType, int i, int bret){
	stringstream destPerfPath, destSolPath;
	
	if(planType == PLANNER_TYPE_ARASTAR){
		if(envType == ENV_TYPE_XYTHETA){
			if(i<10){
				destPerfPath << "./Performance/arastar/xytheta/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xytheta/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/arastar/xytheta/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xytheta/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/arastar/xytheta/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xytheta/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/arastar/xytheta/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xytheta/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV){
			if(i<10){
				destPerfPath << "./Performance/arastar/xythetav/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/arastar/xythetav/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/arastar/xythetav/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/arastar/xythetav/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV_SAFE){
			if(i<10){
				destPerfPath << "./Performance/arastar/xythetav_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/arastar/xythetav_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/arastar/xythetav_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/arastar/xythetav_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetav_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER){
			if(i<10){
				destPerfPath << "./Performance/arastar/xythetavsteer/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/arastar/xythetavsteer/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/arastar/xythetavsteer/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/arastar/xythetavsteer/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER_SAFE){
			if(i<10){
				destPerfPath << "./Performance/arastar/xythetavsteer_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/arastar/xythetavsteer_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/arastar/xythetavsteer_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/arastar/xythetavsteer_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/arastar/xythetavsteer_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else{
			return -1;
		}
		
		rename("arastarperf.txt", destPerfPath.str().c_str());
	}
	else if(planType == PLANNER_TYPE_ANASTARDOUBLE){
		if(envType == ENV_TYPE_XYTHETA){
			if(i<10){
				destPerfPath << "./Performance/anastardouble/xytheta/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xytheta/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/anastardouble/xytheta/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xytheta/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/anastardouble/xytheta/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xytheta/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/anastardouble/xytheta/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xytheta/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV){
			if(i<10){
				destPerfPath << "./Performance/anastardouble/xythetav/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/anastardouble/xythetav/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/anastardouble/xythetav/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/anastardouble/xythetav/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV_SAFE){
			if(i<10){
				destPerfPath << "./Performance/anastardouble/xythetav_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/anastardouble/xythetav_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/anastardouble/xythetav_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/anastardouble/xythetav_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetav_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER){
			if(i<10){
				destPerfPath << "./Performance/anastardouble/xythetavsteer/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/anastardouble/xythetavsteer/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/anastardouble/xythetavsteer/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/anastardouble/xythetavsteer/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER_SAFE){
			if(i<10){
				destPerfPath << "./Performance/anastardouble/xythetavsteer_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/anastardouble/xythetavsteer_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/anastardouble/xythetavsteer_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/anastardouble/xythetavsteer_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/anastardouble/xythetavsteer_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else{
			return -1;
		}
		
		rename("anastardoubleperf.txt", destPerfPath.str().c_str());
	}
	else if(planType == PLANNER_TYPE_ADSTAR){
		if(envType == ENV_TYPE_XYTHETA){
			if(i<10){
				destPerfPath << "./Performance/adstar/xytheta/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xytheta/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/adstar/xytheta/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xytheta/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/adstar/xytheta/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xytheta/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/adstar/xytheta/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xytheta/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV){
			if(i<10){
				destPerfPath << "./Performance/adstar/xythetav/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/adstar/xythetav/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/adstar/xythetav/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/adstar/xythetav/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAV_SAFE){
			if(i<10){
				destPerfPath << "./Performance/adstar/xythetav_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/adstar/xythetav_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/adstar/xythetav_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/adstar/xythetav_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetav_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER){
			if(i<10){
				destPerfPath << "./Performance/adstar/xythetavsteer/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/adstar/xythetavsteer/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/adstar/xythetavsteer/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/adstar/xythetavsteer/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer/" << i << "_Sol" << bret << ".txt";
			}
		}
		else if(envType == ENV_TYPE_XYTHETAVSTEER_SAFE){
			if(i<10){
				destPerfPath << "./Performance/adstar/xythetavsteer_safe/000" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer_safe/000" << i << "_Sol" << bret << ".txt";
			}
			else if(i<100){
				destPerfPath << "./Performance/adstar/xythetavsteer_safe/00" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer_safe/00" << i << "_Sol" << bret << ".txt";
			}
			else if(i<1000){
				destPerfPath << "./Performance/adstar/xythetavsteer_safe/0" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer_safe/0" << i << "_Sol" << bret << ".txt";
			}
			else{
				destPerfPath << "./Performance/adstar/xythetavsteer_safe/" << i << "_Perf" << bret << ".txt";
				destSolPath << "./Performance/adstar/xythetavsteer_safe/" << i << "_Sol" << bret << ".txt";
			}
		}
		else{
			return -1;
		}
		
		rename("adstarperf.txt", destPerfPath.str().c_str());
	}
	else{
		return -1;
	}
	
	rename("solContinuous.txt", destSolPath.str().c_str());
	
	return 0;
}

int executePlan(EnvironmentType envType, PlannerType planType, int width, int height, double * map, int obsthresh, int cost_inscribed_thresh, int cost_possibly_circumscribed_thresh, double nominalvel, double timetoturn, int numtheta, int numV, int numSteers, vector<double> * velocities, vector<sbpl_2Dpt_t> * perimeterptsV, completepose_t * startpose, vector<completepose_t> * goalposes, char * motionFileName, double cellsize_m, bool bforwardsearch, bool bsearchuntilfirstsolution, double initialEpsilon, double allocated_time_secs_foreachplan){
	DiscreteSpaceInformation * environment;
	SBPLPlanner * planner = NULL;
	
	for(int i=0;i<goalposes->size();i++){
		MDPConfig MDPCfg;
		
		environment = initEnvironment(environment, envType, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, velocities, &MDPCfg, perimeterptsV, startpose, &goalposes->at(i), motionFileName, cellsize_m);
		planner = initPlanner(planner, planType, environment, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, &MDPCfg);
		
		vector<int> solution_stateIDs_V;
		int bret = planner->replan(allocated_time_secs_foreachplan, &solution_stateIDs_V);
		
		if(bret){
			printf("Solution found for goal %d...\n", i);
		}
		else{
			printf("Solution not found for goal %d...\n", i);
		}
		
		saveSolution(environment, envType, &solution_stateIDs_V);
		
		delete planner;
		delete environment;
		
		moveFiles(envType, planType, i, bret);
	}
}

vector<sbpl_2Dpt_t> createPerimeters(double halflength, double halfwidth){
	vector<sbpl_2Dpt_t> perimeterptsV;
	sbpl_2Dpt_t pt_m;
	
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
	
	return perimeterptsV;
}

int planTest(char * occupancyFileName, char * startFileName, char * goalsFileName){
	//Set various parameters for various environments
	int obsthresh = 1;
	int cost_inscribed_thresh = 1;
	int cost_possibly_circumscribed_thresh = 0;
	int cellsize_m = 1;
	double nominalvel = 3.0;
	double timetoturn = 1000000;
	int numtheta = 16;
	int numV = 7;
	vector<double> velocities;
	velocities.push_back(-9.0);
	velocities.push_back(-3.0);
	velocities.push_back(-1.5);
	velocities.push_back(0);
	velocities.push_back(1.5);
	velocities.push_back(3.0);
	velocities.push_back(9.0);
	int numSteers = 5;
	
	//Read map used to test
	int width;
	int height;
	double * map;
	FILE * occupancyFile = fopen(occupancyFileName, "r");
	char temp[1024];
	
	fscanf(occupancyFile, "%s", temp);
	fscanf(occupancyFile, "%d", &height);
	fscanf(occupancyFile, "%s", temp);
	fscanf(occupancyFile, "%d", &width);
	map = new double[width * height];
	int i = 0;
	while(!feof(occupancyFile) && i < width*height){
		double c;
		fscanf(occupancyFile, "%lf", &c);
		if(c >= 0 && c <= 1){
			map[i] = c;
			i++;
		}
	}
	
	fclose(occupancyFile);
	
	//Read start state
	completepose_t startpose;
	float dTemp;
	occupancyFile = fopen(startFileName, "r");
	
	fscanf(occupancyFile, "%lf", &(startpose.x));
	fscanf(occupancyFile, "%lf", &(startpose.y));
	fscanf(occupancyFile, "%lf", &(startpose.theta));
	fscanf(occupancyFile, "%lf", &(startpose.v));
	fscanf(occupancyFile, "%lf", &(startpose.steer));
	
	fclose(occupancyFile);
	
	//Read all end states
	vector<completepose_t> goalposes;
	occupancyFile = fopen(goalsFileName, "r");
	int contgoals = 0;
	
	while(!feof(occupancyFile) && contgoals < NUMBER_OF_GOALS){
		completepose_t temppose;
		fscanf(occupancyFile, "%lf", &(temppose.x));
		fscanf(occupancyFile, "%lf", &(temppose.y));
		fscanf(occupancyFile, "%lf", &(temppose.theta));
		fscanf(occupancyFile, "%lf", &(temppose.v));
		fscanf(occupancyFile, "%lf", &(temppose.steer));
		goalposes.push_back(temppose);
		contgoals++;
	}
	
	fclose(occupancyFile);

	for(int contgoals=0;i<goalposes.size();contgoals++)
		goalposes.at(contgoals).v = 0.0;
	
	//Set parameters for planner
	double allocated_time_secs_foreachplan = 1500.0; // in seconds
	double initialEpsilon = 3.0;
	bool bsearchuntilfirstsolution = false;
	bool bforwardsearch = true;
	
	// set the perimeter of the robot
	vector<sbpl_2Dpt_t> perimeterptsV = createPerimeters(1, 1);
	
	//ARA* and xyt
	//executePlan(ENV_TYPE_XYTHETA, PLANNER_TYPE_ARASTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xytheta_prim/xytheta_not_unif_min_forward.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ANA* and xyt
	//executePlan(ENV_TYPE_XYTHETA, PLANNER_TYPE_ANASTARDOUBLE, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xytheta_prim/xytheta_not_unif_min_forward.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//AD* and xyt
	executePlan(ENV_TYPE_XYTHETA, PLANNER_TYPE_ADSTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xytheta_prim/xytheta_not_unif_min_forward.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ARA* and xytv
	//executePlan(ENV_TYPE_XYTHETAV, PLANNER_TYPE_ARASTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/theta_not_reg_rid_pos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ANA* and xytv
	//executePlan(ENV_TYPE_XYTHETAV, PLANNER_TYPE_ANASTARDOUBLE, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/theta_not_reg_rid_pos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);

	//AD* and xytv
	executePlan(ENV_TYPE_XYTHETAV, PLANNER_TYPE_ADSTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/theta_not_reg_rid_pos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ARA* and xytvp
	//executePlan(ENV_TYPE_XYTHETAVSTEER, PLANNER_TYPE_ARASTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetavsteer_prim/intro_steer_reduced_vpos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ANA* and xytvp
	//executePlan(ENV_TYPE_XYTHETAVSTEER, PLANNER_TYPE_ANASTARDOUBLE, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetavsteer_prim/intro_steer_reduced_vpos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//AD* and xytvp
	executePlan(ENV_TYPE_XYTHETAVSTEER, PLANNER_TYPE_ADSTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetavsteer_prim/intro_steer_reduced_vpos.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ARA* and xytv_safe
	//executePlan(ENV_TYPE_XYTHETAV_SAFE, PLANNER_TYPE_ARASTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/xythetav_safe.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//ANA* and xytv_safe
	//executePlan(ENV_TYPE_XYTHETAV_SAFE, PLANNER_TYPE_ANASTARDOUBLE, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/xythetav_safe.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	//AD* and xytv_safe
	velocities.clear();
	velocities.push_back(-9.0);
        velocities.push_back(-3.0);
        velocities.push_back(-1.5);
        velocities.push_back(0.1);
        velocities.push_back(1.5);
        velocities.push_back(3.0);
        velocities.push_back(9.0);

	for(int contgoals=0;i<goalposes.size();contgoals++)
                goalposes.at(contgoals).v = 0.1;

	executePlan(ENV_TYPE_XYTHETAV_SAFE, PLANNER_TYPE_ADSTAR, width, height, map, obsthresh, cost_inscribed_thresh, cost_possibly_circumscribed_thresh, nominalvel, timetoturn, numtheta, numV, numSteers, &velocities, &perimeterptsV, &startpose, &goalposes, "../myprimitives/xythetav_prim/xythetav_safe.txt", cellsize_m, bforwardsearch, bsearchuntilfirstsolution, initialEpsilon, allocated_time_secs_foreachplan);
	
	delete [] map;
	
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
		PrintUsage(argv);
		return MAIN_RESULT_SUCCESS;
	}
	else if (argc != 1) {
		PrintUsage(argv);
		return MAIN_RESULT_TOO_ARGS;
	}

	// Launch the correct example given the planner and an environment file.
	int plannerRes = 0;
	
	plannerRes = planTest("./occupancyFile.txt", "./startFile.txt", "./goalsFile.txt");

	return plannerRes == 1 ? MAIN_RESULT_SUCCESS : MAIN_RESULT_FAILURE;
}
