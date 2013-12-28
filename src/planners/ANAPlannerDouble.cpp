/*
 * This code was used for generating experimental data for the purpose of understanding the performance of
 * the Anytime Nonparametric A* (ANA*) algorithm.  
 * 
 */

#include <cmath>
#include <limits>
#include <sbpl/discrete_space_information/environment.h>
#include <sbpl/planners/ANAPlannerDouble.h>
#include <sbpl/utils/heapdouble.h>
#include <sbpl/utils/list.h>
//#include <sstream>

using namespace std;

//-----------------------------------------------------------------------------------------------------

anaPlannerDouble::anaPlannerDouble(DiscreteSpaceInformation* environment, bool bSearchForward)
{
    bforwardsearch = bSearchForward;

    environment_ = environment;

    bsearchuntilfirstsolution = false;
    finitial_eps = MAX_DOUBLE;
    searchexpands = 0;
    MaxMemoryCounter = 0;

    fDeb = fopen("debug.txt", "w");

    pSearchStateSpace_ = new anaDoubleSearchStateSpace_t;

    //create the ana planner
    if (CreateSearchStateSpace(pSearchStateSpace_) != 1) {
        printf("ERROR: failed to create statespace\n");
        return;
    }

    //set the start and goal states
    if (InitializeSearchStateSpace(pSearchStateSpace_) != 1) {
        printf("ERROR: failed to create statespace\n");
        return;
    }
}

anaPlannerDouble::~anaPlannerDouble()
{
    if (pSearchStateSpace_ != NULL) {
        //delete the statespace
        DeleteSearchStateSpace( pSearchStateSpace_);
        delete pSearchStateSpace_;
    }
    fclose( fDeb);
}

void anaPlannerDouble::Initialize_searchinfo(CMDPSTATE* state, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    anaDoubleState* searchstateinfo = (anaDoubleState*)state->PlannerSpecificData;

    searchstateinfo->MDPstate = state;
    InitializeSearchStateInfo(searchstateinfo, pSearchStateSpace);
}

CMDPSTATE* anaPlannerDouble::CreateState(int stateID, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* state = NULL;

#if DEBUG
    if(environment_->StateID2IndexMapping[stateID][anaMDP_STATEID2IND] != -1)
    {
        printf("ERROR in CreateState: state already created\n");
        exit(1);
    }
#endif

    //adds to the tail a state
    state = pSearchStateSpace->searchMDP.AddState(stateID);

    //remember the index of the state
    environment_->StateID2IndexMapping[stateID][anaMDP_STATEID2IND] =
            pSearchStateSpace->searchMDP.StateArray.size() - 1;

#if DEBUG
    if (state != 
        pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][anaMDP_STATEID2IND]])
    {
        printf("ERROR in CreateState: invalid state index\n");
        exit(1);
    }
#endif

    //create search specific info
    state->PlannerSpecificData = (anaDoubleState*)malloc(sizeof(anaDoubleState));
    Initialize_searchinfo(state, pSearchStateSpace);
    MaxMemoryCounter += sizeof(anaDoubleState);

    return state;
}

CMDPSTATE* anaPlannerDouble::GetState(int stateID, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    if (stateID >= (int)environment_->StateID2IndexMapping.size()) {
        SBPL_ERROR("ERROR in GetState: stateID %d is invalid\n", stateID);
        throw new SBPL_Exception();
    }

    if (environment_->StateID2IndexMapping[stateID][anaMDP_STATEID2IND] == -1)
        return CreateState(stateID, pSearchStateSpace);
    else
        return pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][anaMDP_STATEID2IND]];
}

//-----------------------------------------------------------------------------------------------------

int anaPlannerDouble::ComputeHeuristic(CMDPSTATE* MDPstate, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    //compute heuristic for search

    if (bforwardsearch) {

#if MEM_CHECK == 1
        //int WasEn = DisableMemCheck();
#endif

        //forward search: heur = distance from state to searchgoal which is Goal anaState
        int retv = environment_->GetGoalHeuristic(MDPstate->StateID);

#if MEM_CHECK == 1
        //if (WasEn)
        //	EnableMemCheck();
#endif

        return retv;

    }
    else {
        //backward search: heur = distance from searchgoal to state
        return environment_->GetStartHeuristic(MDPstate->StateID);
    }
}

//initialization of a state
void anaPlannerDouble::InitializeSearchStateInfo(anaDoubleState* state, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    state->g = std::numeric_limits<unsigned int>::max();
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->heapindex = 0;
    state->listelem[ana_INCONS_LIST_ID] = 0;
    state->numofexpands = 0;

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR
    if(pSearchStateSpace->searchgoalstate != NULL)
		state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    else
		state->h = 0;
#else
    state->h = 0;
#endif
}

//re-initialization of a state
void anaPlannerDouble::ReInitializeSearchStateInfo(anaDoubleState* state, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    state->g = std::numeric_limits<unsigned int>::max();
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->heapindex = 0;
    state->listelem[ana_INCONS_LIST_ID] = 0;
    state->numofexpands = 0;

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR

    if (pSearchStateSpace->searchgoalstate != NULL) 
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    else
		state->h = 0;
#else

    state->h = 0;

#endif

}

void anaPlannerDouble::DeleteSearchStateData(anaDoubleState* state)
{
    //no memory was allocated
    MaxMemoryCounter = 0;
    return;
}

double anaPlannerDouble::get_e_value(anaDoubleSearchStateSpace_t* pSearchStateSpace, int stateID)
{

    CMDPSTATE* MDPstate = GetState(stateID, pSearchStateSpace);
    anaDoubleState* searchstateinfo = (anaDoubleState*)MDPstate->PlannerSpecificData;

    if (searchstateinfo->h == 0) {
        if (searchstateinfo->g >= pSearchStateSpace->G) {
            return 0.0;
        }
        else {
            return std::numeric_limits<double>::max();
        }
    }
    else {
        return ((double)(pSearchStateSpace->G - searchstateinfo->g)) / (double)searchstateinfo->h;
    }
}

//used for backward search
void anaPlannerDouble::UpdatePreds(anaDoubleState* state, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    vector<int> PredIDV;
    vector<int> CostV;
    CKeyDouble key;
    anaDoubleState *p;

    environment_->GetPreds(state->MDPstate->StateID, &PredIDV, &CostV);

    //iterate through predecessors of s
    for (int pind = 0; pind < (int)PredIDV.size(); pind++) {
        CMDPSTATE* PredMDPState = GetState(PredIDV[pind], pSearchStateSpace);
        p = (anaDoubleState*)(PredMDPState->PlannerSpecificData);
        if (p->callnumberaccessed != pSearchStateSpace->callnumber) ReInitializeSearchStateInfo(p, pSearchStateSpace);

        //see if we can improve the value of p

        if (p->g > state->g + CostV[pind]) {
            p->g = state->g + CostV[pind];
            p->bestnextstate = state->MDPstate;

			if(p->g + p->h < pSearchStateSpace->G){
				key.key[0] = -get_e_value(pSearchStateSpace, p->MDPstate->StateID);
				if (pSearchStateSpace->heap->inheap(p)) {
					pSearchStateSpace->heap->updateheap(p, key);
				}
				else {
					pSearchStateSpace->heap->insertheap(p, key);
				}
			}
        }
    }
}

//used for forward search
void anaPlannerDouble::UpdateSuccs(anaDoubleState* state, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    CKeyDouble key;
    anaDoubleState *n;

    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);

    //iterate through predecessors of s
    for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
        CMDPSTATE* SuccMDPState = GetState(SuccIDV[sind], pSearchStateSpace);
        int cost = CostV[sind];

        n = (anaDoubleState*)(SuccMDPState->PlannerSpecificData);
        if (n->callnumberaccessed != pSearchStateSpace->callnumber) ReInitializeSearchStateInfo(n, pSearchStateSpace);

        //see if we can improve the value of n
        //taking into account the cost of action
        if (n->g > state->g + cost) {
            n->g = state->g + cost;
            n->bestpredstate = state->MDPstate;

			if(n->g + n->h < pSearchStateSpace->G){
				key.key[0] = -get_e_value(pSearchStateSpace, n->MDPstate->StateID);
				
				if (pSearchStateSpace->heap->inheap(n)) {
					pSearchStateSpace->heap->updateheap(n, key);
				}
				else {
					pSearchStateSpace->heap->insertheap(n, key);
				}
			}
        }
    }
}

//TODO-debugmax - add obsthresh and other thresholds to other environments in 3dkin
int anaPlannerDouble::GetGVal(int StateID, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* cmdp_state = GetState(StateID, pSearchStateSpace);
    anaDoubleState* state = (anaDoubleState*)cmdp_state->PlannerSpecificData;
    return state->g;
}

//returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int anaPlannerDouble::ImprovePath(anaDoubleSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs)
{
    int expands;
    anaDoubleState *state, *searchgoalstate;
    CKeyDouble key, maxkey;

    expands = 0;

    if (pSearchStateSpace->searchgoalstate == NULL) {
        SBPL_ERROR("ERROR searching: no goal state is set\n");
        throw new SBPL_Exception();
    }

    //goal state
    searchgoalstate = (anaDoubleState*)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
    if (searchgoalstate->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo( searchgoalstate, pSearchStateSpace);
    }

    //expand states until done
    while (!pSearchStateSpace->heap->emptyheap() &&
           (clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC)
    {
		//get the key
		maxkey.key[0] = -(pSearchStateSpace->heap->getminkeyheap().key[0]);
        //get the state
        state = (anaDoubleState*)pSearchStateSpace->heap->deleteminheap();

        if (maxkey.key[0] < pSearchStateSpace->eps) {
            pSearchStateSpace->eps = maxkey.key[0];
        }
		
        if (state->MDPstate->StateID == searchgoalstate->MDPstate->StateID) {
            pSearchStateSpace->G = state->g;
            searchexpands += expands;
            return 1;
        }

#if DEBUG
        //fprintf(fDeb, "expanding state(%d): h=%d g=%u key=%u v=%u iterc=%d callnuma=%d expands=%d (g(goal)=%u)\n",
        //	state->MDPstate->StateID, state->h, state->g, state->g+(int)(pSearchStateSpace->eps*state->h), state->v,
        //	state->iterationclosed, state->callnumberaccessed, state->numofexpands, searchgoalstate->g);
        //fprintf(fDeb, "expanding: ");
        //PrintSearchState(state, fDeb);
        if (state->listelem[ana_INCONS_LIST_ID] != NULL) {
            fprintf(fDeb, "ERROR: expanding a state from inconslist\n");
            printf("ERROR: expanding a state from inconslist\n");
            exit(1);
        }
        //fflush(fDeb);
#endif
		
        state->iterationclosed = pSearchStateSpace->searchiteration;

        //new expand
        expands++;
        state->numofexpands++;

        if (bforwardsearch == false)
            UpdatePreds(state, pSearchStateSpace);
        else
            UpdateSuccs(state, pSearchStateSpace);
    }

    int retv = 3;
    if (searchgoalstate->g == std::numeric_limits<unsigned int>::max() && pSearchStateSpace->heap->emptyheap()) {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("solution does not exist: search exited because heap is empty\n");
#endif
        retv = 0;
    }
    else if (!pSearchStateSpace->heap->emptyheap() && 0 < maxkey.key[0]) {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("search exited because it ran out of time\n");
#endif
        retv = 2;
    }
    else if (searchgoalstate->g == std::numeric_limits<unsigned int>::max() && !pSearchStateSpace->heap->emptyheap()) {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("solution does not exist: search exited because all "
               "candidates for expansion have infinite heuristics\n");
#endif
        retv = 0;
    }

    searchexpands += expands;

    return retv;
}

void anaPlannerDouble::Reevaluatefvals(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    CKeyDouble key;
    int i;
    CHeapDouble* pheap = pSearchStateSpace->heap;

    //recompute priorities for states in OPEN and reorder it
    for (i = 1; i <= pheap->currentsize; ++i) {
        //anaState* state = (anaState*)pheap->heap[i].heapstate;

        // CHANGED - cast removed

        pheap->heap[i].key.key[0] = -get_e_value(pSearchStateSpace,
                                                       ((anaDoubleState*)pheap->heap[i].heapstate)->MDPstate->StateID);

        //pheap->heap[i].key.key[1] = state->h;
    }
    pheap->makeheap();

    pSearchStateSpace->bReevaluatefvals = false;
}

//creates (allocates memory) search state space
//does not initialize search statespace
int anaPlannerDouble::CreateSearchStateSpace(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    //create a heap
    pSearchStateSpace->heap = new CHeapDouble;
    //pSearchStateSpace->inconslist = new CHeap;
    MaxMemoryCounter += sizeof(CHeapDouble);
    MaxMemoryCounter += sizeof(CList);

    pSearchStateSpace->searchgoalstate = NULL;
    pSearchStateSpace->searchstartstate = NULL;

    searchexpands = 0;

    pSearchStateSpace->bReinitializeSearchStateSpace = false;

    return 1;
}

//deallocates memory used by SearchStateSpace
void anaPlannerDouble::DeleteSearchStateSpace(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap != NULL) {
        pSearchStateSpace->heap->makeemptyheap();
        delete pSearchStateSpace->heap;
        pSearchStateSpace->heap = NULL;
    }

    /*
    if(pSearchStateSpace->inconslist != NULL)
    {
        pSearchStateSpace->inconslist->makeemptyheap();
        delete pSearchStateSpace->inconslist;
        pSearchStateSpace->inconslist = NULL;
    }
    */

    //delete the states themselves
    int iend = (int)pSearchStateSpace->searchMDP.StateArray.size();
    for (int i = 0; i < iend; i++) {
        CMDPSTATE* state = pSearchStateSpace->searchMDP.StateArray[i];
        if (state != NULL && state->PlannerSpecificData != NULL) {
            DeleteSearchStateData((anaDoubleState*)state->PlannerSpecificData);
            free((anaDoubleState*)state->PlannerSpecificData);
            state->PlannerSpecificData = NULL;
        }
    }
    pSearchStateSpace->searchMDP.Delete();
}

//reset properly search state space
//needs to be done before deleting states
int anaPlannerDouble::ResetSearchStateSpace(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    pSearchStateSpace->heap->makeemptyheap();
    //	pSearchStateSpace->inconslist->makeemptyheap();

    return 1;
}

//initialization before each search
void anaPlannerDouble::ReInitializeSearchStateSpace(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    CKeyDouble key;

    //increase callnumber
    pSearchStateSpace->callnumber++;

    //reset iteration
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;
    pSearchStateSpace->G = std::numeric_limits<unsigned int>::max();

#if DEBUG
    fprintf(fDeb, "reinitializing search state-space (new call number=%d search iter=%d)\n",
        pSearchStateSpace->callnumber,pSearchStateSpace->searchiteration );
#endif

    pSearchStateSpace->heap->makeemptyheap();
    //pSearchStateSpace->inconslist->makeemptyheap();
    //reset
    pSearchStateSpace->eps = this->finitial_eps;

    //initialize start state
    anaDoubleState* startstateinfo = (anaDoubleState*)(pSearchStateSpace->searchstartstate->PlannerSpecificData);
    if (startstateinfo->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo( startstateinfo, pSearchStateSpace);
    }

    startstateinfo->g = 0;

    //insert start state into the heap

    // CHANGED - long int cast removed
    key.key[0] = -get_e_value(pSearchStateSpace, startstateinfo->MDPstate->StateID); 

    //key.key[1] = startstateinfo->h;
    pSearchStateSpace->heap->insertheap(startstateinfo, key);

    pSearchStateSpace->bReinitializeSearchStateSpace = false;
    pSearchStateSpace->bReevaluatefvals = false;
}

//very first initialization
int anaPlannerDouble::InitializeSearchStateSpace(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap->currentsize != 0) {
        SBPL_ERROR("ERROR in InitializeSearchStateSpace: heap or list is not empty\n");
        throw new SBPL_Exception();
    }

    pSearchStateSpace->eps = this->finitial_eps;
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;
    pSearchStateSpace->callnumber = 0;
    pSearchStateSpace->bReevaluatefvals = false;

    pSearchStateSpace->G = std::numeric_limits<unsigned int>::max();

    //create and set the search start state
    pSearchStateSpace->searchgoalstate = NULL;
    //pSearchStateSpace->searchstartstate = GetState(SearchStartStateID, pSearchStateSpace);
    pSearchStateSpace->searchstartstate = NULL;

    pSearchStateSpace->bReinitializeSearchStateSpace = true;

    return 1;
}

int anaPlannerDouble::SetSearchGoalState(int SearchGoalStateID, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->searchgoalstate == NULL ||
        pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID)
    {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

        //should be new search iteration
        pSearchStateSpace->bNewSearchIteration = true;
        pSearchStateSpace_->eps = this->finitial_eps;

        //recompute heuristic for the heap if heuristics is used
#if USE_HEUR
        for(int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++)
        {
            CMDPSTATE* MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
            anaDoubleState* state = (anaDoubleState*)MDPstate->PlannerSpecificData;
            state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
        }

        pSearchStateSpace->bReevaluatefvals = true;
#endif
    }

    return 1;
}

int anaPlannerDouble::SetSearchStartState(int SearchStartStateID, anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* MDPstate = GetState(SearchStartStateID, pSearchStateSpace);

    if (MDPstate != pSearchStateSpace->searchstartstate) {
        pSearchStateSpace->searchstartstate = MDPstate;
        pSearchStateSpace->bReinitializeSearchStateSpace = true;
    }

    return 1;
}

int anaPlannerDouble::ReconstructPath(anaDoubleSearchStateSpace_t* pSearchStateSpace)
{
    if (bforwardsearch) //nothing to do, if search is backward
    {
        CMDPSTATE* MDPstate = pSearchStateSpace->searchgoalstate;
        CMDPSTATE* PredMDPstate;
        anaDoubleState *predstateinfo, *stateinfo;

#if DEBUG
        fprintf(fDeb, "reconstructing a path:\n");
#endif

        while (MDPstate != pSearchStateSpace->searchstartstate) {
            stateinfo = (anaDoubleState*)MDPstate->PlannerSpecificData;

#if DEBUG
            PrintSearchState(stateinfo, fDeb);
#endif
            if (stateinfo->g == std::numeric_limits<unsigned int>::max()) {
                return -1;
            }

            if (stateinfo->bestpredstate == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                throw new SBPL_Exception();
            }

            //get the parent state
            PredMDPstate = stateinfo->bestpredstate;
            predstateinfo = (anaDoubleState*)PredMDPstate->PlannerSpecificData;

            //set its best next info
            predstateinfo->bestnextstate = MDPstate;

            //transition back
            MDPstate = PredMDPstate;
        }
    }

    return 1;
}

void anaPlannerDouble::PrintSearchPath(anaDoubleSearchStateSpace_t* pSearchStateSpace, FILE* fOut)
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    anaDoubleState* searchstateinfo;
    CMDPSTATE* state;
    int goalID;
    int PathCost;

    if (bforwardsearch) {
        state = pSearchStateSpace->searchstartstate;
        goalID = pSearchStateSpace->searchgoalstate->StateID;
    }
    else {
        state = pSearchStateSpace->searchgoalstate;
        goalID = pSearchStateSpace->searchstartstate->StateID;
    }
    if (fOut == NULL) fOut = stdout;

    PathCost = ((anaDoubleState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;

    fprintf(fOut, "Printing a path from state %d to the goal state %d\n", state->StateID,
            pSearchStateSpace->searchgoalstate->StateID);
    fprintf(fOut, "Path cost = %d:\n", PathCost);

    environment_->PrintState(state->StateID, false, fOut);

    int costFromStart = 0;
    while (state->StateID != goalID) {
        fprintf(fOut, "state %d ", state->StateID);

        if (state->PlannerSpecificData == NULL) {
            fprintf(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (anaDoubleState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            fprintf(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == std::numeric_limits<unsigned int>::max()) {
            fprintf(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        int costToGoal = PathCost - costFromStart;
        int transcost = searchstateinfo->g - ((anaDoubleState*)(searchstateinfo->bestnextstate->PlannerSpecificData))->g;
        if (bforwardsearch) transcost = -transcost;

        costFromStart += transcost;

        fprintf(fOut, "g=%d-->state %d, h = %d ctg = %d  ", searchstateinfo->g,
                searchstateinfo->bestnextstate->StateID, searchstateinfo->h, costToGoal);

        state = searchstateinfo->bestnextstate;

        environment_->PrintState(state->StateID, false, fOut);
    }
#endif
}

void anaPlannerDouble::PrintSearchState(anaDoubleState* state, FILE* fOut)
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    fprintf(fOut, "state %d: h=%d g=%u iterc=%d callnuma=%d expands=%d heapind=%d inconslist=%d\n",
            state->MDPstate->StateID, state->h, state->g, state->iterationclosed, state->callnumberaccessed,
            state->numofexpands, state->heapindex, state->listelem[ana_INCONS_LIST_ID] ? 1 : 0);
    environment_->PrintState(state->MDPstate->StateID, true, fOut);
#endif
}

int anaPlannerDouble::getHeurValue(anaDoubleSearchStateSpace_t* pSearchStateSpace, int StateID)
{
    CMDPSTATE* MDPstate = GetState(StateID, pSearchStateSpace);
    anaDoubleState* searchstateinfo = (anaDoubleState*)MDPstate->PlannerSpecificData;
    return searchstateinfo->h;
}

vector<int> anaPlannerDouble::GetSearchPath(anaDoubleSearchStateSpace_t* pSearchStateSpace, int& solcost)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    vector<int> wholePathIds;
    anaDoubleState* searchstateinfo;
    CMDPSTATE* state = NULL;
    CMDPSTATE* goalstate = NULL;
    CMDPSTATE* startstate = NULL;

    if (bforwardsearch) {
        startstate = pSearchStateSpace->searchstartstate;
        goalstate = pSearchStateSpace->searchgoalstate;

        //reconstruct the path by setting bestnextstate pointers appropriately
        ReconstructPath(pSearchStateSpace);
    }
    else {
        startstate = pSearchStateSpace->searchgoalstate;
        goalstate = pSearchStateSpace->searchstartstate;
    }

    state = startstate;

    wholePathIds.push_back(state->StateID);
    solcost = 0;

    FILE* fOut = stdout;
    while (state->StateID != goalstate->StateID) {
        if (state->PlannerSpecificData == NULL) {
            fprintf(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (anaDoubleState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            fprintf(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == INFINITECOST) {
            fprintf(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        environment_->GetSuccs(state->StateID, &SuccIDV, &CostV);
        int actioncost = std::numeric_limits<int>::max();
        for (int i = 0; i < (int)SuccIDV.size(); i++) {

            if (SuccIDV.at(i) == searchstateinfo->bestnextstate->StateID) actioncost = CostV.at(i);

        }
        if (actioncost == INFINITECOST) printf("WARNING: actioncost = %d\n", actioncost);

        solcost += actioncost;

#if DEBUG
        anaDoubleState* nextstateinfo = (anaDoubleState*)(searchstateinfo->bestnextstate->PlannerSpecificData);
        if(actioncost != abs((int)(searchstateinfo->g - nextstateinfo->g)))
        {
            fprintf(fDeb, "ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            printf("ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                actioncost,abs((int)(searchstateinfo->g - nextstateinfo->g)));
            PrintSearchState(searchstateinfo, fDeb);
            PrintSearchState(nextstateinfo, fDeb);
        }
#endif

        state = searchstateinfo->bestnextstate;

        wholePathIds.push_back(state->StateID);
    }

    return wholePathIds;
}

bool anaPlannerDouble::Search(anaDoubleSearchStateSpace_t* pSearchStateSpace, vector<int>& pathIds, int & PathCost,
                        bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs)
{
    CKeyDouble key;
    TimeStarted = clock();
    searchexpands = 0;
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 1
	vector<PlannerStatsExtended> planstats;
#endif

#if DEBUG
    fprintf(fDeb, "new search call (call number=%d)\n", pSearchStateSpace->callnumber);
#endif

    if (pSearchStateSpace->bReinitializeSearchStateSpace == true) {
        //re-initialize state space 
        ReInitializeSearchStateSpace(pSearchStateSpace);
    }

    /*if (bOptimalSolution) {
        pSearchStateSpace->eps = 1;
        MaxNumofSecs = INFINITECOST;
    }
    else */if (bFirstSolution) {
        MaxNumofSecs = INFINITECOST;
    }

    //ensure heuristics are up-to-date
    environment_->EnsureHeuristicsUpdated((bforwardsearch == true));

    //the main loop of ana*
    int prevexpands = 0;
    clock_t loop_time;

    // CHANGE MADE TO WHILE LOOP to account for open.empty() == FALSE
    while (!pSearchStateSpace->heap->emptyheap() &&
           (clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC)
    {
        loop_time = clock();

        pSearchStateSpace->searchiteration++;
        pSearchStateSpace->bNewSearchIteration = false;

        //improve or compute path
        int retVal = ImprovePath(pSearchStateSpace, MaxNumofSecs);
        anaDoubleState* state;
        CKeyDouble key;
        CHeapDouble* open = pSearchStateSpace->heap;

		int j=1;
        while(j <= open->currentsize) {
            state = (anaDoubleState*)open->heap[j].heapstate;
            
            if (state->g + state->h >= pSearchStateSpace->G) {
                open->deleteheap_unsafe(state);
            }
            else {
                key.key[0] = -get_e_value(pSearchStateSpace, state->MDPstate->StateID);

                open->updateheap_unsafe(state, key);
                ++j;
            }
        }
        open->makeheap();
		
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        //print the solution cost and eps bound
        if (retVal == 1) {
            //printf("suboptimality=%f expands=%d g(searchgoal)=%d loop_time=%.3f time_elapsed=%.3f memoryCounter=%d\n", pSearchStateSpace->eps_satisfied, searchexpands - prevexpands, ((anaState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,double(clock()-loop_time)/CLOCKS_PER_SEC, double(clock() - TimeStarted)/CLOCKS_PER_SEC, MaxMemoryCounter);

            printf("suboptimality=%f g(searchgoal)=%d time_elapsed=%.3f memoryCounter=%d\n",
                   pSearchStateSpace->eps,
                   ((anaDoubleState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g, double(clock()
                       - TimeStarted) / CLOCKS_PER_SEC, MaxMemoryCounter);
			
			/*int sc;
 			pathIds = GetSearchPath(pSearchStateSpace, sc);
			vector<int> copy(pathIds);
			sols.push_back(copy);*/

            //printf("states expanded: %d\t states considered: %d\t time elapsed: %f\n",searchexpands - prevexpands, pSearchStateSpace->heap->currentsize, double(clock() - TimeStarted)/CLOCKS_PER_SEC);
        }
#elif ANAPLANNERDOUBLE_PERFORMANCE_TEST == 1
		PlannerStatsExtended temp;
		temp.eps = pSearchStateSpace->eps;
		temp.cost = ((anaDoubleState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
		temp.time = double(clock() - TimeStarted) / CLOCKS_PER_SEC;
		temp.memcount = MaxMemoryCounter;
		temp.inctime = double(clock()-loop_time) / CLOCKS_PER_SEC;
		temp.expands = searchexpands;
		planstats.push_back(temp);
#endif

#if DEBUG
        fprintf(fDeb, "eps=%f expands=%d g(searchgoal)=%d time=%.3f\n", pSearchStateSpace->eps, searchexpands - prevexpands,
            ((anaState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,double(clock()-loop_time)/CLOCKS_PER_SEC);
        PrintSearchState((anaState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData, fDeb);
#endif
        prevexpands = searchexpands;

        //if just the first solution then we are done
        if (bFirstSolution) break;

        //no solution exists
        if (((anaDoubleState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g == std::numeric_limits<unsigned int>::max()) break;
    }

#if DEBUG
    fflush(fDeb);
#endif

#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("Suboptimality = %f\n", pSearchStateSpace->eps);
#endif

    PathCost = ((anaDoubleState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
    MaxMemoryCounter += environment_->StateID2IndexMapping.size() * sizeof(int);

#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("MaxMemoryCounter = %d\n", MaxMemoryCounter);
#endif

    int solcost = std::numeric_limits<int>::max();
    bool ret = false;
    if (PathCost == std::numeric_limits<int>::max()) {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("could not find a solution\n");
#endif
        ret = false;
    }
    else {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("solution is found\n");
#elif ANAPLANNERDOUBLE_PERFORMANCE_TEST == 1
		FILE* perfFile;
		
		if(perfFile=fopen("anastardoubleperf.txt", "w")){
			for(int i=0;i<planstats.size();i++){
				PlannerStatsExtended temp = planstats.at(i);
				fprintf(perfFile, "%f %f %f %f %f %f\n", temp.time, temp.inctime, temp.cost, temp.eps, temp.expands, temp.memcount);
			}
			fclose(perfFile);
		}
#endif
        pathIds = GetSearchPath(pSearchStateSpace, solcost);
        ret = true;
    }

#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("total expands this call = %d, planning time = %.3f secs, solution cost=%d\n", searchexpands, (clock()
        - TimeStarted) / ((double)CLOCKS_PER_SEC), solcost);
#endif

    //fprintf(fStat, "%d %d\n", searchexpands, solcost);

    return ret;
}

//-----------------------------Interface function-----------------------------------------------------

//returns 1 if found a solution, and 0 otherwise
int anaPlannerDouble::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V)
{
    int solcost;

    return replan(allocated_time_secs, solution_stateIDs_V, &solcost);
}

//returns 1 if found a solution, and 0 otherwise
int anaPlannerDouble::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V, int* psolcost)
{
    vector<int> pathIds;
    bool bFound = false;
    int PathCost;
    //bool bFirstSolution = true;
    bool bFirstSolution = this->bsearchuntilfirstsolution;
    bool bOptimalSolution = false;
    *psolcost = 0;

#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("planner: replan called (bFirstSol=%d, bOptSol=%d)\n", bFirstSolution, bOptimalSolution);
#endif
    //plan
    if (!(bFound = Search(pSearchStateSpace_, pathIds, PathCost, bFirstSolution, bOptimalSolution,
                          allocated_time_secs)))
    {
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
        printf("failed to find a solution\n");
#endif
    }

    //copy the solution
    *solution_stateIDs_V = pathIds;
    *psolcost = PathCost;

    return (int)bFound;

}

int anaPlannerDouble::set_goal(int goal_stateID)
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("planner: setting goal to %d\n", goal_stateID);
#endif
    environment_->PrintState(goal_stateID, true, stdout);

    if (bforwardsearch) {
        if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            printf("ERROR: failed to set search goal state\n");
            return 0;
        }
    }
    else {
        if (SetSearchStartState(goal_stateID, pSearchStateSpace_) != 1) {
            printf("ERROR: failed to set search start state\n");
            return 0;
        }
    }

    return 1;
}

int anaPlannerDouble::set_start(int start_stateID)
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("planner: setting start to %d\n", start_stateID);
#endif
    environment_->PrintState(start_stateID, true, stdout);

    if (bforwardsearch) {

        if (SetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
            printf("ERROR: failed to set search start state\n");
            return 0;
        }
    }
    else {
        if (SetSearchGoalState(start_stateID, pSearchStateSpace_) != 1) {
            printf("ERROR: failed to set search goal state\n");
            return 0;
        }
    }

    return 1;
}

void anaPlannerDouble::costs_changed(StateChangeQuery const & stateChange)
{
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

void anaPlannerDouble::costs_changed()
{
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

int anaPlannerDouble::force_planning_from_scratch()
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("planner: forceplanfromscratch set\n");
#endif

    pSearchStateSpace_->bReinitializeSearchStateSpace = true;

    return 1;
}

int anaPlannerDouble::set_search_mode(bool bSearchUntilFirstSolution)
{
#if ANAPLANNERDOUBLE_PERFORMANCE_TEST == 0
    printf("planner: search mode set to %d\n", bSearchUntilFirstSolution);
#endif

    bsearchuntilfirstsolution = bSearchUntilFirstSolution;

    return 1;
}

void anaPlannerDouble::print_searchpath(FILE* fOut)
{
    PrintSearchPath(pSearchStateSpace_, fOut);
}

//---------------------------------------------------------------------------------------------------------
