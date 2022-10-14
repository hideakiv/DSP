/*
 * StoModel.h
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */

#ifndef STOMODEL_H_
#define STOMODEL_H_

/** Coin */
#include "SmiScnModel.hpp"
#include "CoinTime.hpp"
/** Dsp */
#include "Utility/DspTypes.h"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"

#include <list>

class StoModel: public SmiScnModel {

	typedef std::map<int,int> ScenNodeMap;
	typedef std::map<int,int> NodeScenMap;

public:

    /** default constructor */
	StoModel();

	/** copy constructor */
	StoModel(const StoModel & rhs);

	/** copy constructor */
	StoModel(const SmiScnModel & rhs);

	/** default destructor */
	virtual ~StoModel();

    /** read SMPS files */
	DSP_RTN_CODE readSmps(const char * filename);

	/** read DRO file */
    // deal with dro later
	//DSP_RTN_CODE readDro(const char * filename);

	/** construct a map that maps variable names to their indices */
	bool mapVarnameIndex(map<string, int> &map_varName_index, const char * corefilename);

	/** read quadratic data file, extending the smps file */
    // deal with quad later
	//DSP_RTN_CODE readQuad(const char * smps, const char * filename);

	void __printData();

	/* print quadratic rows of scenario s, if s == -1, print quadratic rows in core */
    // deal with quad later
	//DSP_RTN_CODE printQuadRows (const int s);
	//DSP_RTN_CODE printQuadRows (const QuadRowData *qc);
public:
	/** get number of scenarios */
    // already exists
	//int getNumScenarios() const {return nscen_;}

	/** get number of stages */
    int getNumStages() {return getCore()->getNumStages();}

	/** get number of nodes */
	int getNumNodes() {return getSmiTree()->wholeTree().size();}

	/** get probability */
    const double ** getProbability() {return prob_;}

	/** get number of rows for a given stage */
	int getNumRows(int stage) {return getCore()->getNumRows(stage);}

	/** get number of columns for a given stage */
	int getNumCols(int stage) {return getCore()->getNumCols(stage);}

	/** get number of integer variables for a given stage */
	int getNumIntegers(int stage) {return (int)getCore()->getIntCols(stage).size();}

	/** get number of rows in core */
	int getNumCoreRows() {return getCore()->getNumRows();}
	/** get number of columns in core */
	int getNumCoreCols() {return getCore()->getNumCols();}
	/** get number of integer variables in core */
	int getNumCoreIntegers() {return getCore()->getBinaryLength() + getCore()->getIntegerLength();}

	// CoinPackedVector * getRowCore(int rowID);	TODO: do we need this?

	/** get row start index for a given stage */
	int getRowStart(int stage) {return getCore()->getRowStart(stage);}
	/** get column start index for a given stage */
	int getColStart(int stage) {return getCore()->getColStart(stage);}

	/** get number of quadratic constraints in core */
    // deal with quad later
	//int getNumCoreQRows() {return qc_row_core_->nqrows;}
	
	/** get number of quadratic constraints of a scenario*/
    // deal with quad later
	//int getNumScenQRows(int scen) {return qc_row_scen_[scen]->nqrows;}
	
	/** get objective function coefficients for a given stage */
	const double * getObjCore(int stage) {return getCore()->getDenseObjCoefficients(stage);}

	const CoinPackedVector * getObjScenario(int scenario, int stage);

	/** get quadratic objective function coefficients for a given stage */
    // deal with quad later
	//const CoinPackedMatrix * getQuadraticObjCore(int stage) {return qobj_core_[stage];}

	//const CoinPackedMatrix * getQuadraticObjScenario(int scenario) {return qobj_scen_[scenario];}

	/** get column type for a given stage */
	const char * getCtypeCore(int stage) {return ctype_core_[stage];}

	/** get initial solutions */
    // not used?
	const Solutions getInitialSolutions() {return init_solutions_;}

	/** get core coefficeints for a given row index */
	const CoinPackedVector * getRowCore(int i) {return rows_core_[i];}

	/** get parameters for quadratic constraints in core*/
    // deal with quad later
	//QuadRowData * getQuaraticsRowCore() const {return qc_row_core_;}

	/** get parameters for quadratic constraints in a scenario */
    // deal with quad later
	//QuadRowData * getQuaraticsRowScenario(int s) const {return qc_row_scen_[s];}

	//bool hasQuadraticRowCore() const {return qc_row_core_ != NULL ? true : false;};
	//bool hasQuadraticRowScenario() const {return qc_row_scen_ != NULL ? true : false;};
	//bool hasQuadraticRow() const {return (hasQuadraticRowCore() || hasQuadraticRowScenario());}

	/** set probability */
	void setProbability(double *probability);

	/** set initial solutions */
	void setSolution(
			int      size,    /**< size of array */
			double * solution /**< solution */);

	/** 
	 * Set the Wasserstein ambiguity set for distributionally robust optimization.
	 * This should be used for stochastic programming models, where the probabilities
	 * are used as the empirical references of the Wasserstein distance.
	 * 
	 * @param lp_norm use the distance norm of p
	 * @param eps maximum distance size
	 * @return DSP_RTN_OK if no error
	 */

    // deal with dro later
	//DSP_RTN_CODE setWassersteinAmbiguitySet(double lp_norm, double eps);

	/** 
	 * Nomalize probability vector
	 */
    // deal with dro later
	//void normalizeProbability();

#if 0
	/** add branching object */
	void addBranchingHyperplane(int nzcnt, int * indices, double * values, int priority);
#endif

public:

	/** split core matrix row for a given stage */
	CoinPackedVector * splitCoreRowVec(
			int i,  /**< row index */
			int stg /**< stage */);

	/** copy core column lower bounds */
	void copyCoreColLower(double * clbd, int stg);

	/** copy core column upper bounds */
	void copyCoreColUpper(double * cubd, int stg);

	/** copy core objective coefficients */
	void copyCoreObjective(double * obj, int stg);

	/** copy core quadratic objective coefficients 
	 *  for coupling terms, do not shift indices, start from x,
	 *  for non-coupling terms, shift indices, start from y,
	 *  only for two-stage problem
	*/
	// void copyCoreQuadraticObjective(
	// 	CoinPackedMatrix *&qobj_coupling,
	// 	CoinPackedMatrix *&qobj_ncoupling,
	// 	int stg);

	/** copy core column types */
	void copyCoreColType(char * ctype, int stg);

	/** copy core row lower bounds */
	void copyCoreRowLower(double * rlbd, int stg);

	/** copy core row upper bounds */
	void copyCoreRowUpper(double * rubd, int stg);

	/** combine random matrix row for a given scenario */
	void combineRandRowVec(
			CoinPackedVector * row, /**< core row vector */
			int i,                  /**< row index */
			int node                /**< node index */);

	/** combine random matrix row for a given stage and scenario */
	void combineRandRowVec(
			CoinPackedVector * row, /**< core row vector */
			int i,                  /**< row index */
			int stg,                /**< stage index */
			int node                /**< node index */);

	/** combine random column lower bounds */
	void combineRandColLower(double * clbd, SmiScnNode* node);

	/** combine random column upper bounds */
	void combineRandColUpper(double * cubd, SmiScnNode* node);

	/** combine random objective coefficients */
	void combineRandObjective(double * obj, SmiScnNode* node, bool adjustProbability = true);

	/** combine random quadratic objective coefficients */
	// void combineRandQuadraticObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg, int scen, bool adjustProbability = true);

	/** combine random row lower bounds */
	void combineRandRowLower(double * rlbd, SmiScnNode* node);

	/** combine random row upper bounds */
	void combineRandRowUpper(double * rubd, SmiScnNode* node);

	/** shift vector indices by offset */
	void shiftVecIndices(
			CoinPackedVector * vec, /**< vector to shift indices */
			int offset,             /**< offset by which indices are shifted */
			int start = 0           /**< index only after which indices are shifted */);

	/** shift vector indices by offset */
	void shiftVecIndices(
			int size,     /**< size of vecind */
			int * vecind, /**< vector indices to shift */
			int offset,   /**< offset by which indices are shifted */
			int start = 0 /**< index only after which indices are shifted */);

	// The following functions are for distributionally robust variant.
	// TODO: Better to create a new inhereted class?
	// virtual void setDro(bool yes) { isdro_ = yes; }
	// virtual bool isDro() {return isdro_;}
	// virtual int getNumReferences() {return nrefs_;}
	// virtual double getWassersteinSize() {return wass_eps_;}
	// virtual double getWassersteinDist(int i, int j);
	// virtual double getReferenceProbability(int i);

protected:

	// /*
	//  * Core data.
	//  */
	std::list<int> * binCols;		/**< indices of binary columns */
	std::list<int> * intCols;		/**< indices of integer columns */
	char **   ctype_core_;          /**< column types for each stage */
	CoinPackedVector ** rows_core_; /**< rows in core matrix */
	// QuadRowData * qc_row_core_;		/**< parameters for quadratic rows in core: current version only accept noncoupling quadratic rows */
	

	//StageData* stage_data_;
    //NodeData* node_data_;
	
	ScenNodeMap scen2node_; /** map from scenario to leafnode */

	Solutions init_solutions_; /**< initial solutions */

	bool fromSMPS_; /**< problem was read from SMPS files? */

	// bool isdro_;				/**< is this distributionally robust? */
	// int nrefs_;                 /**< number of reference scenarios for DRO */
	// double wass_eps_;           /**< size of the Wasserstein ball */
	// double ** wass_dist_;       /**< Wasserstein distances between two realizations */
	// double * refs_probability_; /** probability vector of references */

public:

#if 0
	vector<pair<CoinPackedVector*,int> > branchingHyperplanes_; /**< branching hyper-planes */
#endif

};



#endif /** STOMODEL_H_ */