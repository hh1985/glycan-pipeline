/*
 * =====================================================================================
 *
 *       Filename:  EnvelopFinder.h
 *
 *    Description:  Class for generating list of envelops.
 *
 *        Version:  1.0
 *        Created:  9/6/2012 4:26:45 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu(hh1985@bu.edu), 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef GAG_ENVELOPFINDER_H
#define GAG_ENVELOPFINDER_H

#include <boost/shared_ptr.hpp>
#include <tuple>
#include <GAGPL/SPECTRUM/Envelop.h>
#include <GAGPL/SPECTRUM/EnvelopReference.h>
#include <GAGPL/MISC/Param.h>

namespace gag
{
	// Atom and coefficient.
	typedef std::map<std::string, double> AveragineFormulae;

	typedef std::vector<RichList> IslandList; 

	// Shift and corresponding peak information. Class ProtoEnv can be considered as prototype of class Envelop.
	typedef std::multimap<int, std::pair<RichPeakPtr, InfoPeakPtr> > ProtoEnv;

	// Overall fitting score and corresponding envelop set.
	typedef std::multimap<double, std::set<EnvelopPtr> > EnvelopSetCollection;

	class PeakSet
	{
	public:
		PeakSet(EnvelopBoundary& boundary)
			: bound(boundary) {}
		inline void addPeak(int shift, RichPeakPtr pk, InfoPeakPtr info)
		{
			peak_information.insert(std::make_pair(shift, std::make_pair(pk, info)));
		}
		inline EnvelopBoundary& getBoundary()
		{
			return bound;
		}
		inline ProtoEnv& getPeakInformation()
		{
			return peak_information;
		}
		inline size_t getSize()
		{
			return peak_information.size();
		}
	private:
		ProtoEnv peak_information;
		EnvelopBoundary bound;
	};

	using namespace param;

	class EnvelopFinder
	{
	private:
		Param& param;
		
		// TBD: This information should finally be retrieved from parameter file.
		// No sulfate.
		AveragineFormulae ave1;
		// High sulfate
		AveragineFormulae ave2;

		// Envelop reference table: from peak id to envelop (smart pointer).
		// Envelop will be managed only by EnvelopReference.
		EnvelopReference env_ref;
		
		// Envelop env_pool;
		std::vector<EnvelopPtr> env_pool;

		// Raw data.
		RichList& spectrum;

	private:
		// Interface for executing the process of finding envelops.
		void run();

		// Convert the raw peak list into several independent peak islands.
		IslandList getPeakIslands(RichList& spec);

		// Start from current candidate base peak, find matching envelop and estimate corresponding lambda values.
		void findNextEnvelops(RichList& island, RichPeakListByResolution::iterator iter_area, unsigned int charge, std::set<EnvelopPtr>& env_set);

		// Based on the lambda value stored in the envelop object, the method will detect if there is any events (peak split, peak overlapping, etc) and register the events in EnvelopReference.
		// Relative shifts and corresponding peak sets. This type is used to deal with peak split due to the varied sulfate contents.
		
		// Peakset and maximum number of envelops.
		std::pair<PeakSet, size_t> extendPeakSet(RichList& island, RichPeakListByMZ::iterator iter_mz, unsigned int charge);

		// Consider about the possibilities of using template for two direction expanding.
		size_t extendPeak(PeakSet& pk_set, RichPeakListByMZ& pks_mz, RichPeakListByMZ::iterator iter_mz, unsigned int charge, int direction);

		// Notice: the boundary has to be corrected for given charge.
		EnvelopBoundary createBoundary(const double mass, unsigned int charge);

		// If the PeakSet meet the basic threshold of the envelop.

		std::vector<EnvelopPtr> convertPeakSet(PeakSet& pk_set);

		// Estimate the global lambda value for each envelop. The corresponding estimated_abundance in each InfoPeak will also be updated 
		void updateEnvelopParameter(std::vector<EnvelopPtr> env_vec);
		void updateEnvelopParameter(EnvelopPtr env);

		// Once lambda is set, update the splitting score for envelops based on current base peak abundance. Calculate the dot product between theoretical envelop and the real one.
		void updateSplittingPotential(std::vector<EnvEntry> entry_vec);
		void updateSplittingPotential(EnvelopPtr env);
		void updateSplittingPotential(std::set<EnvelopPtr>& env_set);

		// Recursively optimize the fitting score by adjusting Lambda and Delta.
		// void optimizeEnvelopSet(RichList& pk_list);

		// Find the lambda value that give you the maximum fitting score.
		void adjustLambdaFromProfile(EnvelopPtr env);

		// Split the peak based on splitting potential.
		void splitPeak(RichPeakPtr pk);

		/* Latest methods */
		// Find all candidate envelops.
		void identifyEnvelops();
		// Remove redundant candidate envelops.
		void selectEnvelops();
		// Envelop set optimization using MS-Deconv.
		void optimizeEnvelopSet(std::set<EnvelopPtr>& env_set, double last_score = 0.0);

	public:
		EnvelopFinder(RichList& spec);
		
		inline std::vector<EnvelopPtr> getEnvelops()
		{ return env_pool; }

		double updateFittingScore(std::vector<EnvEntry> entry_vec);
		// If a theoretical peak is missing, count the exp peak as one with intensity 0.
		double updateFittingScore(EnvelopPtr env);
		double updateFittingScore(std::set<EnvelopPtr>& env_set);

		void printEnvelop(EnvelopPtr env);

	};

}


#endif /* GAG_ENVELOPFINDER_H */
