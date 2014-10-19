/*
 * =====================================================================================
 *
 *       Filename:  EnvelopFinder.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  9/6/2012 4:48:20 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (hh1985@bu.edu), 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <algorithm>
#include <GAGPL/MISC/Param.h>
#include <GAGPL/SPECTRUM/EnvelopFinder.h>
#include <GAGPL/SPECTRUM/EnvelopReference.h>
#include <GAGPL/SPECTRUM/IsotopicDistribution.h>
#include <boost/utility.hpp>
//#include <boost/math/distributions/normal.hpp>

namespace gag
{
	EnvelopFinder::EnvelopFinder(RichList& spec)
		: spectrum(spec), param(Param::instance())	
	{
		/* 
		Guess the range of average formula. Usually we assume the difference between composition estimated from the base peak and the one estimated from monoisotopic peak can be safely ignored.
			No Sulfur (100D)		High Sulfur (100D) 
			C: 3.7238523        1.9961864
			H: 5.4425534        2.9942796
			O: 2.8645018        3.1938983
			N: 0.2864502        0.1996186
			S: 0                0.5988559
		*/
		ave1.insert(std::make_pair("C", 3.7238523));
		ave1.insert(std::make_pair("H", 5.4425534));
		ave1.insert(std::make_pair("O", 2.8645018));
		ave1.insert(std::make_pair("N", 0.2864502));
		ave1.insert(std::make_pair("S", 0.0));

		ave2.insert(std::make_pair("C", 1.9961864));
		ave2.insert(std::make_pair("H", 2.9942796));
		ave2.insert(std::make_pair("O", 3.1938983));
		ave2.insert(std::make_pair("N", 0.1996186));
		ave2.insert(std::make_pair("S", 0.5988559));

		run();
	}

	void EnvelopFinder::run()
	{
		// Get/Set parameters.
		unsigned int precursor_charge = param.getParameter<unsigned int>("precursor_charge").first;

		// Divide the spectrum into several independent islands.
		IslandList is_list = this->getPeakIslands(spectrum);

		// Iterate over all independent islands.
		for(IslandList::iterator iter = is_list.begin(); iter != is_list.end(); iter++)
		{
			std::set<EnvelopPtr> temp_env_set;

			RichList island_pks = *iter;
			std::cout << "Found island: " << std::endl;
			island_pks.printPeakList<peak_mz>();
			// Starts from the peak with highest intensity. Assuming it is the base peak.
			RichPeakListByResolution& pks_area = island_pks.getPeakListByType<peak_resolution>();
			
			RichPeakListByResolution::iterator iter_area = pks_area.begin();
			RichPeakListByResolution::iterator end_iter = pks_area.end();
			// The last peak cannot be considered as a base peak.
			end_iter--;

			/* Step 1. Envelop identification */
			for(; iter_area != end_iter; iter_area++)
			{
				// Iterate over all possible charge states.
				// Notice that the sign of the charge is not used here.
				for(unsigned int z=precursor_charge; z>=1; z--)
				{
					// All candidate envelops should be stored for records.
					this->findNextEnvelops(island_pks, iter_area, z, temp_env_set);
				}
				
			}	

			std::cout << "Before optimizing, output all the envelop information..." << std::endl;
			BOOST_FOREACH(EnvelopPtr env, temp_env_set)
				env_ref.printEnvelopInformation(env);

			std::cout << std::endl;

			// Adjust the fitting scores and splitting scores for envelop set.
			this->optimizeEnvelopSet(temp_env_set);

		} 
	}

	IslandList EnvelopFinder::getPeakIslands(RichList& pks)
	{
		double island_space = param.getParameter<double>("island_space").first;

		IslandList is_list;

		RichList* temp_spec = new RichList();

		RichPeakListByMZ& pks_mz = pks.getPeakListByType<peak_mz>();

		for(RichPeakListByMZ::const_iterator iter = pks_mz.begin(); iter != pks_mz.end(); iter++)
		{
			if(iter == pks_mz.begin()) {
				temp_spec->addPeak(*iter);
				continue;
			}

			if((*iter)->mz - (*boost::prior(iter))->mz < island_space) {
				temp_spec->addPeak(*iter);
			} else { // 
				is_list.push_back(*temp_spec);
				
				// Delete the container. Be careful of this section.
				delete temp_spec;

				temp_spec = new RichList();
				temp_spec->addPeak(*iter);
			}
		}
		is_list.push_back(*temp_spec);
		delete temp_spec;
		
		return is_list;
	}

	void EnvelopFinder::findNextEnvelops(RichList& island, RichPeakListByResolution::iterator iter_area, unsigned int charge, std::set<EnvelopPtr>& env_set)
	{
		// 1. Find candidate envelops which extend from current peak.
		
		// Project the area iterator into corresponding mz iterator.
		RichPeakListByMZ::iterator iter_mz = island.getPeakContainer().project<peak_mz>(iter_area);

		std::cout << "Identify base peak: " << (*iter_mz)->mz << " Charge: " << charge << std::endl;
		
		int mode = param.getParameter<int>("mode").first;
		double pre_mz = param.getParameter<double>("precursor_mz").first;
		int pre_z = param.getParameter<int>("precursor_charge").first;
		
		// Convert m/z value to corresponding mass. 
		double pre_mass = calculateMass(pre_mz, pre_z * mode);

		// If the calculated mass is larger than precursor mass, there is no need to continue identification.
		double mass = calculateMass((*iter_mz)->mz, charge * mode);
		
		if(mass > pre_mass)
			return;

		std::pair<PeakSet,size_t> pk_set = this->extendPeakSet(island, iter_mz, charge);
		
		if(pk_set.first.getSize()==1)
			return;

    // Container for temporary storage of envelops.
		std::vector<EnvelopPtr> env_store;
		for(size_t i=0; i < pk_set.second; i++) {
			EnvelopPtr env = createEnvelop(charge, pk_set.first.getBoundary());
			env_store.push_back(env);
		}

		// 2. TBD: Estimate the validity of the pk_set.
		// e.g. the overall shifts and the continuity of the shifts.

		// 2.a If the base peak has been reported as a distal peak of other envelops.
		// 2.b If envelops which share the base peak have the same charge 
		// 2.c No event report for the base peak.
		// 2.d no contribution of new peaks.
		// 2.e fitting score is worse than expected.

		// Get the maximum number of envelops.

		// Dynamically convert pk_set into vector of envelops.
		ProtoEnv& peak_information = pk_set.first.getPeakInformation();

		ProtoEnv::iterator pk_iter = peak_information.begin(); 
		// Iterate over all shift, decide if all envelop will share.
		for(; pk_iter != peak_information.end(); )
		{
			int count = peak_information.count(pk_iter->first);
			std::pair<ProtoEnv::iterator, ProtoEnv::iterator> shift_pair = peak_information.equal_range(pk_iter->first);
			ProtoEnv::iterator shift_iter = shift_pair.first;
			ProtoEnv::iterator last_iter = shift_pair.second;
			last_iter--;

			for(size_t i=0; i<env_store.size(); i++)
			{	
				// Clean the definition of the internal structure.
				int shift = shift_iter->first;
				RichPeakPtr pk = shift_iter->second.first;			
				InfoPeakPtr info = shift_iter->second.second;
				// Register the (peak, env) pair, shift and peak info into the envelop reference table.
				InfoPeakPtr info_copy(new InfoPeak(*info));
				std::cout << "Register peak " << pk->mz << "for envelop " << env_store[i]->id << std::endl;
				env_ref.addDictionaryReference(pk, env_store.at(i), shift, info_copy);

				// Non_unique_key.
				if(shift_iter != last_iter)
					shift_iter++;		

			}

			pk_iter = shift_pair.second;
		}

		// Update envelop lambda information for each envelop.
		this->updateEnvelopParameter(env_store);

		// Store the envelops.
		env_set.insert(env_store.begin(), env_store.end());
	}
	
std::pair<PeakSet, size_t> EnvelopFinder::extendPeakSet(RichList& island, RichPeakListByMZ::iterator iter_mz, unsigned int charge)
	{
		
		// Candidate base peak. Notice that this might not be the monoisotopic peak for molecules with high mass values.
		RichPeakPtr expr_base_pk = *iter_mz;
		int mode = param.getParameter<int>("mode").first;
		double mass = calculateMass(expr_base_pk->mz, charge * mode);
		EnvelopBoundary bound = this->createBoundary(mass, charge);

		// The reason of creating an object called envelop is because during the extension of peakset towards different directions, there might be multiple matching for given shift. 
		PeakSet pk_set(bound);

		InfoPeakPtr pk_infor = boost::make_shared<InfoPeak>(expr_base_pk->resolution, 0.0, 0.5, 0.0);
		
		pk_set.addPeak(0, expr_base_pk, pk_infor);

		RichPeakListByMZ& pks_mz = island.getPeakListByType<peak_mz>();

		// Shift towards right along the mz axis.
		size_t max_env = this->extendPeak(pk_set, pks_mz, iter_mz, charge, 1);
		// TBD: Shift towards left along the mz axis.
		//this->extendPeak(pk_set, pks_mz, iter_mz, charge, -1);

		return std::make_pair(pk_set,max_env);
	}


	EnvelopBoundary EnvelopFinder::createBoundary(const double mass, unsigned int charge)
	{
		// 1. Estimate the rounded version of compositions.
		double coef = mass / 100.0;
		Composition compo1, compo2;

		AveragineFormulae::iterator iter1 = ave1.begin();
		for(; iter1 != ave1.end(); iter1++)
			compo1.addElement(iter1->first, (int)floor(iter1->second * coef + 0.5));

		AveragineFormulae::iterator iter2 = ave2.begin();
		for(; iter2 != ave2.end(); iter2++)
			compo2.addElement(iter2->first, (int)floor(iter2->second * coef + 0.5));
		
		int mode = param.getParameter<int>("mode").first;

		int z = mode > 0 ? charge : (-1)*charge;
		// 2. Get the range of the isotopic distribution
		IsotopicDistribution iso_up(compo1);
		AggregatedIsotopicVariants peakset1 = iso_up.getAggregatedIsotopicVariants(z);
		IsotopicDistribution iso_down(compo2);
		AggregatedIsotopicVariants peakset2 = iso_down.getAggregatedIsotopicVariants(z);

		return EnvelopBoundary(peakset1, peakset2);
	}

	size_t EnvelopFinder::extendPeak(PeakSet& pk_set, RichPeakListByMZ& pks_mz,RichPeakListByMZ::iterator iter_mz, unsigned int charge, int direction)
	{
		//int dist = direction > 0 ? std::distance(iter_mz, pks_mz.end())-1 : std::distance(pks_mz.begin(), iter_mz);
		
		unsigned int shift = 1;
		//int prev_dist = 0;

		// Experimental candidate base peak.
		RichPeakPtr expr_base_pk = *iter_mz;

		int mode = param.getParameter<int>("mode").first;
		// For confidence estimation, using internal_accuracy.
		double internal_accuracy = param.getParameter<double>("internal_accuracy").first;

		// For peak finding and lambda estimation, using external_accuracy.
		double external_accuracy = param.getParameter<double>("external_accuracy").first;

		// Calculate the mass of the candidate base peak.
		// Notice that there will be mechanism where there is electron capture.
		double mass = calculateMass(expr_base_pk->mz, charge * mode);

		// Create theoretical boundary. The boundary has been adjusted to fit the charge.
		//EnvelopBoundary env_bound = this->createBoundary(mass, charge);
		EnvelopBoundary& env_bound = pk_set.getBoundary();

		//std::cout << "Up Boundary: " << std::endl;
		//env_bound.first.printPeakList<peak_mz>();
		//std::cout << "Down Boundary: " << std::endl;
		//env_bound.second.printPeakList<peak_mz>();

		// Get a copy of the base peak iterator for further move operation.
		//RichPeakListByMZ::iterator iter = iter_mz;
		int stop_flag = 1;

		// Dynamically convert pk_set into vector of envelops.
		std::pair<unsigned int, size_t> shift_count(0,1);
		size_t max_env = 0;

		//std::advance(iter_mz, direction);
		if(direction < 0 && iter_mz == pks_mz.begin())
			return 0;
		std::advance(iter_mz, direction);
		if(direction > 0 && iter_mz == pks_mz.end())
			return 0;

		// TBD: the program should efficiently location the position of most likely peak. instead of sequentially search for it.
		while(1)
		{


			std::cout << "Shift: " << shift << std::endl;

			// Get the boundary peaks at current shift.
			PeakPtr theo_pk1 = env_bound.getUpperBound().getPeakByShift<peak_intensity>(direction * shift);
			PeakPtr theo_pk2 = env_bound.getLowerBound().getPeakByShift<peak_intensity>(direction * shift);

			// No extension any more.
			if(theo_pk1->mz == 0.0 && theo_pk2->mz == 0.0)
				break;

			RichPeakPtr current_pk = *iter_mz;
			std::cout << "Examine peak: " << current_pk->mz << std::endl;

			// Experimental distance from current peak to base peak.
			// Notice that the value might be negative.
			double expr_dist = abs(current_pk->mz - expr_base_pk->mz);

			// Theoretical distance from current peak to base peak. Adjusted for charged peak.
			double theo_dist1 = abs(env_bound.getUpperBound().getMassDifferenceByShift<peak_intensity>(0, shift));
			double theo_dist2 = abs(env_bound.getLowerBound().getMassDifferenceByShift<peak_intensity>(0, shift));

			double max_dist = theo_dist1;
			double min_dist = theo_dist2;

			if(max_dist < min_dist) {
				swap(min_dist, max_dist);
			}

			double error1 = 1e6 * (expr_dist - min_dist)/expr_base_pk->mz;
			double error2 = 1e6 * (expr_dist - max_dist)/expr_base_pk->mz;
			
			if(error1 < -1 * external_accuracy) {
				// Before the shift region: keep moving the iterator.
				// std::advance(iter_mz, direction); continue;
			} else if(error2 > external_accuracy) { 
				// Beyond the shift region: update shift.
				if(stop_flag == 0) {
					shift++; stop_flag = 1; continue;
				} else {
					// break means no missing peak is allowed in the middle. This might be controlled by some parameter.
					stop_flag = 1;
					break;
				}
			} else {
				// A matching peak. Notice the effect of charge state.
				std::cout << "A matching peak!" << std::endl;
				
				if(shift == shift_count.first) {
					shift_count.second++;
				} else {
					// Reset shift_count.
					shift_count = std::make_pair(shift, 1);
				}

				// Update maximum number of envelops.
				if(shift_count.second > max_env)
					max_env = shift_count.second;

				stop_flag = 0;
				// Estimate the confidence of the mz in terms of estimating lambda.
				double diff = max_dist - min_dist;

				// Error window is used for confidence estimation.
				double err_win = internal_accuracy * expr_base_pk->mz * 1e-6;
				double prob = diff/(diff + 2*err_win);

				// Indicate the sign of the error window.
				int err_sign = (theo_dist1 < theo_dist2 ? -1 : 1);

				//: lambda_mz should always be in (0, 1). If lambda is closed to 1, the pattern is closed to no-sulfate boundary, and 0 for high-sulfate boundary.
				double lambda_mz = (expr_dist - theo_dist2+err_sign*err_win)/(theo_dist1 - theo_dist2 + 2*err_sign*err_win);

				// Round the abnormal lambda_mz to the boundary.
				if(lambda_mz > 1) {
					lambda_mz = 1.0;
					prob = 1.0;
				} else if(lambda_mz < 0) {
					lambda_mz = 0.0;
					prob = 1.0;
				}

				// TBD: estimate lambda_abd (within 5%) and delta_lambda.
				// Notice: lambda_abd can be < 0 or > 1.
				double lambda_abd = (current_pk->resolution/expr_base_pk->resolution - theo_pk2->intensity)/(theo_pk1->intensity - theo_pk2->intensity);

				InfoPeakPtr pk_infor = boost::make_shared<InfoPeak>(current_pk->resolution, lambda_mz, lambda_abd, prob);
				// Notice: The abundance information should be initialized during the establish of global lambda value.
				//pk_infor.adjusted_abundance = env_bound.getTheoreticalPeak(lambda_abd, shift).intensity;
				
				pk_set.addPeak(shift, current_pk, pk_infor);
				
			}
			if(direction < 0 && iter_mz == pks_mz.begin())
				break;
			std::advance(iter_mz, direction);
			if(direction > 0 && iter_mz == pks_mz.end())
				break;
		}
		return max_env;
	}

	void EnvelopFinder::updateEnvelopParameter(EnvelopPtr env)
	{
		EnvEntry base_entry = env_ref.getEntryByShift(env, 0);

		// 1.a Find lambda value for each of the peak.
		std::vector<EnvEntry> peak_entries = env_ref.getEntryByEnvelop(env);

		double current_prob = 0.0;
		BOOST_FOREACH(EnvEntry& entry, peak_entries)
		{
			if(entry.info->prob > current_prob) {
				current_prob = entry.info->prob;
				env->lambda_mz = entry.info->lambda_mz;
			}

		}
		//// 1.b Normalize the confidence value.
		//double sum_confidence = 0.0;
		//for(std::vector<EnvEntry>::iterator iter = peak_entries.begin(); 
		//	iter != peak_entries.end(); iter++)
		//	sum_confidence += pow(iter->info->prob,2);

		//double lambda_mz = 0.0;
		//for(std::vector<EnvEntry>::iterator iter = peak_entries.begin(); 
		//	iter != peak_entries.end(); iter++) {
		//		// Pass base peak.
		//		if(iter->getShift() == 0) 
		//			continue;

		//		lambda_mz += pow(iter->info->prob,2)/sum_confidence *  iter->info->lambda_mz;
		//}
		//env->lambda_mz = lambda_mz;

		// 2. Update theoretical abundance from the lambda.
		for(std::vector<EnvEntry>::iterator iter = peak_entries.begin(); 
			iter != peak_entries.end(); iter++) {
			// The adjusted abundance is initialized as the exp data.
			iter->info->adjusted_abundance = iter->pk->area;
		}
	}

	void EnvelopFinder::updateEnvelopParameter( std::vector<EnvelopPtr> env_list )
	{
		BOOST_FOREACH(EnvelopPtr env, env_list) {
			this->updateEnvelopParameter(env);

			//this->updateSplittingPotential(env);
		}
	}

	void EnvelopFinder::optimizeEnvelopSet(std::set<EnvelopPtr>& env_set, double last_score)
	{
		/* New version */
		this->updateSplittingPotential(env_set);
		this->updateFittingScore(env_set);

		// Get all the peaks. 
		RichList pk_list = env_ref.getRichPeakList(env_set);

		RichPeakListByMZ& pks_mz = pk_list.getPeakListByType<peak_mz>();

		for(RichPeakListByMZ::iterator iter = pks_mz.begin(); iter != pks_mz.end(); iter++)
		{
			// Split the peak based on fitting score & estimated abundance.
			// Notice that the identification of FP will also happen here.
			this->splitPeak(*iter);
		}
		// Calculate the overall fitting scores.
		double next_score = this->updateFittingScore(env_set);
		std::cout << "Old score: " << last_score << "\t" << "New score: " << next_score << std::endl;
		if(abs(next_score - last_score)/last_score < 1e-5)
			return;
		else
			this->optimizeEnvelopSet(env_set, next_score);
		/* Old version */

		//std::set<RichPeak> processed_peak_set;
		//std::set<RichPeak>::iterator peak_iter;

		//// Generate fitting score for each envelop.
		//// Estimate scaling factor for each envelop.
		//this->updateSplittingPotential(env_pool);
		//this->updateFittingScore(env_pool);

		//sort(env_pool.begin(), env_pool.end(), Envelop::scoreLarger);

		//// Iterate over all envelops with intensity from high to low.
		//BOOST_FOREACH(EnvelopPtr& env, env_pool)
		//{
		//	// If the base peak has been processed.
		//	EnvEntry base_entry = env_ref.getBasePeakForEnvelop(env);
		//	
		//	peak_iter = processed_peak_set.find(base_entry->pk);
		//	if(peak_iter != processed_peak_set.end())
		//		continue;

		//	// Split the base peak.
		//	bool status = this->splitPeak(base_entry->pk);
		//	if(status) {
		//		this->updateSplittingPotential(env);
		//		this->updateFittingScore(env);
		//	}

		//	// Get all belonging peaks.
		//	std::vector<EnvEntry> entry_vec = env_ref.getEntryByEnvelop(env);

		//	// Iterate over all peaks with intensity from high to low.
		//	BOOST_FOREACH(EnvEntry& entry, entry_vec)
		//	{
		//		if(entry->getShift()==0)
		//			continue;
		//		// Split the peak based on 1. splitting potential and 2. scaling factor.
		//		this->splitPeak(entry->pk);
		//		// Update the records.
		//		processed_peak_set.insert(entry->pk);
		//	}
			
		//}		
	}

	void EnvelopFinder::updateSplittingPotential(EnvelopPtr env)
	{
		
		//double scaling_factor = env->scaling_factor;

		// Entry for base peak.
		EnvEntry& base_entry = env_ref.getEntryByShift(env, 0);
		std::cout << "The base peak is: " << base_entry.pk->mz << std::endl;

		// Get all peaks included by the envelop.
		std::vector<EnvEntry> entry_vec = env_ref.getEntryByEnvelop(env);

		// Update the scaling factor. The value for base peak should be 1.
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			double expected_abundance = env_ref.getTheoreticalAbundance(entry);
			// Adjusted abunance / Theoretical abundance.
			double scaling_factor = entry.info->adjusted_abundance/expected_abundance;
			if(scaling_factor < env->scaling_factor)
				env->scaling_factor = scaling_factor;
		}

		// Update the adjusted abundance.
		// 1. Update base entry.
		base_entry.info->adjusted_abundance *= env->scaling_factor;
		// 2. Update all the rest of the peaks.
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			if(entry.getShift()==0)
				continue;

			entry.info->adjusted_abundance = env_ref.getTheoreticalAbundance(entry);
		}
		


		//// Iterate over all theoretical peaks (except base peak), calculate their relative shift.
		//double total_share = 0.0;
		//std::map<int, PeakPtr> pk_map = env->getTheoreticalPeaks();
		//for(std::map<int, PeakPtr>::iterator iter = pk_map.begin(); 
		//	iter != pk_map.end(); iter++) 
		//{
		//	total_share += iter->second->intensity;
		//}
		//	
		//double base_adjusted_abd = 0.0;
		//for(std::map<int, PeakPtr>::iterator iter = pk_map.begin(); 
		//	iter != pk_map.end(); iter++)
		//{
		//	int shift = iter->first; PeakPtr pk = iter->second;

		//	// Get entry corresponding to the shift. Notice the problem if shift does not exist for the envelop..
		//	EnvEntry entry = env_ref.getEntryByShift(env, shift);
	
		//	// If EnvEntry is empty. Which means the experimental peak corresponds to the shift doesn't show up.
		//	if(entry.isEmpty()) continue;

		//	entry.info->relative_share = pk->intensity / total_share;

		//	double expected_abundance = env_ref.getTheoreticalAbundance(entry);

		//	if(expected_abundance > entry.info->adjusted_abundance) {
		//		// Potential overlapping event for base peak. Scale the base peak.
		//		base_adjusted_abd += base_entry.info->adjusted_abundance * (entry.info->adjusted_abundance / expected_abundance) * entry.info->relative_share;
		//	}
		//	// Update adjusted abundance.
		//	entry.info->adjusted_abundance = expected_abundance;
		//}

		//// Update base peak information.
		//base_entry.info->adjusted_abundance = base_adjusted_abd;
	}

	void EnvelopFinder::updateSplittingPotential( std::vector<EnvEntry> entry_vec )
	{ 
		BOOST_FOREACH(EnvEntry entry, entry_vec)
			this->updateSplittingPotential(entry.env);
	}
	void EnvelopFinder::updateSplittingPotential(std::set<EnvelopPtr>& env_set)
	{
		BOOST_FOREACH(EnvelopPtr env, env_set)
			this->updateSplittingPotential(env);
	}

	double EnvelopFinder::updateFittingScore(EnvelopPtr env)
	{
		std::vector<EnvEntry> entry_vec = env_ref.getEntryByEnvelop(env);
		//RichPeakPtr base_pk = env_ref.getBasePeakForEnvelop(env);
		EnvEntry base_entry = env_ref.getEntryByShift(env, 0);
		
		double fitting_score = 0.0;
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			// For base peak, intensity_shift and mz_shift is always 1.
			if(entry.getShift()==0) {			
				fitting_score += sqrt(base_entry.info->adjusted_abundance);
				continue;
			}

			double expected_abundance = env_ref.getTheoreticalAbundance(entry);
			// The fitting score comes from 3 parts: 
			// 1. the shift of intensity: abs(E-T)/T.
			double intensity_shift = 1.0 - abs(expected_abundance - entry.info->adjusted_abundance)/expected_abundance;
			if(intensity_shift < 0)
				intensity_shift = 0.0;
			//std::cout << "Intensity shift: " << intensity_shift << std::endl;
				
			// 2. the relative shift of mz.
			double internal_accuracy = param.getParameter<double>("internal_accuracy").first;
			
			double mz_shift = 1.0 - abs(entry.pk->mz - base_entry.pk->mz - entry.env->getTheoreticalMZDistance(0, entry.getShift()))/(entry.env->getWindowSizeByShift(entry.getShift()));
			//std::cout << "MZ shift: " << mz_shift << std::endl;

			// 3. the square root of expected_abundance.
			double weight = sqrt(expected_abundance);
			//std::cout << "Weight: " << weight << std::endl;
			fitting_score += weight * intensity_shift * mz_shift;
		}
		//std::cout << "Fitting score: " << fitting_score << std::endl;
		env->fitting_score = fitting_score;
		return fitting_score;
	} 

	double EnvelopFinder::updateFittingScore( std::vector<EnvEntry> entry_vec)
	{
		double sum_score = 0.0;
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
			sum_score += this->updateFittingScore(entry.env);
		return sum_score;
	}

	double EnvelopFinder::updateFittingScore(std::set<EnvelopPtr>& env_set)
	{
		double sum_score = 0.0;
		BOOST_FOREACH(EnvelopPtr env, env_set)
			sum_score += this->updateFittingScore(env);
		return sum_score;
	}

	void EnvelopFinder::adjustLambdaFromProfile( EnvelopPtr env )
	{
		double last_score = this->updateFittingScore(env);
		// The value of lambda_bound has to be in between 0.0 and 1.0
		double lambda_bound = 0.0;

		while(1)
		{
			// Adjust lambda. Notice the direction.
			env->lambda_mz = (lambda_bound + env->lambda_mz)/2;

			// Update splitting score.
			this->updateSplittingPotential(env);

			double new_score = this->updateFittingScore(env);

			double score_diff = new_score - last_score;
			
			// Converge.
			if(abs(score_diff) < 1e-6)
				break;
			
			if(score_diff > 0.0) {
					env->lambda_mz = (lambda_bound + env->lambda_mz)/2;
			} else {
					// Change the bound.
					lambda_bound = 1.0-lambda_bound;
					env->lambda_mz = (lambda_bound + env->lambda_mz)/2;
			}

		}

	}

	void EnvelopFinder::splitPeak(RichPeakPtr pk)
	{
		// Get all connecting envelops.
		std::vector<EnvEntry> entry_vec = env_ref.getEntryByPeak(pk);

		if(entry_vec.size() == 1) // No need to split for single envelop.
			return;
		else if(entry_vec.size() == 0)
			throw std::runtime_error("Undefined peak in the dictionary");
		else
			;

		double splitting_score = 0.0;
		
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			// splitting_score += entry.info->adjusted_abundance;
			splitting_score += sqrt(entry.env->fitting_score * entry.info->adjusted_abundance);
		}
			
		// Split the peak.
		BOOST_FOREACH(EnvEntry& entry, entry_vec)
		{
			double scale_coef = sqrt(entry.env->fitting_score * entry.info->adjusted_abundance)/ splitting_score;

			// Update the estimated abundance.
			entry.info->adjusted_abundance = pk->resolution * scale_coef;

			env_ref.printEnvelopInformation(entry.env);
		}

	}

	void EnvelopFinder::printEnvelop( EnvelopPtr env )
	{
		std::cout << "Envelop ID: " << env->id << std::endl;
		std::cout << "Charge: " << env->charge_state << std::endl;
		std::cout << "Lambda: " << env->lambda_mz << std::endl;
		std::cout << "Score: " << env->fitting_score << std::endl;
		std::cout << std::endl;
		// Peak information.
		std::cout << "Covered peaks:" << std::endl;
		std::cout << "Shift\tMZ\tEXP_ABD\tFIT_ABD" << std::endl;

		// Get all peaks.
		std::vector<EnvEntry> pk_entries = env_ref.getEntryByEnvelop(env);
		BOOST_FOREACH(EnvEntry& pk_entry, pk_entries)
		{
			std::cout << pk_entry.getShift() << "\t" << pk_entry.pk->mz << "\t"
				<< pk_entry.pk->resolution << "\t" << pk_entry.info->adjusted_abundance
				<< std::endl;
		}
		std::cout << std::endl;
	}

	//void EnvelopFinder::optimizeEnvelopSet(RichList& pk_list)
	//{
	//	// 1. Iterate over peaks with m/z from low to high
	//	RichPeakListByMZ& pks_by_mz = pk_list.getPeakListSortedByType<peak_mz>();

	//	RichPeakListByMZ::iterator iter_mz = pks_by_mz.begin();

	//	// Envelop set and corresponding fitting score.
	//	EnvelopSetCollection last_frontier;

	//	// TBD: how to decide the cluster group will end.
	//	for(; iter_mz != pks_by_mz.end(); iter_mz++)
	//	{
	//		// 2. For each peak, find all connecting envelops.
	//		std::vector<EnvEntry> entry_vec = env_ref.getEntryByPeak(*iter_mz);

	//		// Empty information. Consider storing the information with a separate container.
	//		if(entry_vec.size() == 0)
	//			continue;

	//		std::vector<std::set<EnvelopPtr> > env_set_collection;
	//		for(size_t i = 0; i < entry_vec.size()-1; i++)
	//			for(size_t j=i; j<entry_vec.size(); j++)
	//			{
	//				std::set<EnvelopPtr> temp_env;
	//				temp_env.insert(entry_vec[i]->env);

	//				if(i != j) // Not just pairs.
	//					temp_env.insert(entry_vec[j]->env);

	//				env_set_collection.insert(temp_env);
	//			}

	//		// New frontier. 
	//		if(last_frontier.size() == 0)
	//		{
	//			std::vector<std::set<EnvelopPtr> >::iterator iter = env_set_collection.begin();
	//			for(; iter != env_set_collection.end(); iter++)
	//			{
	//				double temp_score = this->calculateFittingScore(*iter);
	//				last_frontier.insert(std::make_pair(temp_score, *iter));
	//			}
	//			continue;
	//		}

	//		EnvelopSetCollection next_frontier;
	//		// Update current frontier.
	//		for(EnvelopSetCollection::iterator fron_iter = last_frontier.begin(); fron_iter != last_frontier.end(); fron_iter++)
	//		{
	//			BOOST_FOREACH(std::set<EnvelopPtr>& env_set, env_set_collection)
	//			{
	//				std::set<EnvelopPtr> temp_set(fron_iter->second);
	//				// Append the set.
	//				set_union(temp_set.begin(), temp_set.end(), 
	//					env_set.begin(), env_set.end(),
	//					inserter(temp_set, temp_set.begin()));

	//				// Compare the fitting score between std::vector<EnvelopPtr> with and without new current envelops.
	//				double score1 = fron_iter->first;
	//				double score2 = this->calculateFittingScore(new_set);

	//				if(score2 > score1)
	//					next_frontier.insert(std::make_pair(score2, temp_set));
	//				else
	//					next_frontier.insert(std::make_pair(score1, fron_iter->second));
	//			}
	//		}

	//		last_frontier = next_frontier;
	//	
	//	}

	//	// Select the set with maximum score. 
	//	std::set<EnvelopPtr>& winner_set = last_frontier.rbegin()->second;
	//	
	//	env_pool.insert(env_pool.end(), winner_set.begin(), winner_set.end());
	//	
	//}

	//double EnvelopFinder::calculateFittingScore( std::set<EnvelopPtr> env_set )
	//{
	//	/* Estimate the abundance based on current set selection. */
	//	


	//	double total_score = 0.0;
	//	BOOST_FOREACH(EnvelopPtr& env, env_set)
	//	{
	//		total_score += this->calculateFittingScore(env);
	//	}
	//	return total_score;
	//}

	//double EnvelopFinder::calculateFittingScore(EnvelopPtr env)
	//{
	//	std::vector<EnvEntry> entry_vec = env_ref.getEntryByEnvelop(env);
	//	
	//	if(entry_vec.size()==0)
	//		return 0.0;

	//	// Get base entry.
	//	EnvEntry base_entry = env_ref.getEntryByShift(env, 0);

	//	double fitting_score = 0.0;
	//	BOOST_FOREACH(EnvEntry& entry, entry_vec)
	//	{
	//		// For base peak, intensity_shift and mz_shift is always 1.
	//		if(entry.getShift()==0) {
	//			fitting_score += sqrt(base_entry.info->adjusted_abundance);
	//			continue;
	//		}

	//		double expected_abundance = env_ref.getTheoreticalAbundance(entry);
	//		
	//		// The fitting score comes from 3 parts: 
	//		// 1. the shift of intensity: abs(E-T)/T.
	//		double intensity_shift = 1.0 - abs(expected_abundance - entry.info->adjusted_abundance)/expected_abundance;
	//		if(intensity_shift < 0)
	//			intensity_shift = 0.0;

	//		// 2. the relative shift of mz. Here the accuracy depends on the position of the peak in the envelop.
	//		double mz_shift = 1.0 - abs(entry.pk->mz - base_entry.pk->mz - entry.env->getTheoreticalMZDistance(0, entry.getShift()))/(entry.env->getWindowSizeByShift(entry.getShift()));

	//		// 3. the square root of expected_abundance.
	//		double weight = sqrt(expected_abundance);

	//		fitting_score += weight * intensity_shift * mz_shift;
	//	}

	//	env->fitting_score = fitting_score;
	//	return fitting_score;
	//}
}