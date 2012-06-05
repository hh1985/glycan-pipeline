#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <iostream>
#include <cmath>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace gag;

void exploreInternalSite(InternalSite& site)
{
	cout << "Site Composition: " << site.getCompositionString() << endl;
	std::multiset<FunctionalGroup>& fg = site.getFunctionalGroups();
	for(std::multiset<FunctionalGroup>::iterator iter = fg.begin(); iter != fg.end(); iter++)
	{
		cout << (*iter).getCompositionString() << endl;
	}
}

void exploreMonosaccharide(Monosaccharide& mono)
{
	cout << "Monosaccharide Composition: " << mono.getCompositionString() << endl;
	std::vector<InternalSite>& sites = mono.getInternalSites();

	for(std::vector<InternalSite>::iterator iter = sites.begin(); iter != sites.end(); iter++)
	{
		exploreInternalSite(*iter);
	}
}

void exploreBranch(Branch& bc)
{
	cout << "Branch: " << bc.getBranchID()  << endl;
	cout << "Composition: " << bc.getCompositionString()<< endl;
	std::vector<Monosaccharide>& mono_chain = bc.getGlycanChainUnits();

	for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin(); iter != mono_chain.end(); iter++)
	{
		exploreMonosaccharide(*iter);
	}
	
}
void printCompositionShift(const CompositionSigned& cs)
{
	cout << "Sign: " << cs.first << endl;
	cout << "Composition: " << cs.second.getCompositionString() << endl;
}



void processMassLoss(MassLossWindow& mlw, size_t cur, Composition& cp)
{	
	MassLoss& ml = mlw.at(cur);
	// Go over all number of mass loss.
	// For SO3, the upper limit is decided by the number of sulfate group.
	int max = abs(ml.upper);
	if(ml.loss_compo == Composition("SO3")) {
		std::map<string, int>::const_iterator iter = cp.get().find("S");
		if(iter != cp.get().end())
			max = iter->second;
		else
			max = 0;
	}
	for(int i = ml.lower; i <= max; i++)
	{
		Composition tmp_compo = cp;
		// e.g. 2H2O
		if(i != 0)
			cout << ml.loss_compo.getCompositionString() << ":" << i << " ";
		
		if(i < 0){
			for(int j = 0; j < abs(i); j++) {					
				tmp_compo.add(ml.loss_compo);
			}
		} else if(i > 0) {
			if((tmp_compo.getMass() - i * ml.loss_compo.getMass()) < 0.0)
				break;
			for(int j = 0; j < i; j++) {
				tmp_compo.deduct(ml.loss_compo);
			}
		}
		if(cur == mlw.size()-1) {
			// The end. Print out the composition and mass.
			cout << tmp_compo.getCompositionString() << " ";
			cout << setprecision(6) << tmp_compo.getMass() << "\n  ";
		
		} else {
			processMassLoss(mlw, cur+1,tmp_compo);
		} 

	}
	
	return;
}

void printFragment(Fragment& fg)
{
	FragmentationTable& ft = FragmentationTable::Instance();
	ft.load();

	MassLossWindow mlw = ft.getMassLoss();

	const CleavageCollection& cc = fg.getCleavages();
	
	for(CleavageCollection::const_iterator const_iter = cc.begin(); const_iter != cc.end(); const_iter++ )
	{
		for(std::vector<FragmentPosition>::const_iterator const_iter2 = const_iter->second.begin(); const_iter2 != const_iter->second.end(); const_iter2++)
		{
			cout << " " << const_iter->first << (*const_iter2).mono_id;
			if(const_iter->first == "X" || const_iter->first == "A")
				cout << ":" << const_iter2->xring_first << "-" << const_iter2->xring_second;
		}
	}
	cout << "\t" << fg.getCompositionString() << "\t" << setprecision(10) << fg.getMass() << endl;
	if(fg.getMass() > 200.0) {
		Composition tmp_compo = fg.getComposition();
		processMassLoss(mlw, 0, tmp_compo);	
	}

	cout << endl;
}



void addNRECleavage(Fragment& fg, int re_cl, size_t mono_id)
{
	for(int nre_cl = NRE_INTACT; nre_cl < NRE_N; nre_cl++)
	{
		// No cleavage.
		if(nre_cl == NRE_INTACT) {
			printFragment(fg);
			continue;
		}

		for(size_t k = 0; k < 5; k++)
		{
			//if(k > mono_id || (mono_id == k && re_cl != RE_INTACT))
			if(re_cl != RE_INTACT && (k > mono_id))
				continue;
			if(k == mono_id && (re_cl == A || nre_cl != X))
				continue;

			if(nre_cl == X) {
				for(size_t i = 0; i <= 3; i++)
				{
					for(size_t j = i+2; j<=5; j++)
					{
						Fragment fg1 = fg;
						FragmentPosition fp1 = {0,k,i,j};
						fg1.setFragmentation("X", fp1);
						printFragment(fg1);
					}
				}
			} else {
				Fragment fg2 = fg;
				FragmentPosition fp2 = {0,k,0,0};
				const char c = 'X' + (nre_cl-1);
				std::string str(&c,1);
				
				fg2.setFragmentation(str,fp2);
				printFragment(fg2);
			}
			// Print fragment.
			
		}
	}
}

void readSequence(std::vector<GlycanSequence>& gag_list)
{
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	// Load Monosaccharide unit Table.
	MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	mut.load();

	// Load Fragmentation parameters.
	FragmentationTable& ft = FragmentationTable::Instance();
	ft.load();

	ifstream seqfile;
	string ROW, flag;
	seqfile.open("../../test/input_test.txt");

	using namespace boost::xpressive;
	GlycanSequence* gs_pt(new GlycanSequence);
	Branch* bc_pt(new Branch);

	sregex res_ori = '>' >> (s1= +_w);

	smatch what;
	
	if(seqfile.is_open()) {
		while(!seqfile.eof()) {
			getline(seqfile, ROW);
			//Define the regex.
			//sregex rex = sregex::compile("^>(\w+)\s?(\d*)");


			if(regex_search(ROW, what, res_ori)) {
				flag = what[1];
				//cout << flag << endl;
				if(boost::lexical_cast<std::string>(what[1]) == "END"){
					gs_pt->update();
					gag_list.push_back(*gs_pt);
					delete gs_pt;
				} else if(boost::lexical_cast<std::string>(what[1]) == "SEQUENCE" && !gag_list.empty()){
					gs_pt = new GlycanSequence();				
				} else if(what[1] == "MODIFICATION") {
					// Add branch.
					gs_pt->addBranch(*bc_pt);
					delete bc_pt;				
				}
					
			} else {
				if(flag == "SEQUENCE") { // Dealing with branches.				
				//bc_pt.release();
				
				// Split the sequence into modules.
					sregex res = +_s;
					sregex_token_iterator cur(ROW.begin(), ROW.end(), res, -1), end;
					//Branch bc(*begin);

					size_t mono_id = 0;
					bc_pt = new Branch(boost::lexical_cast<size_t>(*cur));
					
					for(++cur; cur != end; cur++)
					{	
						std::string str = *cur;
						res = bos >> +_w >> eos;
						if(regex_match(str, what, res)) {
							Monosaccharide ms = mut.getMonosaccharideBySymbol(str);
							bc_pt->addUnit(ms);				
						} else {
							// If linkage.
							res = bos >> (s1= +_w) >> (s2= +_d) >> '-' >> (s3= +_d) >> eos;
							if(regex_match(str, what, res)){
								Linkage lk(mono_id, boost::lexical_cast<size_t>(what[2]), boost::lexical_cast<size_t>(what[3]), boost::lexical_cast<std::string>(what[1]));
								bc_pt->addLinkage(lk);
								mono_id++;
							}
						}
					}
				} else if(flag == "MODIFICATION"){
					
					
					sregex res = bos >> (s1= +_d) >> '\t' >> (s2= +_d) >> '\t' >> (s3= +_d) >> (s4= _w) >> (s5= +_w); 
					if(regex_match(ROW, what, res)){
						//bc_hs.addModification(0,2,"NH2","H2SO4","H2O");

						std::string mod;
						if(what[5] == "S")
							mod = "H2SO4";
						else if(what[5] == "Me")
							mod = "CH3OH";

						if(what[4] == "O"){
							gs_pt->getBranchByID(boost::lexical_cast<size_t>(what[1])).addModification(boost::lexical_cast<size_t>(what[2]),boost::lexical_cast<size_t>(what[3]), "OH", mod, "H2O");
						} else if(what[4] == "N") {
							gs_pt->getBranchByID(boost::lexical_cast<size_t>(what[1])).addModification(boost::lexical_cast<size_t>(what[2]),boost::lexical_cast<size_t>(what[3]), "NH2", mod, "H2O");
						}
					}

				} 
			
			} 
			
		}
	}
	seqfile.close();
}

void addCleavage(Fragment& fg)
{
	for(int re_cl = RE_INTACT; re_cl < RE_N; re_cl++)
	{
		if(re_cl == RE_INTACT)
			addNRECleavage(fg, re_cl, 0);
		else {
			for(size_t k = 0; k < 5; k++)
			{
				if(re_cl == A){
					for(size_t i = 0; i <= 3; i++)
					{
						for(size_t j = i+2; j<=5; j++)
						{	
							if(i == 0 && j == 5) // 0-5 cleavage is not allowed.
								continue;
							Fragment fg1 = fg;
							FragmentPosition fp1 = {0,k,i,j};
							// The internal cleavage doesn't work for A type.
							if(fg1.isInternalCleavage(fp1)) {
								continue;
							}
							fg1.setFragmentation("A", fp1);
							addNRECleavage(fg1, re_cl, k);
						}
					}
				} else {
					Fragment fg2 = fg;
					FragmentPosition fp2 = {0,k,0,0};
					const char c = 'A' + (re_cl-1);
					std::string str(&c,1);
					fg2.setFragmentation(str, fp2);
					addNRECleavage(fg2, re_cl, k);
				}
			}
		}		
	}
}

void Cleavages(GlycanSequence& seq)
{
	Fragment fg(seq);
	addCleavage(fg);
}


int main()
{
	//load();

	///* Test Fragmentation parameters. */
	//// X cleavage.
	//cout << "X cleavage" << endl;
	//FragmentationParams& fpx = ft.getFragmentationParams("X");
	//CompositionShift csx = ft.getCleavageShift("X");
	//CompositionShift csx_vec = ft.getCleavageShift("X","CID");
	//for(size_t i=0; i < csx.size(); i++)
	//	printCompositionShift(csx.at(i));
	//for(size_t i=0; i < csx_vec.size(); i++)
	//	printCompositionShift(csx_vec.at(i));

	//// Y cleavage.
	//cout << "Y cleavage" << endl;
	//FragmentationParams& fpy = ft.getFragmentationParams("Y");
	//CompositionShift csy = ft.getCleavageShift("Y");
	//CompositionShift csy_vec = ft.getCleavageShift("Y","CID");
	//for(size_t i=0; i < csy.size(); i++)
	//	printCompositionShift(csy.at(i));
	//for(size_t i=0; i < csy_vec.size(); i++)
	//	printCompositionShift(csy_vec.at(i));

	//// Z cleavage.
	//cout << "Z cleavage" << endl;
	//FragmentationParams& fpz = ft.getFragmentationParams("Z");
	//CompositionShift csz = ft.getCleavageShift("Z");
	//CompositionShift csz_vec = ft.getCleavageShift("Z","CID");
	//for(size_t i=0; i < csz.size(); i++)
	//	printCompositionShift(csz.at(i));
	//for(size_t i=0; i < csz_vec.size(); i++)
	//	printCompositionShift(csz_vec.at(i));

	//// A cleavage.
	//cout << "A cleavage" << endl;
	//FragmentationParams& fpa = ft.getFragmentationParams("A");
	//CompositionShift csa = ft.getCleavageShift("A");
	//CompositionShift csa_vec = ft.getCleavageShift("A","CID");
	//for(size_t i=0; i < csa.size(); i++)
	//	printCompositionShift(csa.at(i));
	//for(size_t i=0; i < csa_vec.size(); i++)
	//	printCompositionShift(csa_vec.at(i));

	//// B cleavage.
	//cout << "B cleavage" << endl;
	//FragmentationParams& fpb = ft.getFragmentationParams("B");
	//CompositionShift csb = ft.getCleavageShift("B");
	//CompositionShift csb_vec = ft.getCleavageShift("B","CID");
	//
	//for(size_t i=0; i < csb.size(); i++)
	//	printCompositionShift(csb.at(i));
	//for(size_t i=0; i < csb_vec.size(); i++)		
	//	printCompositionShift(csb_vec.at(i));

	//// C cleavage.
	//cout << "C cleavage" << endl;
	//FragmentationParams& fpc = ft.getFragmentationParams("C");
	//CompositionShift csc = ft.getCleavageShift("C");
	//CompositionShift csc_vec = ft.getCleavageShift("C","CID");
	//for(size_t i=0; i < csc.size(); i++)
	//	printCompositionShift(csc.at(i));
	//for(size_t i=0; i < csc_vec.size(); i++)
	//	printCompositionShift(csc_vec.at(i));

	///* Construct a glycan sequence */
	//GlycanSequence gs;

	//// Build up a Branch of only two units.
	//Branch bc1(0);
	//Monosaccharide ms1 = mut.getMonosaccharideBySymbol("Man");
	//// Unsaturated unit. This part should be realized by a class called
	//// Reaction in the future.
	//bc1.addUnit(ms1);	
	//Linkage lk1(0,1,6,"alpha");
	//bc1.addLinkage(lk1);
	//gs.addBranch(bc1);

	//Branch bc2(1);
	//Monosaccharide ms2 = mut.getMonosaccharideBySymbol("Man");
	//bc2.addUnit(ms2);
	//Linkage lk2(0,1,3,"alpha");
	//bc2.addLinkage(lk2);
	//gs.addBranch(bc2);

	//Branch bc3(2);
	//Monosaccharide ms3 = mut.getMonosaccharideBySymbol("Man");
	//bc3.addUnit(ms3);
	//Linkage lk3(0,1,6,"alpha");
	//bc3.addLinkage(lk3);
	//gs.addBranch(bc3);
	//// Be careful that the modification occurs on the unit on the glycan sequence, not on the single branch.
	//gs.addBranchLink(0,2);
	//gs.addBranchLink(1,2);
	//gs.updateChildrenIDs(2);

	//Branch bc4(3);
	//Monosaccharide ms4 = mut.getMonosaccharideBySymbol("GlcNAc");
	//bc4.addUnit(ms4);
	//Linkage lk4(0,1,4,"belta");
	//bc4.addLinkage(lk4);
	//gs.addBranch(bc4);

	//Branch bc5(4);
	//Monosaccharide ms5 = mut.getMonosaccharideBySymbol("Man");
	//bc5.addUnit(ms5);
	//Linkage lk5(0,1,3,"alpha");
	//bc5.addLinkage(lk5);
	//gs.addBranch(bc5);

	//Branch bc6(5);
	//Monosaccharide ms6 = mut.getMonosaccharideBySymbol("Man");
	//bc6.addUnit(ms6);
	//Linkage lk6(0,1,4,"belta");
	//bc6.addLinkage(lk6);
	//Monosaccharide ms7 = mut.getMonosaccharideBySymbol("GlcNAc");
	//bc6.addUnit(ms7);
	//Linkage lk7(1,1,4,"belta");
	//bc6.addLinkage(lk7);
	//Monosaccharide ms8 = mut.getMonosaccharideBySymbol("GlcNAc");
	//bc6.addUnit(ms8);
	//gs.addBranch(bc6);
	//gs.addBranchLink(2,5);
	//gs.addBranchLink(3,5);
	//gs.addBranchLink(4,5);
	//gs.updateChildrenIDs(5);
	////gs.update();

	//for(std::vector<Branch>::iterator iter = gs.getBranches().begin(); iter != gs.getBranches().end(); iter++)
	//{
	//	exploreBranch(*iter);
	//	cout << "Children: " << gs.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	//}

	//cout << "Cleavage: 2,4A(GlcNAc)/1,3X(Man)" << endl;
	//Fragment frag1(gs);
	//FragmentPosition fp11 = {5,2,2,4};
	//frag1.setFragmentation("A", fp11);
	//FragmentPosition fp12 = {2,0,1,3};
	//frag1.setFragmentation("X", fp12);
	//std::cout << "Composition: " << frag1.getCompositionString() << " " << "Mass: " << frag1.getMass() << std::endl;

	//cout << "Cleavage: Z 1,5X(Man)" << endl;
	//Fragment frag2(gs);
	//FragmentPosition fp21 = {3,0,0,0};
	//FragmentPosition fp22 = {0,0,1,5};
	//frag2.setFragmentation("Z", fp21);
	//frag2.setFragmentation("X", fp22);
	//cout << "Composition: " << frag2.getCompositionString() << " " << "Mass: " << frag2.getMass() << endl;
	//
	//cout << "Cleavage: 1,3X(GlcNAc) 1,4X(Man)" << endl;
	//Fragment frag3(gs);
	//FragmentPosition fp31 = {5,1,1,3};
	//FragmentPosition fp32 = {2,0,1,4};
	//frag3.setFragmentation("X", fp31);
	//frag3.setFragmentation("X", fp32);
	//cout << "Composition: " << frag3.getCompositionString() << " " << "Mass: " << frag3.getMass() << endl;
	//
	//cout << "Cleavage: 0,3X(GlcNAc) 1,4X(Man)" << endl;
	//Fragment frag4(gs);
	//FragmentPosition fp41 = {3,0,0,3};
	//FragmentPosition fp42 = {2,0,1,4};
	//frag4.setFragmentation("X", fp41);
	//frag4.setFragmentation("X", fp42);
	//cout << "Composition: " << frag4.getCompositionString() << " " << "Mass: " << frag4.getMass() << endl;	
	
	// GAG sequence
	//Branch bc_hs(0);
	//Monosaccharide ms_hs0 = mut.getMonosaccharideBySymbol("GlcN");
	//bc_hs.addUnit(ms_hs0);
	//Linkage lk_hs0(0,1,4,"alpha");
	//bc_hs.addLinkage(lk_hs0);
	//Monosaccharide ms_hs1 = mut.getMonosaccharideBySymbol("GlcA");
	//bc_hs.addUnit(ms_hs1);
	//Linkage lk_hs1(1,1,4,"belta");
	//bc_hs.addLinkage(lk_hs1);
	//Monosaccharide ms_hs2 = mut.getMonosaccharideBySymbol("GlcN");
	//bc_hs.addUnit(ms_hs2);
	//Linkage lk_hs2(2,1,4,"alpha");
	//bc_hs.addLinkage(lk_hs2);
	//Monosaccharide ms_hs3 = mut.getMonosaccharideBySymbol("GlcA");
	//bc_hs.addUnit(ms_hs3);
	//Linkage lk_hs3(3,1,4,"belta");
	//bc_hs.addLinkage(lk_hs3);
	//Monosaccharide ms_hs4 = mut.getMonosaccharideBySymbol("GlcN");
	//bc_hs.addUnit(ms_hs4);
	//// Add modification.
	//bc_hs.addModification(0,2,"NH2","H2SO4","H2O");
	//bc_hs.addModification(0,6,"OH","H2SO4","H2O");
	//bc_hs.addModification(2,2,"NH2","H2SO4","H2O");
	//bc_hs.addModification(2,3,"OH","H2SO4","H2O");
	//bc_hs.addModification(2,6,"OH","H2SO4","H2O");
	//bc_hs.addModification(3,2,"OH","H2SO4","H2O");
	//bc_hs.addModification(4,2,"NH2","H2SO4","H2O");
	//bc_hs.addModification(4,6,"OH","H2SO4","H2O");
	//bc_hs.addModification(4,1,"OH","CH3OH","H2O");

	//GlycanSequence gag;
	//gag.addBranch(bc_hs);
	// Read peak list file.
	std::vector<GlycanSequence> gag_list;
	readSequence(gag_list);
	for(std::vector<GlycanSequence>::iterator iter = gag_list.begin(); iter != gag_list.end(); iter++) {
		cout << "Get GAG composition: " << iter->getCompositionString() << endl;
		cout << "Get GAG mass: " << iter->getMass() << endl;

		//Fragment frag_hs(gag);
		Cleavages(*iter);	
	}

	// Match against the library.
	return EXIT_SUCCESS;
}