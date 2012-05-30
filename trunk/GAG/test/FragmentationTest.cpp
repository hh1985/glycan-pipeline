#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <iostream>

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

int main()
{
// Load Functional Group Table.
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	// Load Monosaccharide unit Table.
	MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	mut.load();

	// Load Fragmentation parameters.
	FragmentationTable& ft = FragmentationTable::Instance();
	ft.load();

	/* Test Fragmentation parameters. */
	// X cleavage.
	cout << "X cleavage" << endl;
	FragmentationParams& fpx = ft.getFragmentationParams("X");
	CompositionShift csx = ft.getCleavageShift("X");
	CompositionShift csx_vec = ft.getCleavageShift("X","CID");
	for(size_t i=0; i < csx.size(); i++)
		printCompositionShift(csx.at(i));
	for(size_t i=0; i < csx_vec.size(); i++)
		printCompositionShift(csx_vec.at(i));

	// Y cleavage.
	cout << "Y cleavage" << endl;
	FragmentationParams& fpy = ft.getFragmentationParams("Y");
	CompositionShift csy = ft.getCleavageShift("Y");
	CompositionShift csy_vec = ft.getCleavageShift("Y","CID");
	for(size_t i=0; i < csy.size(); i++)
		printCompositionShift(csy.at(i));
	for(size_t i=0; i < csy_vec.size(); i++)
		printCompositionShift(csy_vec.at(i));

	// Z cleavage.
	cout << "Z cleavage" << endl;
	FragmentationParams& fpz = ft.getFragmentationParams("Z");
	CompositionShift csz = ft.getCleavageShift("Z");
	CompositionShift csz_vec = ft.getCleavageShift("Z","CID");
	for(size_t i=0; i < csz.size(); i++)
		printCompositionShift(csz.at(i));
	for(size_t i=0; i < csz_vec.size(); i++)
		printCompositionShift(csz_vec.at(i));

	// A cleavage.
	cout << "A cleavage" << endl;
	FragmentationParams& fpa = ft.getFragmentationParams("A");
	CompositionShift csa = ft.getCleavageShift("A");
	CompositionShift csa_vec = ft.getCleavageShift("A","CID");
	for(size_t i=0; i < csa.size(); i++)
		printCompositionShift(csa.at(i));
	for(size_t i=0; i < csa_vec.size(); i++)
		printCompositionShift(csa_vec.at(i));

	// B cleavage.
	cout << "B cleavage" << endl;
	FragmentationParams& fpb = ft.getFragmentationParams("B");
	CompositionShift csb = ft.getCleavageShift("B");
	CompositionShift csb_vec = ft.getCleavageShift("B","CID");
	
	for(size_t i=0; i < csb.size(); i++)
		printCompositionShift(csb.at(i));
	for(size_t i=0; i < csb_vec.size(); i++)		
		printCompositionShift(csb_vec.at(i));

	// C cleavage.
	cout << "C cleavage" << endl;
	FragmentationParams& fpc = ft.getFragmentationParams("C");
	CompositionShift csc = ft.getCleavageShift("C");
	CompositionShift csc_vec = ft.getCleavageShift("C","CID");
	for(size_t i=0; i < csc.size(); i++)
		printCompositionShift(csc.at(i));
	for(size_t i=0; i < csc_vec.size(); i++)
		printCompositionShift(csc_vec.at(i));

	/* Construct a glycan sequence */
	GlycanSequence gs;

	// Build up a Branch of only two units.
	Branch bc1(0);
	Monosaccharide ms1 = mut.getMonosaccharideBySymbol("Man");
	// Unsaturated unit. This part should be realized by a class called
	// Reaction in the future.
	bc1.addUnit(ms1);	
	Linkage lk1(0,1,6,"alpha");
	bc1.addLinkage(lk1);
	gs.addBranch(bc1);

	Branch bc2(1);
	Monosaccharide ms2 = mut.getMonosaccharideBySymbol("Man");
	bc2.addUnit(ms2);
	Linkage lk2(0,1,3,"alpha");
	bc2.addLinkage(lk2);
	gs.addBranch(bc2);

	Branch bc3(2);
	Monosaccharide ms3 = mut.getMonosaccharideBySymbol("Man");
	bc3.addUnit(ms3);
	Linkage lk3(0,1,6,"alpha");
	bc3.addLinkage(lk3);
	gs.addBranch(bc3);
	// Be careful that the modification occurs on the unit on the glycan sequence, not on the single branch.
	gs.addBranchLink(0,2);
	gs.addBranchLink(1,2);
	gs.updateChildrenIDs(2);

	Branch bc4(3);
	Monosaccharide ms4 = mut.getMonosaccharideBySymbol("GlcNAc");
	bc4.addUnit(ms4);
	Linkage lk4(0,1,4,"belta");
	bc4.addLinkage(lk4);
	gs.addBranch(bc4);

	Branch bc5(4);
	Monosaccharide ms5 = mut.getMonosaccharideBySymbol("Man");
	bc5.addUnit(ms5);
	Linkage lk5(0,1,3,"alpha");
	bc5.addLinkage(lk5);
	gs.addBranch(bc5);

	Branch bc6(5);
	Monosaccharide ms6 = mut.getMonosaccharideBySymbol("Man");
	bc6.addUnit(ms6);
	Linkage lk6(0,1,4,"belta");
	bc6.addLinkage(lk6);
	Monosaccharide ms7 = mut.getMonosaccharideBySymbol("GlcNAc");
	bc6.addUnit(ms7);
	Linkage lk7(1,1,4,"belta");
	bc6.addLinkage(lk7);
	Monosaccharide ms8 = mut.getMonosaccharideBySymbol("GlcNAc");
	bc6.addUnit(ms8);
	gs.addBranch(bc6);
	gs.addBranchLink(2,5);
	gs.addBranchLink(3,5);
	gs.addBranchLink(4,5);
	gs.updateChildrenIDs(5);
	//gs.update();

	for(std::vector<Branch>::iterator iter = gs.getBranches().begin(); iter != gs.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
		cout << "Children: " << gs.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	}

	cout << "Cleavage: 2,4A(GlcNAc)/1,3X(Man)" << endl;
	Fragment frag1(gs);
	FragmentPosition fp11 = {5,2,2,4};
	frag1.setFragmentation("A", fp11);
	FragmentPosition fp12 = {2,0,1,3};
	frag1.setFragmentation("X", fp12);
	std::cout << "Composition: " << frag1.getCompositionString() << " " << "Mass: " << frag1.getMass() << std::endl;

	cout << "Cleavage: Z 1,5X(Man)" << endl;
	Fragment frag2(gs);
	FragmentPosition fp21 = {3,0,0,5};
	FragmentPosition fp22 = {0,0,1,5};
	frag2.setFragmentation("Z", fp21);
	frag2.setFragmentation("X", fp22);
	cout << "Composition: " << frag2.getCompositionString() << " " << "Mass: " << frag2.getMass() << endl;
	
	cout << "Cleavage: 1,3X(GlcNAc) 1,4X(Man)" << endl;
	Fragment frag3(gs);
	FragmentPosition fp31 = {5,1,1,3};
	FragmentPosition fp32 = {2,0,1,4};
	frag3.setFragmentation("X", fp31);
	frag3.setFragmentation("X", fp32);
	cout << "Composition: " << frag3.getCompositionString() << " " << "Mass: " << frag3.getMass() << endl;
	
	cout << "Cleavage: 0,3X(GlcNAc) 1,4X(Man)" << endl;
	Fragment frag4(gs);
	FragmentPosition fp41 = {3,0,0,3};
	FragmentPosition fp42 = {2,0,1,4};
	frag3.setFragmentation("X", fp41);
	frag3.setFragmentation("X", fp42);
	cout << "Composition: " << frag4.getCompositionString() << " " << "Mass: " << frag4.getMass() << endl;	
	
	
	
	
	return EXIT_SUCCESS;
}