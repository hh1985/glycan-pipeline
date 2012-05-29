#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
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

int main()
{


	// Load Functional Group Table.
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	// Load Monosaccharide unit Table.
	MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	mut.load();

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

	// Print results.
	for(std::vector<Branch>::iterator iter = gs.getBranches().begin(); iter != gs.getBranches().end(); iter++)
	{
		exploreBranch(*iter);
		cout << "Children: " << gs.getDescendantBranchIDs(iter->getBranchID()).first << "-" << gs.getDescendantBranchIDs(iter->getBranchID()).second << endl;
	}

	return EXIT_SUCCESS;
}