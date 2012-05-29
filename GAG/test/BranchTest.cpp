#include <GAGPL/GLYCAN/Branch.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>

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
	cout << "Monosaccharide Composition: " << endl;
	std::vector<InternalSite>& sites = mono.getInternalSites();

	for(std::vector<InternalSite>::iterator iter = sites.begin(); iter != sites.end(); iter++)
	{
		exploreInternalSite(*iter);
	}
}



int main ( int argc, char *argv[] )
{


	// Load Functional Group Table.
	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	// Load Monosaccharide unit Table.
	MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	mut.load();

	// Build up a Branch of only two units.
	Branch bc(0);
	Monosaccharide ms = mut.getMonosaccharideBySymbol("GlcA");
	// Unsaturated unit. This part should be realized by a class called
	// Reaction in the future.
	std::map<size_t, FunctionalGroup> substraction;
	substraction.insert(std::make_pair(4, fgt.getFunctionalGroupBySymbol("OH")));
	substraction.insert(std::make_pair(5, fgt.getFunctionalGroupBySymbol("H")));
	ms.remove(substraction);
	bc.addUnit(ms);
	
	Linkage lk0(0,1,4,"alpha");
	bc.addLinkage(lk0);
	bc.addUnit(mut.getMonosaccharideBySymbol("GlcN"));
	Linkage lk1(1,1,4,"belta");
	bc.addLinkage(lk1);
	bc.addUnit(mut.getMonosaccharideBySymbol("GlcA"));
	Linkage lk2(2,1,4,"alpha");
	bc.addLinkage(lk2);
	bc.addUnit(mut.getMonosaccharideBySymbol("GlcNAc"));

	// Add Modification.

	bc.addModification(1,2,"NH2","HSO4","H2O");
	bc.addModification(1,6,"OH","HSO4","H2O");
	bc.addModification(2,2,"OH","HSO4","H2O");
	bc.addModification(3,6,"OH","HSO4","H2O");
	// Recalculate the overall composition.
	bc.update();

	// Print results.
	for(size_t i = 0; i < bc.getGlycanChainUnits().size(); i++)
	{
		exploreMonosaccharide(bc.getGlycanChainUnits().at(i));
	}

	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */