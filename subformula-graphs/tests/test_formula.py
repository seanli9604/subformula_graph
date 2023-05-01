import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from formula import Formula, FormulaGenerator

def test_formula_parsing():
    f = Formula.from_string("C4H6N2O")
    assert(str(f) == "C4H6N2O")

    f = Formula.from_string("CH6N2OClFSi")
    assert(str(f) == "CH6N2OFSiCl")

    blank = Formula.from_string("")
    assert(str(blank) == "")
    
    f = Formula.from_string("C4H9NOClSBr")
    assert(f.latex_formula() == "$\mathrm{C}_{4}\mathrm{H}_{9}\mathrm{N}\mathrm{O}\mathrm{S}\mathrm{Cl}\mathrm{Br}$")


def test_formula_arithmetic():
    f = Formula.from_string("C4H6N2O") + Formula.from_string("C4H6N2O")
    assert(str(f) == "C8H12N4O2")

    f = Formula.from_string("C4H6N2O") - Formula.from_string("C3H5N2O")
    assert(str(f) == "CH")

    f1 = Formula.from_string("C4H6N2O")
    f2 = Formula.from_string("CH") + Formula.from_string("C3H5N2O")
    assert(f1 == f2)

    f1 = Formula.from_string("C4H6N2O")
    f2 = Formula.from_string("C5H7N3O")
    assert(f1.is_subformula(f2))


def test_formula_restrictions():
    f = Formula.from_string("C19H4N3")
    assert(f.dbe() == 19.5)
    assert(not f.dbe().is_integer())


def test_mw_calcs():
    assert(Formula.from_string("").nom_mass == 0)
    assert(Formula.from_string("").exact_mass == 0)


    f = Formula.from_string("C11H16Cl2N3OPS3Si")
    assert(f.nom_mass == 431)
    assert(abs(f.exact_mass - 430.93394) < 0.00001)

    f = Formula.from_string("C5H7N3O")
    assert(f.dbe() == 4)
    assert(round(f.mass_deviation(125.058),2) == 7.29)

    f = Formula.from_string("C12H16Cl2N3OPSi4")
    assert(round(f.mass_deviation(430.933),2) == 35.99)


def test_elemental_restrictions():
    generator = FormulaGenerator(["C", "H", "N", "O"], 100)
    default = generator.get_default_restrictions()
    assert(default["C"] == (2, 8))
    assert(default["H"] == (0, 18))

    f = Formula.from_string("C5H7N3O")
    subformula = generator.get_parental_restrictions(f)
    assert(subformula["C"] == (0, 5))
    assert(subformula["H"] == (0, 7))
    assert(subformula["N"] == (0, 3))
    assert(subformula["O"] == (0, 1))

    gfs = generator.get_generating_functions()
    assert(gfs[0] == [24, 36, 48, 60, 72, 84, 96])
    assert(gfs[1] == [i for i in range(0, 19)])


def test_formula_enumeration():
    f = Formula.from_string("C5H7N3O")
    generator = FormulaGenerator(["C", "H", "N", "O"], 100)
    gfs = generator.get_generating_functions()
    nom_mass = [12, 1, 14, 16]
    formulae = generator.enumerate_compomers(gfs)
    assert(len(formulae) == 29)
    
    gf_subformula = generator.get_generating_functions(formula=f)
    formulae = generator.enumerate_compomers(gf_subformula)
    assert(str(formulae[0]) == "C3H6N3O")
    assert(len(formulae) == 1)


def test_formula_filtering():
    generator = FormulaGenerator(["C", "H", "N", "O"], 100)
    formulae = generator.get_formula_list(100.0613, 120)
    assert(len(formulae) == 4)
    assert(str(formulae[0]) == "C2H6N5")
    assert(all([f.mass_deviation(100.0613) < 120 for f in formulae]))
    formulae = generator.get_formula_list(100.0613, 90, DBE_restriction=lambda x: x.dbe() >= 0 and x.dbe().is_integer())
    assert(len(formulae) == 2)

    f = Formula.from_string("C5H7N3O")
    formulae = generator.get_formula_list(100.0613, 90, parent_formula=f)
    assert(len(formulae) == 0)
    formulae = generator.get_formula_list(100.0613, 110, parent_formula=f)
    assert(str(formulae[0]) == "C3H6N3O")

    generator = FormulaGenerator(["C","H", "N", "O", "P", "S","Cl"], 405)
    formulae = generator.get_formula_list(404.95549, 0.2)
    assert(str(formulae[0]) == "C11H16N3OP2SCl3")
    assert(len(formulae) == 11)

    generator = FormulaGenerator(["C","H", "N", "O","S"], 274)
    corrected_mass = 275.0484 - 1.00782503207 + 0.00054858
    formulae = generator.get_formula_list(corrected_mass, 3, DBE_restriction=lambda x: x.dbe() >= 0 and x.dbe().is_integer())
    assert(str(formulae[0]) == "C13H10N2O3S")
    assert(len(formulae) == 6)


def test_fragment_annotation():
    masses = [124.01557645207001, 126.01867645207, 138.01907645207, 154.01407645207001, 184.01207645207]
    f = Formula.from_string("C6H4N2O5")
    generator = FormulaGenerator(["C","H", "N", "O"], 184)
    annotation = generator.compute_all_fragment_formulae(masses, 4, f)
    assert(str(annotation) == "[(C6H4O3,), (C5H4NO3,), (C6H4NO3,), (C6H4NO4,), (C6H4N2O5,)]")





    


    
















