import formula as mf
import frag_tree as ft
import netcdf
import networkx as nx 

##########################################helpful functions##########################################
def print_green_text(s):
    print(f"\033[1;32;40m{s}\033[0;37;40m")
def print_red_text(s):
    print(f"\033[1;31;40m{s}\033[0;37;40m")
def print_blue_text(s):
    print(f"\033[1;34;40m{s}\033[0;37;40m")

def run(exp, ans):
    try:
        assert(exp == ans)
    except AssertionError:
        print_red_text("Expected Output: {} \n Actual Output: {}".format(ans, exp))
        raise AssertionError

def map_mw(comp_lst, alpha, tol, exp_mass):
    masses = [mf.mw(i, alpha) for i in comp_lst]
    return sum([abs(mass - exp_mass) > tol for mass in masses])

##########################################individual tests##########################################

def test_mf_read(): #tests function which reads in molecular formulae and outputs a compomer
    run(mf.mf_to_compomer("CH4",["H","C","N","O"]), [4,1,0,0])
    run(mf.mf_to_compomer("C2H6NO2",["H","C","N","O"]), [6,2,1,2])
    run(mf.mf_to_compomer("H2SO4",["H","C","N","O","P","S"]), [2,0,0,4,0,1])
    run(mf.mf_to_compomer("C12H20",["H","C"]), [20,12])
    run(mf.mf_to_compomer("C2H6Cl3",["C","H","Cl"]), [2,6,3])

def test_enum_mf(): #tests generating function algorithm
    run(mf.enum_mf([12],24,use_heuristics=False), [[2]])
    run(len(mf.enum_mf([1,12,14,16],60,use_heuristics=False)),34)
    run(len(mf.enum_mf([1,12,14,16],168,use_heuristics=False)),424)
    run(len(mf.enum_mf([12,14,16,29,35],100,use_heuristics=False)),10)
    run(len(mf.enum_mf([1,12,16],28)),2)
    run(len(mf.enum_mf([1,12,14,16],60)),11)
    #tests if formulae generated have the right mw
    res = [round(mf.mw(i,["H","C","N","O"]),0) for i in mf.enum_mf([1,12,14,16],128)]
    for i in res:
        run(i, 128)
    res = [round(mf.mw(i,["H","C","N","O"]),0) for i in mf.enum_mf([1,12,14,16],60,use_heuristics=False)]
    for i in res:
        run(i, 60)


def test_hr_ms(): #tests HR-MS formula generator
    
    #tests if all formulas generated are within tol
    comp_lst = mf.HR_CF(["H","C","N","O"], 60.002, 0.1, as_vec=True)
    assert(map_mw(comp_lst,["H","C","N","O"], 60.002, 0.1) == 0)
    comp_lst = mf.HR_CF(["H","C","N","O","P","S"], 160.5, 0.1, as_vec=True)
    assert(map_mw(comp_lst,["H","C","N","O","P","S"], 160.5, 0.1) == 0)
    
    #test if formula in list generated
    assert("H4C" in mf.HR_CF(["H","C","N","O"], 16.03, 0.01))
    assert("H11C9NO2" in mf.HR_CF(["H","C","N","O"], 165.0722, 0.01))
    
    #test if parental constraint works
    assert(mf.HR_CF(["H","C","N","O"], 165.0722, 0.01,parent=[2,0,0,1]) == [])
    assert("H4C" in mf.HR_CF(["H","C","N","O"], 16.03, 0.01,parent=[6,2,0,0]))


def test_graph_analysis(): #evaluates the subformula graph algorithms
    #configuration info: vertex = 20ppm, edge = 10ppm
    MS_test = [65.03823, 66.0448, 77.03829, 78.04662, 79.05343, 91.05465, 94.041275, 105.033936, 106.041565, 122.03469, 123.04333, 135.04373, 137.05801, 138.07011, 168.07547]
    graphs = ft.analyse_MS(MS_test)
    run(len(graphs), 3) 
    run(ft.score_graph(graphs[0][0]),6.2) ##edge count!
    run(list(nx.nodes(graphs[0][0]))[-1], 'C9H12O3')

def test_ms_fetch(): #evaluates the netcdf module
    return 1 

##########################################running all tests##########################################


def test_all(tests):
    num_tests = len(tests)
    num_succeeded = 0
    num_failed = 0
    num_unimplemented = 0
    for test in tests:
        print_blue_text(f"Testing {test.__name__} ...")
        try:
            test()
            print_green_text("\tTest passed!")
            num_succeeded += 1
        except AssertionError:
            print_red_text("\tTest failed due to failed assertion!")
            num_failed += 1
        except NotImplementedError:
            print_red_text("\tTest failed due to unimplemented test/code!")
            num_unimplemented += 1

    print("-"*64)
    print_blue_text(f"Total number of tests: {num_tests}")
    print_green_text(f"Total number of succeeded tests: {num_succeeded}")
    print_green_text(f"Success rate: {num_succeeded / num_tests * 100}%")
    print_green_text(f"Implementation rate: {(num_succeeded + num_failed)/num_tests * 100}%")
    print_red_text(f"Total number of failed tests: {num_failed}")
    print_red_text(f"Total number of missing implementations: {num_unimplemented}")
    print_red_text(f"Failure rate: {num_failed / num_tests * 100}%")

if __name__ == "__main__":

    # put the name of the test in here:
    tests = [
        test_mf_read,
        test_enum_mf,
        test_hr_ms,
        test_graph_analysis
    ]

test_all(tests)