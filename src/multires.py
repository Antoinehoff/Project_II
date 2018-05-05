from common import ProblemType
from solver import Solver
from gui import Gui
import images
import connection_table as ct
import numpy


def multires(nelx, nely, params, bc):
    # Allocate design variables for the first level
    x = None
    x_comp = None

    # Dynamic parameters
    downsampling = 2**(params.numLevels - 1)
    params.exemplarDownsampling *= downsampling

    #To record Compliance and Appearance evolution
    aggregated_hist={}

    # Multires synthesis
    for level in range(params.numLevels):
        print("*** Level " + str(level))
        if x is not None:
            # Upsample previous solution
            x = images.upsample(x, nelx, nely)
            if params.problemType == ProblemType.AppearanceWithMaxCompliance or params.problemType==ProblemType.AppearanceWithMaxComplianceAndSymmetry:
                x_comp = images.upsample(x_comp, nelx, nely)
        gui = None
        if params.hasGui:
            gui = Gui(nelx, nely)
        if params.problemType==ProblemType.AppearanceWithMaxCompliance:
            params.complianceMax = 0
            solver = Solver(nelx, nely, params, ProblemType.Compliance, bc, gui)
            x_comp = solver.optimize(x_comp)
            min_compliance = solver.last_optimum_value()
            params.complianceMax = min_compliance * params.complianceMaxFactor
            print("")

        ###Added by Antoine Hoffmann EPFL 2018
        if params.problemType == ProblemType.AppearanceWithMaxComplianceAndSymmetry \
        or params.problemType == ProblemType.ComplianceWithSymmetry:
            print("-> Constructing connection table...")
            a_array = numpy.array(params.a)
            c_array = numpy.array(params.c)
            scale_matrix = numpy.array([[nelx,0],[0,nely]])
            c_array = [numpy.dot(scale_matrix,c_array[i]) for i in range(len(c_array))]
            a_array = [numpy.array(a/numpy.sqrt(numpy.dot(a,a))) for a in a_array]
            print("a_array : " + str(a_array))
            print("c_array : " + str(c_array))
            connection_table = ct.construct_connection_table(a_array,c_array,nelx,nely)
            mapping_vector = ct.construct_mapping_vector(connection_table)
            if params.problemType==ProblemType.AppearanceWithMaxComplianceAndSymmetry:
                params.complianceMax = 0
                solver = Solver(nelx, nely, params, ProblemType.ComplianceWithSymmetry,\
                                bc, gui, mapping_vector)
                x_comp = solver.optimize(x_comp)
                min_compliance = solver.last_optimum_value()
                params.complianceMax = min_compliance * params.complianceMaxFactor
                print("")
            # Solve problem
            solver = Solver(nelx, nely, params, params.problemType,\
                            bc, gui, mapping_vector)
        ###
        else:
            solver = Solver(nelx, nely, params, params.problemType, bc, gui)

        x = solver.optimize(x, enforce_constraints=(level > 0))
        if params.hasGui:
            solver.filtering.filter_variables(x, solver.x_phys)
            gui.update(solver.x_phys)

        ### Antoine Hoffmann 2018 EPFL
        #Store histories of the compliance and appearance evolution
        if params.record_histories==True:
            aggregated_hist['level'+str(level)]=solver.get_histories()
        ###

        # Go to next level
        if level < params.numLevels - 1:
            nelx *= 2
            nely *= 2
            params.exemplarDownsampling /= 2.0
            params.maxSolverStep //= 2
            params.lengthSquare /= 2.0
        print("")

    # Filter last result to obtain physical variables
    solver.filtering.filter_variables(x, solver.x_phys)
    results = {
        "last_optimum": solver.last_optimum_value(),
        "volume": sum(solver.x_phys) / len(solver.x_phys)}
    if params.problemType == ProblemType.AppearanceWithMaxCompliance  or params.problemType==ProblemType.AppearanceWithMaxComplianceAndSymmetry:
        results["compliance_factor"] = solver.compliance_max / min_compliance
    return (solver.x_phys, nelx, nely, results, aggregated_hist)
