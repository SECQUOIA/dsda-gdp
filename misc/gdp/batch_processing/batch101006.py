from pyomo.environ import *
from pyomo.gdp import *
import os
import random

def build_model():
    model = AbstractModel()

    # TODO: it looks like they set a bigM for each j. Which I need to look up how to do...
    model.BigM = Suffix(direction=Suffix.LOCAL)
    model.BigM[None] = 1000


    ## Constants from GAMS
    StorageTankSizeFactor = 2*5 # btw, I know 2*5 is 10... I don't know why it's written this way in GAMS?
    StorageTankSizeFactorByProd = 3
    MinFlow = -log(10000)
    VolumeLB = log(300)
    VolumeUB = log(3500)
    StorageTankSizeLB = log(100)
    StorageTankSizeUB = log(15000)
    UnitsInPhaseUB = log(6)
    UnitsOutOfPhaseUB = log(6)
    # TODO: YOU ARE HERE. YOU HAVEN'T ACTUALLY MADE THESE THE BOUNDS YET, NOR HAVE YOU FIGURED OUT WHOSE
    # BOUNDS THEY ARE. AND THERE ARE MORE IN GAMS.


    ##########
    # Sets
    ##########

    model.PRODUCTS = Set()
    model.STAGES = Set(ordered=True)
    model.PARALLELUNITS = Set(ordered=True)

    # TODO: this seems like an over-complicated way to accomplish this task...
    def filter_out_last(model, j):
        return j != model.STAGES.last()
    model.STAGESExceptLast = Set(initialize=model.STAGES, filter=filter_out_last)


    # TODO: these aren't in the formulation??
    #model.STORAGE_TANKS = Set()


    ###############
    # Parameters
    ###############

    model.HorizonTime = Param()
    model.Alpha1 = Param()
    model.Alpha2 = Param()
    model.Beta1 = Param()
    model.Beta2 = Param()

    model.ProductionAmount = Param(model.PRODUCTS)
    model.ProductSizeFactor = Param(model.PRODUCTS, model.STAGES)
    model.ProcessingTime = Param(model.PRODUCTS, model.STAGES)

    # These are hard-coded in the GAMS file, hence the defaults
    model.StorageTankSizeFactor = Param(model.STAGES, default=StorageTankSizeFactor)
    model.StorageTankSizeFactorByProd = Param(model.PRODUCTS, model.STAGES,
                                              default=StorageTankSizeFactorByProd)

    # TODO: bonmin wasn't happy and I think it might have something to do with this?
    # or maybe issues with convexity or a lack thereof... I don't know yet.
    # I made PRODUCTS ordered so I could do this... Is that bad? And it does index
    # from 1, right?
    def get_log_coeffs(model, k):
        return log(model.PARALLELUNITS.ord(k))

    model.LogCoeffs = Param(model.PARALLELUNITS, initialize=get_log_coeffs)

    # bounds
    model.volumeLB = Param(model.STAGES, default=VolumeLB)
    model.volumeUB = Param(model.STAGES, default=VolumeUB)
    model.storageTankSizeLB = Param(model.STAGES, default=StorageTankSizeLB)
    model.storageTankSizeUB = Param(model.STAGES, default=StorageTankSizeUB)
    model.unitsInPhaseUB = Param(model.STAGES, default=UnitsInPhaseUB)
    model.unitsOutOfPhaseUB = Param(model.STAGES, default=UnitsOutOfPhaseUB)


    ################
    # Variables
    ################

    # TODO: right now these match the formulation. There are more in GAMS...

    # unit size of stage j
    # model.volume = Var(model.STAGES)
    # # TODO: GAMS has a batch size indexed just by products that isn't in the formulation... I'm going
    # # to try to avoid it for the moment...
    # # batch size of product i at stage j
    # model.batchSize = Var(model.PRODUCTS, model.STAGES)
    # # TODO: this is different in GAMS... They index by stages too?
    # # cycle time of product i divided by batch size of product i
    # model.cycleTime = Var(model.PRODUCTS)
    # # number of units in parallel out-of-phase (or in phase) at stage j
    # model.unitsOutOfPhase = Var(model.STAGES)
    # model.unitsInPhase = Var(model.STAGES)
    # # TODO: what are we going to do as a boundary condition here? For that last stage?
    # # size of intermediate storage tank between stage j and j+1
    # model.storageTankSize = Var(model.STAGES)

    # variables for convexified problem
    # TODO: I am beginning to think these are my only variables actually.
    # GAMS never un-logs them, I don't think. And I think the GAMs ones
    # must be the log ones.
    def get_volume_bounds(model, j):
        return (model.volumeLB[j], model.volumeUB[j])
    model.volume_log = Var(model.STAGES, bounds=get_volume_bounds)
    model.batchSize_log = Var(model.PRODUCTS, model.STAGES)
    model.cycleTime_log = Var(model.PRODUCTS)
    def get_unitsOutOfPhase_bounds(model, j):
        return (0, model.unitsOutOfPhaseUB[j])
    model.unitsOutOfPhase_log = Var(model.STAGES, bounds=get_unitsOutOfPhase_bounds)
    def get_unitsInPhase_bounds(model, j):
        return (0, model.unitsInPhaseUB[j])
    model.unitsInPhase_log = Var(model.STAGES, bounds=get_unitsInPhase_bounds)
    def get_storageTankSize_bounds(model, j):
        return (model.storageTankSizeLB[j], model.storageTankSizeUB[j])
    # TODO: these bounds make it infeasible...
    model.storageTankSize_log = Var(model.STAGES, bounds=get_storageTankSize_bounds)

    # binary variables for deciding number of parallel units in and out of phase
    model.outOfPhase = Var(model.STAGES, model.PARALLELUNITS, within=Binary)
    model.inPhase = Var(model.STAGES, model.PARALLELUNITS, within=Binary)

    ###############
    # Objective
    ###############

    def get_cost_rule(model):
        return model.Alpha1 * sum(exp(model.unitsInPhase_log[j] + model.unitsOutOfPhase_log[j] + \
                                              model.Beta1 * model.volume_log[j]) for j in model.STAGES) +\
            model.Alpha2 * sum(exp(model.Beta2 * model.storageTankSize_log[j]) for j in model.STAGESExceptLast)
    model.min_cost = Objective(rule=get_cost_rule)


    ##############
    # Constraints
    ##############

    def processing_capacity_rule(model, j, i):
        return model.volume_log[j] >= log(model.ProductSizeFactor[i, j]) + model.batchSize_log[i, j] - \
            model.unitsInPhase_log[j]
    model.processing_capacity = Constraint(model.STAGES, model.PRODUCTS, rule=processing_capacity_rule)

    def processing_time_rule(model, j, i):
        return model.cycleTime_log[i] >= log(model.ProcessingTime[i, j]) - model.batchSize_log[i, j] - \
            model.unitsOutOfPhase_log[j]
    model.processing_time = Constraint(model.STAGES, model.PRODUCTS, rule=processing_time_rule)

    def finish_in_time_rule(model):
        return model.HorizonTime >= sum(model.ProductionAmount[i]*exp(model.cycleTime_log[i]) \
                                        for i in model.PRODUCTS)
    model.finish_in_time = Constraint(rule=finish_in_time_rule)


    ###############
    # Disjunctions
    ###############

    def storage_tank_selection_disjunct_rule(disjunct, selectStorageTank, j):
        model = disjunct.model()
        def volume_stage_j_rule(disjunct, i):
            return model.storageTankSize_log[j] >= log(model.StorageTankSizeFactor[j]) + \
                model.batchSize_log[i, j]
        def volume_stage_jPlus1_rule(disjunct, i):
            return model.storageTankSize_log[j] >= log(model.StorageTankSizeFactor[j]) + \
                model.batchSize_log[i, j+1]
        def batch_size_rule(disjunct, i):
            return inequality(-log(model.StorageTankSizeFactorByProd[i,j]),
                              model.batchSize_log[i,j] - model.batchSize_log[i, j+1],
                              log(model.StorageTankSizeFactorByProd[i,j]))
        def no_batch_rule(disjunct, i):
            return model.batchSize_log[i,j] - model.batchSize_log[i,j+1] == 0

        if selectStorageTank:
            disjunct.volume_stage_j = Constraint(model.PRODUCTS, rule=volume_stage_j_rule)
            disjunct.volume_stage_jPlus1 = Constraint(model.PRODUCTS,
                                                      rule=volume_stage_jPlus1_rule)
            disjunct.batch_size = Constraint(model.PRODUCTS, rule=batch_size_rule)
        else:
            # The formulation says 0, but GAMS has this constant.
            # 04/04: Francisco says volume should be free:
            # disjunct.no_volume = Constraint(expr=model.storageTankSize_log[j] == MinFlow)
            disjunct.no_batch = Constraint(model.PRODUCTS, rule=no_batch_rule)
    model.storage_tank_selection_disjunct = Disjunct([0,1], model.STAGESExceptLast,
                                           rule=storage_tank_selection_disjunct_rule)

    def select_storage_tanks_rule(model, j):
        return [model.storage_tank_selection_disjunct[selectTank, j] for selectTank in [0,1]]
    model.select_storage_tanks = Disjunction(model.STAGESExceptLast, rule=select_storage_tanks_rule)

    # though this is a disjunction in the GAMs model, it is more efficiently formulated this way:
    # TODO: what on earth is k?
    def units_out_of_phase_rule(model, j):
        return model.unitsOutOfPhase_log[j] == sum(model.LogCoeffs[k] * model.outOfPhase[j,k] \
                                                   for k in model.PARALLELUNITS)
    model.units_out_of_phase = Constraint(model.STAGES, rule=units_out_of_phase_rule)

    def units_in_phase_rule(model, j):
        return model.unitsInPhase_log[j] == sum(model.LogCoeffs[k] * model.inPhase[j,k] \
                                                for k in model.PARALLELUNITS)
    model.units_in_phase = Constraint(model.STAGES, rule=units_in_phase_rule)

    # and since I didn't do the disjunction as a disjunction, we need the XORs:
    def units_out_of_phase_xor_rule(model, j):
        return sum(model.outOfPhase[j,k] for k in model.PARALLELUNITS) == 1
    model.units_out_of_phase_xor = Constraint(model.STAGES, rule=units_out_of_phase_xor_rule)

    def units_in_phase_xor_rule(model, j):
        return sum(model.inPhase[j,k] for k in model.PARALLELUNITS) == 1
    model.units_in_phase_xor = Constraint(model.STAGES, rule=units_in_phase_xor_rule)

    return model


def solve_with_big_m(m):
    TransformationFactory('gdp.bigm').apply_to(m)
    # SOLVE
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    #opt = SolverFactory('ipopt')
    #results = opt.solve(m)
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='baron')
    results = opt.solve(m, tee=True,
                        # Uncomment the following lines if you want to save GAMS models
                        #keepfiles=True,
                        #tmpdir=gams_path,
                        #symbolic_solver_labels=True,
                        add_options=[
                            'option reslim = 120;'
                            'option optcr = 0.0;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            # 'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            # '$onecho > baron.opt \n'
                            # 'CompIIS 1 \n'
                            # '$offecho'
                            # 'display(execError);'
                        ])
    return results

def solve_with_hull(m):
    TransformationFactory('gdp.hull').apply_to(m)
    # SOLVE
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    #opt = SolverFactory('ipopt')
    #results = opt.solve(m)
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='baron')
    results = opt.solve(m, tee=True,
                        # Uncomment the following lines if you want to save GAMS models
                        #keepfiles=True,
                        #tmpdir=gams_path,
                        #symbolic_solver_labels=True,
                        add_options=[
                            'option reslim = 120;'
                            'option optcr = 0.0;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            # 'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            # '$onecho > baron.opt \n'
                            # 'CompIIS 1 \n'
                            # '$offecho'
                            # 'display(execError);'
                        ])
    return results

def solve_with_gdpopt(m):
    # SOLVE
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    #opt = SolverFactory('ipopt')
    #results = opt.solve(m)
    solvername = 'gdpopt'
    opt = SolverFactory(solvername)
    results = opt.solve(m, tee=True,
                    strategy='LOA',
                    # strategy='GLOA',
                    time_limit=3600,
                    mip_solver='gams',
                    mip_solver_args=dict(solver='cplex', warmstart=True),
                    nlp_solver='gams',
                    nlp_solver_args=dict(solver='ipopth', warmstart=True,),
                    minlp_solver='gams',
                    minlp_solver_args=dict(solver='dicopt', warmstart=True),
                    subproblem_presolve=False,
                    # init_strategy='no_init',
                    set_cover_iterlim=20,
                    # calc_disjunctive_bounds=True
                    )
    return results



def external_ref(m, x, logic_expr = None):
    ext_var={}
    p=0
    for stage in m.STAGES:
        ext_var[stage,'outPhase']=x[p]
        p=p+1

    for stage in m.STAGES:
        ext_var[stage,'inPhase']=x[p]
        p=p+1

    for stage in m.STAGES:
        for parallel in m.PARALLELUNITS:
            if parallel==ext_var[stage,'outPhase']:
                m.outOfPhase[stage,parallel].fix(True)
            else:
                m.outOfPhase[stage,parallel].fix(False)

    for stage in m.STAGES:
        for parallel in m.PARALLELUNITS:
            if parallel==ext_var[stage,'inPhase']:
                m.inPhase[stage,parallel].fix(True)
            else:
                m.inPhase[stage,parallel].fix(False)      

    TransformationFactory('core.logical_to_linear').apply_to(m)
    #TransformationFactory('gdp.fix_disjuncts').apply_to(m)         #NOT WORKING HERE!!
    TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(m, tmp=False, ignore_infeasible=True)
    
    return m          

    
if __name__ == "__main__":
    m = build_model().create_instance('data_101006.dat')
    # results = solve_with_big_m(m)
    # results = solve_with_hull(m)
    # results = solve_with_gdpopt(m)
    # print(results)
    # m.min_cost.display()


    #External reformulation test (this can be deleted)
    inputval=[random.randrange(1, 7) for i in range(20)]
    external_ref(m, inputval, logic_expr = None)
    p=0
    for stage in m.STAGES:
        print('\n External variable '+str(p+1)+'='+str(inputval[p]))
        p=p+1
        for parallel in m.PARALLELUNITS:
            print(str( m.outOfPhase[stage,parallel])+'='+str(m.outOfPhase[stage,parallel].value))
    for stage in m.STAGES:
        print('\n External variable '+str(p+1)+'='+str(inputval[p]))
        p=p+1
        for parallel in m.PARALLELUNITS:
            print(str( m.inPhase[stage,parallel])+'='+str(m.inPhase[stage,parallel].value))