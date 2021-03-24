# Import reformulation functios
from reformulation_functions import get_external_information, external_ref
# Import models
from small_batch_model import build_small_batch
from column_model import build_column
from cstr_model import build_cstrs
# Import pyomo
import pyomo.environ as pe

######TEST 1: small batch ##################################################
print('\n***BATCH\n')
# 1:User must declare the model
m = build_small_batch()

# 2:User must specify this information: Boolean variable and rest where the reformulation occurs
Ext_Ref = {m.Y: m.k}

# 3: Now we can get information from the model
# It is enough to run this function once (it depends on the model only)
get_info = get_external_information(m, Ext_Ref)
dict_extvar = get_info[0]  # dictionary required by external_ref function
num_extvar = get_info[1]  # number of external variables
lower = get_info[2]  # lower bounds
upper = get_info[3]  # upper bounds

# 4: User must specify equivalences between disjunctions and Boolean variables in the model
# List of list Value 1: Independent Boolean var or expression of independent boolean vars, Value 2: Equivalent indicator var or boolean var
# [[A,b],[C,d],[E,f]], where A<->b, C<->d, E<->f. A,C,E are Independent Boolean variables or expressions and b,d, and f are indicator or boolean vars
logic_expr = []
for k in m.k:
    for j in m.j:
        logic_expr.append([m.Y[k, j], m.Y_exists[k, j].indicator_var])
        logic_expr.append([~m.Y[k, j], m.Y_not_exists[k, j].indicator_var])

# 5: Now we can apply the external variables reformulation

external_ref(m, [1, 2, 1], dict_extvar, logic_expr)

######TEST 2: cstr network #################################################
print('\n***CSTR\n')
# 1:User must declare the model
m = build_cstrs()

# 2:User must specify this information: Boolean variable and rest where the reformulation occurs
Ext_Ref = {m.YF: m.N, m.YR: m.N}

# 3: Now we can get information from the model
# It is enough to run this function once (it depends on the model only)
get_info = get_external_information(m, Ext_Ref)
dict_extvar = get_info[0]
num_extvar = get_info[1]
lower = get_info[2]
upper = get_info[3]


# 4: User must specify equivalences between disjunctions and Boolean variables in the model
# List of list Value 1: Independent Boolean var or expression of independent boolean vars, Value 2: Equivalent indicator var or boolean var
# [[A,b],[C,d],[E,f]], where A<->b, C<->d, E<->f. A,C,E are Independent Boolean variables or expressions and b,d, and f are indicator or boolean vars
logic_expr = []

for n in m.N:
    logic_expr.append([m.YR[n], m.YR_is_recycle[n].indicator_var])
    logic_expr.append([~m.YR[n], m.YR_is_not_recycle[n].indicator_var])
    logic_expr.append([pe.lor(pe.land(~m.YF[n2] for n2 in range(
        1, n)), m.YF[n]), m.YP_is_cstr[n].indicator_var])
    logic_expr.append([~pe.lor(pe.land(~m.YF[n2] for n2 in range(
        1, n)), m.YF[n]), m.YP_is_bypass[n].indicator_var])

# 5: Now we can apply the external variables reformulation

external_ref(m, [1, 1], dict_extvar, logic_expr)

######TEST 3: column #######################################################
print('\n***COLUMN\n')

# 1:User must declare the model
m = build_column(8, 17, 0.95, 0.95)

# 2:User must specify this information: Boolean variable and rest where the reformulation occurs
Ext_Ref = {m.YB: m.intTrays, m.YR: m.intTrays}

# 3: Now we can get information from the model
# It is enough to run this function once (it depends on the model only)
get_info = get_external_information(m, Ext_Ref)
dict_extvar = get_info[0]
num_extvar = get_info[1]
lower = get_info[2]
upper = get_info[3]


# 4: User must specify equivalences between disjunctions and Boolean variables in the model
# List of list Value 1: Independent Boolean var or expression of independent boolean vars, Value 2: Equivalent indicator var or boolean var
# [[A,b],[C,d],[E,f]], where A<->b, C<->d, E<->f. A,C,E are Independent Boolean variables or expressions and b,d, and f are indicator or boolean vars
logic_expr = []
for n in m.intTrays:
    logic_expr.append([pe.land(~m.YR[n] for n in range(
        m.reboil_tray+1, m.feed_tray)), m.YR_is_down])
    logic_expr.append([pe.land(~m.YB[n]
                      for n in range(m.feed_tray+1, m.max_trays)), m.YB_is_up])
for n in m.conditional_trays:
    logic_expr.append([pe.land(pe.lor(m.YR[j] for j in range(n, m.max_trays)), pe.lor(
        pe.land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])), m.tray[n].indicator_var])
    logic_expr.append([~pe.land(pe.lor(m.YR[j] for j in range(n, m.max_trays)), pe.lor(
        pe.land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])), m.no_tray[n].indicator_var])


# 5: Now we can apply the external variables reformulation
external_ref(m, [13, 2], dict_extvar, logic_expr)
