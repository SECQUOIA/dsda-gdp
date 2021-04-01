from pyomo.environ import *

m = ConcreteModel()
m.x = Var()

m.N = RangeSet(1, 5, doc='Set of units in the superstructure')
m.YF = BooleanVar(m.N)
m.YP = BooleanVar(m.N)

# create a Pyomo expression
e1 = m.x + 5

# create another Pyomo expression
# e1 is copied when generating e2
e2 = e1 + m.x

m.I = Set(initialize=['A', 'B'], doc='Set of components')

# Create list of expressions
logic_expr = []
for n in m.N:
    logic_expr.append(lor(land(~m.YF[n2] for n2 in range(1, n)), m.YF[n]))

# Create list of logical constraints
m.logic_list = LogicalConstraintList()

# Populate constraint list with expressions
for n in m.N:
    m.logic_list.add(m.YP[n].equivalent_to(logic_expr[n-1]))


for n in m.N:
    if n == 4:
        m.YF[n].fix(True)
    else:
        m.YF[n].fix(False)


print(logic_expr)
for i in logic_expr:
    print(value(i))
# pe.value(pe.lor(pe.land(~m.YF[n2]
#                         for n2 in range(1, n)), m.YF[n]))
m.conlist = ConstraintList()
m.con = Constraint(expr=m.x**2 == e2)

for i in m.I:
    m.conlist.add(m.x**2 == e2)

m.x.fix(1)
print(value(e2))

m.pprint()
