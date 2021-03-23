from gdp_small_batch import build_small_batch
import pyomo.environ as pe
"""
This function retrieves informaton from a pyomo GDP model that 

"""

# def get_info(m,ExtRef):

#     sete=list(ExtRef['x1'][1])
#     AllVars2=list(ExtRef['x1'][0])
#     newdict={}
#     k=0
#     for i in range(9):
#         for j in  range(3):
#             if AllVars2[i][0] == sete[j]:
#                 newdict[k]==AllVars2[i][0]
#                 k=k+1
#     print(sete)
#     print(AllVars2)
#     print(newdict)
#     return AllVars2


#DECLARE YOUR MODEL
m=build_small_batch()
#m = k.create_instance()

#instance.display()
pe.TransformationFactory('core.logical_to_linear').apply_to(m)

for v in m.component_objects(pe.BooleanVar, descend_into=True):
    print("FOUND BOOLEAN:" + v.name)
    print(v)
    print(v.index_set())
    print(v.index_set()._sets[1])
    print(len(v.index_set()._sets))
    v.pprint()
    

for v_data in m.component_data_objects(pe.BooleanVar, descend_into=True):
    print("Found: "+v_data.name+", value = "+str(v_data.value))
    print(v_data.index()[1])
    #print(dir(v_data))
    print(v_data.is_fixed())
    v_data.pprint()
for c in m.component_data_objects(pe.LogicalConstraint, descend_into=True):
#    print(type(c.index()))
    print(c.index())
#    print(str(c.body.getname()))
#    print(c.is_indexed())
    print(c.body.args[3]._component()._name)
    print(c.body.args[3].index()[1])
    print(c.body.nargs())
    c.pprint()

#for c in m.component_data_objects(pe.LogicalConstraint, descend_into=True):
#    print(c)
#    c.pprint()  
#    print(c.body.getname())



for j in m.component_data_objects(pe.LogicalConstraint, descend_into=True):






#SPECIFY EXTERNAL VARIABLES
# ExtRef={'x1':[m.Y,m.k]}

#GET INFORMATION FROM MODEL
# get_info(ExtRef)




#PRINT INFORMATION

##print(AllVars)
#print(k.value for k in AllVars)

#_VarData.has_lb(m.Y)
