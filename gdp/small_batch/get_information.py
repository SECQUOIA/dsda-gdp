from gdp_small_batch import build_small_batch
import pyomo.environ as pe


def get_external_information(m,Ext_Ref):
    # We work with the linearized logical model to extract the information
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    # Identify Boolean variables introduced by the user and sets where they are defined
    # original_Boolean_vars={}
    # Boolean_without_ref={}
    # for i in Ext_Ref:
    #     original_Boolean_vars[i]=[]
    #     Boolean_without_ref[i]=[]
    #     partial_index=0
    #     for sets in i.index_set()._sets:
    #         original_Boolean_vars[i].append(sets)
    #         if sets.name != Ext_Ref[i].name:
    #             Boolean_without_ref[i].append(sets)




    # for i in Ext_Ref:
    #     original_Boolean_vars[i.name]=[]
    #     Boolean_without_ref[i.name]=[]
    #     for sets in i.index_set()._sets:
    #         original_Boolean_vars[i.name].append(sets.name)
    #         if sets.name != Ext_Ref[i].name:
    #             Boolean_without_ref[i.name].append(sets.name)



    # print(original_Boolean_vars) #Variables that are going to be reformulated
    # print(Ext_Ref) #Set with respect the reformualtio occurs
    # print(Boolean_without_ref)        #Sets where the reformulation does not occur



    ref_index={}   #index of the set where reformultion can be applied for a given boolean variable
    no_ref_index={} #index of the sets where the reformulation cannot be applied for a given boolean variable
    for i in Ext_Ref:
        ref_index[i]=[]
        no_ref_index[i]=[]
        for index_set in range(len(i.index_set()._sets)):
            if i.index_set()._sets[index_set].name ==Ext_Ref[i].name:
                ref_index[i].append(index_set)
            else:
                no_ref_index[i].append(index_set)
    print(ref_index)
    print(no_ref_index)    





        

    #Identify the variables that can be reformualted by performing a loop over logical constraints
    # For the moment we will work with exactly 1 type constraints only
    exactly_dict={}
    for c in m.component_data_objects(pe.LogicalConstraint, descend_into=True):
        if c.body.getname()=='exactly':
            exactly_number=c.body.args[0]
            print(exactly_number)
            for possible_Boolean in Ext_Ref:
            
                expected_Boolean=possible_Boolean.name #expected boolean variable where the reformualtion is going to be applied
                print(expected_Boolean)

                Boolean_name_list=[]
                Boolean_name_list=Boolean_name_list+[c.body.args[1:][k]._component()._name for k in range(len(c.body.args[1:]))]
                print(Boolean_name_list)
                if all(x==expected_Boolean for x in Boolean_name_list):
                    expected_ordered_set_index=ref_index[possible_Boolean] #expected ordered set index where the reformulation is going to be applied
                    print(expected_ordered_set_index)
                    index_of_other_sets=no_ref_index[possible_Boolean] #index of sets where the reformulation is not applied
                    print(index_of_other_sets)
                    if len(index_of_other_sets)>=1:   #If there are other indexes
                        Other_Sets_listOFlists=[]
                        verification_Other_Sets_listOFlists=[]
                        for j in index_of_other_sets:
                            Other_Sets_listOFlists.append([c.body.args[1:][k].index()[j] for k in range(len(c.body.args[1:]))])
                            if all(c.body.args[1:][x].index()[j]==c.body.args[1:][0].index()[j] for x in range(len(c.body.args[1:]))):
                                verification_Other_Sets_listOFlists.append(True)
                            else:
                                verification_Other_Sets_listOFlists.append(False)
                        print(verification_Other_Sets_listOFlists)
                        print(Other_Sets_listOFlists)
                        if all(verification_Other_Sets_listOFlists): #If we get to this point and it is true, it means that we can apply the reformulation for this combination of boolean var and exactly variable
                            print('we can apply the reformulation')


                    else:  #If there is only one index, then we can apply the reformulation at this point
                        print('we can apply the reformulation')
                    
                


                 
            

     


#1:Declare your model
m=build_small_batch()

#2:Specify information to be reformulated

#a1:  for the moment we will deal with Boolean variables only
Ext_Ref={m.Y:m.k}


get_external_information(m,Ext_Ref)








