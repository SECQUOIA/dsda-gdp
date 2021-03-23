import pyomo.environ as pe

def get_external_information(m,Ext_Ref):
    # We work with the linearized logical model to extract the information
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)


    try:
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

    except:
        ref_index={}   #index of the set where reformultion can be applied for a given boolean variable
        no_ref_index={} #index of the sets where the reformulation cannot be applied for a given boolean variable
        for i in Ext_Ref:
            ref_index[i]=[]
            no_ref_index[i]=[]
            if i.index_set().name ==Ext_Ref[i].name:
                ref_index[i].append(0)
            else:
                no_ref_index[i].append(0)

  
    #Identify the variables that can be reformualted by performing a loop over logical constraints
    # For the moment we will work with exactly 1 type constraints only
    count=1;
    reformulation_dict={} #dict of dicts: it contains information from the exactly variables that can be reformualted into external variables.
    for c in m.component_data_objects(pe.LogicalConstraint, descend_into=True):
        if c.body.getname()=='exactly':
            exactly_number=c.body.args[0]
            for possible_Boolean in Ext_Ref:
            
                expected_Boolean=possible_Boolean.name #expected boolean variable where the reformualtion is going to be applied
                Boolean_name_list=[]
                Boolean_name_list=Boolean_name_list+[c.body.args[1:][k]._component()._name for k in range(len(c.body.args[1:]))]
                if all(x==expected_Boolean for x in Boolean_name_list):
                    expected_ordered_set_index=ref_index[possible_Boolean] #expected ordered set index where the reformulation is going to be applied
                    index_of_other_sets=no_ref_index[possible_Boolean] #index of sets where the reformulation is not applied
                    if len(index_of_other_sets)>=1:   #If there are other indexes
                        Other_Sets_listOFlists=[]
                        verification_Other_Sets_listOFlists=[]
                        for j in index_of_other_sets:
                            Other_Sets_listOFlists.append([c.body.args[1:][k].index()[j] for k in range(len(c.body.args[1:]))])
                            if all(c.body.args[1:][x].index()[j]==c.body.args[1:][0].index()[j] for x in range(len(c.body.args[1:]))):
                                verification_Other_Sets_listOFlists.append(True)
                            else:
                                verification_Other_Sets_listOFlists.append(False)
                        if all(verification_Other_Sets_listOFlists): #If we get to this point and it is true, it means that we can apply the reformulation for this combination of boolean var and exactly variable
                            reformulation_dict[count]={}
                            reformulation_dict[count]['exactly_number']=exactly_number
                            sorted_args=sorted(c.body.args[1:], key=lambda x:x.index()[expected_ordered_set_index[0]])# rearange boolean vars in cosntranit
                            #Now work with the ordered version sorted_args instead of c.body.args[1:]  
                            reformulation_dict[count]['Boolean_vars']=sorted_args
                            reformulation_dict[count]['Boolean_vars_names']=[sorted_args[k].name for k in range(len(sorted_args))]
                            reformulation_dict[count]['Boolean_vars_ordered_index']=[sorted_args[k].index()[expected_ordered_set_index[0]] for k in range(len(sorted_args))]
                            reformulation_dict[count]['Ext_var_lower_bound']=1
                            reformulation_dict[count]['Ext_var_upper_bound']=len(sorted_args)

                            count=count+1

                    else:  #If there is only one index, then we can apply the reformulation at this point
                        reformulation_dict[count]={}
                        reformulation_dict[count]['exactly_number']=exactly_number      
                        sorted_args=sorted(c.body.args[1:], key=lambda x:x.index())# rearange boolean vars in cosntranit
                            #Now work with the ordered version sorted_args instead of c.body.args[1:]  
                        reformulation_dict[count]['Boolean_vars']=sorted_args
                        reformulation_dict[count]['Boolean_vars_names']=[sorted_args[k].name for k in range(len(sorted_args))]
                        reformulation_dict[count]['Boolean_vars_ordered_index']=[sorted_args[k].index() for k in range(len(sorted_args))]
                        reformulation_dict[count]['Ext_var_lower_bound']=1
                        reformulation_dict[count]['Ext_var_upper_bound']=len(sorted_args)

                        count=count+1
                    
    number_of_external_variables=sum(reformulation_dict[j]['exactly_number'] for j in reformulation_dict)


    print('\n------------------------Reformulation Summary---------------------\n')
    exvar_num=0
    for i in reformulation_dict:
        for j in range(reformulation_dict[i]['exactly_number']):
            print('External variable x['+str(exvar_num)+'] '+' is associated to '+str(reformulation_dict[i]['Boolean_vars_names'])+' and it must be within '+str(reformulation_dict[i]['Ext_var_lower_bound'])+' and '+str(reformulation_dict[i]['Ext_var_upper_bound'])+'.')
            exvar_num=exvar_num+1

    print('\nThere are '+str(number_of_external_variables)+' external variables in total')

    return reformulation_dict,number_of_external_variables


def external_ref(m, x, dict_extvar,logic_expr):

    ext_var_position=0
    for i in dict_extvar:
        for j in range(dict_extvar[i]['exactly_number']):
            for k in range(1,len(dict_extvar[i]['Boolean_vars'])+1):
                if x[ext_var_position]==k:
                    dict_extvar[i]['Boolean_vars'][k-1].fix(True)  #fix True variables
            ext_var_position=ext_var_position+1 

    for i in dict_extvar:
        for j in range(dict_extvar[i]['exactly_number']):
            for k in range(1,len(dict_extvar[i]['Boolean_vars'])+1):
                if dict_extvar[i]['Boolean_vars'][k-1].is_fixed()==False:
                    dict_extvar[i]['Boolean_vars'][k-1].fix(False) #fix False variables

    for i in logic_expr:
        i[1].fix(pe.value(i[0]))

    ###pe.TransformationFactory('core.logical_to_linear').apply_to(m) ESTA EN get_external_information
    pe.TransformationFactory('gdp.fix_disjuncts').apply_to(m)
    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(m, tmp=False, ignore_infeasible=True)

    print('\n------------------------Fixed variables at current iteration---------------------\n')
    print('\n Independent Boolean variables\n')
    for i in dict_extvar:
        for k in range(1,len(dict_extvar[i]['Boolean_vars'])+1):
            print(dict_extvar[i]['Boolean_vars_names'][k-1]+'='+str(dict_extvar[i]['Boolean_vars'][k-1].value))

    print('\n Dependent Boolean variables and disjunctions\n')
    for i in logic_expr:
        print(i[1].name+'='+str(i[1].value))
    
    return m