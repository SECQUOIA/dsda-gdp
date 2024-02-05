Run main.gms

main.gms runs 2 codes:

The first (DISTILLATION.gms) generates de GDP catalytic distillation model. 
It also contains the required information to perform the reformulation with external variables.
Inputs related to the reformulation (ordered sets, independent Boolean variables, dependent Boolean variables and neighborhood) are introduced in the last part of this code.



The second (D_SDA.gms) runs the D-SDA, using the information provided in DISTILLATION.gms.