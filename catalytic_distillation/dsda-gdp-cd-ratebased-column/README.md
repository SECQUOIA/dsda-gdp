Run main.gms

main.gms runs 2 codes:

The first (DISTILLATION.gms) generates de GDP rate-based catalytic distillation model.
It also contains the required information to perform the reformulation with external variables.
Inputs related to the reformualtion (ordered sets, independent Boolean variables, dependent Boolean variables and beighborhood) are itroduced in the last part of this code.



The second (D_SDA.gms) runs the D-SDA, using the information provided in DISTILLATION.gms.
