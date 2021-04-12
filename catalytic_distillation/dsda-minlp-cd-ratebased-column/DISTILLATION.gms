$ontext
================================================================================
                Dise�o de una columna de destilacion catalitica para
                la producci�n de ETBE con un modelo de no equilibrio

                               CODIGO NLP


          David Esteban Bernal Neira-David Alejandro Li�an Romero
                                 2020

================================================================================
$offtext
*-------------------------------------------------------------------------------
*                                Secci�n 1
*           Conjuntos para definir operacion en estado estable
*-------------------------------------------------------------------------------
*Se usa solo el elemento inicial con el punto de colocacion inicial: estado estable
*Esto es equivalente a quitar j y N de la formualcion
sets j "1 punto de colocacion (0)" /1/
     N "1 elemento finito" /1/;
*-------------------------------------------------------------------------------
*                                Secci�n 2
*       Conjuntos, variables, par�metros y ecuaciones principales del sistema
*-------------------------------------------------------------------------------
*Conjuntos
set comp "lista de componentes que intervienen en el sistema" /iButene, Ethanol, nButene, ETBE/;
set compR(comp) "componentes para calculos matriciales" /iButene, Ethanol, nButene/;
set compREF(comp) "componente de referencia" /ETBE/;
alias(comp,comp1,comp2,comp3);
alias(compR,compR1,compR2);
parameter ident(comp) "identificador de componente"
/
iButene 1
Ethanol 2
nButene 3
ETBE 4
/;


set Net "Todas las etapas de la columna reales y sobrantes incluyendo con y rev" /1*45/;
alias(Net,Net1);

set F "Conjunto de alimentaciones" /1,2/;

*Variables principales
Variables
Nbonita(N,j,comp,net)  "Transferencia de masa interfacial por componente [mol/min]"
Ebonita(N,j,net)  "Transferencia de energia interfacial [kJ/min]"
ktransVap(N,j,net,comp,comp1) "Coeficiente de transferencia de masa multicomponente en el vapor [m/min]"
ktransLiq(N,j,net,comp,comp1) "Coeficiente de transferencia de masa molticomponente en el liquido [m/min]"
;
positive variables
L(N,j,Net)  "Flujo de l�quido [mol/min]"
V(N,j,Net)  "Flujo de vapor [mol/min]"
x(N,j,comp,Net)     "Porcentaje molar en el l�qudio (bulk) [%]"
y(N,j,comp,Net)     "Porcentaje molar en el vapor (bulk) [%]"
xI(N,j,comp,Net)     "Porcentaje molar en el l�qudio (Interface) [%]"
yI(N,j,comp,Net)     "Porcentaje molar en el vapor (Interface) [%]"
TempL(N,j,Net)       "Temperatura de operaci�n (Liquido) [K]"
TempV(N,j,Net)       "Temperatura de operaci�n (Vapor) [K]"
TempI(N,j,net)       "Temperatura de operacion (Interface) [K]"
P(N,j,Net)  "Presi�n por etapa [bar]"
Z(N,j,Net)  "Coeficiente de compresibilidad [-]"
RR(N,j)      "Relaci�n molar de reflujo [-]"
Qc(N,j)      "Carga t�rmica del condensador [kJ/min]"
Qr(N,j)      "Carga t�rmica del rehervidor [kJ/min]"
BR(N,j)       "Boil up [-]"
;

*Par�metros hidr�ulicos
parameter aI(N,j,net) "Area interfacial de transferencia por etapa (Resultados del modelo son independientes de este valor) [m^2]" ;
aI(N,j,net)=1;
parameter
da      "Di�metro de los agujeros [m]"  /2E-3/
ep      "Espesor del plato [m]" /0.002/
pitch   "Distancia entre agujeros [m]"  /0.009/
Sfactor "Factor de seguridad altura de la columna [-]" /0.15/
poro  "Porosidad del plato [-]"
K0    "Coeficiente de orificio [-]"
;
poro=0.907*sqr(da/pitch);
K0=(880.6-(67.7*da/ep)+(7.32*(sqr(da/ep)))-(0.338*(power(da/ep,3))))*1E-3;
*Variables hidraulicas

positive variable D "Di�metro de la columna [m]";
positive variable hw      "Weir height [m]";
positive variable  HS      "Altura de cada plato [m]";
positive variable Htotal "Altura total de la columna [m]";
positive variable At "Area activa [m2]";
positive variable Ad  "Area de derramadero [m2]";
positive variable Lw "Weir length [m]";
positive variable A0 "Area agujerada [m2]";

*Ecuaciones hidraulicas
Equation EqHwmin;
EqHwmin.. hw=g=0.05*HS;
Equation EqHwmax;
EqHwmax.. hw=l=HS/3;
equation EqAt;
EqAt.. At=e=sqr(D/2)*(pi-(1.854590-0.96));
equation EqAd;
EqAd.. Ad=e=sqr(D/2)*(0.5*(1.854590-0.96));
equation EqLw;
EqLw.. Lw=e=0.8*D;
equation EqA0;
EqA0.. A0=e=At*poro;

*Alimentacion 1 (butenos en el caso de ETBE)

parameter FB "Flujo de alimentaci�n de butenos [mol/min]" /5.774/
parameter zb(N,j,comp) "Porcentaje molar en la alimentaci�n de butenos";
zb(N,j,'iButene')=30;
zb(N,j,'nButene')=100-zb(N,j,'iButene');
zb(N,j,'Ethanol')=0;
zb(N,j,'ETBE')=0;

*Alimentacion 2 (Etanol en el caso de ETBE)

parameter
FE  "Flujo de alimentaci�n de etanol [mol-h]" /1.7118/
ze(comp)        "Porcentaje molar en la alimentaci�n de etanol"
/
iButene 0
Ethanol 100
nButene 0
ETBE 0
/
;
*Parametros de operacion

parameter
Pop     "Presi�n de operaci�n condensador [bar]"        /9.5/
TaliB   "Temperatura de alimentaci�n de butenos [K]"    /323/
TaliE   "Temperatura de alimentaci�n etanol [K]"        /342.38/
xBetbe  "Composici�n molar de ETBE en fondos deseada"   /83/
MCR     "Retenci�n constante en rehervidor y condensador [mol]" /1/
cR      "Constante de los gases [m3*bar/K*mol]" /0.00008314/
;

*-------------------------------------------------------------------------------
*                                Secci�n 3
*                    Parametro de conversion de unidades
*-------------------------------------------------------------------------------
parameter hora    "Si estamos en an�lisis por minutos u hora [s]" /60/;
*-------------------------------------------------------------------------------
*                                Secci�n 4
*                          Restricciones de pureza
*-------------------------------------------------------------------------------

equations pureza0(Net);
pureza0(Net)$(ord(Net) eq card(Net)).. x('1','1','ETBE',Net)=g=xBetbe;

*-------------------------------------------------------------------------------
*                                Secci�n 5
*    C�lculo de presiones de saturaci�n por medio de la ecuaci�n de Antoine
*-------------------------------------------------------------------------------
*Constantes de la ecuaci�n de Antoine expandida

parameters
C1a(comp)
/
iButene 66.4970745
Ethanol 61.7910745
nButene 40.3230745
ETBE    52.67507454
/
C2a(comp)
/
iButene -4634.1
Ethanol -7122.3
nButene -4019.2
ETBE    -5820.2
/
C3a(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C4a(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C5a(comp)
/
iButene -8.9575
Ethanol -7.1424
nButene -4.5229
ETBE    -6.1343
/
C6a(comp)
/
iButene 1.3413E-5
Ethanol 2.8853E-6
nButene 4.8833E-17
ETBE    2.1405E-17
/
C7a(comp)
/
iButene 2
Ethanol 2
nButene 6
ETBE    6
/
;

positive variables PsatI(N,j,comp,Net) presi�n de saturaci�n interfacial (bar);
equations EqPsatI(N,j,comp,Net);
EqPsatI(N,j,comp,Net).. PsatI(N,j,comp,Net)=e=exp( C1a(comp) + (C2a(comp)/(TempI(N,j,Net)+C3a(comp))) + (C4a(comp)*TempI(N,j,Net)) + (C5a(comp)*log(TempI(N,j,Net)) + (C6a(comp)*power(TempI(N,j,Net),C7a(comp)))) );
*-------------------------------------------------------------------------------
*                                Secci�n 6
*     C�lculo de densidades de l�quido por medio de la ecuaci�n IK-CAPI
*     C�lculo de densidades de l�quido por medio de la ecuaci�n DIPPR cr�tica
*     C�lculo de densidades de gas por medio de ecuaci�n de gas ideal corregida
*-------------------------------------------------------------------------------
*Constantes de la ecuaci�n DIPPR

parameters
MW(comp) "Peso molecular [kg/kmol]"
/
iButene 56.10752
Ethanol 46.06904
nButene 56.10752
ETBE    102.17656
/
Tcrit(comp) "Temperatura cr�tica [K]"
/
iButene 417.9
Ethanol 516.2
nButene 419.6
ETBE    509.4
/
Pcrit(comp) "Presi�n cr�tica [bar]"
/
iButene 38.98675
Ethanol 60.35675
nButene 39.18675
ETBE    28.32675
/
C1rh(comp)
/
iButene        8.9711123119
Ethanol        -2.932961888E-2
nButene        5.956235579
ETBE           -1.323678817E-1
/
C2rh(comp)
/
iButene        0
Ethanol        6.9361857406E-4
nButene        0
ETBE           2.1486345729E-3
/
C3rh(comp)
/
iButene 0
Ethanol -1.962897037E-6
nButene 0
ETBE    -6.092181735E-6
/
C4rh(comp)
/
iButene 0
Ethanol 2.089632106E-9
nButene 0
ETBE    6.4627035532E-9
/
C5rh(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C6rh(comp)
/
iButene -1.4666609E-10
Ethanol 0
nButene -9.3717935E-11
ETBE    0
/
C7rh(comp)
/
iButene 1.286186216E-12
Ethanol 0
nButene 8.150339357E-13
ETBE    0
/
C8rh(comp)
/
iButene -4.33826109E-15
Ethanol 0
nButene -2.72421122E-15
ETBE    0
/
C9rh(comp)
/
iButene 6.619652613E-18
Ethanol 0
nButene 4.115761136E-18
ETBE    0
/
C10rh(comp)
/
iButene -3.8362103001E-21
Ethanol 0
nButene -2.3593237507E-21
ETBE    0
/
C1r(comp)
/
iButene        1.1446
Ethanol        1.6288
nButene        1.0877
ETBE        0.66333
/
C2r(comp)
/
iButene        0.2724
Ethanol        0.27469
nButene        2.6454E-01
ETBE        2.6135E-01
/
C3r(comp)
/
iButene        0.28172
Ethanol        0.23178
nButene        0.2843
ETBE        0.28571
/
C4r(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
;

positive variable Tcritm(N,j,Net);

equation EqTcritm(N,j,Net);
EqTcritm(N,j,Net).. Tcritm(N,j,Net) =e= (sqr(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/sqrt(Pcrit(comp)))))/(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/Pcrit(comp)));
positive variables rho(N,j,comp,Net) "Densidad molar por componente de l�quido [mol/m^3]";
equation Eqrho(N,j,comp,Net);
Eqrho(N,j,comp,Net).. rho(N,j,comp,Net)=e=( C1r(comp)/(C2r(comp)**(1+((1-(TempL(N,j,Net)/Tcritm(N,j,Net)))**C4r(comp)))) )*1000;
equation DomError1(N,j,Net);
DomError1(N,j,Net)..TempL(N,j,Net)/Tcritm(N,j,Net)=l=1-(1e-3);


positive variable rhoV(N,j,Net) "Densidad molar de vapor [mol/m^3]";
equation EqurhoV(N,j,Net);
EqurhoV(N,j,Net).. rhoV(N,j,Net)=e=P(N,j,Net)/(0.00008314*TempV(N,j,Net)*(Z(N,j,Net)));


positive variable rhoL(N,j,Net) "Densidad molar de liquido [mol/m^3]";
equation EqurhoL(N,j,Net);
EqurhoL(N,j,Net)..rhoL(N,j,Net)=e=sum(comp,(x(N,j,comp,Net)/100)*rho(N,j,comp,Net));

*-------------------------------------------------------------------------------
*                                Secci�n 7
*     C�lculo de tensi�n superficial por medio de la ecuaci�n DIPPR cr�tica
*-------------------------------------------------------------------------------

*Constantes de la ecuaci�n DIPPR
parameters
C1sig(comp)
/
iButene        0.05544
Ethanol        0.03764
nButene        0.055945
ETBE        0.071885
/
C2sig(comp)
/
iButene        1.2453
Ethanol        -2.157E-5
nButene        1.2402
ETBE        2.1204
/
C3sig(comp)
/
iButene        0.0
Ethanol        1.025E-7
nButene        0
ETBE        -1.5583
/
C4sig(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0.76657
/
;
positive variables sigma(N,j,Net) "Tensi�n superficial l�quido vapor [N/m]";
equation Eqsigma(N,j,Net);
Eqsigma(N,j,Net)..  sigma(N,j,Net)=e=sum(comp,(x(N,j,comp,Net)/100)*C1sig(comp)*(1-(TempL(N,j,Net)/Tcritm(N,j,Net)))**(C2sig(comp)+C3sig(comp)*(Templ(N,j,Net)/Tcritm(N,j,Net))+C4sig(comp)*(sqr(Templ(N,j,Net)/Tcritm(N,j,Net)))));
*-------------------------------------------------------------------------------
*                                Secci�n 8
*          C�lculo de coeficientes de actividad por medio del modelo NRTL
*-------------------------------------------------------------------------------
*Parametros a b c
table a_nrtl(comp,comp) Par�metro a de NRTL
                       iButene            Ethanol            nButene            ETBE
iButene                0.0                0.0                0.0                0.0
Ethanol                0.0                0.0                0.0                0.0
nButene                0.0                0.0                0.0                0.0
ETBE                   0.0                0.0                0.0                0.0
;

table b_nrtl(comp,comp) Par�metro b de NRTL
                iButene                Ethanol            nButene            ETBE
iButene         0.0                    623.5810010        107.526499         219.73407
Ethanol         141.9632130            0.0                164.57256          187.104064
nButene         -93.24546420           595.5299820        0.0                226.373398
ETBE            -172.59152             344.481315         -177.88565         0.0
;

table c_nrtl(comp,comp) Par�metro c de NRTL
                       iButene            Ethanol            nButene            ETBE
iButene                0.0                0.3                0.3                0.3
Ethanol                0.3                0.0                0.3                0.3
nButene                0.3                0.3                0.0                0.3
ETBE                   0.3                0.3                0.3                0.0
;
parameter alfa_nrtl(comp,comp);
alfa_nrtl(comp,comp1)$(ord(comp) ne ord(comp1))=c_nrtl(comp,comp1);

*Par�metros G y Tao equilibrio
variables tao_nrtlI(N,j,comp,comp1,Net);
equations Eq_tao_nrtlI(N,j,comp,comp1,Net),Eq_tao_nrt2I(N,j,comp,comp1,Net);
Eq_tao_nrtlI(N,j,comp,comp1,Net)$(ord(comp) ne ord(comp1)).. tao_nrtlI(N,j,comp,comp1,Net)=e=a_nrtl(comp,comp1) + (b_nrtl(comp,comp1)/TempI(N,j,Net));
Eq_tao_nrt2I(N,j,comp,comp1,Net)$(ord(comp) eq ord(comp1)).. tao_nrtlI(N,j,comp,comp1,Net)=e=0;

variables g_nrtlI(N,j,comp,comp1,Net);
equations Eq_g_nrtlI(N,j,comp,comp1,Net),Eq_g_nrt2I(N,j,comp,comp1,Net);
Eq_g_nrtlI(N,j,comp,comp1,Net)$(ord(comp) ne ord(comp1)).. g_nrtlI(N,j,comp,comp1,Net)=e=exp( -alfa_nrtl(comp,comp1)*tao_nrtlI(N,j,comp,comp1,Net));
Eq_g_nrt2I(N,j,comp,comp1,Net)$(ord(comp) eq ord(comp1)).. g_nrtlI(N,j,comp,comp1,Net)=e=1;

*Coeficiente de actividad equilibrio(gammaI)
variables gammaI(N,j,comp,Net);
equations EqgammaI(N,j,comp,Net);
EqgammaI(N,j,comp,Net).. gammaI(N,j,comp,Net)=e=
        exp(sum(comp1,xI(N,j,comp1,Net)*tao_nrtlI(N,j,comp1,comp,Net)*
        g_nrtlI(N,j,comp1,comp,Net))/sum(comp1,xI(N,j,comp1,Net)*
        g_nrtlI(N,j,comp1,comp,Net))+sum(comp1,xI(N,j,comp1,Net)*
        g_nrtlI(N,j,comp,comp1,Net)/sum(comp2,xI(N,j,comp2,Net)*
        g_nrtlI(N,j,comp2,comp1,Net))*(tao_nrtlI(N,j,comp,comp1,Net)-
        sum(comp2,xI(N,j,comp2,Net)*tao_nrtlI(N,j,comp2,comp1,Net)*
        g_nrtlI(N,j,comp2,comp1,Net))/sum(comp3,xI(N,j,comp3,Net)*
        g_nrtlI(N,j,comp3,comp1,Net)))));

*Factor termodinamico en la interfase
variable termofacLiq(N,j,net,comp,comp1) "factor termodinamico liquido";
parameter Kronecker(comp,comp1) "matriz identidad";
Kronecker(comp,comp1)$(ord(comp) eq ord(comp1))=1;

equation deftermofacLiq(N,j,net,comp,comp1,compREF);
deftermofacLiq(N,j,net,comp,comp1,compREF)$(ord(Net)>1 and ord(Net)<card(Net))..termofacLiq(N,j,net,comp,comp1)=e=Kronecker(comp,comp1)+xI(N,j,comp,Net)*(((((g_nrtlI(N,j,comp,comp1,Net))*((tao_nrtlI(N,j,comp,comp1,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp1,Net)*tao_nrtlI(N,j,comp3,comp1,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp1,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp1,Net))))+(((g_nrtlI(N,j,comp1,comp,Net))*((tao_nrtlI(N,j,comp1,comp,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net)*tao_nrtlI(N,j,comp3,comp,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net)))) -sum(comp2,(xI(N,j,comp2,Net)*((g_nrtlI(N,j,comp,comp2,Net)*(((g_nrtlI(N,j,comp1,comp2,Net))*((tao_nrtlI(N,j,comp1,comp2,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)*tao_nrtlI(N,j,comp3,comp2,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))+(g_nrtlI(N,j,comp1,comp2,Net)*((g_nrtlI(N,j,comp,comp2,Net))*((tao_nrtlI(N,j,comp,comp2,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)*tao_nrtlI(N,j,comp3,comp2,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net))))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))-((((g_nrtlI(N,j,comp,compREF,Net))*((tao_nrtlI(N,j,comp,compREF,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,compREF,Net)*tao_nrtlI(N,j,comp3,compREF,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,compREF,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,compREF,Net))))+(((g_nrtlI(N,j,compREF,comp,Net))*((tao_nrtlI(N,j,compREF,comp,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net)*tao_nrtlI(N,j,comp3,comp,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp,Net))))-sum(comp2,(xI(N,j,comp2,Net)*((g_nrtlI(N,j,comp,comp2,Net)*(((g_nrtlI(N,j,compREF,comp2,Net))*((tao_nrtlI(N,j,compREF,comp2,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)*tao_nrtlI(N,j,comp3,comp2,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))+(g_nrtlI(N,j,compREF,comp2,Net)*((g_nrtlI(N,j,comp,comp2,Net))*((tao_nrtlI(N,j,comp,comp2,Net))-(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)*tao_nrtlI(N,j,comp3,comp2,Net))/sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net)))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net))))))/(sum(comp3,xI(N,j,comp3,Net)*g_nrtlI(N,j,comp3,comp2,Net))))));
*-------------------------------------------------------------------------------
*                                Secci�n 9
*                           C�lculo de reacci�n qu�mica
*-------------------------------------------------------------------------------
Parameter
Nu(comp) "Coeficientes estequiom�tricos en la reacci�n"
/
iButene -1
Ethanol -1
nButene 0
ETBE    1
/
mcat  "Masa del catalizador"     /0.4/
;
variable Ketbe(N,j,Net) "Constante de equilibrio [-]";
equation EqKetbe(N,j,Net);
EqKetbe(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. Ketbe(N,j,Net) =e=
                exp(10.387+4060.59/(Templ(N,j,Net))
                -2.89055*log(Templ(N,j,Net))
                -0.01915144*Templ(N,j,Net)
                +0.0000528586*power(Templ(N,j,Net),2)
                -0.0000000532977*power(Templ(N,j,Net),3));

positive variable Krate(N,j,Net) "Tasa de avance de reacci�n [mol/(kg_cat.min)]";
equation EqKrate(N,j,Net);
EqKrate(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Krate(N,j,Net) =e= 7.41816E15*exp(-60400.0/(8.314*Templ(N,j,Net)))*hora/3600;
positive variable Ka(N,j,Net) "Tasa de adsorci�n";
equation EqKa(N,j,Net);
EqKa(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Ka(N,j,Net) =e= exp(-1.0707+1323.1/Templ(N,j,Net));
variable Rx(N,j,Net) "Tasa de reacci�n [mol/(kg_cat.min)]";
equation EqRx(N,j,Net);
EqRx(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Rx(N,j,Net)*(power(1+Ka(N,j,Net)*gammaI(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100,3))*Ketbe(N,j,Net) =e=
                        (Krate(N,j,Net)*(gammaI(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100))
                        *((Ketbe(N,j,Net)*gammaI(N,j,'iButene',Net)*x(N,j,'iButene',Net)/100*gammaI(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100)
                        -(gammaI(N,j,'ETBE',Net)*x(N,j,'ETBE',Net)/100));

*-------------------------------------------------------------------------------
*                                Secci�n 10
*                           Ecuaci�n de estado (calculo de phi)
*-------------------------------------------------------------------------------
parameter
Omega(comp) "Factor ac�ntrico [-]"
/
iButene 0.19484
Ethanol 0.643558
nButene 0.184495
ETBE    0.316231
/
TcritSRK(comp) "Temperatura cr�tica de Soave-Redlich-Kwong [K]"
/
iButene 417.9
Ethanol 514
nButene 419.5
ETBE    509.4
/
mEOS(comp) "Parameter m in EOS"
biEOS(comp) "Parameter bi in EOS"
;
mEOS(comp)=0.48508+1.55171*Omega(comp)-0.15613*sqr(Omega(comp));
biEOS(comp)=0.08664*0.00008314*TcritSRK(comp)/Pcrit(comp);
positive variable alphaEOS(N,j,comp,Net);
equation EqAlphaEOS(N,j,comp,Net);
EqAlphaEOS(N,j,comp,Net).. alphaEOS(N,j,comp,Net) =e= sqr(1+mEOS(comp)*(1-sqrt(TempV(N,j,Net)/Tcritm(N,j,Net))));

positive variable aiEOS(N,j,comp,Net);
equation EqaiEOS(N,j,comp,Net);
EqaiEOS(N,j,comp,Net).. aiEOS(N,j,comp,Net) =e= alphaEOS(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);


positive variable bEOS(N,j,Net);
equation EqbEOS(N,j,Net);
EqbEOS(N,j,Net).. bEOS(N,j,Net) =e= sum(comp,(y(N,j,comp,Net)/100)*biEOS(comp));


positive variable aEOS(N,j,Net);
equation EqaEOS(N,j,Net);
EqaEOS(N,j,Net).. aEOS(N,j,Net) =e= sum(comp,sum(comp1, (y(N,j,comp,Net)/100)*(y(N,j,comp1,Net)/100)*sqrt(aiEOS(N,j,comp,Net)*aiEOS(N,j,comp1,Net))));


equation VaporZ(N,j,Net);
VaporZ(N,j,Net).. power(Z(N,j,Net),3)-sqr(Z(N,j,Net))+(Z(N,j,Net))
                *((aEOS(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*TempV(N,j,Net))))
                -(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*TempV(N,j,Net)))
                -sqr(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*TempV(N,j,Net))))
                -((aEOS(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*TempV(N,j,Net)))))
                *(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*TempV(N,j,Net))) =e= 0;



positive variable phiI(N,j,comp,Net);
equation EqPhi(N,j,comp,Net);
EqPhi(N,j,comp,Net).. phiI(N,j,comp,Net) =e= exp(((Z(N,j,Net))-1)*biEOS(comp)/bEOS(N,j,Net)
                                        -log((Z(N,j,Net))-bEOS(N,j,Net))
                                        -(aEOS(N,j,Net)/bEOS(N,j,Net))
                                        *(2*(sqrt(aiEOS(N,j,comp,Net)/aEOS(N,j,Net)))
                                        -biEOS(comp)/bEOS(N,j,Net))*log(((Z(N,j,Net))
                                        -bEOS(N,j,Net))/(Z(N,j,Net))));

*-------------------------------------------------------------------------------
*                                Secci�n 11
*                           C�lculo de entalp�as
*-------------------------------------------------------------------------------
*Constantes de Cp (kJ/mol.K) gas ideal

parameters
C1c(comp)
/
iButene        0.016052191
Ethanol        0.00901418
nButene        -0.00299356
ETBE        -0.014651654
/
C2c(comp)
/
iButene        0.000280432
Ethanol        0.000214071
nButene        0.000353198
ETBE        0.000698631
/
C3c(comp)
/
iButene        -0.00000010914988
Ethanol        -0.000000083903472
nButene        -0.00000019904047
ETBE           -0.00000044791741
/
C4c(comp)
/
iButene        0.0000000000090979164
Ethanol        0.0000000000013732704
nButene        0.000000000044631288
ETBE           0.00000000011636811
/
C5c(comp)
/
iButene        0
Ethanol        0
nButene        0
ETBE           0
/
C6c(comp)
/
iButene        0
Ethanol        0
nButene        0
ETBE           0
/
;

*Entalp�a de formaci�n del gas ideal y temperatura de referencia

parameter Tref "Temperatura de referencia [K]" /298.15/;
parameter Hform(comp) "Entalp�a de formaci�n (kJ/mol)"
/
iButene -16.9147
Ethanol -234.963
nButene -0.125604
ETBE        -313.9
/;

*Temperatura de ebullici�n a la presi�n de referencia

parameter Tb(comp)      "Temperatura de ebullici�n de los componentes a P=9.5bar [K]"
/
iButene 341.7
Ethanol 421.9
nButene 342.6
ETBE    438.8
/;

*Entalp�a de la fase vapor (kJ/mol) -- Int(CpdT)
variable HVi(N,j,comp,Net),HV(N,j,Net);
equations EqHVi(N,j,comp,Net),EqHV(N,j,Net);
EqHVi(N,j,comp,Net).. HVi(N,j,comp,Net)=e=( (C1c(comp)*(TempV(N,j,Net)-Tref)) + ((C2c(comp)/2)*(sqr(TempV(N,j,Net))-(sqr(Tref))))
                                   + ((C3c(comp)/3)*(power(TempV(N,j,Net),3)-(power(Tref,3)))) + ((C4c(comp)/4)*(power(TempV(N,j,Net),4)-(power(Tref,4))))
                                   + ((C5c(comp)/5)*(power(TempV(N,j,Net),5)-(power(Tref,5)))) + ((C6c(comp)/6)*(power(TempV(N,j,Net),6)-(power(Tref,6)))) + Hform(comp)
                                   + (8.314/1000)*TempV(N,j,Net)*(Z(N,j,Net)-1)+(1+mEOS(comp))*(sqrt(aEOS(N,j,Net))/bEOS(N,j,Net))*log(Z(N,j,Net)/(Z(N,j,Net)+(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*TempV(N,j,Net))))));


EqHV(N,j,Net).. HV(N,j,Net)=e=sum(comp,HVi(N,j,comp,Net)*y(N,j,comp,Net)/100);
*Constantes de entalp�a de vaporizaci�n (kJ/mol)

parameter
C1v(comp)
/iButene        32.614
Ethanol        55.789
nButene        33.774
ETBE        45.29
/
C2v(comp)
/iButene        0.38073
Ethanol        0.31245
nButene        0.5107
ETBE        0.27343
/
C3v(comp)
/iButene        0
Ethanol        0
nButene        -0.17304
ETBE        0.21645
/
C4v(comp)
/iButene        0
Ethanol        0
nButene        0.05181
ETBE        -0.11756
/
C5v(comp)
/iButene        0
Ethanol        0
nButene        0
ETBE        0
/
;

*Temperaturas reducidas

parameter Tred(comp);
Tred(comp)=Tb(comp)/Tcrit(comp);


parameter alphaEOSb(comp), aiEOSb(comp);
alphaEOSb(comp)=sqr(1+mEOS(comp)*(1-sqrt(Tb(comp)/Tcrit(comp))));
aiEOSb(comp)=alphaEOSb(comp)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);
positive variable Zboil(N,j,comp,Net);
equation VaporZb(N,j,comp,Net);
VaporZb(N,j,comp,Net).. power(Zboil(N,j,comp,Net),3)-sqr(Zboil(N,j,comp,Net))+(Zboil(N,j,comp,Net))
                        *((aiEOSb(comp)*P(N,j,Net)/(sqr(0.00008314*Tb(comp))))
                        -(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp)))
                        -sqr(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp))))
                        -((aiEOSb(comp)*P(N,j,Net)/(sqr(0.00008314*Tb(comp)))))
                        *(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp))) =e= 0;

*Entalp�a de vaporizaci�n (kJ/mol)

parameter DHvap(comp), Hvib(comp);
DHVap(comp)=( C1v(comp)*( (1-Tred(comp))**( C2v(comp) + (C3v(comp)*Tred(comp)) + (C4v(comp)*(sqr(Tred(comp)))) + (C5v(comp)*(power(Tred(comp),3)) ) ) ) );
HVib(comp)=( (C1c(comp)*(Tb(comp)-Tref)) + ((C2c(comp)/2)*((sqr(Tb(comp))-sqr(Tref)))) + ((C3c(comp)/3)*((power(Tb(comp),3))-(power(Tref,3)))) + ((C4c(comp)/4)*((power(Tb(comp),4))-(power(Tref,4)))) + ((C5c(comp)/5)*((power(Tb(comp),5))-(power(Tref,5)))) + ((C6c(comp)/6)*((power(Tb(comp),6))-(power(Tref,6)))) + Hform(comp));
variable depHvib(N,j,comp,Net);
equation EqdepHvib(N,j,comp,Net);
EqdepHvib(N,j,comp,Net).. depHvib(N,j,comp,Net) =e= (8.314/1000)*Tb(comp)*(Zboil(N,j,comp,Net)-1)
                                                +(1+mEOS(comp))*(sqrt(aiEOSb(comp))/biEOS(comp))
                                                *log(Zboil(N,j,comp,Net)/(Zboil(N,j,comp,Net)+(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp)))));

*Constantes de Cp (kJ/mol.K) de l�quido

parameter
C1l(comp)
/
iButene         0.08768
Ethanol         0.10264
nButene         0.18205
ETBE            0.11096
/
C2l(comp)
/iButene        0.0002171
Ethanol         -0.00013963
nButene         -0.001611
ETBE            0.00031422
/
C3l(comp)
/iButene        -9.15300E-07
Ethanol         -3.03410E-08
nButene         1.19630E-05
ETBE            1.74800E-07
/
C4l(comp)
/iButene        2.2660E-09
Ethanol         2.0386E-09
nButene         -3.7454E-08
ETBE            0
/
C5l(comp)
/iButene        0
Ethanol         0
nButene         4.5027E-11
ETBE            0
/
;

*Entalp�a de la fase liquida (kJ/mol)
variable HLi(N,j,comp,Net),HL(N,j,Net);
equation EqHLi(N,j,comp,Net),EqHL(N,j,Net);
EqHLi(N,j,comp,Net).. HLi(N,j,comp,Net)=e=HVib(comp)-DHVap(comp)
                        +((C1l(comp)*(TempL(N,j,Net)-Tb(comp))) + ((C2l(comp)/2)*(sqr(TempL(N,j,Net))-(sqr(Tb(comp)))))
                        +((C3l(comp)/3)*(power(TempL(N,j,Net),3)-(power(Tb(comp),3)))) + ((C4l(comp)/4)*(power(TempL(N,j,Net),4)-(power(Tb(comp),4))))
                        +((C5l(comp)/5)*(power(TempL(N,j,Net),5)-(power(Tb(comp),5)))))+depHvib(N,j,comp,Net);

EqHL(N,j,Net).. HL(N,j,Net)=e=sum(comp,HLi(N,j,comp,Net)*x(N,j,comp,Net)/100);

*-------------------------------------------------------------------------------
*                                Secci�n 12
*                  C�lculo de entalp�a de alimentaci�n
*-------------------------------------------------------------------------------
*Entalp�a de la alimentaci�n 1 (Butenos en el caso de ETBE)

parameter HV_b(comp)    "Entalp�a de vapor de la alimentaci�n [kJ/mol]"
          Tred_b(comp)  "Temperatura reducida alimentaci�n [-]"
          DHVap_b(comp) "Entalp�a de vaporizaci�n alimentaci�n [kJ/mol]"
          HL_b(comp)    "Entalp�a de l�quido de la alimentaci�n [kJ/mol]";
HV_b(comp)=( (C1c(comp)*(TaliB-Tref)) + ((C2c(comp)/2)*((TaliB**2)-(Tref**2))) + ((C3c(comp)/3)*((TaliB**3)-(Tref**3))) + ((C4c(comp)/4)*((TaliB**4)-(Tref**4))) + ((C5c(comp)/5)*((TaliB**5)-(Tref**5))) + ((C6c(comp)/6)*((TaliB**6)-(Tref**6))) + Hform(comp));
Tred_b(comp)=TaliB/Tcrit(comp);
DHVap_b(comp)=( C1v(comp)*( (1-Tred_b(comp))**( C2v(comp) + (C3v(comp)*Tred_b(comp)) + (C4v(comp)*(Tred_b(comp)**2)) + (C5v(comp)*(Tred_b(comp)**3)) ) ) );
HL_b(comp)=HV_b(comp)-DHVap_b(comp);
parameter alphaEOSbut(comp), aiEOSbut(comp), aEOSbut(N,j), bEOSbut(N,j);
alphaEOSbut(comp)=(1+mEOS(comp)*(1-(TaliB/Tcrit(comp))**(1/2)))**2;
aiEOSbut(comp)=alphaEOSbut(comp)*0.42747*((0.00008314*TcritSRK(comp))**2)/Pcrit(comp);
bEOSbut(N,j)=sum(comp,(zb(N,j,comp)/100)*biEOS(comp));
aEOSbut(N,j)=sum(comp,sum(comp1, (zb(N,j,comp)/100)*(zb(N,j,comp1)/100)*(aiEOSbut(comp)*aiEOSbut(comp1))**0.5));

*Z alimentaci�n 1 se calcula para todas las etapas internas

positive variable Zbut(N,j,Net);
equation VaporZbut(N,j,Net);
VaporZbut(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. power(Zbut(N,j,Net),3)-sqr(Zbut(N,j,Net))+(Zbut(N,j,Net))
                        *((aEOSbut(N,j)*P(N,j,Net)/(sqr(0.00008314*TaliB)))
                        -(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB))
                        -sqr(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB)))
                        -((aEOSbut(N,j)*P(N,j,Net)/(sqr(0.00008314*TaliB))))
                        *(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB)) =e= 0;

*entalpia alimentacion 1 se calcula para todas las etapas internas

variable HFB(N,j,Net) "Entalp�a de la alimentaci�n de butenos";
equation EqHFB(N,j,Net);
EqHFB(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFB(N,j,Net) =e= sum(comp,(zb(N,j,comp)/100)*(HL_b(comp)+(8.314/1000)*TaliB*(Zbut(N,j,Net)-1)
                        +(1+mEOS(comp))*(sqrt(aEOSbut(N,j))/bEOSbut(N,j))
                        *log(Zbut(N,j,Net)/(Zbut(N,j,Net)+(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB))))));

*Entalpia de la alimentaci�n 2 (Etanol en el caso de ETBE)

parameter HV_e(comp)    "Entalp�a de vapor de la alimentaci�n [kJ/mol]"
          Tred_e(comp)  "Temperatura reducida alimentaci�n [K]"
          DHVap_e(comp) "Entalp�a de vaporizaci�n alimentaci�n [kJ/mol]"
          HL_e(comp)    "Entalp�a de l�quido de la alimentaci�n [kJ/mol]";
HV_e(comp)=( (C1c(comp)*(TaliE-Tref)) + ((C2c(comp)/2)*((TaliE**2)-(Tref**2))) + ((C3c(comp)/3)*((TaliE**3)-(Tref**3))) + ((C4c(comp)/4)*((TaliE**4)-(Tref**4))) + ((C5c(comp)/5)*((TaliE**5)-(Tref**5))) + ((C6c(comp)/6)*((TaliE**6)-(Tref**6))) + Hform(comp));
Tred_e(comp)=TaliE/Tcrit(comp);
DHVap_e(comp)=( C1v(comp)*( (1-Tred_e(comp))**( C2v(comp) + (C3v(comp)*Tred_e(comp)) + (C4v(comp)*(Tred_e(comp)**2)) + (C5v(comp)*(Tred_e(comp)**3)) ) ) );
HL_e(comp)=HV_e(comp)-DHVap_e(comp);
parameter alphaEOSeth(comp), aiEOSeth(comp), aEOSeth, bEOSeth;
alphaEOSeth(comp)=(1+mEOS(comp)*(1-(TaliE/Tcrit(comp))**(1/2)))**2;
aiEOSeth(comp)=alphaEOSeth(comp)*0.42747*((0.00008314*TcritSRK(comp))**2)/Pcrit(comp);
bEOSeth=sum(comp,(ze(comp)/100)*biEOS(comp));
aEOSeth=sum(comp,sum(comp1, (ze(comp)/100)*(ze(comp1)/100)*(aiEOSeth(comp)*aiEOSeth(comp1))**0.5));

*Z alimentaci�n 2 se calcula para todas las etapas internas

positive variable Zeth(N,j,Net);
equation VaporZeth(N,j,Net);
VaporZeth(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. (Zeth(N,j,Net))**3-(Zeth(N,j,Net))**2+(Zeth(N,j,Net))
                        *((aEOSeth*P(N,j,Net)/((0.00008314*TaliE)**2))
                        -(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))
                        -(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))**2)
                        -((aEOSeth*P(N,j,Net)/((0.00008314*TaliE)**2)))
                        *(bEOSeth*P(N,j,Net)/(0.00008314*TaliE)) =e= 0;

*entalpia alimentacion 2 se calcula para todas las etapas internas

variable  HFE(N,j,Net)   "Entalp�a de la alimentaci�n de etanol [kJ/mol]";
equation EqHFE(N,j,Net);
EqHFE(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFE(N,j,Net) =e= sum(comp,(ze(comp)/100)*(HL_e(comp)+(8.314/1000)*TaliE*(Zeth(N,j,Net)-1)
                        +(1+mEOS(comp))*(sqrt(aEOSeth)/bEOSeth)
                        *log(Zeth(N,j,Net)/(Zeth(N,j,Net)+(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))))));
*-------------------------------------------------------------------------------
*                                Secci�n 13
*          Definicion de parametros, restricciones y variables binarias
*-------------------------------------------------------------------------------

*Parametro que determina si las etapas de rxn deben considerarse como etapas de equilibrio
scalar CASE "0 indica que en las etapas de rxn si hay equilibrio"/0/;
*Existencia de catalizador
parameter yc(Net) "1 indica que en la etapa si hay catalizador";

*Existencia de reflujo
parameter yr(Net) "1 indica que en la etapa si hay reflujo";
yr(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=0;
yr('2')=1;
*Existencia de boil up
parameter yb(Net) "1 indica que en la etapa si hay boil up";

*Permite saber si la etapa es real o sobrante
parameter par(Net) "1 indica que la etapa es real fisicamente";
par('1')=1;
par(Net)$(ord(Net) eq card(Net))=1;
*Existencia de relaciones de equilibrio
parameter ye(Net) "1 indica que la etapa es de equilibrio";

*Existencia de alimentacion
parameter yf1(Net) "1 indica que en la etapa hay alimentacion de F";
parameter yf2(Net) "1 indica que en la etapa hay alimentacion de F";
*-------------------------------------------------------------------------------
*                                Secci�n 14
*                           Restricciones logicas
*-------------------------------------------------------------------------------
scalar cmej /1/;
scalar NCmax "numero maximo de etapas reactivas" /3/;

equation logic1(Net) "The boil up stage is below the reflux stage";
logic1(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*(yb(Net));

equation logic2 "There is one reflux stage";
logic2..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yr(Net1)))=e=cmej*1;

equation logic3 "There is one boil up stage";
logic3..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yb(Net1)))=e=cmej*1;

equation logic4a,logic4b "There is one feed stage of EtOH and there is one feed stage of butenes";
logic4a..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yf1(Net1)))=e=cmej*1;
logic4b..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yf2(Net1)))=e=cmej*1;

equation logic6 "There is a maximum number of  catalytic stages";
logic6.. cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yc(Net1))) =e=cmej*NCmax;

equation logic7a(Net),logic7b(Net) "Both feed stages are below the reflux";
logic7a(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*yf1(Net);
logic7b(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*yf2(Net);

equation logic8a(Net),logic8b(Net) "The boil up stage is below the feed stages";
logic8a(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf1(Net1)))=g=cmej*yb(Net);
logic8b(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf2(Net1)))=g=cmej*yb(Net);

equation logic9(Net)  "The EtOH feed is above the butenes feed";
logic9(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf1(Net1)))=g=cmej*yf2(Net);

equation logic10(Net) "The catalytic stages are below the EtOH feed stage";
logic10(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf1(Net1)))=g=cmej*yc(Net);

equation logic11(Net) "The catalytic stages are above the butenes  feed stage";
logic11(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*((sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf2(Net1)))-(yf2(Net)))=l=cmej*(1-yc(Net));

equation logic12(Net) "The catalytic stages are below the reflux stage";
logic12(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*yc(Net);

equation logic13(Net) "The catalytic stages are above the boil up stage";
logic13(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*((sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)))-(yb(Net)))=l=cmej*(1-yc(Net));

*-------------------------------------------------------------------------------
*                                Secci�n 15
*                         Ecuaciones del condensador
*-------------------------------------------------------------------------------
*Condiciones iniciales (operaci�n en estado estable)
equation BalMasaC0,BalMasaParcialC0(comp),SumaC0,EquilibrioC0(comp),BalEnergiaC0;
BalMasaC0.. 0=e=V('1','1','2')-V('1','1','1')*(1+RR('1','1'));
BalMasaParcialC0(comp).. 0=e=V('1','1','2')*y('1','1',comp,'2')-V('1','1','1')*x('1','1',comp,'1')*(1+RR('1','1'));
SumaC0.. sum(comp,y('1','1',comp,'1')-x('1','1',comp,'1'))=e=0;
EquilibrioC0(comp).. yI('1','1',comp,'1')*P('1','1','1')*phiI('1','1',comp,'1')=e=PsatI('1','1',comp,'1')*gammaI('1','1',comp,'1')*xI('1','1',comp,'1');
BalEnergiaC0.. 0=e=V('1','1','2')*HV('1','1','2')-V('1','1','1')*(1+RR('1','1'))*HL('1','1','1')-QC('1','1');

*Flujo de liquido fijo
equation fixedL(N,j);
fixedL(N,j)..L(N,j,'1')=e=0;

*Condicion de equilibrio en el condensador
equation EquilC1(comp),EquilC2(comp),EquilC3,EquilC4;
EquilC1(comp)..yI('1','1',comp,'1')=e=y('1','1',comp,'1');
EquilC2(comp)..xI('1','1',comp,'1')=e=x('1','1',comp,'1');
EquilC3..tempV('1','1','1')=e=tempI('1','1','1');
EquilC4..tempL('1','1','1')=e=tempI('1','1','1');
*-------------------------------------------------------------------------------
*                                Secci�n 16
*                    Ecuaciones de la columna- Modelo Rate-based
*-------------------------------------------------------------------------------

*parameters
parameter termofacVap(N,j,net,comp,comp1) "factor termodinamico vapor";
termofacVap(N,j,net,comp,comp1)$(ord(comp) eq ord(comp1))=1;

*Variables
variable matRvap(N,j,net,comp,comp1) "Elementos de matriz R para el vapor [min/m]";
variable matRliq(N,j,net,comp,comp1) "Elementos de matriz R para el liquido [min/m]";

*positive variable
positive variable kbinVap(N,j,net,comp,comp1) "Coef transferencia de masa binario en el vapor [m/min]";
positive variable kbinLiq(N,j,net,comp,comp1) "Coef transferencia de masa binario en el liquido [m/min]";

*RATE-BASED MODEL: Ecuaciones de definicion

*Definicion de matriz R para comp-1 componentes
equation defmatRvap(N,j,net,compR,compR1,compREF),defmatRvapDiag(N,j,net,compR,compR1,compREF);
defmatRvap(N,j,net,compR,compR1,compREF)$( (ord(Net)>1 and ord(Net)<card(Net)) and (ord(compR) ne ord(compR1)) )..0=e=100*matRvap(N,j,net,compR,compR1)+y(N,j,compR,net)*((1/kbinVap(N,j,net,compR,compR1))-(1/kbinVap(N,j,net,compR,compREF)));
defmatRvapDiag(N,j,net,compR,compR1,compREF)$( (ord(Net)>1 and ord(Net)<card(Net)) and (ord(compR) eq ord(compR1)) )..0=e=100*matRvap(N,j,net,compR,compR1)-(((y(N,j,compR,net))/(kbinVap(N,j,net,compR,compREF)))+sum(comp$(ident(comp) ne ident(compR)),(y(N,j,comp,net))/(kbinVap(N,j,net,compR,comp))));

equation defmatRliq(N,j,net,compR,compR1,compREF),defmatRliqDiag(N,j,net,compR,compR1,compREF);
defmatRliq(N,j,net,compR,compR1,compREF)$( (ord(Net)>1 and ord(Net)<card(Net)) and (ord(compR) ne ord(compR1)) )..0=e=100*matRliq(N,j,net,compR,compR1)+x(N,j,compR,net)*((1/kbinliq(N,j,net,compR,compR1))-(1/kbinliq(N,j,net,compR,compREF)));
defmatRliqDiag(N,j,net,compR,compR1,compREF)$( (ord(Net)>1 and ord(Net)<card(Net)) and (ord(compR) eq ord(compR1)) )..0=e=100*matRliq(N,j,net,compR,compR1)-(((x(N,j,compR,net))/(kbinliq(N,j,net,compR,compREF)))+sum(comp$(ident(comp) ne ident(compR)),(x(N,j,comp,net))/(kbinliq(N,j,net,compR,comp))));

*Definicion de coeficientes de transferencia de masa multicomponente para comp-1 componentes
equation defktransVap(N,j,net,compR,compR1),defktransLiq(N,j,net,compR,compR1);
defktransVap(N,j,net,compR,compR1)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=sum(compR2,ktransVap(N,j,net,compR,compR2)*matRvap(N,j,net,compR2,compR1))-termofacVap(N,j,net,compR,compR1);
defktransLiq(N,j,net,compR,compR1)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=sum(compR2,ktransLiq(N,j,net,compR,compR2)*matRliq(N,j,net,compR2,compR1))-termofacliq(N,j,net,compR,compR1);

*Viscosidad de componentes puros
positive variable muL(N,j,comp,net) "Viscosidad liquido puro [Pa s]";
positive variable muV(N,j,comp,net)"Viscosidad vapor puro a baja presion [Pa s]";

parameter C1MLiq(comp)
/
ETBE        -1.065700E+01
Ethanol        7.875000E+00
nButene        -1.077300E+01
iButene        -1.038500E+01
/;
parameter C2MLiq(comp)
/
ETBE        8.687700E+02
Ethanol        7.819800E+02
nButene        5.916100E+02
iButene        5.995900E+02
/;
parameter C3MLiq(comp)
/
ETBE        0.000000E+00
Ethanol        -3.041800E+00
nButene        0.000000E+00
iButene        -4.608800E-02
/;


parameter C1MVap(comp)
/
ETBE        1.444900E-07
Ethanol        1.061300E-07
nButene        6.974400E-07
iButene        9.098100E-07
/;
parameter C2MVap(comp)
/
ETBE        7.340800E-01
Ethanol        8.066000E-01
nButene        5.462000E-01
iButene        4.928800E-01
/;
parameter C3MVap(comp)
/
ETBE        1.146100E+02
Ethanol        5.270000E+01
nButene        3.052500E+02
iButene        2.600800E+02
/;
equation defmuL(N,j,comp,net),defmuV(N,j,comp,net);
defmuL(N,j,comp,net)$(ord(Net)>1 and ord(Net)<card(Net))..muL(N,j,comp,net)=e=exp(C1MLiq(comp)+(C2MLiq(comp)/TempL(N,j,net))+C3MLiq(comp)*log(TempL(N,j,net)));
defmuV(N,j,comp,net)$(ord(Net)>1 and ord(Net)<card(Net))..muV(N,j,comp,net)=e=(C1MVap(comp)*(TempV(N,j,net)**C2MVap(comp)))/(1+(C3MVap(comp)/TempV(N,j,net)));

*Viscosidad de la mezcla (liquido)
positive variable mumixL(N,j,net) "Viscosidad de la mezcla liquida [Pa s]";
equation defmumixL(N,j,net) "Andrade mixing rule";
defmumixL(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..mumixL(N,j,net)=e=exp(sum(comp,(x(N,j,comp,net)/100)*(log(muL(N,j,comp,net)))));

*Parametros para Difusividades binarias (liquido y vapor)

parameter phiwilke(comp) "Factor de asociacion"
/
ETBE 1
Ethanol 1.5
nButene 1
iButene 1
/;
parameter Vb(comp) "Liquid molar volumen at normal boiling point [m^3/kmol]"
/
ETBE 0.149303
Ethanol 0.0626953
nButene 0.0896971
iButene 0.0895495
/;
parameter Vfuller(comp) "Volumen de difusion molecular [m^3/mol]"
/
ETBE 133.85
Ethanol 51.77
nButene 82.08
iButene 82.08
/;

*Variables y ecuaciones relacionadas con la saturacion del vapor y el liquido

positive variables
y_Bubble(N,j,comp,Net)     "Porcentaje molar en el vapor para el punto de burbuja [%]"
Z_Bubble(N,j,Net) "Factor z en punto de burbuja"
Temp_Bubble(N,j,Net)       "Temperatura del punto de burbuja [K]"
x_dew(N,j,comp,Net)     "Porcentaje molar en el l�qudio para el punto de rocio [%]"
Z_Dew(N,j,Net) "Factor z en punto de rocio"
Temp_Dew(N,j,Net)       "Temperatura del punto de rocio [K]"
;
y_Bubble.up(N,j,comp,Net) = 100;
x_dew.up(N,j,comp,Net) = 100;

positive variables Psat_Bubble(N,j,comp,Net) presi�n de saturaci�n Burbuja (bar);
equations EqPsatBubble(N,j,comp,Net);
EqPsatBubble(N,j,comp,Net).. Psat_Bubble(N,j,comp,Net)=e=exp( C1a(comp) + (C2a(comp)/(Temp_Bubble(N,j,Net)+C3a(comp))) + (C4a(comp)*Temp_Bubble(N,j,Net)) + (C5a(comp)*log(Temp_Bubble(N,j,Net)) + (C6a(comp)*(power(Temp_Bubble(N,j,Net),C7a(comp))))) );

positive variables Psat_Dew(N,j,comp,Net) presi�n de saturaci�n rocio(bar);
equations EqPsatDew(N,j,comp,Net);
EqPsatDew(N,j,comp,Net).. Psat_Dew(N,j,comp,Net)=e=exp( C1a(comp) + (C2a(comp)/(Temp_Dew(N,j,Net)+C3a(comp))) + (C4a(comp)*Temp_Dew(N,j,Net)) + (C5a(comp)*log(Temp_Dew(N,j,Net)) + (C6a(comp)*(power(Temp_Dew(N,j,Net),C7a(comp))))) );

positive variable Tcritm_Dew(N,j,Net);

equation EqTcritmDew(N,j,Net);
EqTcritmDew(N,j,Net).. Tcritm_Dew(N,j,Net) =e= (sqr(sum(comp,(x_Dew(N,j,comp,Net)/100)*Tcrit(comp)/sqrt(Pcrit(comp)))))/(sum(comp,(x_Dew(N,j,comp,Net)/100)*Tcrit(comp)/Pcrit(comp)));

variables tao_nrtl_Bubble(N,j,comp,comp1,Net);
equations Eq_tao_nrtlBubble(N,j,comp,comp1,Net);
Eq_tao_nrtlBubble(N,j,comp,comp1,Net).. tao_nrtl_Bubble(N,j,comp,comp1,Net)=e=a_nrtl(comp,comp1) + (b_nrtl(comp,comp1)/Temp_Bubble(N,j,Net));

variables g_nrtl_Bubble(N,j,comp,comp1,Net);
equations Eq_g_nrtlBubble(N,j,comp,comp1,Net);
Eq_g_nrtlBubble(N,j,comp,comp1,Net).. g_nrtl_Bubble(N,j,comp,comp1,Net)=e=exp( -alfa_nrtl(comp,comp1)*tao_nrtl_Bubble(N,j,comp,comp1,Net));

variables gamma_Bubble(N,j,comp,Net);
equations EqgammaBubble(N,j,comp,Net);
EqgammaBubble(N,j,comp,Net).. gamma_Bubble(N,j,comp,Net)=e=
        exp(sum(comp1,x(N,j,comp1,Net)*tao_nrtl_Bubble(N,j,comp1,comp,Net)*
        g_nrtl_Bubble(N,j,comp1,comp,Net))/sum(comp1,x(N,j,comp1,Net)*
        g_nrtl_Bubble(N,j,comp1,comp,Net))+sum(comp1,x(N,j,comp1,Net)*
        g_nrtl_Bubble(N,j,comp,comp1,Net)/sum(comp2,x(N,j,comp2,Net)*
        g_nrtl_Bubble(N,j,comp2,comp1,Net))*(tao_nrtl_Bubble(N,j,comp,comp1,Net)-
        sum(comp2,x(N,j,comp2,Net)*tao_nrtl_Bubble(N,j,comp2,comp1,Net)*
        g_nrtl_Bubble(N,j,comp2,comp1,Net))/sum(comp3,x(N,j,comp3,Net)*
        g_nrtl_Bubble(N,j,comp3,comp1,Net)))));

variables tao_nrtl_Dew(N,j,comp,comp1,Net);
equations Eq_tao_nrtlDew(N,j,comp,comp1,Net);
Eq_tao_nrtlDew(N,j,comp,comp1,Net).. tao_nrtl_Dew(N,j,comp,comp1,Net)=e=a_nrtl(comp,comp1) + (b_nrtl(comp,comp1)/Temp_Dew(N,j,Net));

variables g_nrtl_Dew(N,j,comp,comp1,Net);
equations Eq_g_nrtlDew(N,j,comp,comp1,Net);
Eq_g_nrtlDew(N,j,comp,comp1,Net).. g_nrtl_Dew(N,j,comp,comp1,Net)=e=exp( -alfa_nrtl(comp,comp1)*tao_nrtl_Dew(N,j,comp,comp1,Net));

variables gamma_Dew(N,j,comp,Net);
equations EqgammaDew(N,j,comp,Net);
EqgammaDew(N,j,comp,Net).. gamma_Dew(N,j,comp,Net)=e=
        exp(sum(comp1,x_Dew(N,j,comp1,Net)*tao_nrtl_Dew(N,j,comp1,comp,Net)*
        g_nrtl_Dew(N,j,comp1,comp,Net))/sum(comp1,x_Dew(N,j,comp1,Net)*
        g_nrtl_Dew(N,j,comp1,comp,Net))+sum(comp1,x_Dew(N,j,comp1,Net)*
        g_nrtl_Dew(N,j,comp,comp1,Net)/sum(comp2,x_Dew(N,j,comp2,Net)*
        g_nrtl_Dew(N,j,comp2,comp1,Net))*(tao_nrtl_Dew(N,j,comp,comp1,Net)-
        sum(comp2,x_Dew(N,j,comp2,Net)*tao_nrtl_Dew(N,j,comp2,comp1,Net)*
        g_nrtl_Dew(N,j,comp2,comp1,Net))/sum(comp3,x_Dew(N,j,comp3,Net)*
        g_nrtl_Dew(N,j,comp3,comp1,Net)))));


positive variable alphaEOS_Bubble(N,j,comp,Net);
equation EqAlphaEOSBubble(N,j,comp,Net);
EqAlphaEOSBubble(N,j,comp,Net).. alphaEOS_Bubble(N,j,comp,Net) =e= sqr(1+mEOS(comp)*(1-sqrt(Temp_Bubble(N,j,Net)/Tcritm(N,j,Net))));


positive variable aiEOS_Bubble(N,j,comp,Net);
equation EqaiEOSBubble(N,j,comp,Net);
EqaiEOSBubble(N,j,comp,Net).. aiEOS_Bubble(N,j,comp,Net) =e= alphaEOS_Bubble(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);


positive variable bEOS_Bubble(N,j,Net);
equation EqbEOSBubble(N,j,Net);
EqbEOSBubble(N,j,Net).. bEOS_Bubble(N,j,Net) =e= sum(comp,(y_Bubble(N,j,comp,Net)/100)*biEOS(comp));


positive variable aEOS_Bubble(N,j,Net);
equation EqaEOSBubble(N,j,Net);
EqaEOSBubble(N,j,Net).. aEOS_Bubble(N,j,Net) =e= sum(comp,sum(comp1, (y_Bubble(N,j,comp,Net)/100)*(y_Bubble(N,j,comp1,Net)/100)*sqrt(aiEOS_Bubble(N,j,comp,Net)*aiEOS_Bubble(N,j,comp1,Net))));


equation VaporZBubble(N,j,Net);
VaporZBubble(N,j,Net).. power(Z_Bubble(N,j,Net),3)-sqr(Z_Bubble(N,j,Net))+(Z_Bubble(N,j,Net))
                *((aEOS_Bubble(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*Temp_Bubble(N,j,Net))))
                -(bEOS_Bubble(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Bubble(N,j,Net)))
                -sqr(bEOS_Bubble(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Bubble(N,j,Net))))
                -((aEOS_Bubble(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*Temp_Bubble(N,j,Net)))))
                *(bEOS_Bubble(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Bubble(N,j,Net))) =e= 0;

positive variable phi_Bubble(N,j,comp,Net);
equation EqPhiBubble(N,j,comp,Net);
EqPhiBubble(N,j,comp,Net).. phi_Bubble(N,j,comp,Net) =e= exp(((Z_Bubble(N,j,Net))-1)*biEOS(comp)/bEOS_Bubble(N,j,Net)
                                        -log((Z_Bubble(N,j,Net))-bEOS_Bubble(N,j,Net))
                                        -(aEOS_Bubble(N,j,Net)/bEOS_Bubble(N,j,Net))
                                        *(2*(sqrt(aiEOS_Bubble(N,j,comp,Net)/aEOS_Bubble(N,j,Net)))
                                        -biEOS(comp)/bEOS_Bubble(N,j,Net))*log(((Z_Bubble(N,j,Net))
                                        -bEOS_Bubble(N,j,Net))/(Z_Bubble(N,j,Net))));

positive variable alphaEOS_Dew(N,j,comp,Net);
equation EqAlphaEOSDew(N,j,comp,Net);
EqAlphaEOSDew(N,j,comp,Net).. alphaEOS_Dew(N,j,comp,Net) =e= sqr(1+mEOS(comp)*(1-sqrt(Temp_Dew(N,j,Net)/Tcritm_Dew(N,j,Net))));


positive variable aiEOS_Dew(N,j,comp,Net);
equation EqaiEOSDew(N,j,comp,Net);
EqaiEOSDew(N,j,comp,Net).. aiEOS_Dew(N,j,comp,Net) =e= alphaEOS_Dew(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);


positive variable aEOS_Dew(N,j,Net);
equation EqaEOSDew(N,j,Net);
EqaEOSDew(N,j,Net).. aEOS_Dew(N,j,Net) =e= sum(comp,sum(comp1, (y(N,j,comp,Net)/100)*(y(N,j,comp1,Net)/100)*sqrt(aiEOS_Dew(N,j,comp,Net)*aiEOS_Dew(N,j,comp1,Net))));

equation VaporZDew(N,j,Net);
VaporZDew(N,j,Net).. power(Z_Dew(N,j,Net),3)-sqr(Z_dEW(N,j,Net))+(Z_Dew(N,j,Net))
                *((aEOS_Dew(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*Temp_Dew(N,j,Net))))
                -(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Dew(N,j,Net)))
                -sqr(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Dew(N,j,Net))))
                -((aEOS_Dew(N,j,Net)*P(N,j,Net)/(sqr(0.00008314*Temp_Dew(N,j,Net)))))
                *(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp_Dew(N,j,Net))) =e= 0;


positive variable phi_Dew(N,j,comp,Net);
equation EqPhiDew(N,j,comp,Net);
EqPhiDew(N,j,comp,Net).. phi_Dew(N,j,comp,Net) =e= exp(((Z_Dew(N,j,Net))-1)*biEOS(comp)/bEOS(N,j,Net)
                                        -log((Z_Dew(N,j,Net))-bEOS(N,j,Net))
                                        -(aEOS_Dew(N,j,Net)/bEOS(N,j,Net))
                                        *(2*(sqrt(aiEOS_Dew(N,j,comp,Net)/aEOS_Dew(N,j,Net)))
                                        -biEOS(comp)/bEOS(N,j,Net))*log(((Z_Dew(N,j,Net))
                                        -bEOS(N,j,Net))/(Z_Dew(N,j,Net))));


*RATE-BASED MODEL: mass balances
equations ValMassParcialVapor0(comp,net,net1), ValmassParcialLiquid0(comp,net),Sumvapor0(net),Sumliquid0(net);
ValMassParcialVapor0(comp,net,net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=BR('1','1')*L('1','1',Net1)*yb(Net)*y('1','1',comp,Net1)+V('1','1',Net+1)*y('1','1',comp,Net+1)-V('1','1',Net)*y('1','1',comp,Net)-100*Nbonita('1','1',comp,net);
ValmassParcialLiquid0(comp,net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=yf1(Net)*FE*ze(comp)+yf2(Net)*FB*zb('1','1',comp)+RR('1','1')*V('1','1','1')*yr(Net)*x('1','1',comp,'1')+L('1','1',Net-1)*x('1','1',comp,Net-1)-L('1','1',Net)*x('1','1',comp,Net)+100*Nbonita('1','1',comp,net)+100*yc(Net)*(Nu(comp)*mcat*Rx('1','1',Net));
Sumvapor0(net)$(ord(Net)>1 and ord(Net)<card(Net))..sum(comp,y('1','1',comp,Net))=e=100;
Sumliquid0(net)$(ord(Net)>1 and ord(Net)<card(Net))..sum(comp,x('1','1',comp,Net))=e=100;

*RATE-BASED MODEL: Mecanismo de transferencia de masa
equations Notrans(comp,net),Trans1Vap(CompR,net),Trans2Liq(CompR,net),Trans3Equil(comp,net),Trans4SumVap(net),Trans5SumLiq(net);

Notrans(comp,net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=(1-ye(net))*Nbonita('1','1',comp,net);
Trans1Vap(compR,net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*(  100*Nbonita('1','1',compR,net)-rhoV('1','1',net)*aI('1','1',net)*sum(compR1,ktransVap('1','1',net,compR,compR1)*(y('1','1',compR1,Net)-yI('1','1',compR1,Net)))-y('1','1',compR,Net)*sum(comp,Nbonita('1','1',comp,net))   );
Trans2Liq(CompR,net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*(  100*Nbonita('1','1',compR,net)-rhoL('1','1',net)*aI('1','1',net)*sum(compR1,ktransLiq('1','1',net,compR,compR1)*(xI('1','1',compR1,Net)-x('1','1',compR1,Net)))-x('1','1',compR,Net)*sum(comp,Nbonita('1','1',comp,net))   );
Trans3Equil(comp,net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(Net)*((yI('1','1',comp,Net)*P('1','1',Net)*phiI('1','1',comp,Net))-(PsatI('1','1',comp,Net)*gammaI('1','1',comp,Net)*xI('1','1',comp,Net)));
Trans4SumVap(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*(sum(comp,yI('1','1',comp,Net))-100);
Trans5SumLiq(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*(sum(comp,xI('1','1',comp,Net))-100);

*RATE-BASED MODEL: Balances de energia
equations EnergyVapor0(net,net1),EnergyLiquid0(net);
EnergyVapor0(net,net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=BR('1','1')*L('1','1',Net1)*yb(Net)*HV('1','1',Net1)+V('1','1',Net+1)*HV('1','1',Net+1)-V('1','1',Net)*HV('1','1',Net)-Ebonita('1','1',net);
EnergyLiquid0(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=yf1(Net)*FE*HFE('1','1',Net)+yf2(Net)*FB*HFB('1','1',Net)+RR('1','1')*V('1','1','1')*yr(Net)*HL('1','1','1')+L('1','1',Net-1)*HL('1','1',Net-1)-L('1','1',Net)*HL('1','1',Net)+ Ebonita('1','1',net);

*RATE-BASED MODEL: Mecanismo de transferencia de energia
equation ENotrans(net),Etrans1(net),Etrans2(net);
ENotrans(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=(1-ye(net))*Ebonita('1','1',net);
Etrans1(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*( TempV('1','1',Net)-Temp_Dew('1','1',Net)  );
Etrans2(net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(net)*( Temp_Bubble('1','1',Net)-TempL('1','1',Net)  );



equations EquilibrioBubble0(comp,Net),EquilibrioDew0(comp,Net) ,SumaDew0(net),SumaBubble0(net);
SumaBubble0(net)..0=e=ye(Net)*(sum(comp,y_Bubble('1','1',comp,Net))-100);
SumaDew0(net)..0=e=ye(Net)*(sum(comp,x_Dew('1','1',comp,Net))-100);
EquilibrioBubble0(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(Net)*((y_Bubble('1','1',comp,Net)*P('1','1',Net)*phi_Bubble('1','1',comp,Net))-(Psat_Bubble('1','1',comp,Net)*gamma_Bubble('1','1',comp,Net)*x('1','1',comp,Net)));
EquilibrioDew0(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(Net)*((y('1','1',comp,Net)*P('1','1',Net)*phi_Dew('1','1',comp,Net))-(Psat_Dew('1','1',comp,Net)*gamma_Dew('1','1',comp,Net)*x_Dew('1','1',comp,Net)));
*-------------------------------------------------------------------------------
*                                Secci�n 17
*                        Ecuaciones del rehervidor
*-------------------------------------------------------------------------------

*Condiciones iniciales (operaci�n en estado estable)
equation BalMasaR0(Net),BalMasaParcialR0(comp,Net),SumaR0(Net),EquilibrioR0(comp,Net),BalEnergiaR0(Net);
BalMasaR0(Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)-L('1','1',Net)*(1+BR('1','1'));
BalMasaParcialR0(comp,Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)*x('1','1',comp,Net-1)-L('1','1',Net)*(x('1','1',comp,Net)+BR('1','1')*y('1','1',comp,Net));
SumaR0(Net)$(ord(Net) eq card(Net)).. sum(comp,y('1','1',comp,Net)-x('1','1',comp,Net))=e=0;
EquilibrioR0(comp,Net)$(ord(Net) eq card(Net))..yI('1','1',comp,Net)*P('1','1',Net)*phiI('1','1',comp,Net)=e=PsatI('1','1',comp,Net)*gammaI('1','1',comp,Net)*xI('1','1',comp,Net);
BalEnergiaR0(Net)$(ord(Net) eq card(Net))..0=e=QR('1','1')+L('1','1',Net-1)*HL('1','1',Net-1)-L('1','1',Net)*HL('1','1',Net)-BR('1','1')*L('1','1',Net)*HV('1','1',Net);

*Variable fija del flujo de vapor en la ultima etapa
equation fixedV(N,j,Net);
fixedV(N,j,Net)$(ord(Net) eq card(Net))..V(N,j,Net)=e=0;

*Condicion de equilibrio en el rehervidor
equation EquilR1(comp,net),EquilR2(comp,net),EquilR3(net),EquilR4(net);
EquilR1(comp,net)$(ord(Net) eq card(Net))..yI('1','1',comp,net)=e=y('1','1',comp,net);
EquilR2(comp,net)$(ord(Net) eq card(Net))..xI('1','1',comp,net)=e=x('1','1',comp,net);
EquilR3(net)$(ord(Net) eq card(Net))..tempV('1','1',net)=e=tempI('1','1',net);
EquilR4(net)$(ord(Net) eq card(Net))..tempL('1','1',net)=e=tempI('1','1',net);

*-------------------------------------------------------------------------------
*                                Secci�n 18
*               Relaciones hidr�ulicas para todas las etapas internas y
*               coeficientes de transferencia de masa binarios
*-------------------------------------------------------------------------------
*Caracteristicas del catalizador
scalar fracvol /0.3/;
scalar fracEnvelop /0.5/;

*Definici�n de velocidad de vapor
positive variables far(N,j,Net) "Factor de areaci�n [-]";
equations Eqfa(N,j,Net);
Eqfa(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(far(N,j,Net))=e=par(Net)*(0.981*exp(-0.411*((V(N,j,Net)/(rhoV(N,j,Net))/hora)*(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)**(0.5))/At));

positive variable hD(N,j,Net)   "Altura del l�quido por encima del divisor [m]";
equations EqhD(N,j,Net);
EqhD(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (hD(N,j,Net))=e=(0.6*(((((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Lw))**(2/3)));

positive variable uhv(N,j,Net) "Velocidad del vapor por los agujeros [m/s]";
equations Equhv(N,j,Net);
Equhv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(uhv(N,j,Net))=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))/hora)/A0);

positive variable unv(N,j,Net) "Velocidad del vapor por el plato [m/s]";
equations Equnv(N,j,Net);
Equnv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*unv(N,j,Net)=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))/hora)/At);

*Definicion de velocidad del liquido
positive variable ul(N,j,Net) "Velocidad del l�quido en el derramadero [m/s]";
equations Equl(N,j,Net);
Equl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*ul(N,j,Net)=e=par(Net)*((L(N,j,Net)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Ad);

*Carga de liquido
positive variable hcl(N,j,Net)  "Altura del l�quido libre en r�gimen de spray [m]"
equation Eqhcl(N,j,Net);
scalar consmach /1e-20/;
Eqhcl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*hcl(N,j,Net)=e=par(Net)*((0.157*(poro**(-0.791))/(1+1.04E-4*(((((L(N,j,Net)+consmach)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Lw)**(-0.59))
                                                        *(poro**(-1.791))))*(da**0.833)
                                                        *(996/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**(0.5*(1-0.91*da/poro)));
positive variable Csbf(N,j,Net);
equation EqCsbf(N,j,Net);
EqCsbf(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(Csbf(N,j,Net))=e=par(Net)*(0.37*(((sqr(da)*sigma(N,j,Net)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)))**0.125)
                                                        *((((rhoV(N,j,Net))*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**0.1)
                                                        *(sqrt(HS/hcl(N,j,Net))));
*Carga de liquido en etapas cataliticas
positive variable Lload(N,j,Net) "carga de liquido en etapas cataliticas [m_s]";
equation eqLload(N,j,Net);
eqLload(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..(1-fracvol)*((3.1415926/4)*(sqr(D)))*Lload(N,j,Net)=e=ul(N,j,net)*Ad ;

*Factor de flujo de vapor en etapas cataliticas
positive variable Ffactor(N,j,Net) "Factor de flujo de vapor para etapas cataliticas [Pa**0.5]"
equation eqFfactor(N,j,Net);
eqFfactor(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (1-fracvol)*(3.1415926/4)*(sqr(D))*(sqrt((rhov(N,j,net))*(sum(comp,(y(N,j,comp,Net)/100)*MW(comp)))*(1/1000)))*(Ffactor(N,j,Net))=e=(V(N,j,net)*(1/60))*(sum(comp,(y(N,j,comp,Net)/100)*MW(comp)*(1/1000)));

*Caida de presion
positive variables DPL(N,j,Net) "Ca�da de presi�n por la presencia de l�quido [bar]";
equations EqDPL(N,j,Net);
EqDPL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPL(N,j,Net))=e=((far(N,j,Net)*9.81*(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)*(hD(N,j,Net)+hw))/100000);

positive variables DPS(N,j,Net) "Ca�da de presi�n debido a la presencia de los agujeros - seco [bar]";
equations EqDPS(N,j,Net);
EqDPS(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPS(N,j,Net))=e=((1/(2*sqr(K0)))*( (((sqr(V(N,j,Net)/(rhoV(N,j,Net))/hora)/A0)) )*((rhoV(N,j,Net))*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)*(1-sqr(poro)))/100000);

positive variable DPq(N,j,Net)      "Ca�da de presi�n en el derramadero [bar]";
equations EqDPq(N,j,Net);
EqDPq(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. DPq(N,j,Net)=e=(1/(100000))*1.62*((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))/(sqr(Lw*hw))*(sqr((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)+sqr((V(N,j,Net)/(rhoV(N,j,Net))/hora)));

positive variables DP(N,j,Net)  "Ca�da de presi�n total [bar]";
positive variable dPcat(N,j,net)    "caida de presion por catalizador en etapas cataliticas [bar]";
equations EqDP(N,j,Net),EqDPR(N,j,Net),EqdPcat(N,j,net),EqP(N,j,Net),EqPC(N,j,Net),EqPR(N,j,Net) "Definici�n de presi�n por etapa [bar]";

EqDPR(N,j,Net)$(ord(Net) eq card(Net)).. DP(N,j,Net)=e=DP(N,j,Net-1);
EqdPcat(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..dPcat(N,j,net)=e=hs*fracEnvelop*(0.001)*(   (5.69228924748553E-06)*((Lload(N,j,Net)*60*60)**3.05308055949085)*((Ffactor(N,j,Net))**7.851695947) + 1.367015225*((Ffactor(N,j,Net))**1.764157687)    );
EqDP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DP(N,j,Net))=e=(DPS(N,j,Net)+DPL(N,j,Net))+yc(Net)*dPcat(N,j,net);
EqP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. P(N,j,Net)=e=P(N,j,Net-1)+par(Net)*DP(N,j,Net);
EqPC(N,j,Net)$(ord(Net) eq 1).. P(N,j,Net)=e=Pop;
EqPR(N,j,Net)$(ord(Net) eq card(Net)).. P(N,j,Net)=e=P(N,j,Net-1);

*Efectos indeseados en la columna
*Downflow flooding (inundaci�n en los derramaderos)
equation DownFlood(N,j,Net);
DownFlood(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=g=((HD(N,j,Net)+((DP(N,j,Net)+DPq(N,j,Net))*100000)
                                                        /(9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))))-(HS))*par(Net);
*Entrainment flooding (inundaci�n por arrastre de l�quido)
equation EntrainFloodV(N,j,Net), EntrainFloodL(N,j,Net);
EntrainFloodV(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..par(Net)*((unv(N,j,Net))-(Csbf(N,j,Net)*sqrt(((((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)))
                                                        /(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))))=l=0;

EntrainFloodL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..par(Net)*((ul(N,j,Net))-((sigma(N,j,Net)*9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))
                                                        /(sqr(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)))**(1/4)))=l=0 ;

*Weeping (lloriqueo)
equation Weep(N,j,Net);
Weep(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. 0=g=(((0.68-0.12)/(sqrt((rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)
                                                                /((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)
                                                                *9.81*far(N,j,Net)*(hw+hd(N,j,Net))))))-(uhv(N,j,Net)))*par(Net);

*Catalyst flooding (inundaci�n del empaque del catalizador)
equation catflood(N,j,net);
catflood(N,j,net)..yc(net)*(dPcat(N,j,net)-(12E-3)*hs*fracEnvelop)=l=0;

*Construcci�n de la columna
equation Size "Tama�o del equipo";
Size.. 1*Htotal =e= 1*((1+Sfactor)*sum(Net$(ord(Net)>1 and ord(Net)<card(Net)),HS*par(Net)));

equation  ammountcat "Espacio disponible para el catalizador";
ammountcat..mcat=l=(fracvol)*((3.1415926/4)*(sqr(D)))*(hs*fracEnvelop)*770;

equation DtoLratio "Relacion entre el diametro y la altura";
DtoLratio..htotal/D=l=20;

*Ecuaciones para el calculo de coeficietnes de transferencia de masa binarios
positive variable unvmax(N,j,Net) "Velocidad superficial maxima para el vapor [m/s]";
equation maxVapVel(N,j,net);
maxVapVel(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..unvmax(N,j,Net)=e=(Csbf(N,j,Net)*sqrt(((((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))-(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)))/(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)));
positive variable Ffa(N,j,net) "Fraccion de flooding";
equation defFfa(N,j,net);
defFfa(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..Ffa(N,j,net)=e=(unv(N,j,Net))/(unvmax(N,j,Net));
positive variable alphaBennett(N,j,net) "Densidad relativa de la espuma";
equation defalphaBennett(N,j,net);
defalphaBennett(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..alphaBennett(N,j,net)=e=exp(-12.55*((unv(N,j,Net)*(sqrt((rhoV(N,j,net)*sum(comp,MW(comp)*(y(N,j,comp,net)/100))*(1/1000))/((rhoL(N,j,net)*sum(comp,MW(comp)*(x(N,j,comp,net)/100))*(1/1000))-(rhoV(N,j,net)*sum(comp,MW(comp)*(y(N,j,comp,net)/100))*(1/1000))))))**0.91));
positive variable hBennett(N,j,net) "Altura de liquido libre de espuma [m]";
equation defhBennett(N,j,net);
defhBennett(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..hBennett(N,j,net)=e=hD(N,j,Net)+hw;
positive variable Fsa(N,j,net) "Factor superficial [kg^0.5/s*m^0.5]";
equation defFsa(N,j,net);
defFsa(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..Fsa(N,j,net)=e=unv(N,j,Net)*(sqrt(rhoV(N,j,net)*sum(comp,MW(comp)*(y(N,j,comp,net)/100))*(1/1000)));
positive variable hfroth(N,j,net) "Altura de region de espuma [m]";
equation defhfroth(N,j,net);
defhfroth(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..hfroth(N,j,net)=e=hBennett(N,j,net)/alphaBennett(N,j,net);

*Definicion de coeficientes de transferencia de masa binarios (kbin) para todos los componentes
equation defkbinvap(N,j,net,comp,comp1),defkbinliq(N,j,net,comp,comp1);
defkbinvap(N,j,net,comp,comp1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(comp) ne ord(comp1)))..0=e= par(net)*(kbinVap(N,j,net,comp,comp1)-((    (0.776+4.57*hw-0.238*(Fsa(N,j,net))+104.8*((L(N,j,Net)*(1/rhoL(N,j,net))*(1/60))/(lw)))/(sqrt((sum(comp3,muV(N,j,comp3,net)*(y(N,j,comp3,net)/100)))/((rhoV(N,j,net)*sum(comp3,MW(comp3)*(y(N,j,comp3,net)/100))*(1/1000))*( (1.013*(10**(-2)))*(TempV.l(N,j,net)**1.75)*((sqrt((MW(comp)+MW(comp1))/(MW(comp)*MW(comp1))))/((P.l(N,j,Net)*(10**5))*(sqr(((Vfuller(comp))**(1/3))+((Vfuller(comp1))**(1/3))))))  ))))      )*((V(N,j,net)*(1/rhoV(N,j,net))*(1/60))/aI(N,j,net))*60));
defkbinliq(N,j,net,comp,comp1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(comp) ne ord(comp1)))..kbinLiq(N,j,net,comp,comp1)=e=19700*(sqrt(((((((1.17282e-16)*(sqrt(phiwilke(comp1)*MW(comp1)))*tempL(N,j,net))/(muL(N,j,comp1,net)*(VB(comp)**0.6)))*(muL(N,j,comp1,net)))**((x(N,j,comp1,net))/(x(N,j,comp,net)+x(N,j,comp1,net))))*(((((1.17282e-16)*(sqrt(phiwilke(comp)*MW(comp)))*tempL(N,j,net))/(muL(N,j,comp,net)*(VB(comp1)**0.6)))*(muL(N,j,comp,net)))**((x(N,j,comp,net))/(x(N,j,comp,net)+x(N,j,comp1,net)))))/(mumixL(N,j,net))))*(0.4*Fsa(N,j,net)+0.17)*((alphaBennett(N,j,net)*At*hfroth(N,j,net))/aI(N,j,net))*60;

*-------------------------------------------------------------------------------
*                                Secci�n 19
*                            Funci�n objetivo
*-------------------------------------------------------------------------------
parameters alfa1 "Peso para pureza de ETBE" /1e4/
           alfa2 "Peso para carga t�rmica del rehervidor" /500/
           alfa3 "Peso para relaci�n de reflujo" /100/
           CostQr "Costo de la carga t�rmica del rehervidor [$/yr]"
           CostQc "Costo de la carga t�rmica del condensador [$/yr]"
           CostB  "Ganacia del ETBE en fondos [$/yr]"
           CostEth"Costo de alimentaci�n de etanol [$/yr]"
           CostBut"Costo de alimentaci�n de butanos [$/yr]"
           year   "Operational hours per year [hr/yr]" /8000/
           CostCat"Costo del catalizador [$/kg]" /7.7/
           AF     "Factor de anualizaci�n (5 a�os, 5% de tasa de inter�s) [1/yr]"
           MS     "Marshall & Swift coefficient" /1050/
           FM     "Material factor (carbon steel)" /1/
           FPres  "Pressure factor (up to 200psi=13.78bar)" /1.15/
           C0     "Costo de inversi�n inicial AF(Cr1+Cc1) [$]" /10000/
           CT     "Costo de inversi�n para etapas [$]"
           Csh    "Costo de inversi�n para coraza [$]"
           Fcol   "Factor de costo de la columna [-]";
Fcol=FM*FPres;
AF=(0.05/(1-1/(1+0.05)**5));
CT=AF*(MS/280)*4.7*Fcol;
Csh=AF*(MS/280)*101.9*(2.18+Fcol);
CostQr=146.8/(hora);
CostQc=24.5/(hora);
CostB=25.3*3600*year/(1000*hora);
CostEth=15*3600*year/(1000*hora);
CostBut=8.25*3600*year/(1000*hora);

variables zobj;

equation Fobj(Net);


Fobj(Net)$(ord(Net) eq card(Net)).. zobj=e=1*((((CostEth*FE+CostBut*FB+(CostQr*Qr('1','1'))+(CostQc*Qc('1','1')))/year)*year)
+(C0+AF*(mcat*CostCat*sum(Net1,yc(Net1))+CT*((D/0.3048)**1.55)*(sum(Net1,HS*par(Net1))/0.3048)+Csh*((D/0.3048)**1.066)*((Htotal/0.3048)**0.802)))
-((((CostB*L('1','1',Net)))/year)*year));


*-------------------------------------------------------------------------------
*                                Secci�n 20
*                            Cotas en las variables
*-------------------------------------------------------------------------------
*bounds
D.up=0.3;
D.lo=0.1;

hw.up=1;
hw.lo=0.0001;

hs.up=1;
hs.lo=0.1;

htotal.up=10;

at.lo=1e-6;
at.up=0.1;

ad.lo=0.0001;
ad.up=0.01;

lw.up=1;

A0.lo=1e-12;
A0.up=0.01;

RR.lo(N,j)=1;
RR.up(N,j)=500;

Qr.up(N,j)=14000;
Qr.lo(N,j)=100;

Qc.up(N,j)=15000;

L.up(N,j,Net)=200;
V.up(N,j,Net)=200;

BR.up(N,j)=50;

P.up(N,j,Net)= 10;
P.lo(N,j,Net)=Pop;

zboil.up(N,j,comp,Net)=1.5;
zboil.lo(N,j,comp,Net)=0.7;

Zbut.up(N,j,Net)= 1.5;
Zbut.lo(N,j,Net)= 0.5;

Zeth.up(N,j,Net)= 1.5;
Zeth.lo(N,j,Net)= 0.5;

PsatI.up(N,j,comp,Net)=100;
Psat_Bubble.up(N,j,comp,Net)=100;
Psat_Dew.up(N,j,comp,Net)=100;

tao_nrtlI.up(N,j,comp,comp1,Net)    =5;
tao_nrtlI.lo(N,j,comp,comp1,Net)    =-5;

g_nrtlI.up(N,j,comp,comp1,Net)=2;
g_nrtlI.lo(N,j,comp,comp1,Net)=0;

gammaI.up(N,j,comp,Net)=50;
gammaI.lo(N,j,comp,Net)=0;

Ketbe.lo(N,j,Net)=0;
Ketbe.up(N,j,Net)=100;

Krate.up(N,j,Net)=1000000000;
Krate.lo(N,j,Net)=-10;
Krate.scale(N,j,Net)=10000000;

Ka.up(N,j,Net)=100;

Rx.up(N,j,Net)=100;
Rx.lo(N,j,Net)=-100000;
Rx.scale(N,j,Net)=100000;


alphaEOS.lo(N,j,comp,Net) = 0.1;
alphaEOS.up(N,j,comp,Net)=5;

alphaEOS_Bubble.lo(N,j,comp,Net) = 0.1;
alphaEOS_Bubble.up(N,j,comp,Net) = 5;
alphaEOS_Dew.lo(N,j,comp,Net) = 0.1;
alphaEOS_Dew.up(N,j,comp,Net) = 5;

aiEOS.up(N,j,comp,Net)=1e-3;
aiEOS_Bubble.up(N,j,comp,Net) = 1e-3;
aiEOS_Dew.up(N,j,comp,Net) = 1e-3;


aiEOS.lo(N,j,comp,Net) = alphaEOS.lo(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);
aiEOS_Bubble.lo(N,j,comp,Net) = alphaEOS_Bubble.lo(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);
aiEOS_Dew.lo(N,j,comp,Net) = alphaEOS_Dew.lo(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);



bEOS.lo(N,j,Net) = smin(comp,biEOS(comp));
bEOS.up(N,j,Net)=0.01;
bEOS_Bubble.lo(N,j,Net) = smin(comp,biEOS(comp));




aEOS.up(N,j,Net)=1e-3;
aEOS_Bubble.up(N,j,Net) = 1e-3;
aEOS_Dew.up(N,j,Net) = 1e-3;

aEOS.lo(N,j,Net) = smin(comp,aiEOS.lo(N,j,comp,Net));
aEOS_Bubble.lo(N,j,Net) = smin(comp,aiEOS_Bubble.lo(N,j,comp,Net));
aEOS_Dew.lo(N,j,Net) = smin(comp,aiEOS_Dew.lo(N,j,comp,Net));


phiI.up(N,j,comp,Net)=2;
phiI.lo(N,j,comp,Net)=0.01;
phi_dew.up(N,j,comp,Net)=2;
phi_dew.lo(N,j,comp,Net)=0.01;
phi_bubble.up(N,j,comp,Net)=2;
phi_bubble.lo(N,j,comp,Net)=0.01;

HVi.up(N,j,comp,Net)=1000;
HVi.lo(N,j,comp,Net)=-1000;

HV.up(N,j,Net)=1000;
HV.lo(N,j,Net)=-1000;

depHvib.lo(N,j,comp,Net)=-10;
depHvib.up(N,j,comp,Net)=10;

HLi.lo(N,j,comp,Net)=-1000;
HLi.up(N,j,comp,Net)=1000;

HL.lo(N,j,Net)=-1000;
HL.up(N,j,Net)=1000;

HFB.lo(N,j,Net)=-50;
HFB.up(N,j,Net)=1;

HFE.lo(N,j,Net)=-500;
HFE.up(N,j,Net)=0;

far.up(N,j,Net)=2;
far.lo(N,j,Net)=0.1;

hD.up(N,j,Net)=0.1;
hD.lo(N,j,Net)=0.0001;

DPL.up(N,j,Net)=0.1;

DPS.up(N,j,Net)=0.01;

DP.up(N,j,Net)=0.1;

DPq.up(N,j,Net)=0.1;

uhv.up(N,j,Net)=10;
uhv.lo(N,j,Net)=0.4;

hcl.up(N,j,Net)=0.1;
hcl.lo(N,j,Net)=1e-6;

Lload.up(N,j,Net)=0.008333;
Ffactor.up(N,j,Net)=3.575;
dPcat.up(N,j,net)=10;

x.lo(N,j,comp,Net)= 1e-11;
x.up(N,j,comp,Net)= 100;

y.lo(N,j,comp,Net)= 0;
y.up(N,j,comp,Net)= 100;

xI.lo(N,j,comp,Net)= 0;
xI.up(N,j,comp,Net)= 100;

yI.lo(N,j,comp,Net)= 0;
yI.up(N,j,comp,Net)= 100;

z.lo(N,j,Net)= 0.5 ;
z.up(N,j,Net)=  1.3  ;

z_bubble.lo(N,j,Net)= 0.5 ;
z_bubble.up(N,j,Net)=  1.3  ;

z_dew.lo(N,j,Net)= 0.5 ;
z_dew.up(N,j,Net)=  1.3  ;

TempI.lo(N,j,Net)= 200;
TempI.up(N,j,Net)= 700;

TempV.lo(N,j,Net)= 200;
TempV.up(N,j,Net)= 700;

TempL.lo(N,j,Net)= 200;
TempL.up(N,j,Net)= 700;

Temp_Bubble.lo(N,j,Net)= 200;
Temp_Bubble.up(N,j,Net)= 700;

Temp_Dew.lo(N,j,Net)= 200;
Temp_Dew.up(N,j,Net)= 700;

Tcritm.up(N,j,Net)=600;
Tcritm.lo(N,j,Net)=smin(comp, Tcrit(comp));
Tcritm_dew.lo(N,j,Net)=smin(comp, Tcrit(comp));

*Bounded flooding variables cause problems with DICOPT:
unv.lo(N,j,Net)=0.01;
unv.up(N,j,Net)=0.89;

ul.up(N,j,Net)=30;
ul.lo(N,j,Net)=0.001;

rho.up(N,j,comp,Net)=25000;
rho.lo(N,j,comp,Net)=8000;
rhoL.lo(N,j,Net)=smin(comp,rho.lo(N,j,comp,Net));
rhoL.up(N,j,Net)=smax(comp,rho.up(N,j,comp,Net));

rhoV.up(N,j,Net)=700;
rhoV.lo(N,j,Net)=Pop/(0.00008314*TempV.up(N,j,Net)*(Z.up(N,j,Net)));

Csbf.up(N,j,Net)=1;
Csbf.lo(N,j,Net)=0.001;

sigma.up(N,j,Net)=0.03;
sigma.lo(N,j,Net)=0.005;

*Nuevas variables (modelo de no equilibrio)
ktransVap.lo(N,j,net,comp,comp1)=-1;
ktransVap.up(N,j,net,comp,comp1)=1;

ktransLiq.lo(N,j,net,comp,comp1)=-1;
ktransLiq.up(N,j,net,comp,comp1)=1;

matRvap.up(N,j,net,compR,compR1)=1000;

matRliq.up(N,j,net,compR,compR1)=1000;

termofacLiq.lo(N,j,net,comp,comp1)=-3;
termofacLiq.up(N,j,net,comp,comp1)=3;

kbinVap.lo(N,j,net,comp,comp1)=0;
kbinVap.up(N,j,net,comp,comp1)=1;

kbinliq.lo(N,j,net,comp,comp1)=0;
kbinliq.up(N,j,net,comp,comp1)=1;

muL.lo(N,j,comp,net)=8e-5;
muL.up(N,j,comp,net)=1.2e-3;


muV.lo(N,j,comp,net)=1e-7;
muV.up(N,j,comp,net)=1e-4;

mumixL.lo(N,j,net)=1e-6;
mumixL.up(N,j,net)=1e-2;

unvmax.lo(N,j,Net)=1e-12;
unvmax.up(N,j,Net)=20;

Ffa.lo(N,j,net)=1e-5;
Ffa.up(N,j,net)=5;

alphaBennett.up(N,j,net)=5;

hBennett.lo(N,j,net)=hw.lo+hD.lo(N,j,Net);
hBennett.up(N,j,net)=hw.up+hD.up(N,j,Net);

Fsa.lo(N,j,net)=1e-7;
Fsa.up(N,j,net)=10;

hfroth.lo(N,j,net)=1e-9;
hfroth.up(N,j,net)=1;

kbinVap.lo(N,j,net,comp,comp1)=1e-8;
kbinVap.up(N,j,net,comp,comp1)=10;

kbinliq.lo(N,j,net,comp,comp1)=1e-6;
kbinliq.up(N,j,net,comp,comp1)=10;


*------INITIAL VALUE OF BINARY TERMS--------------------------------------------
yb(Net)=0;
yc(Net)=0;
yf1(Net)=0;
yf2(Net)=0;

yb('17')=1;
yf1('6')=1;
yf2('15')=1;
yc('6')=1;
yc('11')=1;
yc('15')=1;

par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);

*----------------------MODEL AND SOLUTION---------------------------------------
*the name of the model must be DSDA
model DSDA /all/;
option nlp=conopt;

*--------------------DATA FROM THE MODEL------------------------------------
$embeddedCode Python:
import csv
import os
###User params 1: independent & dependent binar terms---------------------------
Indep={"x1":["yf1",None,"Net"],
"x2":["yf2",None,"Net"],
"x3":["yc",None,"Net"],
"x4":["yb",None,"Net"]}

Dep=["par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));",
"ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);"]

###User params 2: neighborhood--------------------------------------------------
#use: Infinity,Separable,Mflat,Lflat
neighborhood="Infinity"

###User params 3:User defined reformulation,initialization----------------------
###and_ inequality constraints--------------------------------------------------

#use: Automatic,User_defined
Reform_type="User_defined"
Reform=["fex('x1_1')=xvalue_p('x1_1');",
"fex('x2_1')=xvalue_p('x4_1')-xvalue_p('x2_1');",
"fex('x3_1')=xvalue_p('x1_1')+xvalue_p('x3_1'); ",
"fex('x3_2')=xvalue_p('x1_1')+xvalue_p('x3_1')+xvalue_p('x3_2');",
"fex('x3_3')=xvalue_p('x4_1')-xvalue_p('x2_1')-xvalue_p('x3_3');","fex('x4_1')=xvalue_p('x4_1'); "]
Initial=["xvalinit('x1_1')=6;","xvalinit('x2_1')=2;","xvalinit('x3_1')=0; ","xvalinit('x3_2')=5; ","xvalinit('x3_3')=0;  ","xvalinit('x4_1')=17; "]

#use: Automatic,User_defined
Ineq_type="User_defined"
#merge_constraints: 1 if_ user constraints must be merged with automatic constraints.
#merge_constraints: 0 if_ user defined constraints are going to be used alone.
merge_constraints="0"
Ineq=["2-fex('x1_1');","fex('x2_1')-fex('x4_1');","fex('x3_1')+1-fex('x3_2'); ","fex('x3_2')+1-fex('x3_3');","fex('x1_1')-fex('x3_1');","fex('x3_3')-fex('x2_1'); ","fex('x4_1')-21; "]

###user params 3: Decide if_ logical constranits are going to remain in the subProblems
#logic_sub:1 if_ logical constraints appear in the subproblems, gams is allowed to skip the solution stage for_ infeasible neighbors
#logic_sub:0 if_ logical constraints do_ not_ appear in the subproblems, the master problem oversees the logic_
logic_sub="1"

###Create list of external variabs (extvar)-------------------------------------
extvar=list()
for i in Indep:
 Indep[i][1]=str(sum(list(gams.get(Indep[i][0],KeyType.INT,KeyFormat.SKIP))))
 numb=int(float(Indep[i][1]))
 for j in range(1,numb+1):
  extvar.append(i+"_"+str(j))
  Indep[i].append(i+"_"+str(j))

#Total numb of external var
numb=len(extvar)


#print("number of external variables=",numb)
#print("set of external variables=",extvar)


#Total numb of neighbors
if neighborhood=="Infinity":
 num_neigh=(3**numb)-1
elif neighborhood=="Separable":
 num_neigh=2*numb
elif neighborhood=="Mflat":
 num_neigh=numb*(numb+1)
elif neighborhood=="Lflat":
 num_neigh=2**(numb+1)-2

#extvar set_ and_ subsets
Ext="set extvar /"+','.join(extvar)+"/; \n"
Ext_sub=""
for i in Indep:
 Ext_sub=Ext_sub+"set s_"+i+"(extvar) /"+','.join(Indep[i][3:None])+"/; \n"
#print(Ext)
#print(Ext_sub)

neigh="set neigh /d1*d"+str(num_neigh)+"/; \n"
#print(neigh)

#create de csv file_ with the niegihbors
os.chdir('./neighborhood_MATLAB')
file=str(numb)+'_'+neighborhood+'.csv'
with open(file) as csvfile:
 readCSV = csv.reader(csvfile, delimiter=',')
 neigh_values=","+','.join(extvar)+"\n"
 for row in readCSV:
  neigh_values=neigh_values+str(','.join(row))+"\n"
os.chdir('..')

#print(neigh_values)

#get and_ add ordered sets_in list format
for i in Indep:
 Indep[i].append(list(gams.get(Indep[i][2],KeyType.STRING)))

#create aux variables_ for_ reformulation
Ext_aux=""
for i in Indep:
 Ext_aux=Ext_aux+"parameter yaux_"+i+"(extvar,"+Indep[i][2]+");\nyaux_"+i+"(extvar,"+Indep[i][2]+")=0;\n"

#get and_ add value of intependente binary_ terms for_ initializations
for i in Indep:
 Indep[i].append(list(gams.get(Indep[i][0],KeyType.STRING)))

#print("\n",str(Indep))


#Creation of external variables_ functions_ and_ initialization
reform_string=""
initialization_string=""
if Reform_type=="User_defined":
 reform_string='\n'.join(Reform)
 initialization_string='\n'.join(Initial)
if Reform_type=="Automatic":
 cuent=0;
 for i in Indep:
  reform_string=reform_string+"fex(s_"+i+")=xvalue_p(s_"+i+");\n"
  for j in Indep[i][-1]:
   initialization_string=initialization_string+"xvalinit(\'"+extvar[cuent]+"\')="+str(j[0])+";\n"
   cuent=cuent+1

#print(reform_string)
#print(initialization_string)
#Total number of inequality constraints and_ definiton of ineqs
number_user_defined_ineq=len(Ineq)
numb_ineq=numb*2
ineq_string_Automatic=""
ineq_string_User=""
ineq_string=""
cuent1=1;
for i in Indep:
 cuent2=0;
 for j in Indep[i][-1]:
  cuent2=cuent2+1;
  ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1)+"\')=1-fex(\'"+i+"_"+str(cuent2)+"\');\n"
  ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1+1)+"\')=fex(\'"+i+"_"+str(cuent2)+"\')-"+str(len(Indep[i][-2]))+";\n"
  cuent1=cuent1+2
  if cuent2>=2:
   ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1)+"\')=1+fex(\'"+i+"_"+str(cuent2-1)+"\')-fex(\'"+i+"_"+str(cuent2)+"\');\n"
   numb_ineq=numb_ineq+1
   cuent1=cuent1+1

if Ineq_type=="Automatic":
 ineq_string=ineq_string_Automatic
elif Ineq_type=="User_defined" and merge_constraints=="1":
 numb_ineq=numb_ineq+int(number_user_defined_ineq)
 for i in Ineq:
  ineq_string_User=ineq_string_User+"gval(\'g"+str(cuent1)+"\')="+i+"\n"
  cuent1=cuent1+1
 ineq_string=ineq_string_Automatic+ineq_string_User
elif Ineq_type=="User_defined" and merge_constraints=="0":
 numb_ineq=int(number_user_defined_ineq)
 cuent1=1;
 for i in Ineq:
  ineq_string_User=ineq_string_User+"gval(\'g"+str(cuent1)+"\')="+i+"\n"
  cuent1=cuent1+1
 ineq_string=ineq_string_User

Ineq_set="set extineq /g1*g"+str(numb_ineq)+"/; \n"
#print(Ineq_set)
#print(ineq_string)

#Reformulation with external variables_
recalculation_string="\n"
for i in Indep:
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")=0;\n"
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")$(ord("+Indep[i][2]+") eq floor(fex(extvar)))=1-mod(fex(extvar),1);\n"
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")$(ord("+Indep[i][2]+") eq 1+floor(fex(extvar)))=mod(fex(extvar),1);\n\n"
 recalculation_string=recalculation_string+Indep[i][0]+"("+Indep[i][2]+")=sum(s_"+i+",yaux_"+i+"(s_"+i+","+Indep[i][2]+"));\n"
recalculation_string=recalculation_string+'\n'.join(Dep)

#print(recalculation_string)

#treatment of logical constraints
logic_string1=""
logic_string2=""

if logic_sub=="1":
 logic_string1=logic_string1+"maxExecError=100000;"
 logic_string2=logic_string2+"ExecError=0;"
###Writing txts-----------------------------------------------------------------
file = open("neigh_values.csv",'w')
file.write(neigh_values)
file.close()

file=open("param_sets.txt",'w')
file.write(Ext+Ext_sub+neigh+Ineq_set)
file.close()

file=open("param_extaux.txt",'w')
file.write(Ext_aux+initialization_string)
file.close()

file=open("recalculation.txt",'w')
file.write(recalculation_string)
file.close()

file=open("reformulation.txt",'w')
file.write(reform_string)
file.close()

file=open("Inequalities.txt",'w')
file.write(ineq_string)
file.close()

file=open("logics1.txt",'w')
file.write(logic_string1)
file.close()

file=open("logics2.txt",'w')
file.write(logic_string2)
file.close()
$endEmbeddedCode
