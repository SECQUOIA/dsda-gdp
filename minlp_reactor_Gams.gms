
*------------------------SETS---------------------------------------------------
set N "Superstructure units" /1*5/;
set I "Chemical components" /A,B/;
alias(N,N1);

*-----------------------PARAMETERS AND SCALARS----------------------------------
***Kinetic parameters
scalar m1 "Partial order of reaction 1 []" /1/;
scalar m2 "Partial order of reaction 2 []" /1/;
scalar k "kinetic constant [L/(mol*s)]"/2/;
***Feed parameters
scalar QF0 "Inlet volumetric flow [L/s]" /1/;
parameter C0(i) "Initial concentration of reagents [mol/L]"
/
A 0.99
B 0.01
/
;
parameter F0(i) "Inlet molar flow [mol/s]";
F0(i)=C0(i)*QF0;

*-----------------------BINARY VARIABLES AND CONSTRAINTS------------------------
***Independent binary variables
binary variable yf(n) "Location of the last reactor (feed reactor) in the superstructure";
binary variable yr(n) "Location of the recycle flow";
***Dependent binary variables
positive variable yp(n) "Unit operation. 1: Reaction. 0: Simple input-output";
equation defyp(n) "Definition of yp(n) in terms of the independent binary terms";
defyp(n)..yp(n)=e=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));
***Logical constraints
equations logic1 "There is only one 'final' reactor in the superstructure";
equation logic2   "There is only one recycle location";
equation logic3(n) "The recycle can be located at n if there is a reaction operation";
logic1..sum(n,yf(n))=e=1;
logic2..sum(n,yr(n))=e=1;
logic3(n)..yr(n)=l=yp(n);

*-----------------------REAL VARIABLES------------------------------------------
***Network variables
positive variable Q(n) "Outlet flow rate of the superstricture unit [L/s]";
positive variable F(i,n) "Molar flow [mol/s]";
variable  rate(i,n) "Reaction rate [mol/(L*s)]";
positive variable V(n) "Reactor volume [L]";
***Splitter variables
positive variable QR "Recycle flow rate  [L/s]";
positive variable QP "Product flow rate  [L/s]";
positive variable R(i) "Recycle molar flow [mol/s]";
positive variable P(i) "Product molar flow [mol/s]";

*-----------------------SUPERSTRUCTURE CONSTRAINTS------------------------------

***Kinetic constraints
equation net_rate(i,n) "Reactor Network: Reaction rates";
net_rate(i,n)..(rate('A',n)*((Q(n))**m1)*((Q(n))**m2)+k*((F('A',n))**m1)*((F('B',n))**m2))$(ord(i) eq 1)+(rate('B',n)+rate('A',n))$(ord(i) eq 2)=e=0;

***Network constraints
equation net_comp(i,n) "Reactor Network: Component molar balance";
equation net_cont(n) "Reactor Network: continuity equation";
net_comp(i,n)$(ord(n) ne card(n))..F(i,n+1)+yr(n)*R(i)-F(i,n)+yp(n)*rate(i,n)*V(n)=e=0;
net_cont(n)$(ord(n) ne card(n))..Q(n+1)+yr(n)*QR-Q(n)=e=0;


***Feed unit constraints
equation feed_comp(i,n) "Feed unit: Component molar balance";
equation feed_cont(n) "Feed unit: continuity equation";
feed_comp(i,n)$(ord(n) eq card(n))..F0(i)+yr(n)*R(i)-F(i,n)+yp(n)*rate(i,n)*V(n)=e=0;
feed_cont(n)$(ord(n) eq card(n))..QF0+yr(n)*QR-Q(n)=e=0;

***Splitter constraints
equation spl_comp(i) "Splitter: Component molar balance";
equation spl_cont "Splitter: continuity equation";
equation spl_ad(i) "Splitter: additional splitter constraints";
spl_comp(i)..F(i,'1')-P(i)-R(i)=e=0;
spl_cont..Q('1')-QP-QR=e=0;
spl_ad(i)..P(i)*Q('1')-F(i,'1')*QP=e=0;

*---------------------VOLUMEN CONSTRAINT----------------------------------------
equation eqvol(n);
eqvol(n)$(ord(n) ne 1)..V(n)=e=V(n-1);

*---------------------PRODUCT QUALITY CONSTRAINT--------------------------------
equation qual;
qual..QP*0.95=e=P('B');
*----------------------OBJECTIVE FUNCTION---------------------------------------

variables zobj;

equation Fobj;
Fobj..zobj=e=sum(n,V(n)*yp(n));

*-----------------------BOUNDS ON VARIABLES-------------------------------------
Q.up(n)=10;
QR.up=10;
QP.up=10;
R.up(i)=10;
P.up(i)=10;
F.up(i,n)=10;
rate.lo(i,n)=-10;
rate.up(i,n)=10;
V.up(n)=10;
*----------------------VARIABLE INITIALIZATION----------------------------------
yf.l(n)=0;
yr.l(n)=0;
yp.l(n)=0;
Q.l(n)=0;
F.l(i,n)=0;
rate.l(i,n)=0;
V.l(n)=0;
QR.l=0;
QP.l=0;
R.l(i)=0;
P.l(i)=0;
*----------------------SOLUTION-------------------------------------------------
model SUP_CSTR /all/;
SUP_CSTR.optcr=0;
option reslim = 18000;
option minlp=baron;
option nlp = conopt;
option threads=0;
$onecho > baron.opt
EpsR 0.0001
$offecho
yf.l('1') = 1;
yr.l('1') = 1;
SUP_CSTR.OptFile = 1;
option sysout = on ;
solve SUP_CSTR using minlp minimizing zobj;
execute_unload "MINLP"





