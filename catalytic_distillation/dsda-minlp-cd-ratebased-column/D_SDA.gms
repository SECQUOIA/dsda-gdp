***Iterations Sets
*iter: Select the maximum value of iterations allowed. If the polyhedron is known
*to be bounded, this should be as large as possible.
set iter "Iterations" /i1*i100/;
*iner1: This corresponds to the set It=In.
set iner1 "Inner iterations 1" /i1/;
*iner2: Select the maximum value of inner iterations in the line search. If iner2
*=/i1/, the line search procedure won't be executed.
set iner2 "Inner iterations 2" /i1*i100/;
$include param_sets.txt

***Neighborhood
table directions(neigh,extvar)
$ondelim
$include neigh_values.csv
$offdelim
;

***Steps defintion
*hstep: So far, the infinity-norm and 2-norms have been compared with an hstep of
*1. Considering other values of hstep is left as future work
scalar hstep "step" /1/;
parameter hiner1(iner1) "step: inner iterations 1";
loop(iner1,
if(ord(iner1) eq 1,
hiner1(iner1)=0;
else
hiner1(iner1)=hiner1(iner1-1)+(hstep/((card(iner1))-1));
);
);
parameter hiner2(iner2) "step: inner iterations 2";
loop(iner2,
if(ord(iner2) eq 1,
hiner2(iner2)=hstep;
else
hiner2(iner2)=hiner2(iner2-1)+hstep;
);
);
parameter extvarhstep(extvar,neigh) "Neighborhood considering hstep";
extvarhstep(extvar,neigh)=hstep*directions(neigh,extvar);

***Iterative process parameters and scalars
*stop1: stopping criterion in S4 for local optimality
scalar stop1 "Main stopping criterion";
*stop2: Stopping criterion that decides if the objective function value
*of a neighbor was propperly calculated (e.g, if the nlp solver finds an infeasible
*solution, the additional convergence procedure is executed).
scalar stop2 "Stopping criterion 2";
*stop3: Stopping criterion for the line search
scalar stop3 "Stopping criterion 3";
*stop4: Stopping criterion that decides if the objective function value
*in the line search was propperly calculated.
scalar stop4 "Stopping criterion 4";
scalar count;
parameter xvalue(iter,extvar) "value of x at each iteration";
parameter xvalueselect(iter,extvar) "selected value of x for the next iteration";
parameter gval(extineq) "value of inequality constraints at x";
parameter dvs(iter,neigh) "f(neigh)-f(x) at the iteration";
parameter dmin(iter) "Minimum value of dvs";
parameter fvalue(iter) "Objective function value";
parameter fplushvalue(iter,neigh) "Objective function to calculate dvs";
parameter fvalueiner(iter,iner2) "Objective function: inner iterations";
parameter selectd(iter,neigh) "ds: Selected direction for line search";
scalar infeasinit "if this parameter is 1, it means that the initialization was infeasible."
parameter mstatS2(iter)"model status: Feasibility of initialization point";
parameter mstatS4(iter,neigh,iner1) "model status: Feasibility of neighbors";
parameter mstatS5(iter,neigh,iner1) "model status: Feasibility of neighnors with convergence procedure" ;
parameter mstatS7(iter,iner2,iner1) "model status: Feasibility in line search";
parameter mstatS8(iter,iner2,iner1) "model status: Feasibility in line search with convergence procedure";
parameter objval "Objective function value";
parameter xvalinit(extvar) "Initialization of external variables";
parameter CPUtime;
scalar CPUtimeactual;
parameter xvalue_p(extvar) "partial value of the external variable";
parameter fex(extvar)"partial value of the f functions";
scalar cuent_gineq "count for inequality constraint";

CPUtimeactual=0;
infeasinit=0;
xvalue(iter,extvar)=0;
xvalueselect(iter,extvar)=0;
gval(extineq)=0;
dvs(iter,neigh)=0;
dmin(iter)=0;
fvalue(iter)=0;
fplushvalue(iter,neigh)=0;
fvalueiner(iter,iner2)=0;
selectd(iter,neigh)=0;
mstatS2(iter)=0;
mstatS4(iter,neigh,iner1)=0;
mstatS5(iter,neigh,iner1)=0;
mstatS7(iter,iner2,iner1)=0;
mstatS8(iter,iner2,iner1)=0;
xvalue_p(extvar)=0;
fex(extvar)=0;
$include param_extaux.txt

***S1. Initialization of external variables
xvalue('i1',extvar)=xvalinit(extvar);


$include logics1.txt

stop1=0;
loop(iter$((stop1 eq 0) and (ord(iter) ne card(iter))),
***S2. Feasibility of external variables
         xvalue_p(extvar)=xvalue(iter,extvar);
$include reformulation.txt
         if(ord(iter) eq 1,
$include Inequalities.txt
                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0,
                 stop1=1;
                 );
         );
***S1. Initialization if continuous variables
         if(ord(iter) eq 1,
         execute_loadpoint "init";
         else
         execute_loadpoint "DLR";
         );
***S2. Feasibilyty of continuous variables
$include recalculation.txt
$include logics2.txt
solve DSDA using nlp minimizing zobj;
mstatS2(iter)=DSDA.modelstat;
         if((mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19) and ord(iter) eq 1,
         stop1=1;
         infeasinit=1;
*deberia sobrar
*         elseif (mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19),
*         stop1=1;
*         fvalue(iter)=zobj.l;
         else
         execute_unload "DLR";
         fvalue(iter)=zobj.l;
                 if(ord(iter) eq 1,
                 execute_unload "DSDA_solution";
                 );
         );

         loop(neigh$(infeasinit ne 1),
***S3. k-Neighborhood
                  xvalue_p(extvar)=xvalue(iter,extvar)+extvarhstep(extvar,neigh);
$include reformulation.txt
$include Inequalities.txt
                          if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0 ,
                          dvs(iter,neigh)=1;
                          else
                          dvs(iter,neigh)=-1;
                          );

                          if(dvs(iter,neigh) le 0,
                          stop2=0;
   loop(iner1$(stop2 eq 0),

                                         if(ord(iner1) eq 1,
***S4. Local optimality: Computation of the objective function value for the neighbors
                                         execute_loadpoint "DLR";
$include recalculation.txt
$include logics2.txt
                                         solve DSDA using nlp minimizing zobj;
                                         mstatS4(iter,neigh,iner1)=DSDA.modelstat;
                                                 if(mstatS4(iter,neigh,iner1) eq 3 or mstatS4(iter,neigh,iner1) eq 4 or mstatS4(iter,neigh,iner1) eq 5 or mstatS4(iter,neigh,iner1) eq 6 or mstatS4(iter,neigh,iner1) eq 11 or mstatS4(iter,neigh,iner1) eq 12 or mstatS4(iter,neigh,iner1) eq 13 or mstatS4(iter,neigh,iner1) eq 14 or mstatS4(iter,neigh,iner1) eq 18 or mstatS4(iter,neigh,iner1) eq 19,
                                                 stop2=0;
                                                 dvs(iter,neigh)=1;
                                                 else
                                                 stop2=1;
                                                 fplushvalue(iter,neigh)=zobj.l;
                                                 dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                 );

                                         else
***S4. Local optimality: Computation of the objective function value for the neighbors with the convergence procedure

                                                 if(ord(iner1) eq 2,
                                                 execute_loadpoint "DLR";
                                                 else
                                                 execute_loadpoint "DLR1";
                                                 );
                                         xvalue_p(extvar)=xvalue(iter,extvar)+((extvarhstep(extvar,neigh))/(hstep))*(hiner1(iner1));
$include reformulation.txt
$include recalculation.txt
$include logics2.txt
                                         solve DSDA using nlp minimizing zobj;
                                         mstatS5(iter,neigh,iner1)=DSDA.modelstat;
                                         execute_unload "DLR1";

                                                 if(ord(iner1) eq card(iner1),
                                                         if(mstatS5(iter,neigh,iner1) eq 3 or mstatS5(iter,neigh,iner1) eq 4 or mstatS5(iter,neigh,iner1) eq 5 or mstatS5(iter,neigh,iner1) eq 6 or mstatS5(iter,neigh,iner1) eq 11 or mstatS5(iter,neigh,iner1) eq 12 or mstatS5(iter,neigh,iner1) eq 13 or mstatS5(iter,neigh,iner1) eq 14 or mstatS5(iter,neigh,iner1) eq 18 or mstatS5(iter,neigh,iner1) eq 19 ,
                                                         dvs(iter,neigh)=1;
                                                         else
                                                         fplushvalue(iter,neigh)=zobj.l;
                                                         dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                         );
                                                 );

                                         );

                                 );

                         );

         );

***S4. Local optimality: Stopping criterion verification
dvs(iter,neigh)=round(dvs(iter,neigh),8);
dmin(iter)=smin((neigh),dvs(iter,neigh));
         if(dmin(iter)  ge 0 ,
         stop1=1;
         );
***S5. Steepest descent: Selection of ds
count=0;
         loop(neigh,
                         if(dvs(iter,neigh) eq dmin(iter) and count eq 0,
                         selectd(iter,neigh)=1;
                         count=count+1;
                         );
         );

stop3=0;
         loop(iner2$(stop1=0 and stop3=0 and infeasinit ne 1),
***S6. Line search: Feasibility of external variables
         stop4=0;
         xvalue_p(extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);
$include reformulation.txt
$include Inequalities.txt
                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0  ,
                 stop4=1;
                 stop3=1;
                 );

                 loop(iner1$(stop4 eq 0),
                         if(ord(iner1) eq 1 and ord(iner2) eq 1 ,
                         execute_loadpoint "DLR";
                         elseif ord(iner1) eq 2 and ord(iner2) eq 1,
                         execute_loadpoint "DLR";
                         else
                         execute_loadpoint "DLR2";
                         );


                           if(ord(iner1) eq 1 ,
***S6. Line search: Computation of the objective function value
$include recalculation.txt
$include logics2.txt
                           solve DSDA using nlp minimizing zobj;
                           mstatS7(iter,iner2,iner1)=DSDA.modelstat;
                           execute_unload "DLR2";

                                 if( mstatS7(iter,iner2,iner1) eq 3 or mstatS7(iter,iner2,iner1) eq 4 or mstatS7(iter,iner2,iner1) eq 5 or mstatS7(iter,iner2,iner1) eq 6 or mstatS7(iter,iner2,iner1) eq 11 or mstatS7(iter,iner2,iner1) eq 12 or mstatS7(iter,iner2,iner1) eq 13 or mstatS7(iter,iner2,iner1) eq 14 or mstatS7(iter,iner2,iner1) eq 18 or mstatS7(iter,iner2,iner1) eq 19  ,
                                 fvalueiner(iter,iner2)=1E+10;
                                 stop4=0;
                                 else
                                 stop4=1;
                                 fvalueiner(iter,iner2)=zobj.l;
                                 );
***S6. Line search: Computation of the objective function value  with the convergence procedure
                           else
                           xvalue_p(extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*(hiner2(iner2)+hiner1(iner1)-hstep);
$include reformulation.txt
$include recalculation.txt
$include logics2.txt
                           solve DSDA using nlp minimizing zobj;
                           mstatS8(iter,iner2,iner1)=DSDA.modelstat;
                           execute_unload "DLR2";

                                 if(ord(iner1) eq card(iner1),
                                 stop4=1;

                                         if(mstatS8(iter,iner2,iner1) eq 3 or mstatS8(iter,iner2,iner1) eq 4 or mstatS8(iter,iner2,iner1) eq 5 or mstatS8(iter,iner2,iner1) eq 6 or mstatS8(iter,iner2,iner1) eq 11 or mstatS8(iter,iner2,iner1) eq 12 or mstatS8(iter,iner2,iner1) eq 13 or mstatS8(iter,iner2,iner1) eq 14 or mstatS8(iter,iner2,iner1) eq 18 or mstatS8(iter,iner2,iner1) eq 19  ,
                                         fvalueiner(iter,iner2)=1E+10;
                                         else
                                         fvalueiner(iter,iner2)=zobj.l;
                                         );

                                 );

                           );
                 );
***S6. Line search: Stopping criterion
                 if(ord(iner2) ne 1 and  fvalueiner(iter,iner2) ge fvalueiner(iter,iner2-1),
                 stop3=1;
                 );
                 if(stop3=0 and stop1=0,
                 objval=zobj.l;
                 execute_unload "DLR";
                 execute_unload "DSDA_solution";
                         loop(extvar,
                          xvalueselect(iter,extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);
                         );
                 );


         );
***Initial value of x for the next iteration
 if(stop1 =0,
         xvalue(iter+1,extvar)= xvalueselect(iter,extvar);
         else
         xvalue(iter+1,extvar)=xvalue(iter,extvar);
         );



);
CPUtime=timeElapsed;
execute_unload "DSDA_solution_INFO" CPUtime,infeasinit;

embeddedCode Python:
import os
os.remove("neigh_values.csv")
os.remove("param_sets.txt")
os.remove("param_extaux.txt")
os.remove("recalculation.txt")
os.remove("reformulation.txt")
os.remove("Inequalities.txt")
os.remove("logics1.txt")
os.remove("logics2.txt")
os.remove("DLR.gdx")
os.remove("DLR1.gdx")
os.remove("DLR2.gdx")
endEmbeddedCode
