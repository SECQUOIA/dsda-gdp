"""Distillation column model for 2018 PSE conference"""

from __future__ import division

from math import fabs

from pyomo.environ import SolverFactory, Suffix, value
from pyomo.util.infeasible import log_infeasible_constraints

from column import build_column



def main():

#---solve  (write DSDA instead)
    m = build_column(
        min_trays=8,
        max_trays=17,
        xD=0.95,
        xB=0.95,
        # x_input=[16,2], # Original initialization
        x_input=[13,4], # To find optimal solution in PSE paper
        provide_init=False,
        init={})

#---display results
    display_column(m[0])



#-----other functions------------------
def display_column(m):
    print('Objective: %s' % value(m.obj))
    print('Qc: {: >3.0f}kW  DB: {: >3.0f} DT: {: >3.0f} dis: {: >3.0f}'
          .format(value(m.Qc * 1E3),
                  value(m.D['benzene']),
                  value(m.D['toluene']),
                  value(m.dis)))
    for t in reversed(list(m.trays)):
        print('T{: >2.0f}-{:1.0g} T: {: >3.0f} '
              'F: {: >4.0f} '
              'L: {: >4.0f} V: {: >4.0f} '
              'xB: {: >3.0f} xT: {: >3.0f} yB: {: >3.0f} yT: {: >3.0f}'
              .format(t,
                      fabs(value(m.tray[t].indicator_var))
                      if t in m.conditional_trays else 1,
                      value(m.T[t]) - 273.15,
                      value(sum(m.feed[c] for c in m.comps))
                      if t == m.feed_tray else 0,
                      value(m.liq[t]),
                      value(m.vap[t]),
                      value(m.x['benzene', t]) * 100,
                      value(m.x['toluene', t]) * 100,
                      value(m.y['benzene', t]) * 100,
                      value(m.y['toluene', t]) * 100
                      ))
    print('Qb: {: >3.0f}kW  BB: {: > 3.0f} BT: {: >3.0f} bot: {: >3.0f}'
          .format(value(m.Qb * 1E3),
                  value(m.B['benzene']),
                  value(m.B['toluene']),
                  value(m.bot)))
    for t in reversed(list(m.trays)):
        print('T{: >2.0f}-{:1.0g} '
              'FB: {: >3.0f} FT: {: >3.0f} '
              'LB: {: >4.0f} LT: {: >4.0f} VB: {: >4.0f} VT: {: >4.0f}'
              .format(t,
                      fabs(value(m.tray[t].indicator_var))
                      if t in m.conditional_trays else 1,
                      value(m.feed['benzene']) if t == m.feed_tray else 0,
                      value(m.feed['toluene']) if t == m.feed_tray else 0,
                      value(m.L['benzene', t]),
                      value(m.L['toluene', t]),
                      value(m.V['benzene', t]),
                      value(m.V['toluene', t])
                      ))
    print('RF: {: >3.2f} RB: {: >3.2f}'
          .format(value(m.reflux_frac / (1 - m.reflux_frac)),
                  value(m.boilup_frac / (1 - m.boilup_frac))))
# Show value of Boolean variables
    for k in m.conditional_trays:
       print(str(m.tray[k].indicator_var)+'='+str(m.tray[k].indicator_var.value))
          

if __name__ == "__main__":
    m = main()


