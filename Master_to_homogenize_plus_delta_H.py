#

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys

# Intial candidates for fit, per FU: - thus, the E vs V input data has to be per FU
E0_init = -941.510817926696  # -1882.50963222/2.0 
V0_init = 63.54960592453 #125.8532/2.0 
B0_init = 76.3746233515232 #74.49 
B0_prime_init = 4.05340727164527 #4.15

def BM(x, E0, V0, B0, B0_prime):
        return  E0+ (2.293710449E+17)*(1E-21)*( (9.0/16.0)*(V0*B0) * (  (((V0/x)**(2.0/3.0)-1.0)**3.0)*B0_prime  + ((V0/x)**(2.0/3.0)-1)**2  *  (6.0-4.0*(V0/x)**(2.0/3.0))  ))


def P(V, E0, V0, B0, B0_prime):
    f0=(3.0/2.0)*B0
    f1=((V0/V)**(7.0/3.0))-((V0/V)**(5.0/3.0))
    f2=((V0/V)**(2.0/3.0))-1
    pressure= f0*f1*(1+(3.0/4.0)*(B0_prime-4)*f2)
    return pressure

# Calcite I (Red triangles):
V_not_p_f_unit_C_I, E_not_p_f_unit_C_I = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/calcite_I.dat', skiprows = 1).T

# Calcite II - Trapped (Green triangles):
V_not_p_f_unit_C_II, E_not_p_f_unit_C_II = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/calcite_II_trapped.dat', skiprows = 1).T

# 14 (Empty grey triangles):
V_14_not_p_f_unit, E_14_not_p_f_unit = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/calcite_II.dat', skiprows = 1).T

# Calcite II - Trapped 0.98 (Green triangles):
V_not_p_f_unit_C_II_0_98, E_not_p_f_unit_C_II_0_98 = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/calcite_II_trapped_0_98.dat', skiprows = 1).T


# Calcite II - Trapped 0.87 (Green triangles):
V_not_p_f_unit_C_II_0_87, E_not_p_f_unit_C_II_0_87 = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/calcite_II_trapped_0_87.dat', skiprows = 1).T

# If the data is not per f unit, do this:
nFU_C_I = 2.0
nFU_C_II = 4.0
E_C_I = E_not_p_f_unit_C_I/nFU_C_I
V_C_I = V_not_p_f_unit_C_I/nFU_C_I

E_C_II = E_not_p_f_unit_C_II/nFU_C_II
V_C_II = V_not_p_f_unit_C_II/nFU_C_II

E_14 = E_14_not_p_f_unit/nFU_C_II
V_14 = V_14_not_p_f_unit/nFU_C_II

E_C_II_0_98 = E_not_p_f_unit_C_II_0_98/nFU_C_II
V_C_II_0_98 = V_not_p_f_unit_C_II_0_98/nFU_C_II

E_C_II_0_87 = E_not_p_f_unit_C_II_0_87/nFU_C_II
V_C_II_0_87 = V_not_p_f_unit_C_II_0_87/nFU_C_II


######

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_C_I, pcov_C_I = curve_fit(BM, V_C_I, E_C_I, p0=init_vals)
popt_C_II, pcov_C_II = curve_fit(BM, V_C_II, E_C_II, p0=init_vals)
popt_14, pcov_14 = curve_fit(BM, V_14, E_14, p0=init_vals)

# Linspace for plotting the fitting curves:
V_C_I_lin = np.linspace(V_C_I[0], V_C_I[-1], 100)
V_C_II_lin = np.linspace(V_C_II[0], V_C_II[-1], 100)
V_14_lin = np.linspace(V_14[0], V_14[-1], 100)

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6, p11, p12), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II', "$P\overline{1}$; scan from the neg. mode\nobserved in Calcite II at [1/4, 0, 0]", '$P1$;  scan from neg.mode'), prop=fontP)
plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)

pressures_per_F_unit_C_I = P(V_C_I, *popt_C_I)
print 'CI = ', popt_C_I
output_array_2 = np.vstack((E_not_p_f_unit_C_I, V_not_p_f_unit_C_I, E_C_I, V_C_I, pressures_per_F_unit_C_I)).T
np.savetxt('Volumes_and_pressures.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

pressures_per_F_unit_14 = P(V_14, *popt_14)
print 'C14 = ', popt_14
#output_array_2 = np.vstack((E_not_p_f_unit_14, V_not_p_f_unit_14, E_14, V_14, pressures_per_F_unit_14)).T
#np.savetxt('Volumes_and_pressures_14.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")


plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
#plt.legend()
#plt.show()
#plt.closefig()
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot.pdf', bbox_inches='tight')
plt.figure()
# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
#p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
#p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
#p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
#p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6, p11, p12), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II', "$P\overline{1}$; scan from the neg. mode\nobserved in Calcite II at [1/4, 0, 0]", '$P1$;  scan from neg.mode'), prop=fontP)
plt.legend((p1, p2), ("Calcite I", "BM fit Calcite I"), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
#plt.legend()
#plt.show()
#plt.closefig()
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()
plt.figure()
# Plotting the fitting curves:
#p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
#p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
#p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
#p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6, p11, p12), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II', "$P\overline{1}$; scan from the neg. mode\nobserved in Calcite II at [1/4, 0, 0]", '$P1$;  scan from neg.mode'), prop=fontP)
plt.legend((p3, p4), ('Calcite II ("trapped")', 'BM fit Calcite II ("trapped")'), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
#plt.legend()
#plt.show()
#plt.closefig()
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_II_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()

plt.figure()
# Plotting the fitting curves:
#p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
#p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
#p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
#p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
plt.scatter(V_C_II_0_87, E_C_II_0_87, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
#p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6, p11, p12), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II', "$P\overline{1}$; scan from the neg. mode\nobserved in Calcite II at [1/4, 0, 0]", '$P1$;  scan from neg.mode'), prop=fontP)
plt.legend() #(p3, p4), ('Calcite II ("trapped")', 'BM fit Calcite II ("trapped")'), prop=fontP)

plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$", fontsize=10)
#plt.legend()
plt.ticklabel_format(useOffset=False)
#plt.closefig()
plt.ylim(-941.231, -941.217)
#plt.show()
plt.savefig('calcite_II_0_87_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()


plt.figure()
# Plotting the fitting curves:
#p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
#p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
#p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

# Plotting the scattered points: 
#p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
plt.scatter(V_C_II_0_98, E_C_II_0_98, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
#p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

plt.legend() #p3, 'Calcite II ("trapped")' ) #, prop=fontP)
#plt.legend((p3, p4), ('Calcite II ("trapped")', 'BM fit Calcite II ("trapped")'), prop=fontP)
plt.xlabel('V / Formula unit (Angstrom$^{3}$)')
plt.ylabel('E / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$", fontsize=10)
#plt.legend()
#plt.closefig()
plt.ylim(-941.231, -941.226) # adjust
plt.ticklabel_format(useOffset=False)
#plt.show()
plt.savefig('calcite_II_0_98_summary_better_plot.pdf', bbox_inches='tight')
plt.clf()
plt.cla()

# H = E + PV

H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I
H_14 = E_14 + pressures_per_F_unit_14 * V_14

output_array_3 = np.vstack((E_C_I, V_C_I, pressures_per_F_unit_C_I, H_C_I)).T
np.savetxt('E_V_P_H__C_I.dat', output_array_3, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

output_array_4 = np.vstack((E_14, V_14, pressures_per_F_unit_14, H_14)).T
np.savetxt('E_V_P_H__14.dat', output_array_4, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 



# Fitting Delta_H;

# Quadratic fit of T=T(P):
#c, d, f = np.polyfit(P1, T1, 2)
fitting_C_I = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
fit_C_I = np.poly1d(fitting_C_I)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 1)
fit_14 = np.poly1d(fitting_14)


print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit_C_I = """, fit_C_I
print 'fit_14 = ', fit_14

fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
fit = np.poly1d(fitting)

print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit = """, fit

fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
fit = np.poly1d(fitting)

print """
HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
fit = """, fit


# If we want the Regression coefficcient:
# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

All_in_one_1st_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
All_in_one_2nd_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
All_in_one_3rd_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
All_in_one_4th_degree = polyfit(pressures_per_F_unit_C_I, H_C_I, 4)

print 'All_in_one_1st_degree = ',All_in_one_1st_degree
print 'All_in_one_2nd_degree = ',All_in_one_2nd_degree
print 'All_in_one_3rd_degree = ',All_in_one_3rd_degree
print 'All_in_one_4th_degree = ',All_in_one_4th_degree


# Plotting Delta_H:

# Linspace for plotting the fitting curves:
#V_C_I_lin = np.linspace(V_C_I[0], V_C_I[-1], 100)
#V_C_II_lin = np.linspace(V_C_II[0], V_C_II[-1], 100)
#V_14_lin = np.linspace(V_14[0], V_14[-1], 100)

# Plotting the fitting curves:
#p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='grey', label='BM fit Calcite I' )
#p4, = plt.plot(V_C_II_lin, BM(V_C_II_lin, *popt_C_II), 'k--', label='BM fit Calcite II ("trapped")')
#p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='BM fit Calcite II')

#********* 1st degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 1)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()

EnergyCI, VolumeCI, PressureCI, EnthalpyCI  = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/plus_delta_H/E_V_P_H__C_I.dat', skiprows = 1).T

print 'PressureCI[0] = ', PressureCI[0]
print 'PressureCI[-1] = ', PressureCI[-1]

Energy14, Volume14, Pressure14, Enthalpy14  = np.loadtxt('/home/david/Trabajo/structures/Calcite_I_and_II/PBE__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8/plus_delta_H/E_V_P_H__14.dat', skiprows = 1).T

print 'Pressure14[0] = ', Pressure14[0]
print 'Pressure14[-1] = ', Pressure14[-1]

xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
#p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='linear fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='linear fit')

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)
plt.legend((p1, p11, p5, p55), ("Calcite I", 'linear fit', 'Calcite II', 'linear fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)
plt.ticklabel_format(useOffset=False)

#fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 1)
#fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 1)
print 'fitting = ', fitting
print 'fitting_14 = ', fitting_14
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
ax = fig.add_subplot(111)
ax.annotate('Intersection\nP=5.8630061 GPa\nH = -603.16506544 a.u.', xy=(5.8630061, -603.16506544), xytext=(5.8630061+2.7767, -603.16506544-162.27),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_1st_degree.pdf', bbox_inches='tight')


#********* 2nd degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 2)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 2)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
#p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='2nd degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='2nd degree pol. fit')

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)
plt.legend((p1, p11, p5, p55), ("Calcite I", '2nd degree pol. fit', 'Calcite II', '2nd degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
ax = fig.add_subplot(111)
ax.annotate('Intersection\nP= 1.83135388 GPa\nH = -826.74261077 a.u.', xy=(1.83135388, -826.74261077), xytext=(1.83135388+2.7767, -826.74261077-162.27),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_2nd_degree.pdf', bbox_inches='tight')


#********* 3rd degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 3)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 3)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
#p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='3rd degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='3rd degree pol. fit')

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)
plt.legend((p1, p11, p5, p55), ("Calcite I", '3rd degree pol. fit', 'Calcite II', '3rd degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
ax = fig.add_subplot(111)
ax.annotate('Intersection\nP= 0.8433419 GPa\nH = -887.29940088 a.u.', xy=(0.8433419, -887.29940088), xytext=(0.8433419+2.7767, -887.29940088-162.27),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_3rd_degree.pdf', bbox_inches='tight')

#********* 4th degree:
fitting = np.polyfit(pressures_per_F_unit_C_I, H_C_I, 4)
fit = np.poly1d(fitting)

fitting_14 = np.polyfit(pressures_per_F_unit_14, H_14, 4)
fit_14 = np.poly1d(fitting_14)

fig = plt.figure()
xp_C_I = np.linspace(PressureCI[-1], PressureCI[0], 100)
xp_14 = np.linspace(Pressure14[-1], Pressure14[0], 100)


# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
#p3 = plt.scatter(V_C_II, E_C_II, color='lawngreen', marker="^", label='Calcite II ("trapped")', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p11, = plt.plot(xp_C_I, fit(xp_C_I), "black", label='4th degree pol. fit')
p55, = plt.plot(xp_14, fit_14(xp_14), "grey", label='4th degree pol. fit')

#p11 = plt.scatter(V_y, E_y, color='yellow', marker="^", facecolors='none', label='$P2$;  scan from\nneg.mode at [1/4, 0, 0]', s=100)
#p12 = plt.scatter(V_int_P1, E_int_P1, color='orange', marker="^", facecolors='none', label='$P1$; scan from\nneg.mode', s=100)

fontP = FontProperties()
fontP.set_size('small')

#plt.legend((p1, p2, p3, p4, p5, p6), ("Calcite I", "BM fit Calcite I", 'Calcite II ("trapped")', 'BM fit Calcite II ("trapped")', "Calcite II", 'BM fit Calcite II'), prop=fontP)
plt.legend((p1, p11, p5, p55), ("Calcite I", '4th degree pol. fit', 'Calcite II', '4th degree pol. fit'), prop=fontP)

plt.xlabel('P / Formula unit (GPa)')
plt.ylabel('H / Formula unit (a.u.)')
plt.suptitle("PBE, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.98 - 1.08)$V_{eq}$ and (0.87 - 0.98)$V_{eq}$", fontsize=10)

plt.ticklabel_format(useOffset=False)
crossing_x = np.roots(fitting - fitting_14)
crossing_y = fit(crossing_x)
print 'crossing_x = ', crossing_x
print 'crossing_y = ', crossing_y
ax = fig.add_subplot(111)
ax.annotate('Intersection\nP =  [ 20.99234265+15.3046569j\n       20.99234265-15.3046569j\n       -2.39980734 +0.j\n      0.48016630 +0.j  ]\nH =  [  301.97135381+700.65039446j\n       301.97135381-700.65039446j\n       -1101.92633691 +0.j\n      -910.37471868  +0.j ]', xy=(0.17560462, -929.86990214), xytext=(0.17560462+2.7767, -929.86990214-320.27),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='blue'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_4th_degree.pdf', bbox_inches='tight')
plt.show()

