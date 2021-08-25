import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.integrate as integrate


# Helper function for converting scientific notation in imported files to python floats
def float_notation(scientific):
    scientific_temp = list(scientific)
    scientific_temp[scientific_temp.index("E")] = "e"
    scientific = "".join(scientific_temp)
    floatified = float(scientific)
    return floatified


# Helper function for duplicating rows
def dup_rows(array, index, num_dups):
    return np.insert(array, [index+1] * num_dups, array[index], axis=0)


# Helper function for duplicating columns
def dup_cols(array, index, num_dups):
    return np.insert(array, [index+1] * num_dups, array[:,[index]], axis=1)


# Imports spectrum responsivity data at various angles to find EL intensity spectrum
spectrums_list = []
n = [1, 2, 3, 4, 5, 6, 7, 8, 9]
for a in n:
    with open('spectrum' + str(a) + '.IRR', 'r') as file:
        spectrums_list.extend([line.strip() for line in file])
        del spectrums_list[(600 * (a-1)):(600 * (a-1) + 2)]
        del spectrums_list[-1:]

# Separates spectrum responsivity data into lists based on angle, and separates wavelengths into another list
spectrum0_list = spectrums_list[0:600]

wavelength_list = []
for a in spectrum0_list:
    wavelength_list.append(float(a[0:6]))

spectrum1_list = spectrums_list[600:1200]
spectrum2_list = spectrums_list[1200:1800]
spectrum3_list = spectrums_list[1800:2400]
spectrum4_list = spectrums_list[2400:3000]
spectrum5_list = spectrums_list[3000:3600]
spectrum6_list = spectrums_list[3600:4200]
spectrum7_list = spectrums_list[4200:4800]
spectrum8_list = spectrums_list[4800:5400]
spectrumall_list = [spectrum0_list, spectrum1_list, spectrum2_list, spectrum3_list, spectrum4_list, spectrum5_list,
                    spectrum6_list, spectrum7_list, spectrum8_list] # Also called Intensity_lambda_theta in MatLab

# Removes wavelengths from lists and converts data into python floats
for spectrum_list in spectrumall_list:
    for a in spectrum_list:
        spectrum_list[spectrum_list.index(a)] = float_notation(a[6:])

# Finds area under the curve (trapezoidal sum) for each spectrum
trapezoid_list = []
for spectrum_list in spectrumall_list:
    trapezoid_list.append(integrate.trapezoid(spectrum_list, wavelength_list))

# Find rho lambda theta (normalized spectrum so the maximum value is one) for each spectrum
rho_lambda_thetaRAW0 = []
rho_lambda_thetaRAW1 = []
rho_lambda_thetaRAW2 = []
rho_lambda_thetaRAW3 = []
rho_lambda_thetaRAW4 = []
rho_lambda_thetaRAW5 = []
rho_lambda_thetaRAW6 = []
rho_lambda_thetaRAW7 = []
rho_lambda_thetaRAW8 = []
rho_lambda_thetaRAWall = [rho_lambda_thetaRAW0, rho_lambda_thetaRAW1, rho_lambda_thetaRAW2, rho_lambda_thetaRAW3,
                          rho_lambda_thetaRAW4, rho_lambda_thetaRAW5, rho_lambda_thetaRAW6, rho_lambda_thetaRAW7,
                          rho_lambda_thetaRAW8]

# Normalizes each spectrum by its trapezoidal sum
a = 0
for spectrum_list in spectrumall_list:
    for element in spectrum_list:
        rho_lambda_thetaRAWall[a].append(100 * element/trapezoid_list[a])
    a += 1

# Plots every spectrum on top of each other
fig1, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
for rho_lambda_thetaRAW in rho_lambda_thetaRAWall:
    ax1.plot(wavelength_list, rho_lambda_thetaRAW)
ax1.title.set_text("Spectrum Plots")
ax1.legend(["0", "10", "20", "30", "40", "50", "60", "70", "80"])
ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("Rho Lambda Theta Raw")

# Creates angles to be used in calculations
theta = np.ndarray.tolist(np.linspace(0, 90 * (np.pi/180), num=10)) # 0 to 90 degrees, 10 degree interval
theta_interp = np.ndarray.tolist(np.linspace(0, 89 * (np.pi/180), num=100)) # 0 to 89 degrees, 100 total values

# Imports "angle.txt" file to process photodiode current at different angles data
with open('angle.txt', 'r') as file:
    angle_txt_list = [line.strip() for line in file]

# Separate lists for data in angle.txt
time_list = []
voltage_list = []
current_list = []
i_pd_list = []

# Finds device area from angle.txt
for a in angle_txt_list:
    if a.find("Device Area (cm^2)") != -1:
        A_LED_index = angle_txt_list.index(a)
A_LED = float_notation(angle_txt_list[A_LED_index][-11:]) * 0.0001

# Deletes all the initial information
del angle_txt_list[0:(angle_txt_list.index("Time (s)\tVoltage (V)\tCurrent (mA)\tV Det (V)")+1)]

# Separates data in angle.txt into separate lists
for a in angle_txt_list:
    b = a[0:a.index("\t")]
    time_list.append(float_notation(b))
    a = a[a.index("\t") + 1:]
    c = a[0:a.index("\t")]
    voltage_list.append(float_notation(c))
    a = a[a.index("\t") + 1:]
    d = a[0:a.index("\t")]
    current_list.append(float_notation(d))
    a = a[a.index("\t") + 1:]
    i_pd_list.append(float_notation(a))

# Converts current measurement data to averages of 10 measurements at each angle
I_PD_averages = []
counter = 0
running_sum = 0
for a in i_pd_list:
    running_sum += a
    counter += 1
    if counter == 10:
        I_PD_averages.append(running_sum / 10)
        running_sum = 0
        counter = 0
I_PD_averages.append(0) # Appends zero for a 90 degree measurement

# Interpolates and then normalizes I_PD
I_PD_interpolated = np.ndarray.tolist(interp1d(theta, I_PD_averages)(theta_interp))
I_PD_interpolated_max = max(I_PD_interpolated)
for a in I_PD_interpolated:
    I_PD_interpolated[I_PD_interpolated.index(a)] = a/I_PD_interpolated_max

# Plot of interpolated I_PD values, graphed against a cosine function
cosine_function = []
theta_interp_degrees = []
for a in theta_interp:
    cosine_function.append(math.cos(a))
    theta_interp_degrees.append(a * 180/np.pi)
ax3.plot(theta_interp_degrees, I_PD_interpolated)
ax3.plot(theta_interp_degrees, cosine_function)
ax3.title.set_text("I_PD(A)")
ax3.legend(["Measured I_PD", "Lambertian"])
ax3.set_xlabel("Degrees")
ax3.set_ylabel("I_PD")

# Adds 600 zeroes to the spectrum list for a 90 degree measurement
spectrum9_list = [0] * 600
spectrumall_list.append(spectrum9_list)

# Interpolates spectrum data, so that each index of intensity_interpolated corresponds to a new angle, ex. 68 deg)
spectrum_array = np.array(spectrumall_list)
intensity_interpolated = np.ndarray.tolist(np.transpose((interp1d(theta, np.transpose(spectrum_array))
                                                            (theta_interp))))

# Normalizes interpolated data by trapezoidal rule
trapezoidall_list = []
for spectrum_list in intensity_interpolated:
    trapezoidall_list.append(integrate.trapezoid(spectrum_list, wavelength_list))
intensity_interpolated_array = np.array(intensity_interpolated)
rho_lambda_theta = np.ndarray.tolist(np.transpose(intensity_interpolated_array) / trapezoidall_list * 100)

# Normalizes by maximum to find peak normalized forward (EL) spectrum
rho_lambda_theta_max = 0
for spectrum_list in rho_lambda_theta:
    for element in spectrum_list:
        if element > rho_lambda_theta_max:
            rho_lambda_theta_max = element
peak_normalized_forward_spectrum = np.ndarray.tolist(np.transpose(np.array(rho_lambda_theta.copy())))
for spectrum_list in peak_normalized_forward_spectrum:
    for element in spectrum_list:
        spectrum_list[spectrum_list.index(element)] = element / rho_lambda_theta_max

# Plots peak normalized forward (EL) spectrum (question here - which dimensions to use?)
for peak_spectrum in peak_normalized_forward_spectrum:
    ax2.plot(wavelength_list, peak_spectrum)
ax2.title.set_text("Peak Normalized Forward (EL) Spectrum")
ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Rho Lambda Theta")



# Imports and processes DC_data.txt
with open('DC_data.dat', 'r') as file:
    dc_data_list = [line.strip() for line in file]

# Separate lists for data in DC_data.dat
voltage_dc_list = []
j_list = []
I_PD_dc_list = []

# Deletes all the initial information
del dc_data_list[0:(dc_data_list.index('V (V)\tJ_source(mA/cm2)\tI_PD (A)\tLuminascence (cd/m2)\tEQE (%)\tPE (lm/W)')+1)]

# Separates data in angle.txt into separate lists
for a in dc_data_list:
    b = a[0:a.index("\t")]
    voltage_dc_list.append(float_notation(b))
    a = a[a.index("\t") + 1:]
    c = a[0:a.index("\t")]
    j_list.append(float_notation(c))
    a = a[a.index("\t") + 1:]
    d = a[0:a.index("\t")]
    I_PD_dc_list.append(float_notation(d))

# J-V (current density - voltage) calculations and graph (MatLab produces something different - looks wrong)
ax4.plot(voltage_dc_list, j_list)
ax4.title.set_text("Current Density vs. Voltage")
ax4.set_xlabel("Voltage (V)")
ax4.set_ylabel("Current Density (mA/cm^2)")
plt.show()

# Current density is multiplied by 10 for unit conversion: mA/cm^2 to A/m^2
for j in j_list:
    j *= 10

# Total current at a given DC voltage calculations
current_dc_list = []
for j in j_list:
    current_dc_list.append(j * A_LED)

# Interpolating I_PD DC data so that there is data at angles other than normal
np.set_printoptions(threshold=np.inf) # Shows full arrays - useful for debugging
I_PD_dc_array = np.vstack(np.asarray(I_PD_dc_list))
I_PD_dc_array = dup_cols(array=I_PD_dc_array, index=0, num_dups=99)
I_PD_interpolated_array = np.hstack(np.asarray(I_PD_interpolated)).reshape(1, 100)
I_PD_interpolated_array = dup_rows(array=I_PD_interpolated_array, index=0, num_dups=160)
I_PD_V_theta = I_PD_dc_array * I_PD_interpolated_array

# Imports the eye's spectrum sensitivity
with open('R2.csv', newline='') as f:
    reader = csv.reader(f)
    R2 = list(reader)
for line in R2:
    del line[1]
    R2[R2.index(line)] = float(line[0])
R2.remove(R2[600]) # 600th index is at 1000nm, but data only goes up to 999nm
R2_array = np.array(R2)

# Imports the photodiode's spectrum sensitivity
with open('f.csv', newline='') as f:
    reader = csv.reader(f)
    F = list(reader)
for line in F:
    del line[1]
    F[F.index(line)] = float(line[0])
F.remove(F[600]) # 600th index is at 1000nm, but data only goes up to 999nm
F_array = np.array(F)

# Variables for calculations
r = 59.5 * 1e-3
A_PD = 12.96 * (1e-3)**2
# A_LED was previously parsed from angle.txt

# Color scheme list
C = ['b', 'r', 'g', 'y', 'm', 'k', 'c', [.2, .5, 1.0], [.2, .6, .2], [1.0, .6, .3], [1.0, .3, .8], [.3, .1, .91],
     [.6, .99, .05], [.99, .13, .1], [.1, .5, .3]]

# External Quantum Efficiency (EQE) calculations and graph
rho_lambda_theta_array = np.transpose(np.array(rho_lambda_theta))
current_dc_array = np.array(current_dc_list)
EQE_num1 = (I_PD_V_theta * r**2) * np.sin(dup_rows(array=np.hstack(np.asarray(theta_interp)).reshape(1, 100), index=0,
                                               num_dups=160))
EQE_num2 = []
EQE_num3 = []
for a in range(len(theta_interp)):
    EQE_num2.append(integrate.trapezoid(rho_lambda_theta_array[a] * wavelength_list, wavelength_list))
    EQE_num3.append(A_PD * integrate.trapezoid((rho_lambda_theta_array[a] * R2_array), wavelength_list))
EQE_num4 = []
for a in range(len(I_PD_dc_list)):
    EQE_num4.append(2 * np.pi * integrate.trapezoid((EQE_num1[a] * EQE_num2 / EQE_num3), theta_interp))
EQE = []
for a in EQE_num4:
    EQE.append(1/100000000 * a / (1240e-9 * current_dc_array[EQE_num4.index(a)]))
for j in j_list:
    j /= 10 # Converts current density back to mA/cm^2
fig2, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5)
ax1.plot(EQE, j_list)
ax1.title.set_text("EQE vs. Current Density")
ax1.set_xlabel("J (mA/cm^2)")
ax1.set_ylabel("EQE (A/W)")
print("EQE (A/W) is: ")
print(EQE)

# Radiance calculations
Radiance_W = []
Radiance_lm = []
for a in range(len(I_PD_dc_list)):
    Radiance_W.append(100 * (I_PD_dc_list[a] * r**2) / (A_PD * integrate.trapezoid((rho_lambda_theta_array[1] * R2_array),
                                                                             wavelength_list)) / A_LED)
for a in Radiance_W:
    Radiance_lm.append(1/100 * 683 * a * integrate.trapezoid((rho_lambda_theta_array[1] * F_array), wavelength_list))
print("Radiance (W) is: ")
print(Radiance_W)
print("Radiance (lm) is: ")
print(Radiance_lm)

# Luminance calculations and graph
Luminance = []
for a in Radiance_W:
    Luminance.append(1/100 * 683 * a * integrate.trapezoid((rho_lambda_theta_array[1] * F_array), wavelength_list))
ax2.plot(Luminance, voltage_dc_list)
ax2.title.set_text("Luminance vs. Voltage")
ax2.set_xlabel("V (V)")
ax2.set_ylabel("Forward Luminance (cd/m^2)")
print("Forward luminance (cd/m^2) is: ")
print(Luminance)

# Power Efficiency calculations and graph
PE_num1 = []
for a in range(len(theta_interp)):
    PE_num1.append(integrate.trapezoid((rho_lambda_theta_array[a] * F_array), wavelength_list))
PE_num2 = EQE_num1
PE_num3 = EQE_num3
PE_num4 = []
for a in range(len(I_PD_dc_list)):
    PE_num4.append(integrate.trapezoid((PE_num2[a] * PE_num1 / PE_num3), theta_interp))
PE = np.array([])
voltage_dc_list[10] = 0.1 # Converts zero value to a very small value
PE = PE_num4 / (current_dc_array * np.array(voltage_dc_list))
PE_array = PE * (2 * np.pi * 683 / 10)
PE_list = np.ndarray.tolist(PE_array)
ax3.plot(PE_list, Luminance)
ax3.title.set_text("Power Efficiency vs. Luminance")
ax3.set_xlabel("Forward Luminance (cd/m^2)")
ax3.set_ylabel("Power Efficiency (lm/W)")
print("Power efficiency (lm/W) is: ")
print(PE_list)

# Current efficacy calculations and graph
CE = np.ndarray.tolist(A_LED * (np.array(Luminance) / current_dc_array / 10))
ax4.plot(CE, Luminance)
ax4.title.set_text("Current Efficacy vs. Luminance")
ax4.set_xlabel("Forward Luminance (cd/m^2)")
ax4.set_ylabel("Current Efficacy (cd/A)")
print("Current efficacy (cd/A) is: ")
print(CE)

# PD power calculations
PD_Power_V_Theta = []
for a in range(len(theta_interp)):
    PD_Power_V_Theta.append((np.transpose(I_PD_V_theta))[a] / integrate.trapezoid((rho_lambda_theta_array[a] * R2_array), wavelength_list))
PD_Power_V_Theta = np.transpose(np.array(PD_Power_V_Theta))
for a in PD_Power_V_Theta:
    ax5.plot(a, theta_interp_degrees)
ax5.title.set_text("Total PD Power at 6.5V")
ax5.set_xlabel("Degrees")
ax5.set_ylabel("PD Power")

# Other printed lists
print("Voltage (V) is: ")
print(voltage_dc_list)
print("Current density (J) is: ")
print(j_list)
plt.show()

# CODE GRAVEYARD: THE PLACE OLD CODE GOES TO DIE
# # External Quantum Efficiency (EQE) calculations
# def eqe(normalized_el_intensity_spectrum, i_pd, area_pd, radius, current, lambda1, lambda2):
#     step1 = integrate.quad(lambda Lambda: normalized_el_intensity_spectrum * Lambda, lambda1, lambda2)
#     # Below step assumes R(Lambda) = C * Lambda
#     step2 = lambda theta: area_pd * integrate.quad(lambda Lambda: normalized_el_intensity_spectrum * constants.c *
#                                                                   Lambda, lambda1, lambda2)
#     step3 = lambda theta: ((i_pd * (radius**2) * math.sin(theta)) * step1) / step2
#     step4 = (2 * math.pi * (integrate.quad(step3, 0, (math.pi / 2)))) / (1240 * current)
#     return step4
#
#
# # Power efficiency calculations
# def pwr(normalized_el_intensity_spectrum, i_pd, area_pd, radius, current, voltage, lambda1, lambda2):
#     # Below step multiplies normalized EL intensity spectrum with the photopic response of the human eye
#     step1 = lambda Lambda, x: normalized_el_intensity_spectrum * \
#                                   (1.019 * math.e**(-285.4 * ((Lambda/1000) - 0.559))**2)
#     step2 = lambda theta: i_pd * (radius**2) * math.sin(theta) * integrate.quad(step1, lambda1, lambda2)
#     step3 = lambda theta: area_pd * integrate.quad(lambda Lambda: normalized_el_intensity_spectrum * constants.c *
#                                                                   Lambda, lambda1, lambda2)
#     step4 = lambda theta: step2 * (1/step3)
#     step5 = (2 * math.pi * 683 * (integrate.quad(step4, 0, (math.pi / 2)))) / (current * voltage)
#     return step5
#
#
# # Current efficacy calculations
# def cur(forward_light_output, current):
#     return forward_light_output / current
#
#
# # Optimizing the i_PD averages in polar coordinates order to produce a function
# r = i_pd_averages
# degrees = [0, 10, 20, 30, 40, 50, 60, 70, 80]
# thetas = []
# for degree in degrees:
#     thetas.append(math.radians(degree))
# rdata = np.array(r)
# thetadata = np.array(thetas)
#
# # Finding Cartesian Coordinates for plotting
# x = []
# y = []
# for degree in degrees:
#     x.append(r[degrees.index(degree)] * math.cos(thetas[degrees.index(degree)]))
#     y.append(r[degrees.index(degree)] * math.sin(thetas[degrees.index(degree)]))
#
# # Centering x data for polar coordinates
# xoffset = sum(x[1:8]) / 8
# xcentered = []
# for xdatum in x:
#     xcentered.append(xdatum - xoffset)
#
# # Converting centered data back into polar coordinates
# rnew = []
# thetanew = []
# for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
#     rnew.append(np.sqrt(xcentered[i]**2 + y[i]**2))
#     thetanew.append(np.arctan(y[i] / xcentered[i]))
#
# # Converting data to numpy arrays
# xdata = np.array(xcentered)
# ydata = np.array(y)
# rnewdata = np.array(rnew)
# thetanewdata = np.array(thetanew)
#
# # Fit with 3rd order polynomial because there is one inflection point
# # p = np.polyfit(degrees, r, 2)
# # x_fit = np.linspace(0, 80)
# # y_fit = np.polyval(p, x_fit)
# # print("Fit line is: y = " + str(p[0]) + "x^2 + " + str(p[1]) + "x + " + str(p[2]))
# #
# # # Final plotting
# # x_out = y_fit * np.cos(np.radians(x_fit))
# # y_out = y_fit * np.sin(np.radians(x_fit))
# # plt.scatter(x, y)
# # plt.plot(x_out, y_out, '--')
# # plt.grid()
# # plt.title('Photodiode Current (i_PD) Output Fit Line')
# # plt.legend(['Fit', 'Original Data'])
# # plt.show()
#
#
# # Model function for an ellipse, used in curve fitting for I_PD
# def ellipse(theta, a, b):
#     return (a * b) / np.sqrt((b * np.cos(theta))**2 + (a * np.sin(theta))**2)
#
# # Fit with ellipse formula in terms of theta
# plt.plot(xdata, ydata, 'b-', label='Original Data')
# popt, pcov = curve_fit(ellipse, thetanewdata, rnewdata)
# xfitdata = []
# yfitdata = []
# degreeplotdata = np.linspace(0, 180)
# thetaplotdata = []
# for degree in degreeplotdata:
#     thetaplotdata.append(math.radians(degree))
# for theta in thetaplotdata:
#     xfitdata.append((ellipse(theta, popt[0], popt[1])) * np.cos(theta))
#     yfitdata.append((ellipse(theta, popt[0], popt[1])) * np.sin(theta))
# plt.plot(xfitdata, yfitdata, 'r-', label='Fit: a=%5.12f, b=%5.3f' % tuple(popt))
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.legend()
# print("Fit line for photodiode is: r = " + str(popt[0]) + " * " + str(popt[1]) + " / sqrt(" + str(popt[1])
#       + " * cos(theta))^2 + (" + str(popt[0]) + " * sin(theta))^2)")
# plt.show()

# def polyfitter(spectrum_list):
#     wavelengths = []
#     responsivities = []
#     for spectrum in spectrum_list:
#         wavelengths.append(float(spectrum[0:6]))
#         responsivities.append(float_notation(spectrum[-11:]))
#     wavelengthdata = np.array(wavelengths)
#     responsivitydata = np.array(responsivities)
#
#     # Fit with 2nd order polynomial
#     pfit = np.polyfit(wavelengthdata, responsivitydata, 4)
#     p = np.poly1d(pfit)
#     xp = np.linspace(400, 1000)
#     _ = plt.plot(wavelengthdata, responsivitydata, '.', xp, p(xp), '-')
#     plt.show()
#     return("Fit line for spectrum is: y = " + str(p[0]) + "x^2 + " + str(p[1]) + "x + " + str(p[2]))
#
# print(polyfitter(spectrum1_list))
# print(polyfitter(spectrum2_list))
# print(polyfitter(spectrum3_list))
# print(polyfitter(spectrum4_list))
# print(polyfitter(spectrum5_list))
# print(polyfitter(spectrum6_list))
# print(polyfitter(spectrum7_list))
# print(polyfitter(spectrum8_list))
# print(polyfitter(spectrum9_list))

