# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 22:16:39 2022

@author: Reena Shrestha
"""

# Estimate the Evapotranspiration using Penman-Monteith Method
# Datas extracted from https://power.larc.nasa.gov/data-access-viewer

# Import Libraries
import pandas as pd
import os
import math
from matplotlib import pyplot as plt

# Define Constants
exp = 2.7183      # base of natural logarithm
z_B = 5           # Elevation above sea level, m , Beaumont
z_A = 1118        # Elevation above sea level, m , Amarillo
Gsc = 0.0820      # Solar Constant, MJ per square meter per minute
a = 0.23          # albedo or canopy reflection coefficient
L_B = 30.086153   # Latitude of Beaumont
L_A = 35.202614   # Latitude of Amarillo
sig = 4.903e-9    # Stefan-Boltzmann Constant, MJ/K4/m2/day

# Set working directory and read data
path = "C:/Users/Reena Shrestha/Desktop/Magenta"
os.chdir(path)

# Read data from csv file
fname = 'Data.csv'
A = pd.read_csv(fname)

# Function to calculate Evapotranspiration
def ET(t_max,t_min,R,U2,RH,J,L,z):
    '''t_max = Maximum Temperature in degree Celsius,
    t_min = Minimum Temperature in degree Celsius,
    R = Solar Radiation in Watt per meter square per day,
    U2 = Wind Speed at 2 meter above ground in meter per second,
    RH = Relative Humidity in percentage,
    J = Number of the day in the year,
    L = latitude in decimal degree,
    z = Elevation above sea level in meter'''
 
    T_mean = (t_max+t_min)/2    # Mean Temperature
    Rs = R*0.0864    #Convert Solar Radiation to Megajoules/squaremeter/day
    
    # Numerator for calculation of slope of saturation vapour pressure
    Nsvp = 4098*(0.6108*(exp**((17.27*T_mean)/(T_mean+237.3))))
    # Denominator for calculation of slope of saturation vapour pressure
    Dsvp = (T_mean + 237.3)**2
    SVP = Nsvp/Dsvp  #Slope of saturation vapour pressure
    
    P = 101.3*((293-(0.0065*z))/293)**5.26  # Atmospheric Pressure in kPa
    PC =0.000665*P       # Psychrometric Constant, kPa per degree celsius
    DT = SVP/(SVP+(PC*(1+0.34*U2)))     # Delta Term
    PT = PC/(SVP+(PC*(1+0.34*U2)))    # Psi Term
    TT = (900/(T_mean+273))*U2      # Temperature Term
    
    # Saturation vapour pressure at Maximum temperature
    e_max = 0.6108*(exp**((17.27*t_max)/(t_max + 237.3)))
    #Saturation vapour pressure at Minimum temperature
    e_min = 0.6108*(exp**((17.27*t_min)/(t_min + 237.3)))
    es = (e_max+e_min)/2     #Mean Saturation Vapour Pressure
    ea = (RH/100)*es     #Actual Vapour pressure
    
    # Inverese relative distance Earth-Sun
    dr = 1 + 0.033*math.cos((2*math.pi*J)/365) 
    ds = 0.409*math.sin(((2*math.pi*J)/365)-1.39) #Solar declination
    Lr = (math.pi/180)*L     #Conversion of Latitudes in degrees to radians
    ws =(math.cos((-math.tan(Lr)*math.tan(ds))))**-1     #Sunset hour angle
    
    #Breaking the equation for extraterrestrial radiation 
    T1 = ((24*60)/math.pi)*Gsc*dr
    T2 = ws*math.sin(Lr)*math.sin(ds)
    T3 = math.cos(Lr)*math.cos(ds)*math.sin(ws)
    Ra = T1*(T2 + T3)
    
    Rso = (0.75+2e-5*z)*Ra    #Clear sky Solar Radiation
    Rns = (1-a)*Rs # Net Solar radiation,Megajoules per square meter per day
    #Breaking the equation for Net outgoing long wave solar radiation
    E1 = ((t_max+273.16)**4+(t_min+273.16)**4)/2
    E2 = (0.34-0.14*(ea**(1/2))) 
    E3 = (1.35*(Rs/Rso))-0.35
    Rnl = sig*E1*E2*E3     #Net outgoing long wave solar radiation
    
    Rn = Rns - Rnl   #Net Radiation
    Rng = 0.408*Rn   #Expressing Net Radiation  equivalent of evaporation(mm)
    
    ET_rad = DT*Rng  #Radiation Term
    ET_wind = PT*TT*(es-ea)  #Wind Term
    ET_o = ET_wind + ET_rad  #Final Reference Evapotranspiration, mm/day
    ET_o = round(ET_o,4)
    return(ET_o)

J = A[['Day']] # Extract day in J

# For Beaumont
t_maxB = A[['t_max_B']] # Extract Max Temperature of Beaumont in t_maxB
t_minB = A[['t_min_B']] # Extract Min Temperature of Beaumont in t_minB
R_B = A[['R_B']] # Extract Solar radiation of Beaumont in R_B
U2_B = A[['U2_B']] # Extract Wind Speed at 2m of Beaumont in U2_B
RH_B = A[['RH_B']] # Extract Relative Humidity of Beaumont in RH_B

# For Amarillo
t_maxA = A[['t_max_A']] # Extract Max Temperature of Amarillo in t_maxA
t_minA = A[['t_min_A']] # Extract Min Temperature of Amarillo in t_minA
R_A = A[['R_A']] # Extract Solar radiation of Amarillo in R_A
U2_A = A[['U2_A']] # Extract Wind Speed at 2m of Amarillo in U2_A
RH_A = A[['RH_A']] # Extract Relative Humidity of Amarillo in RH_A


ETo_B = []
ETo_A = []

idx = len(J)
idx = range(0,idx,1) # creating a range for the loop
# Loop through all the variables for each city for the given year
for z in idx:
    J1 = float(J.iloc[z])
    
    # For Beaumont
    tmaxB = float(t_maxB.iloc[z])
    tminB = float(t_minB.iloc[z])
    RB = float(R_B.iloc[z])
    UB = float(U2_B.iloc[z])
    RHB = float(RH_B.iloc[z])
    ET1 = ET(tmaxB,tminB,RB,UB,RHB,J1,L_B,z_B)
    ETo_B.append(ET1)
    
    # For Amarillo
    tmaxA = float(t_maxA.iloc[z])
    tminA = float(t_minA.iloc[z])
    RA = float(R_A.iloc[z])
    UA = float(U2_A.iloc[z])
    RHA = float(RH_A.iloc[z])
    ET2 = ET(tmaxA,tminA,RA,UA,RHA,J1,L_A,z_A)
    ETo_A.append(ET2)

ETo_B
ETo_A

# Append the results to the dataframe
A['ETo_B'] = ETo_B
A['ETo_A'] = ETo_A

# Store the results to a csv file
A.to_csv('FinalData.csv')

D = A.Day


# Create a plot 
fig,axs = plt.subplots(2,sharex=True) # Two plots sharing x-axis
axs[0].grid() #Add grid to plot 1
axs[1].grid() #Add grid to plot 2
axs[0].plot(D,ETo_B,color='r',marker='o',linestyle='none') #Plot1
axs[1].plot(D,ETo_A,color='b',marker='o',linestyle='none') #Plot2
axs[0].set_title('For Beaumont')
axs[1].set_title('For Amarillo')
plt.xlabel('Days') # Add X-axis label(Same for both plots)
plt.ylabel('Evapotranspiration(mm/day)') #Add Y-axis label(Same for both plots)
axs[0].legend(['ETo_B'],loc='upper right')
axs[1].legend(['ETo_A'],loc='upper right')
plt.show()


import PySimpleGUI as sg   

sg.theme('DarkBlue4')

# Define the window's contents
layout = [[sg.Text("Enter max Temperature in Celsius:"),sg.Input(key='-IN1-')],
          [sg.Text("Enter min Temperature in Celsius :"),sg.Input(key="-IN2-")],
          [sg.Text("Enter Solar Radiation:"),sg.Input(key="-IN3-")],
          [sg.Text("Enter Wind Speed:"),sg.Input(key="-IN4-")],
          [sg.Text("Enter Relative Humidity:"),sg.Input(key="-IN5-")],
          [sg.Text("Enter Elevation Above Sea Level:"),sg.Input(key="-IN6-")],
          [sg.Text("Enter Latitude:"),sg.Input(key="-IN7-")],
          [sg.Text("Enter Julian:"),sg.Input(key="-IN8-")],            
          [sg.Text("The Evapotranspiration is:",key="-OUT-")],
          [sg.Button("Calculate"), sg.Exit()]]

# Create the window with title and layout
window = sg.Window("Estimate Evapotranspiration", layout )

# Create an event loop
while True:
    event, values = window.read()
    if event == 'Calculate':
        t_max = float(values["-IN1-"])
        t_min = float(values["-IN2-"])
        R = float(values["-IN3-"])
        U2 = float(values["-IN4-"])
        RH = float(values["-IN5-"])
        z = float(values["-IN6-"])
        L = float(values["-IN7-"])
        J = float(values["-IN8-"])
        ETo = ET(t_max,t_min,R,U2,RH,J,L,z)
        strc = "The Evapotranspiration is:" + str(ETo) + "mm/day"
        window["-OUT-"].update(strc)
        print(event, values)
        
# Close the window if user closes window or press the Exit
    if event == sg.WIN_CLOSED or event == "Exit":
        break
 
window.close()
