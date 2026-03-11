# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 19:35:38 2025

Complete Code for Computng Exercise 2. Utilises Fresnel equation to generate 
diffraction patterns in varying conditions using both dblquad and Monte Carlo

@author: james
"""

#first import necessarry libraries
from scipy import integrate #integration function
import numpy as np #numerical python
import matplotlib.pyplot as plt #plotting software
from scipy.constants import epsilon_0 as e0 #import exact value of Permittivity of Free Space
from scipy.constants import c #import exact value of Speed of light
import sys #allows user to exit program if they want to
import time #allows the user to see how much time it took for the integration to run

"""a) Define our functions"""

#First our real component
def Fresnel2dreal (yp, xp, y, x, k, z):
    """ 
    This function inputs general coordinates for the screen (x, y) and the apeture (xp, yp), and combines
    them with the screen aperture seperation (z) and k, a function of the wavelength of light used.
    The function then outputs the kernel of the integral of the real component of the Fresnel equation, 
    which cna be integrated to find the real component of the field (E).
    """
    E0 = 1 #V/m
    Fresnel2dreal = ((k * E0)/(2 * np.pi * z)) * np.cos( (k)/(2 * z) * ( (x - xp)**2 + (y - yp)**2 ) ) #V/m
    
    return Fresnel2dreal

#Now our imaginary component
def Fresnel2dimag (yp, xp, y, x, k, z):
    """
    This function has the same purpose, and uses the same variables as "Fresnel2dreal",
    however it instead calculates the imaginary component of the Electric field (E)
    in terms of y and x. We will add the j component on at the end to show that it's imaginary.
    """
    E0 = 1 #V/m
    Fresnel2dimag = ((k * E0)/(2 * np.pi * z)) * np.sin( (k)/(2 * z) * ( (x - xp)**2 + (y - yp)**2 ) ) #V/m
    
    return Fresnel2dimag

#Functions for y limits when apeture is circular
def yminfunc (xp):
    """
    This function calculates the upper limit of the integral when the aperture is circular by taking
    in the Radius of the apeture (R), as well as a general x coordinate of the apeture (xp).
    """
    return - np.sqrt( R**2 - xp**2 ) #ymin

def ymaxfunc (xp):
    """
    This function calculates the lower limit of the integral when the aperture is circular by taking
    in the Radius of the apeture (R), as well as a general x coordinate of the apeture (xp).
    """
    return np.sqrt( R**2 - xp**2 ) #ymax

#First our General Fresnel Equation
def Fresnel2d_Mc (yp, xp, y, x, k, z):
    """ 
    This function operates similarly to "Fresnel2dreal" and "Fresnel2dimag", but combines them into complex polar form
    since our Monte Carlo function can handle complex numbers.  This is then called in "MC" to get a close estimate for the 
    integral of the function.
    """
    E0 = 1 #V/m
    return ((k * E0)/(2 * np.pi * z)) * np.exp( (1j * k)/(2 * z) * ( (x - xp)**2 + (y - yp)**2 ) ) #V/m

#Now our Monte Carlo Function
def MC (Fresnel2d_Mc, xmin, xmax, ymin, ymax, xp, yp, k, z, N):
    """
    This function will obtain a value for the integral of our Fresnel equation, using f = Fresnel2d_Mc. 
    It does this by multiplying the area we are concerned with by the mean value 
    of the function, that way we get all the values under the curve. Before we
    do this however we must get random values to get an average for the function, this
    cannot be done circilarly so we need to do it with a square where length = width 
    = radius, and then remove any values that lie outside the circle. This is done
    using a by multiplying out array of random points by a Boolean array.
    """
    x_samples = np.random.uniform(low=xmin, high=xmax, size=N) #create our random x values for Monet Carlo
    y_samples = np.random.uniform(low=ymin, high=ymax, size=N) #create our random y values for Monte Carlo
    
    r2 = x_samples**2 + y_samples**2 #Pythragoras' theorem to turn our x,y values into "distance from the origin"
    boolean_array = np.less_equal(r2, R**2) #Covert our array into a bunch of "Boolean values" (True or False: True if inside circular apeture region, False if not)
   
    values = Fresnel2d_Mc(y_samples, x_samples, yp, xp, k, z) #Call our Fresnel2d_Mc to give us our function to be pluuged into Mc
    values = values * boolean_array #Since our Boolean array converts True and False into 1 and 0 respectively, multiplying our values array by this gets rid of coordinates outside region
    
    mean = values.sum() / N #Generate our mean (sum of values divided by number of values)
    meansq = (values * np.conjugate(values)).sum() / N #Generate our mean squared (square of values added then divide by number of values)
    
    results = a * mean #Formula for our integral using the Monte-Carlo method, given in the Exeriese 2 intructions
    results_error = a * np.sqrt( (meansq - mean * mean) / N ) #Formula for error in Monte Carlo given in the Exercise 2 instructions
    return (results, results_error)

"""Give the User a Choice for near field or far field""" #set parameters depending on options picked
while True: #Gives the user the option to go through the script again if needed
    MyInput = '0'
    while MyInput != "q":  #setting a while loop to ask what shape apeture the user wants, giving an option to quit with "q"
        MyInput = input('Enter a choice, "1" for near-field (Fresnel) case, or "2" for far-field (Fraunhofer) case, or "q" tp quit.')
        if MyInput == '1':
            print('You have selected "1": near-field (Fresnel) case\n') #"\n" creates a space in the console text to make it easier to read
            
            #Defining our parameters for Near field (Fresnel diffraction)
            z = 0.005 # /m
            ymin_s, ymax_s = -2.5e-4, 2.5e-4 #y Screen limits /m
            xmin_s, xmax_s = -2.5e-4, 2.5e-4 #x Screen limits /m
            xmin, xmax = -1e-4, 1e-4 #x square apeture limits /m
            ymin, ymax = -1e-4, 1e-4 #y square apeture limits /m
            R = 1e-4 #Fresnel-Limit to the radius of circular apeture
            
            break #break loop so we can continue with program
#=========================================================================================================================                                                    
        elif MyInput == '2':
            print('You have selected Part 2: far-field (Fraunhofer) case\n') #"\n" creates a space in the console text to make it easier to read
            
            #Defining our parameters for Far field (Fraunhofer diffraction)
            z = 0.05 #Screen Distance /m
            ymin_s, ymax_s = -2.5e-3, 2.5e-3 #y Screen limits /m
            xmin_s, xmax_s = -2.5e-3, 2.5e-3 #x Screen limits /m
            xmin, xmax = -1e-5, 1e-5 #x square apeture limits /m
            ymin, ymax = -1e-5, 1e-5 #y square apeture limits /m
            R = 1e-5 #Fraunhofer-Limit to the radius of circular apeture
           
            break #break loop so we can continue with program
#==========================================================================================================================        
        elif MyInput != 'q':
            print('This is not a valid choice\n') #"\n" creates a space in the console text to make it easier to read
#=============================================================================================================================
        elif MyInput == 'q':
            print('Program ending')
            sys.exit() #exit program
#======================================================================================================================
    """Create set parameters"""
    #Let's define our set parameters 
    N = 1000 #Number of random points we want to generate for Monte Carlo
    
    #Apeture area and width
    a = (xmax - xmin) * (ymax - ymin) #Area of our square /m^2
    w = (xmax - xmin) #Width of apeture, displayed on plots
    
    #Number of divisions in axes
    xlengthsteps = 50 #For x
    ylengthsteps = 50 #For y
    
    # Screen dimensions
    y, x = np.linspace(ymin_s, ymax_s, ylengthsteps), np.linspace(xmin_s, xmax_s, xlengthsteps) #Array of x, y coordinates to integrate over /m
    
    #Define k using wavelength for Fresnel2dreal
    k = ( 2 * np.pi ) / (500e-9) #2pi/wavelength /m
    
    #Create our lists
    real_results = [] #For Real Integral
    imag_results = [] #For Imaginary Integral
    realerror_results = [] #Error in the real component 
    imagerror_results = [] #Error in the Imaginary Component
    E = [] #Results list for Monte Carlo Method
    E_error = [] #Error list for Monte Carlo Method


    """Asking for method of choice and Integrating"""
    MyInput = '0'
    while MyInput != "q":  #setting a while loop to ask what shape apeture the user wants, giving an option to quit with "q"
        MyInput = input('Enter a choice, "1" for one dimensional square apeture (dblquad), "2" for two dimensional square apeture (dblquad), "3" for two dimensional circular apeture (dblquad), "4" for two dimensional circular apeture (Monte Carlo), "b" to go back, or "q" to quit\n')
        if MyInput == '1':
            print('You have selected Part 1: One dimensional square apeture (dblquad)\n') #"\n" creates a space in the console text to make it easier to read
            Part = 1 #Keepts track of what part user chose so we know which plot to use
            start_time = time.time() #Start the stopwatch to track how long integral took
                
            #Integrate our Functions using dblquad. dblquad won't take in arrays, so we need to create a for loop to repeatedly compute the integral for each x value in my array
            for xi in x:
                y = 0 #No y dimension when in one dimension (only x)
                realpart, realerror = integrate.dblquad (Fresnel2dreal, xmin ,xmax, ymin, ymax, args=(y, xi, k, z), epsabs=1.49e-08, epsrel=1.49e-08) #integral of the real component
                imagpart, imagerror = integrate.dblquad (Fresnel2dimag, xmin ,xmax, ymin, ymax, args=(y, xi, k, z), epsabs=1.49e-08, epsrel=1.49e-08) # integral of the imaginary component
                real_results.append(realpart) #Append each real integral result to our list of real results
                imag_results.append(imagpart) #Append each imaginary integral result to our list of imaginary results
                realerror_results.append(realerror) #Append each real error to our list of errors in real component integration
                imagerror_results.append(imagerror) #Append each imaginary error to our list of errors in imaginary component integration
            end_time = time.time() #End the stopwatch to track how long integral took
            break #break loop so we can continue with program
#==================================================================================================================================                                                    
        elif MyInput == '2':
            print('You have selected Part 2: Two dimensional square apeture (dblquad)\n') #"\n" creates a space in the console text to make it easier to read
            Part = 2 #Keepts track of what part user chose so we know which plot to use
            start_time = time.time() #Start the stopwatch to track how long integral loop took
            y, x = np.linspace(ymin_s, ymax_s, ylengthsteps), np.linspace(xmin_s, xmax_s, xlengthsteps) #Recreate arrays since bug arises if 1D is selected first due to "y=0" line
            
            #Now let's integrate using the Dblquad method in both dimensions
            for yi in y:
                for xi in x:
                    realpart, realerror = integrate.dblquad (Fresnel2dreal, xmin ,xmax, ymin, ymax, args=(yi, xi, k, z), epsabs=1.49e-03, epsrel=1.49e-03) #integral of the real component
                    imagpart, imagerror = integrate.dblquad (Fresnel2dimag, xmin ,xmax, ymin, ymax, args=(yi, xi, k, z), epsabs=1.49e-03, epsrel=1.49e-03) # integral of the imaginary component
                    real_results.append(realpart) #Append each real integral result to our list of real results
                    imag_results.append(imagpart) #Append each imaginary integral result to our list of imaginary results
                    realerror_results.append(realerror) #Append each real error to our list of errors in real component integration
                    imagerror_results.append(imagerror) #Append each imaginary error to our list of errors in imaginary component integration
            end_time = time.time() #End the stopwatch to track how long integral loop took
            break #break loop so we can continue with program
#===============================================================================================================================================           
        elif MyInput == '3':   
            print('You have selected Part 3: Two dimensional circular apeture (dblquad)\n') #"\n" creates a space in the console text to make it easier to read
            Part = 3 #Keepts track of what part user chose so we know which plot to use
            start_time = time.time() #Start the stopwatch to track how long integral loop took
            xmin, xmax, ymin, ymax = -R, R, -R, R #Change apeture dimensions to align with R /m
            y, x = np.linspace(ymin_s, ymax_s, ylengthsteps), np.linspace(xmin_s, xmax_s, xlengthsteps) #Recreate arrays since bug arises if 1D is selected first due to "y=0" line
            
            #Now let's integrate using the Dblquad Method in both dimensions (with our y limits as functions of x)
            for yi in y:
                for xi in x:
                    realpart, realerror = integrate.dblquad (Fresnel2dreal, xmin ,xmax, yminfunc, ymaxfunc, args=(yi, xi, k, z), epsabs=1.49e-03, epsrel=1.49e-03) #integral of the real component
                    imagpart, imagerror = integrate.dblquad (Fresnel2dimag, xmin ,xmax, yminfunc, ymaxfunc, args=(yi, xi, k, z), epsabs=1.49e-03, epsrel=1.49e-03) # integral of the imaginary component
                    real_results.append(realpart) #Append each real integral result to our list of real results
                    imag_results.append(imagpart) #Append each imaginary integral result to our list of imaginary results
                    realerror_results.append(realerror) #Append each real error to our list of errors in real component integration
                    imagerror_results.append(imagerror) #Append each imaginary error to our list of errors in imaginary component integration
            end_time = time.time() #End the stopwatch to track how long integral loop took
            break #break loop so we can continue with program
#==============================================================================================================================================        
        elif MyInput == '4':
            print('You have selected Part 4: Two dimensional circular apeture (Monte Carlo)\n') #"\n" creates a space in the console text to make it easier to read
            Part = 4 #Keepts track of what part user chose so we know which plot to use
            start_time = time.time() #Start the stopwatch to track how long MC loop took
            xmin, xmax, ymin, ymax = -R, R, -R, R #Change apeture dimensions to align with R /m
            y, x = np.linspace(ymin_s, ymax_s, ylengthsteps), np.linspace(xmin_s, xmax_s, xlengthsteps) #Recreate arrays since bug arises if 1D is selected first due to "y=0" line
            
            #Now let's integrate using the Monte Carlo Method (in both dimensions)
            for y_ in y: #To get our y value for every y_ in y (y_ is a dummy variable)
                for x_ in x:  #To get our x value for every x_ in x (x_ is a dummy variable)
                    results, results_error = MC(Fresnel2d_Mc, xmin, xmax, ymin, ymax, x_, y_, k, z, N) 
                    E.append(results) #store out results into "E" list
                    E_error.append(results_error) #store our error results int our error list
            end_time = time.time() #End the stopwatch to track how long MC loop took
            break #break loop so we can continue with program
#==========================================================================================================================    
        elif MyInput == 'b':
            break
#=========================================================================================================================================
        elif MyInput != 'q':
            print('This is not a valid choice\n') #"\n" creates a space in the console text to make it easier to read
#=============================================================================================================================
        elif MyInput == 'q':
            print('Program ending')
            sys.exit() #exit program
#================================================================================================================================================

    #Convert lists into arrays for vectorisation
    real_results = np.asarray(real_results) #Real Component 
    imag_results = np.asarray(imag_results) #Imaginary Component 
    realerror_results = np.asarray(realerror_results) #Real error 
    imagerror_results = np.asarray(imagerror_results) #Imagnary error 
    E = np.asarray(E)
    E_error = np.asarray(E_error)
    
    #Now let's compute the intensities for a given x value, being careful to do this differently if using Monte Carlo since it uses polar form
    if Part in (1, 2, 3):
        I = e0 * c * ( real_results**2 + imag_results**2 ) #value of Intensity W/m^2
        I_rel = I / I.max() #Relative value of Intensity W/m^2
        I_error = np.sqrt (e0 * c * ( (2 * real_results * realerror_results )**2 + ( 2 * imag_results * imagerror_results )**2 ) ) #Error in I using partial differentials method W/m^2
        I_errormax = np.round(np.absolute(I_error.max()), decimals = 11) #Gives us the max error to call rounded to 11 decimal places, so that in standard form python will call it as a number with only a few decimal places independant of the method (1, 2, or 3) used (so it doesnt take up entre screen)
        
    else:
        #Computing I when using polar form (exclusive for MC)
        I = e0 * c * (E * np.conjugate(E) ).real #value of Intensity, calling our values from the Monte-Carlo method. Using conjugate since E is complex so to square we must do E * conjugate of E (also result must be real)
        I_rel = I / I.max() #Relative value of Intensity W/m^2
        I_error = E_error #W/m^2
        I_errormax = np.round(np.absolute(I_error.max()), decimals = 6) #Gives us the max error to call rounded to 6 decimal places, this must be different to the other parts since error in MC is much lower than dblquad
        
    time_taken = end_time -start_time #calculate how long chosen method took
    
    """c) Let's Plot our Diffraction Patterns"""
    if Part == 1:
        #Plotting the results
        plt.plot(x, I, label='Intensity')
        plt.xlabel("x (m)") #label x axis with units of meters
        plt.ylabel("Intensity (W/m^2)") #label y axis with units of Inensity
        plt.title('1D Dblquad Square-Apeture Diffraction;\n z = {} ,\n Time Elapsed = {:.2f} seconds, Apeture Width = {}' .format(z, time_taken, w)) #Let user know parameters used
        plt.show()

        #Plotting the errors
        plt.plot(x, I_error, color='red') #reating plott with array of the error in the Intensity against array of x values
        plt.xlabel("x (m)") #label x axis with units of meters
        plt.ylabel("Error in the Intensity (W/m^2)") #label y axis with units of Inensity
        plt.title('Error in 1D Dblquad Square-Apeture Diffraction;\n z = {} , Time Elapsed = {:.2f} seconds,\n Apeture Width = {}' .format(z, time_taken, w)) #Let user know parameters used
        plt.show()
#===========================================================================================================================================
    elif Part == 2:
        #We need to do a couple things first so that out "I" data is ready to be put into "plt.imshow".
        #Firstly, our I array looks like "y1x1,y1x2...,y1x9,y1x10, y2x1, y2x2, y2,x3..., y2x9, y2x10 ... y10x10" because of the way the for loop works
        #So all I need to do is make the code into a 2D array that increases in y every 10 values (of x).
        #This is done using "numpy.reshape"
        I_2D = np.reshape(I_rel, (ylengthsteps, xlengthsteps)) #This creates a 2D array for our results described as above to be put into "plt.imshow"
        extents = (x.min(), x.max(), y.min(), y.max()) # Sets the limits for the plot, taking the smallest and largets values from our x,y arrays which correspond to the 4 corners of our square limit

        #Now for our actual plot
        plt.imshow(I_2D,vmin=0.0,vmax=1.0*I_2D.max(),extent=extents, #Plot our variables
        origin="lower",cmap="nipy_spectral_r")  #Colour map
        plt.xlabel("x (m)") #x axis label
        plt.ylabel("y (m)") #y axis label
        plt.title('2D Dblquad Square-Apeture Diffraction;\n Maximun Error = {} W/m^2, Time Elapsed = {:.2f} seconds,\n z = {}, Apeture Width = {}' .format(I_errormax, time_taken, z, w)) #Let user know parameters used
        plt.colorbar() #Include colour bar
        plt.show()
#=====================================================================================================================================  
    elif Part == 3:
        #See input == 3 plot for explanation
        I_2D = np.reshape(I_rel, (ylengthsteps, xlengthsteps)) #This creates a 2D array described as above to be put into "plt.imshow"
        extents = (x.min(), x.max(), y.min(), y.max()) # Sets the limits for the plot, taking the smallest and largets values from our x,y arrays which correspond to the 4 corners of our square limit

        #Now for our actual plot
        plt.imshow(I_2D,vmin=0.0,vmax=1.0*I_2D.max(),extent=extents,\
        origin="lower",cmap="nipy_spectral_r") #Colour map
        plt.xlabel("x (m)") #x label
        plt.ylabel("y (m)") #y label
        plt.title('2D Dblquad Circular-Apeture Diffraction;\n Maximun Error = {} W/m^2, Time Elapsed = {:.2f} seconds,\n z = {}, Apeture Width = {}' .format(I_errormax, time_taken, z, w)) #Let user know parameters used
        plt.colorbar() #Make sure to include colour bar
        plt.show()
#==========================================================================================================================
    elif Part == 4:
        
        #See input == 3 plot for explanation
        I_2D = np.reshape(I_rel, (ylengthsteps, xlengthsteps)).real #This creates a 2D array described as above to be put into "plt.imshow". It is important the the y values are the outside loop for this to work
        extents = (x.min(), x.max(), y.min(), y.max()) # Sets the limits for the plot, taking the smallest and largets values from our x,y arrays which correspond to the 4 corners of our square limit

        #Now for our actual plot
        plt.imshow(I_2D,vmin=0.0,vmax=1.0*I_2D.max(),extent=extents,\
        origin="lower",cmap="nipy_spectral_r")  #Colour map
        plt.xlabel("x (m)") #x axis label
        plt.ylabel("y (m)") #y axis label
        plt.title('2D Monte Carlo Circular-Apeture Diffraction;\n Maximun Error = {} W/m^2, Time Elapsed = {:.2f} seconds,\n z = {}, Apeture Width = {}' .format(I_errormax, time_taken, z, w)) #Let user know parameters used
        plt.colorbar()
        plt.show()
#=======================================================================================================================
    print(f'Chosen Method took {time_taken:.2f} seconds./n') #print how long chosen method took, rounding the time to two decimal places]
    
    #Give the user to break the loop to restart the program, or to exit completley
    while True: #repeats loop until user enters a valid input
        MyInput = input('Type "r" to restart, "q" to quit/n') #needs to be inside loop as to reset MyInput each time
        if MyInput == 'r':
            break #breaks while true loop to go back to the beginning (other while true loop)
        elif MyInput =='q':
            sys.exit()
        else: print('Invalid Input/n')

