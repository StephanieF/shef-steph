import os
import sys
import datetime
from math import *                    # Get all the math functions
from time import sleep
#                           imports from plot
import matplotlib
matplotlib.use('Agg') #COMMENT THIS LINE IF YOU RUN INTO DISPLAY/RENDERING ERRORS

import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.gridspec import GridSpec
#
#**********************************************************************************************
#                                   Globals
#
site_long=-72.9426                  # '-' West of Meridian, Burlington
site_lat=41.76983
#
#site_long=-73.001883               # '-' West of Meridian, New Hartford
#site_lat=41.827747
#
#DST=0
#                               Key Varibles
j2000=2451545.0
h2freq=float(1420405752)                                    # Hydrogen Spin Flip quanta
c=float(299792458)                                          # Speed of Light, m/s
wl=float(c/h2freq)                                          # the wavelength meters
#
RtoD=float(180/3.14159265)
DtoR=float(3.14159265/180)                                  # Back and forth Rad / Degrees
#
#**********************************************************************************************
#

#***********************************************************************************************
#
#                       Check to see if Current time is Start Time
#
def Check_go():
    #                                             
    #           This function will check the time NOW and see if 'we' have
    #           exceeded the requested Start time.0
    #
    FnctDT = datetime.datetime.now()   
    ctDT = FnctDT.strftime('%Y-%m-%d-%H-%M-%S')      #Get Current Time
    print('\n Current Computer time: ' + ctDT)
    
    ct=ctDT.split("-",6)
    ct_yr=int(ct[0])
    ct_mn=int(ct[1])
    ct_dy=int(ct[2])
    ct_hr=int(ct[3])
    ct_mi=int(ct[4])
    
    yr_Ok=0
    mn_Ok=0
    dy_Ok=0
    hr_Ok=0
    mi_Ok=0                                         # Define Flags               

    if ct_yr == st_yr:
        yr_Ok=1
 
    if ct_mn == st_mn:
        mn_Ok=1
   
    if ct_dy == st_dy:
        dy_Ok=1
        
    if ct_hr == st_hr:
        hr_Ok=1
       
    if ct_mi >= st_mi:
        mi_Ok=1

    go_flag=yr_Ok + mn_Ok + dy_Ok + hr_Ok + mi_Ok 

    if go_flag==5:
        go=1
    else:
        go=0                                    #Set as False for Tests

    return go
#
#    
#***********************************************************************************************
#
#                   Calculate the Local Mean Sidereal Time
# 
def LMST(dT):                                           # Calc Julian Date and LMST
    #
    #       This function will calculate the Local Mean Sidereal Time
    #       based on the Time and predefined Longitude.
    #
    FnctJD_UTC = datetime.datetime.utcnow()   
    jdUTC = FnctJD_UTC.strftime('%Y-%m-%d-%H-%M-%S')      #Get Current Time
    #print('\n UTC Start Time ' + jdUTC+' for this data block')
    
    ct=jdUTC.split("-",6)
    jd_yr=int(ct[0])                                    # Year
    jd_mn=int(ct[1])                                    # Month
    jd_dy=int(ct[2])                                    # Day
    jd_hr=int(ct[3])                                    # Hour
    jd_mi=int(ct[4])                                    # Minute
        
    dTmin=int(dT/60)                    # 1/2 the Data Acq time, minutes
                                        # Add to AJD for LMST Calc  
#    jd_yr=int(2020)
#    jd_mn=int(7)
#    jd_dy=int(1)
#    jd_hr=int(8+4)                      # Debug
#    jd_mi=int(30)
#    dTmin=0                             # Debug/Test

    if jd_mn <3:
        jd_yr=jd_yr-1
        jd_mn=jd_mn+12                                  # Year/Month adjustment
    
    A=jd_yr/100                                         # year ALWAYS > 1582
    B=2-A+int(A/4)
    
    jd1=(float(365.25*jd_yr))
    jd2=int(float(30.6001*(jd_mn+1)))
    jd3=1720994.5+B+jd_dy
    
    AJD=float(jd1+jd2+jd3)

    #print('\n The Adjusted Julian Date is :  ' +str(AJD))
    
    T=float((AJD-j2000)/36525)
    theta0=100.4606184+36000.77004*T+.000387933*T*T
    thetaCorr=theta0 % 360
    sidetime_Grnwich=thetaCorr*24/360
    
    #print(' The Sidereal time 0Hrs GMt :   ' +str(sidetime_Grnwich))
                        
    hours=float(jd_hr)+float(jd_mi+dTmin)/60    # advance time to 1/2 data acq time
    
    #print(' UTC time middle of Data Acq:   '+str(hours))
    
    theta1=(thetaCorr+360.98564724*(hours/24))
    thetag=theta1 % 360
    SRg=(thetag*24/360) % 24
    
    #print(' The Sidereal Time at input UTC:' + str(SRg))
    
    theta=thetag+site_long
    SRu=(theta*24/360) % 24
    if SRu<0:
        SRu=SRu+24                                      # Always '+'
     
    #print(' Local Mean Sidereal time, hrs: ' + str(SRu))    
    
    return SRu                                          # Return LMST
#
#***********************************************************************************************
#
#                   Calculate the Beam RightAscension and Declination
# 
def AzEl_to_RaDec(SRang):                              # Calculate the Beam  RA/Declination
    #
    #           This function will calculate the Beams pointing RA and Dec
    #           based on the passed Sidereal time, in Degrees.v The HourAngle
    #           and finally the RA is calculated. The Declination and RA are
    #           returned.
    #
    #           This subroutine issued a 'Math Domain Error' when the HourAngle Calcuation
    #           was being executed. 'Apaprently' The arcos() has issues when the argument
    #           is one, or the diviion does not have suffieient decimal places. The print
    #           debug statement left so the used can SEE whats going on. The IF statement
    #           traps the potential error and defines the minimum HourAngle based on the
    #           compare constnt.
    #
    dec=RtoD*asin((sin(DtoR*site_lat)*sin(DtoR*El)+cos(DtoR*site_lat)*cos(DtoR*El)*cos(DtoR*Az)))
    
    #print(' The Beam Declination, Degrees: '+ str(dec))
    
    num=sin(DtoR*El)-sin(DtoR*dec)*sin(DtoR*site_lat)
    dem=cos(DtoR*dec)*cos(DtoR*site_lat)
    
    #print " Num ",num
    #print " Denm",dem
    #print " El ",El
    #print " Lat ",site_lat
    
    if num/dem<.9999998:                        # Checking to preven Math Domain Error
        Ha=RtoD*acos((sin(DtoR*El)-sin(DtoR*dec)*sin(DtoR*site_lat))/(cos(DtoR*dec)*cos(DtoR*site_lat)))
        #print "Norm Ha "
    else:
        Ha=.0431                                 # Min Calcuable angle, Deg,, found by 'Math Domain Error'

    #print(' The HourAngle, Degrees, is:    '+str(Ha))
    
    if (Az-180)<0:
        Ha=360-Ha                               # HourAngle Correction
    
    auc=SRang-Ha                                # RAW RA, Degrees

    if auc<0:
        auc=360+auc
    
    ra=(auc*24/360)%24                          # The RA, Decimal Hours

    #print(' The Beam Right Ascension, Hrs  '+str(ra))
    
    return dec,ra                               # Beam Looking Ra/Dec
#
#***********************************************************************************************\
#  
#                 Calculate the Beam Looking Galactic Latitude/Longitude   
#    
def RaDec_to_GlatGlong(RA,Dec):
    #
    #           This function calculates the Beams pointing Galactic Latitude
    #           and Longitude based on the passed pointing Ra/DEC.
    #
    #RtoD=float(180/3.14159265)
    #DtoR=float(3.14159265/180)                  # Back and forth Rad / Degrees
   
    RAdeg=float(RA*360/24)                      # Ra Hours to Degrees
    decrad=float(DtoR*Dec)
    rarad=float(DtoR*RAdeg)
    k1=float(DtoR*27.4)
    #k2=float(DtoR*192.25)
    k3=float(DtoR*(192.25-RAdeg))               # defining some Constants
      
    # Meeus Eq 8.8
    lat=RtoD*asin(sin(decrad)*sin(k1)+cos(decrad)*cos(k1)*cos(k3))
    
    # Meeus Eq 8.7
    long=303-RtoD*atan2((sin(k3)),(cos(k3)*sin(k1)-tan(decrad)*cos(k1)))
    Long=long%360
    
#    print(' The Galactic Latitude:         '+str(lat))
#    print(' The Galactic Longitude:        '+str(Long))
    
    return lat,Long                             # Beam Looking Galactic Lat/Long
#
#************************************************************************************
# 
#
print('\n======================================================================================')
sleep(0.5)
print('VIRGO: An easy-to-use spectrometer & radiometer for Radio Astronomy based on GNU Radio')
print('Taken from 0xCoto, Apostolos' )
sleep(0.5)
print('  ')
print('                   Revised, Shef Robotham shefrobotham@comcast.net                   ')
print('  ')
print('             This Python Program calls GNU Radio with User parameters                ')
print('             and writes CSV files for further Processing with a Spread               ')
print('             Sheet. This version requests the Frequency, BW, FFT Bins,               ')
print('             Total Observation time and the length time of the Data Blocks           ')
print('             within the Total Time.                                                  ')
print('  ')
print('             This code calculates the LMST, the Beams Ra/Dec and Galactic            ')
print('             Lat/Long and is included in the CSV file.                               ')
print('  ')
print('             Frequency and Velocity info are written as column headers.              ')
print('             "+" Velocity is defined as the object is coming toward the              ')
print('             the observer.                                                           ')
print('  ')
print('             The data acqusition can begin immediately or at a desired               ')
print('             Start Time. If the time entered, shown for validation, is               ')
print('             incorrect, use Control-C.                                               ')
print('  ')    
print('             Version as of:   12/30/2020 11:21                                       ')
print('  ')
print('             GNU will fail depending on internal register allocations                ')
print('             Thus the integration time is adjusted f (FFT Bins/BW). See              ')
print('             the explanation during execution.                                       ')
print(' ')
print('======================================================================================\n')
sleep(3)
print('[*] Please enter your desired observation parameters...\n')
sleep(0.5)
#
#               Input USER observation parameters as STRINGS

f_center =      str(input('Center frequency [MHz]:                 '))
bandwidth =     str(input('Bandwidth [MHz]:                        '))
fft_chan =      str(input('Number of FFT Bins [FFT size]:          '))
tot_dur=        str(input('Total Observation time, [Min]:          '))
sub_durm =      str(input('Requested Data Block time [Min]:        '))
antAz=          str(input('Enter Antenna Azimuth [Deg, 0 to 360]   '))

Az=float(antAz)
#
#       Antenna Az input checking
#
if Az<0 or Az>360:
    Az_err=1
else:
    Az_err=0

while Az_err==1:
    print('\n ! Az Enter Error, Re Enter Az 0>Az<360' )
    antAz=str(input('   Enter Antenna Azimuth [Deg, 0 to 360]  '))
    print('\n')
    Az=float(antAz)
    if Az<0 or Az>360:
        Az_err=1
    else:
        Az_err=0

antEl=          str(input('Enter Antenna Elevation [Deg,-90 to 90]:'))
El=float(antEl)                                     # String to floats
#
#       Antenna El input checking
# 
if El<-90 or El>90:
    El_err=1
else:
    El_err=0 

while El_err==1:
    print('\n ! El Enter Error, Re Enter El -90>El<90' )
    antEl=str(input('   Enter Antenna Elevation [Deg,-90 to 90]:'))
    print('\n')
    El=float(antEl)
    if El<-90 or El>90:
        El_err=1
    else:
        El_err=0 
#
#       Add Session Notes
#
Session_Notes=raw_input("Enter Session Notes   ")
#print(Session_Notes)        
#
#       Time to GO?
#        
yes = {'y', 'ye', 'yes'}

strt_flag =    raw_input("Start Immediately [y or n]             ") 

if strt_flag.lower() in yes:                            # Delay if Needed                       
    strt_flag = True
else:
    strt_flag = False
#       
#
#       Adjust User Inputs accordingly and keep as String for GNU varible passing

f_center = str(float(f_center)*10**6)
bandwidth = str(float(bandwidth)*10**6)

FFTBins=int(fft_chan)                               # interger Bins                         
BW=float(bandwidth)                                 # BW, Hz
freq=float(f_center)                                # Center Freq, Hz

#                   Convert Min to seconds etc...
tot_dur_s_int=int(tot_dur)*60                       # total Obs time, seconds
sub_durs=int(sub_durm)*60                           # REQUESTED sub block time, seconds
#subdurstr=str(sub_durs)                             # convert to string for GNU
tot_dur_min_int=int(tot_dur)                        # total Obs time, minutes

halfsub_dur=int(sub_durs/2)                         # 1/2 the sub duration, seconds
                                                    # for midway Ra/Dec Calcs
#
#                   Adjust integration time and sub-block times to prevent
#                   GNU Radio 'buffer allocation fault'
# 
tint_max=int(float((133.12*10**6)/BW))              # The MAX integration time for GNURdio
#       IF computed tint_max greater than the requested sub block...
if tint_max> sub_durs:
    tint_max=sub_durs
file_Wrts=int((sub_durs/tint_max)+.5)               # the number of csv file entries
sub_durs_act=int(tint_max*file_Wrts)                # ACTUAL sub length, seconds
sub_durs_act_str=str(sub_durs_act)                  # Convert for GNU Radio,seconds
sub_durm_act_str=str(float(float(sub_durs_act)/60))        # Actual Sub Lenght, Min for Reporting

print('\n  Based on the requested BW and Sub Block Time, the integration Time is adjusted to')
print('  prevent a GNU Buffer Allocation Fault. The maximum allowable integration time is   ')
print('  130 s which ripples thru the actual sub block and total times.')
print('  The adjusted Integration Time, sec:         '+ str(tint_max))
print('  The number of rows written to the csv File: '+ str(file_Wrts))
print('  The ACTUAL Sub Block time, sec:             '+ str(sub_durs_act))
print('  The ACTUAL Sub Block Time Min:              '+ sub_durm_act_str)
 
#       Deciamte factor for GNU Radio Integrator
dec_factr = str(int(float(tint_max) * float(bandwidth)/float(fft_chan)))                    
dec_factr_int= int(dec_factr)               
#t_int_int=int(t_int)                                # make integers 

nbr_gnu_launch = int(tot_dur_s_int/int(sub_durs_act)+.5)
act_tot_time_m= float(float(sub_durs_act)*float(nbr_gnu_launch)/60)
print('  The ACTUAL total time, min:                 '+ str(act_tot_time_m))

print('  The nuber of GNU Radio launches             ' + str(nbr_gnu_launch))

sleep(3)
gnu_launch=nbr_gnu_launch                           # total for user prompting            

#                   Delete pre-existing observation.dat & plot.png files
try:
    os.remove('observation.dat')
    os.remove('plot.png')
except OSError:
    pass

#************************************************************************************
#
#                   When to start the Data Acqusition
#
if strt_flag==False:
#    
#                  Get current datetime
#
    currentDT = datetime.datetime.now()
    obsDT = currentDT.strftime('%Y-%m-%d %H:%M:%S')             # for User presentation
    print('\n Current Computer time: ' + obsDT)
    
    if strt_flag is False:
        st_yr=int(input(' Start Year        :  '))    
        st_mn=int(input(' Start Month [1-12]:  '))
        st_dy=int(input(' Start Day [1-31]:    '))
        st_hr=int(input(' Start Hour [0-23]:   '))
        st_mi=int(input(' Start Minute [0-59]: '))
        print('\n Entered Start Time:')    
        print (st_yr,st_mn,st_dy,st_hr,st_mi)
    
    go_flag=0
    go_flag=Check_go()                                          # Check time for While Loop

    while go_flag==0:
        print(' Waiting')
        sleep(10)
        go_flag=Check_go()            
        
print('\n Starting Data Collection')
#
#
#************************************************************************************
#
currentDT = datetime.datetime.now()
obsDT = currentDT.strftime('%Y-%m-%d %H:%M:%S')             # for User presentation
obsDT2 = currentDT.strftime('%Y-%m-%d-%H-%M-%S')            # Different string Format for CSV filename
print('\n Starting observation at ' + obsDT + ' local computer time.')

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#
#                   Set up for-to loop launching GNU Radio n times
#
k=0                                                         # GNU Launch counter

#//////////////////////////////////////////////////////////////////////////////////////////////////////////

for k in range(nbr_gnu_launch):


    gnu_launch_nbr=str(k)                                   # Convert to String for filename
    gnu_launch_ID=str(k+1)                                  # so Launch # and Tot # agree at end
    
#***********************************************************************************************
#
#                       Calculate Beam Ra and Dec and Galactic Lat/Long   
#    
    userSideTime=LMST(halfsub_dur)                  # Sidereal Time for input parms, hours
                                                    # 1/2 the SubFrame time passed for LMST,sec
    SRang=userSideTime*360/24 % 360                 # Sidereal tiime angel, Degrees

    Dec,RA=AzEl_to_RaDec(SRang)                     # Calc Beam Ra/Declination
    
    gaLat,gaLong=RaDec_to_GlatGlong(RA,Dec)         # Calc Galactic Lat/Long
    
    GaLat=str(gaLat)
    GaLong=str(gaLong)
    
    Strt_UTC = datetime.datetime.utcnow()   
    FrameStartUTC = Strt_UTC.strftime('%Y-%m-%d-%H-%M-%S')  #Get UTC Frame Start Time 
    
    print('\n Frame Start UTC: '+ FrameStartUTC)
    print(' GNU iteration : '+str(gnu_launch_ID))
    print('   Midway Beam Pointing Data:')
    print(' LMST, Hrs     : '+str(userSideTime))
    print(' RA, Hrs       : '+str(RA))
    print(' Dec           : '+str(Dec))
    print(' Galactic Lat  ; '+str(gaLat))
    print(' Galactic Long : '+str(gaLong))
    print('\n')
    
#
#***********************************************************************************************\#

#                   Delete pre-existing observation.dat & plot.png files
    try:
        os.remove('observation.dat')
        os.remove('plot.png')
    except OSError:
        pass

#                           Launch GNU Radio with passed string parameters
#                           Execute top_block.py with parameters
    sys.argv = ['top_block.py', '--c-freq='+f_center, '--samp-rate='+bandwidth, '--nchan='+fft_chan, '--nbin='+dec_factr, '--obs-time='+sub_durs_act_str]
    execfile('top_block.py')

    print('  ')
    print(gnu_launch_ID + ' GNU Observation finished! Saving data...')

#                   Segment of plot.py with CSV file writing

    def decibel(x):
#       return 10.0*np.log10(x)             # Uncomment for dB-scaled Power axes
        return x

    fname = "observation.dat"               # File name from GNU Radio
#

#
#                           Load data from GNU Radio file
    z = np.fromfile(fname, dtype="float32").reshape(-1, FFTBins)/dec_factr_int
    z = z*10000                             #Rescale Power values

#                       Number of sub-integrations in observation data
    nsub = z.shape[0]

#   Now Write a CSV Data file 

    print "\n             Frame DATA"
    print "Observation Start Time: ",obsDT
    print "Session Notes           ", Session_Notes
    print "Center Freq:            ", f_center
    print "BW:                     ", BW
    print "FFT Bins:               ", FFTBins
    print "Total Obs time min:     ", tot_dur 
    print "Actual Total Obs time:  ", act_tot_time_m
    print "Req Sub Block time min  ", sub_durm
    print "Actual Sub Block, min:  ", sub_durm_act_str
    print "Int time sec:           ", tint_max
    print "GNU iteration           ", gnu_launch_ID
    print "Tot GNU Iterations      ", gnu_launch
    print "Antenna Azimuth:        ", antAz
    print "Antenna Elevation:      ", antEl
    print "LMST Midway Acq         ", userSideTime
    print "RA Midway Acq           ", RA
    print "Dec Midway Acq          ", Dec
    print "Galactic Lat Midway     ", GaLat
    print "Galactic Long Midway    ", GaLong

#       Break up date-time to make data file name:
    obs=obsDT2.split("-",5)                     # retrieve time at start
#   print(obs)
    obsYr=obs[0]
    obsMn=obs[1]
    obsDy=obs[2]
    obsHr=obs[3]
    obsMin=obs[4]

    filename=obsMn + obsDy +" " + obsHr + obsMin + " "+ obsYr + " "+ gnu_launch_ID +'.txt'
    
    print('  ')
    print('  ')
    print "The data will be saved to :  ",filename
    print('  ')
    print('  ')

    f = open(filename,'w+')
    
    f.write("Observation Start time:    " + obsDT)
    f.write('\n')
    f.write("Session Notes: {}\n".format(Session_Notes))
    f.write("Center Freq               {}\n".format(f_center))
    f.write("Band Width                {}\n".format(BW))
    f.write("Nbr FFT Bins              {}\n".format(FFTBins))
    f.write("Tot Obs time  Min         {}\n".format(tot_dur))
    f.write("Act Total Obs Min         {}\n".format(act_tot_time_m))
    f.write("Req Block Time Min        {}\n".format(sub_durm))
    f.write("Actual Block Time Min     {}\n".format(sub_durm_act_str))  
    f.write("Integ Time Sec            {}\n".format(tint_max))
    f.write("GNU Iteration #           {}\n".format(gnu_launch_ID))
    f.write("Tot GNU Iterations        {}\n".format(gnu_launch))
    f.write("Site Latitude             {}\n".format(site_lat))
    f.write("Site Longitude            {}\n".format(site_long))
    f.write('\n')
    f.write("This Frame Start Time, UTC "+ FrameStartUTC)
    f.write('\n')
    f.write("Antenna Az, Deg           {}\n".format(antAz))
    f.write("Antenna El, Deg           {}\n".format(antEl))
    f.write("LMST during Acq, Hrs      {}\n".format(userSideTime))
    f.write("RA during Acq, Hrs        {}\n".format(RA))
    f.write("Dec during Acq, Deg       {}\n".format(Dec))
    f.write("Galactic Latitude, Deg    {}\n".format(GaLat))
    f.write("Galactic Longitude, Deg   {}\n".format(GaLong))
    f.write('\n')
    f.write('\n')

#       Create Frequency Header for Excel 
    FStep=BW/FFTBins                                # the dF each FFT Sample, Hz
    Fmin=freq-(FFTBins*FStep)/2                     # the minimum Frequency
    vmin=(h2freq-Fmin)*wl                           # the minimum velocty

    i=0                                             # define and reset counter
    ifreq=Fmin/1000000                              # set var for incrementing, MHz

    for i in range(FFTBins):
        f.write('{},'.format(ifreq))
        ifreq=ifreq+(FStep/1000000)
    f.write('\n')
    
    f.write('\n')                                   # empty row for Exccel Formatting
    
    i=0                                             # Reset
    ifreq=Fmin                                      # set var for incrementing, Hz
    for i in range(FFTBins):
        ivel=((ifreq-h2freq)*wl)/1000              # calc i velocity km/s 
        f.write('{},'.format(ivel))                # '+' Vel, coming toward Obsvr
        ifreq=ifreq+FStep
    f.write('\n')                                   # finish w CR/LF
    
    f.write('\n')                                   # empty row for Exccel Formatting 
    
#       Now Step thru z generating a csv file for Excel
#       the FFT sample will reside in a Column
#       Integrations will reside in Rows

    i=0                                             # Define Counter and reset
    j=0                                             # Define Counter and reset

#   J is the Row counter, Outside, i the FFT sample Counter, Inside

    while j<nsub:
        while i<FFTBins:
            jidata=z[j,i]                           # the j,i data value
        #   print jidata, j, i                      # Debug
            if i==(FFTBins-1):
                f.write('{}\n'.format(jidata))      # end of samples, write CR/LF
            else:
               f.write('{},'.format(jidata))       # write the j,i data to the file
            i=i+1                                   # inc i
        i=0                                         # Reset i Counter
        j=j+1                                       # inc j

    f.close()


#//////////////////////////////////////////////////////////////////////////////////////////////////////////

#
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#
#                   Set up for-to loop launching GNU Radio n times
#
#                   Delete pre-existing observation.dat & plot.png files
try:
    os.remove('observation.dat')
    os.remove('plot.png')
except OSError:
    pass


print('\n======================================================================================')
print('[+] The observation data complete and saved as CSV for Excel Import.')
print('======================================================================================')
