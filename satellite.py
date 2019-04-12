from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import J2_Propagator as j2
import numpy as np
import julian 
import testcase_satellite as tc
import matplotlib.pyplot as plt
import read_ins_2018_tle as rins


#-----Defining a class Satellite-----
class satellite:
    def __init__(self, line1, line2):
        self.line1 = line1 
        self.line2 = line2
        self.jddate = twoline2rv(self.line1, self.line2, wgs72).jdsatepoch #Initializing the Julian Date for given TLE.
        self.gregdate, self.gregtime = self.read_date(self.jddate) #Initializing the Gregorian Date, and Time (24 Hrs. format) for given TLE.
        self.r_init, self.v_init = self.sgpPropagate(self.jddate) #Initializing Inital State Vectors at t = 0 from TLE.
        self.b_star = self.read_b_star(self.line1) #Initializing b_star term for given TLE at epoch.
        
    def read_b_star(self, line):
        """Returns the BSTAR term as a float.
    
           Parameters:
           -----------
           line: string
               Takes in the first line of the TLE, as a string.  
        """ 
        temp = line[54:61].split('-')
        if temp[0] == '00000+0':
            drag = 0
        else:
            drag = (float(temp[0])*10**-5) * (10**(-(float(temp[1]))))
        return drag    
    
    def read_date(self, jd_temp):
        """Returns the Gregorian Date and Time in a string format.
    
           Parameters:
           -----------
           jd_temp: floating-point number
               Takes in the Julian Date as a floating point number.  
        """ 
        a = julian.from_jd(jd_temp)       
        y = a.year
        mont = a.month
        d = a.day
        h = a.hour
        m = a.minute
        s = a.second + a.microsecond*10**-6
        temp_time = str(h) + ':' + str(m) + ':' + str(s) 
        temp_date = str(d) + '/' + str(mont) + '/' + str(y)
        return temp_date, temp_time
        
    def sgpPropagate(self,jd=0):
        """Propogate the State Variables through SGP4 model, and returns the
           propagated R and V state vectors.
    
           Parameters:
           -----------
           jd: floating-point number, optional
               Takes in the time period over which the state vectors have 
               to be propagated, in terms of delta - time, in julian format.
               Note: 1 Julian Time = 24 Hrs.
               
               Calling the default value, gives the initial state vectors of
               the satelite, i.e, at t = 0.
               Default = 0     
        """ 
        if jd == 0:
            jddate = twoline2rv(self.line1, self.line2, wgs72).jdsatepoch
        else:
            jddate = jd
        a = julian.from_jd(jddate)       
        y = a.year
        mont = a.month
        d = a.day
        h = a.hour
        m = a.minute
        s = a.second + a.microsecond*10**-6  
        position, velocity = twoline2rv(self.line1, self.line2,wgs72).propagate(y,mont,d,h,m,s)
        return np.array(position, np.float64), np.array(velocity, np.float64)
        
    def comparej2Vsgp(self,deltatime):
        """Propogate the State Variables through both SGP4 and J2 propagators,
           and returns the R,V - values in SGP4, and J2 in that order, i.e,
           R - SGP4, V - SGP4, R - J2, V - J2.
               
           Parameters:
           -----------
           deltatime: floating-point number, optional
               Takes in the time period over which the state vectors have 
               to be propagated, in terms of delta - time, in hours.  
        """ 
        rsgp, vsgp = self.sgpPropagate(self.jddate + (deltatime/24))
        rj2, vj2 = j2.propagate(self.r_init, self.v_init, deltatime*3600, h_step_size = 1) 
        return rsgp, vsgp, rj2, vj2
    
    
    

def mock_runs(ch):
    """Runs differnet variants of Propagator Analysis, based
        on user-defined input.
    
           Parameters:
           -----------
           ch: integer
               Takes the input of the choice. 

           Notes:
           ------
           1: Block of Code that plots the error in norm of delta-r, delta-v
              w.r.t the propagated state values (The initial states for 
              this obtained from TLE - 1 (sat1) and the predicted state 
              values from a later TLE - 2 (sat2).

           2: Block of code that prints the Propagated R, V state vectors
              from SGP4 and J2 from TLE - 1, along with the Initial state 
              vectors generated from TLE2. It also calculates the error in 
              norm of delta R and delta V for both SGP4 and J2 against the 
              initial state vectors.

           3: Block of code that prints the Pe3rcentage Error in the norm of
              delta R, and V divided by the norm of corresponding R or V 
              of SGP4, over a given time interval.

           4: Block of Code that plots the error in norm of delta-r, delta-v
              w.r.t the propagated state values (The initial states for 
              this obtained from TLE - 1 (sat1) and the predicted state 
              values from a later TLE - 2 (sat2).
              NOTE: This works on the complete data set of TLE's available
    """ 
    
    #----- i takes int-values from 1 to 7, and corresonds to-----
    #            TLEs of older times as i increases .
    i = 1


    if ch == 1:

        erj2, evj2 = [], []
        ers, evs = [], []
        erj2p, evj2p = [], []
        ersp, evsp = [], []
        j = range(1,3)
        for i in j:
            l21, l22 = tc.tleINS[0:2]
            l11, l12 = tc.tleINS[2*i:2*i+2]
            sat1 = satellite(l11,l12)
            sat2 = satellite(l21,l22)
            
            deltaT = (sat2.jddate - sat1.jddate)*24
            rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
            print('Time(Hrs):', np.round(deltaT,2))
            
            R,V = sat2.r_init, sat2.v_init
            erj2.append( (j2.norm(rj - R)) )
            evj2.append( (j2.norm(vj - V )) )
            ers.append( j2.norm(rs - R))
            evs.append( j2.norm(vs - V))
            
            erj2p.append( 100*(j2.norm(rj - R))/j2.norm(R) )
            evj2p.append( 100*(j2.norm(vj - V ))/j2.norm(V) )
            ersp.append( 100*j2.norm(rs - R)/j2.norm(R))
            evsp.append( 100*j2.norm(vs - V)/j2.norm(V))        
        
        plt.plot(j, erj2,'bo')
        plt.plot(j, ers, 'ro')
        plt.legend(["Error Position- J2", "Error Position - SGP4"], fontsize = 'large')
        plt.plot(j, erj2, 'b--')
        plt.plot(j, ers, 'r--')
        plt.title(r'Norm of $\delta r$', fontdict = {'fontsize':20})
        plt.xlabel(r'Iterations $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta r (kms) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, evj2,'bo')
        plt.plot(j, evs,'ro')
        plt.title(r'Norm of $\delta v$', fontdict = {'fontsize':20})
        plt.legend(["Error Velocity - J2", "Error Velocity - SGP4"], fontsize='large')
        plt.plot(j, evj2,'b--')
        plt.plot(j, evs,'r--')
        plt.xlabel(r'Iterations $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta v (km/s) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, erj2p, 'bo')
        plt.plot(j, ersp, 'ro')
        plt.legend(["Error Position- J2", "Error Position - SGP4"], fontsize = 'large')
        plt.title(r'Norm of $\delta r (in  \%)$', fontdict = {'fontsize':24})
        plt.xlabel(r'Iterations  $\rightarrow$', fontdict = {'fontsize':15})
        plt.plot(j, erj2p, 'b--')
        plt.plot(j, ersp, 'r--')
        plt.ylabel(r'$\frac{|\delta r| \times 100}{|R|}$', fontdict = {'fontsize':22} )
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, evj2p,'bo')
        plt.plot(j, evsp,'ro')
        plt.title(r'Norm of $\delta v (in  \%)$', fontdict = {'fontsize':24})
        plt.legend(["Error Velocity - J2", "Error Velocity - SGP4"], fontsize = 'large')
        plt.plot(j, evj2p,'b--')
        plt.plot(j, evsp,'r--')
        plt.xlabel(r'Iterations $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\frac{|\delta v| \times 100} {|V|}$', fontdict = {'fontsize':22})
        plt.grid()
        plt.show()

    elif ch == 2:

        l21, l22 = tc.tleINS[0:2]
        l11, l12 = tc.tleINS[2*i:2*i+2]
        sat1 = satellite(l11,l12)
        sat2 = satellite(l21,l22)
        deltaT = (sat2.jddate - sat1.jddate)*24
        rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
        R,V = sat2.r_init, sat2.v_init
        print('SGP4 ---> ', rs, vs, ' ---- ', sat1.gregdate)
        print('J2 -----> ',rj,vj)
        print('Obs ----> ', R, V, ' ---- ', sat2.gregdate)
        erj2, evj2 = j2.norm(rj - R), j2.norm(vj - V)
        ers, evs = j2.norm(rs - R), j2.norm(vs - V)
        print('Error - R (J2) = ', erj2, '       Error - V(J2) = ',evj2)
        print()
        print('Error - R (SGP4) = ', ers, '       Error - V(SGP4) = ',evs)

    elif ch == 3:

        ERJ, EVJ = [], []
        ERS, EVS = [], []
        l21, l22 = tc.tleINS[0:2]
        l11, l12 = tc.tleINS[2*i:2*i+2]
        sat1 = satellite(l11,l12)
        sat2 = satellite(l21,l22)
        hours = 4 #Set time duration
        scale = np.arange(0,3600*hours,220)
        for deltaT in scale:
            deltaT=deltaT/3600
            rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
            erj2, evj2 = j2.norm(rj - rs)/j2.norm(rs), j2.norm(vj - vs)/j2.norm(vs)
            ERJ.append(erj2*100)
            EVJ.append(evj2*100)
        plt.plot(scale, ERJ,'o')
        plt.plot(scale,EVJ,'o')
        plt.legend(['Error - r', 'Error - v'])
        plt.title('Error - Position, Velocity')
        plt.xlabel('Time (seconds) ---->')
        plt.ylabel('Error (%): (norm(delta vector)/norm(final vector)) * 100:')
        plt.show()

    elif ch == 4:

        erj2, evj2 = [], []
        ers, evs = [], []
        erj2p, evj2p = [], []
        ersp, evsp = [], []
        j = range(0,101,5)
        t_axis = []
        for i in j:
            t = 1
            l11, l12 = tc.tle_2018_complete[(t-1)*2:t*2]
            l21, l22 = tc.tle_2018_complete[2*i:2*i+2]
            sat1 = satellite(l11,l12)
            sat2 = satellite(l21,l22)
            
            deltaT = (sat2.jddate - sat1.jddate)*24
            rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
            print('Time(Hrs):', np.round(deltaT,2))
            t_axis.append(np.round(deltaT,2))
            
            R,V = sat2.r_init, sat2.v_init
            erj2.append( (j2.norm(rj - R)) )
            evj2.append( (j2.norm(vj - V )) )
            ers.append( j2.norm(rs - R))
            evs.append( j2.norm(vs - V))
            
            erj2p.append( 100*(j2.norm(rj - R))/j2.norm(R) )
            evj2p.append( 100*(j2.norm(vj - V ))/j2.norm(V) )
            ersp.append( 100*j2.norm(rs - R)/j2.norm(R))
            evsp.append( 100*j2.norm(vs - V)/j2.norm(V))  
            
            
        
        j = t_axis
        plt.plot(j, erj2,'bo')
        plt.plot(j, ers, 'ro')
        plt.legend(["Error Position- J2", "Error Position - SGP4"], fontsize = 'large')
        plt.plot(j, erj2, 'b--')
        plt.plot(j, ers, 'r--')
        plt.title(r'Norm of $\delta r$', fontdict = {'fontsize':20})
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta r (kms) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, evj2,'bo')
        plt.plot(j, evs,'ro')
        plt.title(r'Norm of $\delta v$', fontdict = {'fontsize':20})
        plt.legend(["Error Velocity - J2", "Error Velocity - SGP4"], fontsize='large')
        plt.plot(j, evj2,'b--')
        plt.plot(j, evs,'r--')
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta v (km/s) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, erj2p, 'bo')
        plt.plot(j, ersp, 'ro')
        plt.legend(["Error Position- J2", "Error Position - SGP4"], fontsize = 'large')
        plt.title(r'Norm of $\delta r (in  \%)$', fontdict = {'fontsize':24})
        plt.xlabel(r'Hours  $\rightarrow$', fontdict = {'fontsize':15})
        plt.plot(j, erj2p, 'b--')
        plt.plot(j, ersp, 'r--')
        plt.ylabel(r'$\frac{|\delta r| \times 100}{|R|}$', fontdict = {'fontsize':22} )
        plt.grid()
        plt.show()
        plt.figure()
        plt.plot(j, evj2p,'bo')
        plt.plot(j, evsp,'ro')
        plt.title(r'Norm of $\delta v (in  \%)$', fontdict = {'fontsize':24})
        plt.legend(["Error Velocity - J2", "Error Velocity - SGP4"], fontsize = 'large')
        plt.plot(j, evj2p,'b--')
        plt.plot(j, evsp,'r--')
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\frac{|\delta v| \times 100} {|V|}$', fontdict = {'fontsize':22})
        plt.grid()
        plt.show()
