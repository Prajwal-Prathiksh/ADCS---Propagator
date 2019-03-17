from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import J2_Propagator as j2
import numpy as np
import julian 
import testcase_satellite as tc
import matplotlib.pyplot as plt


#-----Defining a class Satellite-----
class satellite:
    def __init__(self, line1, line2):
        self.line1 = line1 
        self.line2 = line2
        self.jddate = twoline2rv(self.line1, self.line2, wgs72).jdsatepoch #Initializing the Julian Date for given TLE.
        self.gregdate = julian.from_jd(self.jddate) #Initializing the Gregorian Date for given TLE.
        self.r_init, self.v_init = self.sgpPropagate(self.jddate) #Initializing Inital State Vectors at t = 0 from TLE.
    
    
    def sgpPropagate(self,jd=0):
        """Propogate the State Variables through SGP4 model, and returns the
           propagated R and V state vectors.
    
           Parameters:
           -----------
           jd: floating-point number, optional
               Takes in the time period over which the state vectors have 
               to be propagated, in terms of delta - time, in julian format.
               
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
        
   
#----- i takes int-values from 1 to 7, and corresonds to-----
#            TLEs of older times as i increases .
i = 1


#-----Block of Code that plots the error in norm of delta-r, delta-v-------
#       w.r.t the propagated state values (The initial states for 
#       this obtained from TLE - 1 (sat1) and the predicted state 
#                values from a later TLE - 2 (sat2).
'''
erj2, evj2 = [], []
ers, evs = [], []
j = range(1,8)
for i in j:
    l21, l22 = tc.tleINS[0:2]
    l11, l12 = tc.tleINS[2*i:2*i+2]

    sat1 = satellite(l11,l12)
    sat2 = satellite(l21,l22)
    
    deltaT = (sat2.jddate - sat1.jddate)*24
    rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
    print(np.round(deltaT,2))
    
    R,V = sat2.r_init, sat2.v_init
    erj2.append( (j2.norm(rj - R)) )
    evj2.append( (j2.norm(vj - V ))/j2.norm(V) )
    ers.append( j2.norm(rs - R)/j2.norm(R) )
    evs.append( j2.norm(vs - V)/j2.norm(V) ) 
    
    
plt.plot(j, erj2, 'o')
plt.plot(j, ers, 'o')
plt.legend(["Error Position- J2", "Error Position - SGP4"])
plt.title('Norm of delta-r')
plt.xlabel('Iterations')
plt.ylabel('delta-r (kms)')
plt.show()

plt.figure()
plt.plot(j, evj2,'o')
plt.plot(j, evs,'o')
plt.title('Norm of delta-v')
plt.legend(["Error Velocity - J2", "Error Velocity - SGP4"])
plt.xlabel('Iterations')
plt.ylabel('delta-v (km/s)')
plt.show()
'''


#-----Block of code that prints the Propagated R, V state vectors-----
#     from SGP4 and J2 from TLE - 1, along with the Initial state 
#     vectors generated from TLE2. It also calculates the error in 
#     norm of delta R and delta V for both SGP4 and J2 against the 
#                      initial state vectors.
'''
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
'''


#-----Block of code that prints the Percentage Error in the norm of------
#      delta R, and V divided by the norm of corresponding R or V 
#               of SGP4, over a given time interval.
'''
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
'''
