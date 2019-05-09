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
        
    def comparej2Vsgp(self,deltatime, drag_boolean = False):
        """Propogate the State Variables through both SGP4 and J2 propagators,
           and returns the R,V - values in SGP4, and J2 in that order, i.e,
           R - SGP4, V - SGP4, R - J2, V - J2.
           
           NOTE: BY DEFAULT DRAG == FALSE!!           
               
           Parameters:
           -----------
           deltatime: floating-point number,
               Takes in the time period over which the state vectors have 
               to be propagated, in terms of delta - time, in hours. 
           
           drag_boolean : boolean, optional,
                Sets the drag perturbating factor to either True/False in J2.
                Default = False
        """ 
        rsgp, vsgp = self.sgpPropagate(self.jddate + (deltatime/24))
        rj2, vj2 = j2.propagate(self.r_init, self.v_init, deltatime*3600, h_step_size = 1, drag = drag_boolean) 
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
           Purpose of each integer input :-

           1: Compares propagated J2 & SGP4 state vectors from TLE - 1 (sat1), 
              to the expected state vectors of (t = 0) from TLE -2 (sat2).
              Input:
                - TLE (sat2)    : Can be changed FROM WITHIN the code only, by changing (i).
              Output: 
                - SGP4          : Propagated state vectors
                - J2            : Propagated state vectors
                - TLE - 2       : Predicted state vectors
                - Norm Errors   : B/w J2 - TLE2 and SGP4 - TLE2

           2: Plots the difference between the propagted vectors from SGP4 
              and J2 (sat1), over a time duration.
              Input:
                - Time Duration : In Hrs.
                - TLE (sat1)    : Can be changed FROM WITHIN the code only, by changing (i).
              Output:
                - Plot - 1      : Norm of difference between the position vectors.
                - Plot - 2      : Norm of difference between the velocity vectors.

           3: Plots the norm error in propagated state vectors of SGP4 and 
              J2 (including and exclding drag) from a TLE - 1 (sat1), w.r.t the 
              expected state vectors of (t = 0) from a later TLE - 2 (sat2),
              by keeping a FIXED TLE 1 (sat1), and comparing with successive 
              TLEs (sat2) for a range of TLEs.
              Input:
                - TLE (sat1 & sat2) : Can only be changed WITHIN the code, by changing range of (j).
            
              Output:
                - Plot - 1      : Norm of difference between the position vectors.
                - Plot - 2      : Norm of difference between the velocity vectors.
                - Plot - 3      : Norm difference in the postion vectors, in (%).
                - Plot - 4      : Norm difference in the velocity vectors, in (%).

           4: Plots the norm error in propagated state vectors of SGP4 ONLY, 
              between that of a FIXED TLE - 1 to the expected state vectors 
              of (t = 0) from TLE -2 (satf), by creating multiple instances
              of the same TLE - 1 (sat2, sat3, ... , sat5), but with different
              BSTAR values, to check for the variation in SGP4 propagation 
              w.r.t the BSTAR term.
              Input:
                - Iterations - Range    : Can be changed from FROM WITHIN the code only,
                                          by changing (iterations).
                - Initial TLE (sat1, ..., sat5)   : Can be changed FROM WITHIN the code only,
                                                    by changing (tle_num).
              Output:
                - Plot - 1      : Norm difference in the postion vectors, in (%), of sat1, ... , sat5.
                - Plot - 2      : Norm difference in the velocity vectors, in (%), of sat1, ... , sat5.
    """ 
    

    if ch == 1:
        #----- i takes int-values from 1 to 1374, and -----
        #     corresonds to newer TLEs as i increases.
        i = 2

        l21, l22 = tc.tle_2018_complete[0:2]
        l11, l12 = tc.tle_2018_complete[2*i:2*i+2]
        sat1 = satellite(l11,l12)
        sat2 = satellite(l21,l22)

        deltaT = (sat2.jddate - sat1.jddate)*24
        rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
        R,V = sat2.r_init, sat2.v_init

        print('SGP4 (Propagated R & V) ---> ', rs, vs, ' ---- ', sat1.gregdate)
        print('J2 (Propagated R & V) -----> ',rj,vj)
        print('Obs (TLE - 2 R & V) ----> ', R, V, ' ---- ', sat2.gregdate)

        erj2, evj2 = j2.norm(rj - R), j2.norm(vj - V)
        ers, evs = j2.norm(rs - R), j2.norm(vs - V)

        print('Error - R (J2) = ', erj2, '       Error - V(J2) = ',evj2)
        print()
        print('Error - R (SGP4) = ', ers, '       Error - V(SGP4) = ',evs)

    elif ch == 2:
        #----- i takes int-values from 0 to 1374, and -----
        #     corresonds to newer TLEs as i increases.
        i = 0

        hours = float(input('Enter Time Duration (in Hrs.): '))
        #Accepting Time Duration Input.

        ERJ, EVJ = [], []
        ERS, EVS = [], []
        l11, l12 = tc.tle_2018_complete[2*i:2*i+2]
        sat1 = satellite(l11,l12)
        
        scale = np.arange(0,3600*hours,220)
        for deltaT in scale:
            deltaT=deltaT/3600
            rs,vs,rj,vj = sat1.comparej2Vsgp(deltaT)
            erj2, evj2 = j2.norm(rj - rs), j2.norm(vj - vs)
            ERJ.append(erj2)
            EVJ.append(evj2)
        plt.plot(scale, ERJ,'o')
        plt.legend(['Difference - r'])
        plt.title('Difference Norm - Position Vector', fontsize = 30)
        plt.xlabel(r'Time (seconds) $\rightarrow$', fontsize = 25)
        plt.ylabel(r'$ | \vec{r_{SGP4}} - \vec{r_{J2}}|$', fontsize = 25)
        plt.show()
        plt.grid()
        
        plt.figure()
        plt.plot(scale,EVJ,'o')
        plt.title('Difference Norm - Velocity Vector', fontsize = 30)
        plt.legend(['Difference - v'])
        plt.xlabel(r'Time (seconds) $\rightarrow$', fontsize = 25)
        plt.ylabel(r'$ | \vec{v_{SGP4}} - \vec{v_{J2}}|$', fontsize = 25)
        plt.show()
        plt.grid()

    elif ch == 3:
        #----- (j) represents the range of TLEs over which the  -----
        #    program will compare the propagated state vectors from
        #     TLE - j[0]  with the TLEs represented by succesive 
        #                      elements of (j).
        j = range(0,8)

        # Declaring lists that will store the error values.
        erj2_nd, evj2_nd = [], []
        erj2_d, evj2_d = [], []
        ers, evs = [], []
        erj2p_nd, evj2p_nd = [], []
        erj2p_d,evj2p_d = [], [] 
        ersp, evsp = [], []
        
        t_axis = []

        t = j[0]
        l11, l12 = tc.pratham_tle_complete[t*2 : t*2+2]
        sat1 = satellite(l11,l12)

        for i in j:            

            l21, l22 = tc.pratham_tle_complete[2*i : 2*i+2]            
            sat2 = satellite(l21,l22)
            
            deltaT = (sat2.jddate - sat1.jddate)*24
            rs,vs,rjd,vjd = sat1.comparej2Vsgp(deltaT, drag_boolean=True)
            rs,vs,rjnd,vjnd = sat1.comparej2Vsgp(deltaT, drag_boolean=False)
            print('Time(Hrs):', np.round(deltaT,2))
            t_axis.append(np.round(deltaT,2))
            
            R,V = sat2.r_init, sat2.v_init
            erj2_nd.append( (j2.norm(rjnd - R)) )
            evj2_nd.append( (j2.norm(vjnd - V )) )
            evj2_d.append( (j2.norm(vjd - V )) )
            erj2_d.append( (j2.norm(rjd - R)) )
            ers.append( j2.norm(rs - R))
            evs.append( j2.norm(vs - V))
            
            erj2p_nd.append( 100*(j2.norm(rjnd - R))/j2.norm(R) )
            evj2p_nd.append( 100*(j2.norm(vjnd - V ))/j2.norm(V) )
            erj2p_d.append( 100*(j2.norm(rjd - R))/j2.norm(R) )
            evj2p_d.append( 100*(j2.norm(vjd - V ))/j2.norm(V) )
            ersp.append( 100*j2.norm(rs - R)/j2.norm(R))
            evsp.append( 100*j2.norm(vs - V)/j2.norm(V))  
            
            
        
        j = t_axis
        
        plt.plot(j, erj2_nd,'bo', label = 'Error Position - J2 (No Drag)')
        plt.plot(j, erj2_d, 'go', label = 'Error Position - J2 (Drag)')
        plt.plot(j, ers, 'ro', label = 'Error Position - SGP4') 
        plt.legend(fontsize = 'large')
        plt.plot(j, erj2_nd, 'b--')
        plt.plot(j, erj2_d, 'g--')
        plt.plot(j, ers, 'r--')
        plt.title(r'Norm of $\delta r$', fontdict = {'fontsize':20})
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta r (kms) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        
        plt.figure()
        plt.plot(j, evj2_nd,'bo', label = "Error Velocity - J2 (No Drag)")
        plt.plot(j, evj2_d,'go', label = "Error Velocity - J2 (Drag)")
        plt.plot(j, evs,'ro', label = "Error Velocity - SGP4")
        plt.title(r'Norm of $\delta v$', fontdict = {'fontsize':20})
        plt.legend(fontsize='large')
        plt.plot(j, evj2_nd,'b--')
        plt.plot(j, evj2_d,'g--')
        plt.plot(j, evs,'r--')
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\delta v (km/s) \rightarrow$', fontdict = {'fontsize':15})
        plt.grid()
        plt.show()
        
        plt.figure()
        plt.plot(j, erj2p_nd, 'bo', label = 'Error Position - J2 (No Drag)')
        plt.plot(j, erj2p_d, 'go', label = 'Error Position - J2 (Drag)')
        plt.plot(j, ersp, 'ro', label = 'Error Position - SGP4')
        plt.legend(fontsize = 'large')
        plt.title(r'Norm of $\delta r (in  \%)$', fontdict = {'fontsize':24})
        plt.xlabel(r'Hours  $\rightarrow$', fontdict = {'fontsize':15})
        plt.plot(j, erj2p_nd, 'b--')
        plt.plot(j, erj2p_d, 'g--')
        plt.plot(j, ersp, 'r--')
        plt.ylabel(r'$\frac{|\delta r| \times 100}{|R|} (in \%) \rightarrow$', fontdict = {'fontsize':22} )
        plt.grid()
        plt.show()
        
        plt.figure()
        plt.plot(j, evj2p_nd,'bo', label = "Error Velocity - J2 (No Drag)")
        plt.plot(j, evj2p_d,'go', label = "Error Velocity - J2 (Drag)")
        plt.plot(j, evsp,'ro', label = "Error Velocity - SGP4")
        plt.title(r'Norm of $\delta v (in  \%)$', fontdict = {'fontsize':24})
        plt.legend(fontsize = 'large')
        plt.plot(j, evj2p_nd,'b--')
        plt.plot(j, evj2p_d,'g--')
        plt.plot(j, evsp,'r--')
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\frac{|\delta v| \times 100} {|V|} (in \%) \rightarrow$', fontdict = {'fontsize':22})
        plt.grid()
        plt.show()

    elif ch == 4:
        #----- (iterations) represents the range of TLEs over which the  -----
        #         program will compare the propagated state vectors.
        iterations = range(0,30,2)

        # Declaring lists that will store the error values.
        er1, ev1, erp1, evp1 = [], [], [], []
        er2, ev2, erp2, evp2 = [], [], [], []
        er3, ev3, erp3, evp3 = [], [], [], []
        er4, ev4, erp4, evp4 = [], [], [], []
        er5, ev5, erp5, evp5 = [], [], [], []

        tle_num= 460
        #tle_num = 260       

        t_axis = []
        for i in iterations:

            temp = tc.pisat_tle_complete[tle_num:tle_num+2]
            l11_0, l12 = temp
            l21, l22 = tc.pisat_tle_complete[(tle_num + 2*i): (tle_num + 2*i+2)]

            sat1 = satellite(l11_0,l12)
            satf = satellite(l21,l22)

            init_drag = l11_0[54:61]
            new_drag1 = '25100-4'
            new_drag2 = '27000-4' 
            new_drag3 = '46000-4'
            new_drag4 = '78000-3'

            l11_1 = l11_0.replace(init_drag, new_drag1)
            l11_2 = l11_0.replace(init_drag, new_drag2)
            l11_3 = l11_0.replace(init_drag, new_drag3)
            l11_4 = l11_0.replace(init_drag, new_drag4)


            sat2 = satellite(l11_1, l12)
            sat3 = satellite(l11_2, l12)
            sat4 = satellite(l11_3, l12)
            sat5 = satellite(l11_4, l12)
            
            deltaT = (satf.jddate - sat1.jddate)*24
            t_axis.append(deltaT)

            r1, v1 = sat1.sgpPropagate(satf.jddate)
            r2, v2 = sat2.sgpPropagate(satf.jddate)
            r3, v3 = sat3.sgpPropagate(satf.jddate)
            r4, v4 = sat4.sgpPropagate(satf.jddate)
            r5, v5 = sat5.sgpPropagate(satf.jddate)

            R_predict, V_predict = satf.r_init, satf.v_init

            er1.append( j2.norm(R_predict - r1))
            ev1.append( j2.norm(V_predict - v1))

            er2.append( j2.norm(R_predict - r2))
            ev2.append( j2.norm(V_predict - v2))

            er3.append( j2.norm(R_predict - r3))
            ev3.append( j2.norm(V_predict - v3))

            er4.append( j2.norm(R_predict - r4))
            ev4.append( j2.norm(V_predict - v4))

            er5.append( j2.norm(R_predict - r1))
            ev5.append( j2.norm(V_predict - v1))
            
            erp1.append( 100*(j2.norm(R_predict - r1))/j2.norm(R_predict) )
            evp1.append( 100*(j2.norm(V_predict - v1))/j2.norm(V_predict) )

            erp2.append( 100*(j2.norm(R_predict - r2))/j2.norm(R_predict) )
            evp2.append( 100*(j2.norm(V_predict - v2))/j2.norm(V_predict) )

            erp3.append( 100*(j2.norm(R_predict - r3))/j2.norm(R_predict) )
            evp3.append( 100*(j2.norm(V_predict - v3))/j2.norm(V_predict) )

            erp4.append( 100*(j2.norm(R_predict - r4))/j2.norm(R_predict) )
            evp4.append( 100*(j2.norm(V_predict - v4))/j2.norm(V_predict) )

            erp5.append( 100*(j2.norm(R_predict - r5))/j2.norm(R_predict) )
            evp5.append( 100*(j2.norm(V_predict - v5))/j2.norm(V_predict) )

        j = t_axis
        o_bstar = 'Actual BSTAR --> ' + init_drag
        plt.figure()
        plt.plot(j, erp1, 'o')
        plt.plot(j, erp2, 'o')
        plt.plot(j, erp3, 'o')
        plt.plot(j, erp4, 'o')
        plt.plot(j, erp5, 'o')
        plt.legend([o_bstar, new_drag1, new_drag2, new_drag3, new_drag4], fontsize = 'large')
        plt.plot(j, erp1, 'b-')
        plt.title(r'Norm of $\delta r (in  \%)$', fontdict = {'fontsize':24})
        plt.xlabel(r'Hours  $\rightarrow$', fontdict = {'fontsize':15})        
        plt.ylabel(r'$\frac{|\delta r| \times 100}{|R|} (in \%) \rightarrow$', fontdict = {'fontsize':22} )
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(j, evp1,'o')
        plt.plot(j, evp2,'o')
        plt.plot(j, evp3,'o')
        plt.plot(j, evp4,'o')
        plt.plot(j, evp5,'o')
        plt.title(r'Norm of $\delta v (in  \%)$', fontdict = {'fontsize':24})
        plt.legend([o_bstar, new_drag1, new_drag2, new_drag3, new_drag4], fontsize = 'large')
        plt.plot(j, evp1,'b-')        
        plt.xlabel(r'Hours $\rightarrow$', fontdict = {'fontsize':15})
        plt.ylabel(r'$\frac{|\delta v| \times 100} {|V|} (in \%) \rightarrow$', fontdict = {'fontsize':22})
        plt.grid()
        plt.show()
