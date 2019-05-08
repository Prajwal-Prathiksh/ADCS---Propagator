import numpy as np


#-----Differential Equation - Function------
#        d(x1) = x2
#        -----           ,    x1 = r,   x2 = v
#         dt        

def f1(r, vel):   
    post = vel
    return post



#-----Differential Equation - Function-s----
#        d(x2) =  -u . x1       P             P 
#        -----    -------   +     oblate  +     drag      , x1 = r, x2 = v
#         dt       (x1)^3   
def f2(r,vel):   
    if dragBoolean == True:
        p_o    = oblate_pertubations(r)
        p_d    = drag(r,vel)
        accltn =(((-mu/((norm(r))**3)))*r) + p_o  + p_d 
        return accltn 
    else:
        p_o    = oblate_pertubations(r)
        accltn =(((-mu/((norm(r))**3)))*r) + p_o 
        return accltn   
    


#-----Calculates the Pertubations from the Effects of Oblateness of Earth-----
#			 through J2 Model
    
def oblate_pertubations(post):   
    '''
        J2 Perturbations: Pg 664 HD Curtis
        
    '''
    r = norm(post)
    const_J2 = (1.5*J2*mu*Radius_Earth**2) / (r**4)
    z = post[2]
    x, y, z = post
    k = ((  5*((z/r)**2)  - 1 )/r)
    
    p = np.array([x*k,  y*k,  z*((  5*((z/r)**2)  - 3 )/r)])
    p = p*const_J2
    return p



#-----Calculates the Pertubations from the Effects of Atmospheric Drag-----    
def drag(r,vel):  
    '''
        Assuming velocity of atmosphere at a particular point,
        is appx. equal to cross product of angular velocity of Earth 
        and position vector of that point.
    '''
    vel = np.array(vel, dtype = np.float64)*1000 #Converting to m/s
    r = np.array(r, dtype = np.float64) * 1000 #Converting to m
    
    v_atm = np.cross(w_angular_vel,r)
    v_rel = vel - v_atm #Relative velocity being calculated.
    
    denst = density(r/1000) #Requires argument in Km
    p = (-0.5) * denst * norm(v_rel) * B_coeff * v_rel #Accleration in m/s
    return (p/1000) #Returning Accleration in km/s



#-----Calulates the atomosphereic density from a model based on 
#       exponential decay, where the input is the height
#           of the satellite from surface of the Earth,
#      and the output is density of the atmosphere in kg/m^-3
#  NOTE: The model is made specifically for LEO, i.e, h ranging from 200km to 1000km.

def density(r):
    h = norm(r) - Radius_Earth #Appx height of satellite from surface of Earth.
    
    I, alpha1, alpha2, alpha3, alpha4, beta = (-54.23161450332475,  22123.13823106815, -6906047.242456764, 
                                                           1057730365.397583, -61951792129.727295, 0.004345148162663792)
    
    rho = np.exp(I + beta*h + alpha1/h + alpha2/(h**2) + alpha3/(h**3) + alpha4/(h**4) )
    return rho #In kg/m^3



#-----The Propogator involved in the Orbit Modelling that uses the-----
#         Initial State Vectors and RK - 4 to approximate the 
#		subsequent system of State Variables.
#         RK - 4: X(n+1) = X(n) + h      (a + 2*b + 2*c + d)
#                                ---  *  
#                                 6
def propagate(r,vel, time, h_step_size = 1, drag=False):
    """Propogate the State Variables.
    
       Parameters:
       -----------
       r: list of length (3)
            This initializes the position of the satellite.
          
       vel: list of length (3)
            This initializes the velocity of the satellite.
            
       time: floating-point number
            This initializes the time interval over which the propagator
            runs, in seconds.
            
       h_step_size: floating-point number, optional
            This sets the value of the step size, for RK-4.
            Default = 1
        
       drag: boolean, optional
            Takes Drag into consideration as a perturbation as well.  
            Default = False      
    """    
    x,y = np.array(r,np.float64),np.array(vel,np.float64)
    steps = int(time/h_step_size)
    h = h_step_size
    global dragBoolean
    dragBoolean = drag
    
    for i in range(steps):
        ax,ay = RK4_a(x,y,h)
        bx,by = RK4_b(x,y,ax,ay,h)
        cx,cy = RK4_c(x,y,bx,by,h)
        dx,dy = RK4_d(x,y,cx,cy,h)
        x = x + ((h)*(ax + 2*(bx + cx) + dx))/6
        y = y + ((h)*(ay + 2*(by + cy) + dy))/6

    return x,y

    
def RK4_a(r,vel,h):
    """ To Calculate Value of (a) of RK-4:
        a = f(x).
    """
    ax,ay = np.zeros(3), np.zeros(3)
    ax = f1(r,vel)
    ay = f2(r,vel)
    return ax,ay
    
def RK4_b(r,vel,ax,ay,h):
    """ To Calculate Value of (b) of RK-4:
        b = f(x + (h/2)a)
    """
    l,m = np.zeros(3), np.zeros(3)
    l =  r + ((h/2)*ax)
    m =  vel + ((h/2)*ay)
    bx,by = f1(l,m),f2(l,m)
    return bx,by
    
def RK4_c(r,vel,bx,by,h):
    """ To Calculate Value of (c) of RK-4:
        c = f(x + (h/2)b)
    """
    l,m = np.zeros(3), np.zeros(3)
    l =  r + ((h/2)*bx)
    m =  vel + ((h/2)*by)
    cx,cy = f1(l,m),f2(l,m)
    return cx,cy
    
def RK4_d(r,vel,cx,cy,h):
    """ To Calculate Value of (d) of RK-4:
        d = f(x + (h)c)
    """
    l,m = np.zeros(3), np.zeros(3)
    l =  r + ((h)*cx)
    m =  vel + ((h)*cy)
    dx,dy = f1(l,m),f2(l,m)
    return dx,dy
    
def norm(arr):
    ''' Calculates 2-norm of a vector.
    '''
    return (arr[0]**2 + arr[1]**2 + arr[2]**2)**0.5
    



#-----Values of Orbital Constants being assigned-----
mu = 398600.4415
J2 = 1.082635854e-3 
w_angular_vel = np.array([0,0,7.2921156e-5]) #Radians/sec
c_drag = 2
area = 0.01 #In m^3
mass = 0.9 #In Kg
Radius_Earth = 6378.1363          #In Km.
B_coeff = (c_drag*area)/mass
