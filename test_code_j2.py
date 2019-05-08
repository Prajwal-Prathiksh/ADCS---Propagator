import J2_Propagator as j2
import numpy as np

#-----Test Cases-----

circular_orbit = dict( r = [7000.0, 0.0, 0], v = [0,7.546,0], T = 5828.0) 
pratham_orbit = dict( r = [7088.31342, 55.24511, 0.30484], v = [-0.0074, -1.0677, 7.4116], T = 5906.8)

error_percentage = 0.05


#-----Test Code-----

#Circular - J2, 3 Hrs
r_out, v_out = j2.propagate(circular_orbit['r'], circular_orbit['v'],(3*3600), drag=False)
r_check, v_check = np.array([4403.45981728984, -5436.756085742211, 0.1557151272128546]), np.array([5.872054960160551, 4.745602423244825, -0.0001268547458108781])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Circular Orbit, Time : 3 Hrs, J2'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Circular Orbit, Time : 3 Hrs, J2'

#Circular - J2 + Drag, 3 Hrs
r_out, v_out = j2.propagate(circular_orbit['r'], circular_orbit['v'],(3*3600), drag=True)
r_check, v_check = np.array([4403.475840562699, -5436.73915685563, 0.1557146629250469]), np.array([5.872041201881828, 4.745622736401809, -0.0001268555239062031])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Circular Orbit, Time : 3 Hrs, J2 + Drag'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Circular Orbit, Time : 3 Hrs, J2 + Drag'

#Circular - J2, 3 Orbits
r_out, v_out = j2.propagate(circular_orbit['r'], circular_orbit['v'], 3 * circular_orbit['T'], drag=False)
r_check, v_check = np.array([6991.386254963272, 347.1014884331252, -0.0239416448049494]), np.array([-0.3744231916738881, 7.536708053887089, -0.0003559105360160409])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Circular Orbit, Time : 3 Orbits, J2'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Circular Orbit, Time : 3 Orbits, J2'

#Circular - J2 + Drag, 3 Orbits
r_out, v_out = j2.propagate(circular_orbit['r'], circular_orbit['v'], 3 * circular_orbit['T'], drag=True)
r_check, v_check = np.array([6991.379769463948, 347.1590873904567, -0.02394452090905394]), np.array([-0.3744842584502502,7.536706975549463,-0.0003559109161584875])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Circular Orbit, Time : 3 Orbits, J2 + Drag'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Circular Orbit, Time : 3 Orbits, J2 + Drag'


#Pratham - J2, 3 Hrs
r_out, v_out = j2.propagate(pratham_orbit['r'], pratham_orbit['v'], (3*3600), drag=False)
r_check, v_check = np.array([3353.28587029678, 924.09143819921, -6174.42472030125]), np.array([6.599485418801074, -0.4382185522862251, 3.498972757839065])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Pratham Orbit, Time : 3 Hrs, J2'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Pratham Orbit, Time : 3 Hrs, J2'

#Pratham - J2 + Drag, 3 Hrs
r_out, v_out = j2.propagate(pratham_orbit['r'], pratham_orbit['v'], (3*3600), drag=True)
r_check, v_check = np.array([3353.291963072318, 924.0908544630976, -6174.420323733792]), np.array([6.599482278341525, -0.4382196079299935, 3.498979935728481])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Pratham Orbit, Time : 3 Hrs, J2 + Drag'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Pratham Orbit, Time : 3 Hrs, J2 + Drag'

#Pratham - J2, 3 Orbits
r_out, v_out = j2.propagate(pratham_orbit['r'], pratham_orbit['v'], 3*pratham_orbit['T'], drag=False)
r_check, v_check = np.array([7087.580295708112, 90.72198955101115, -70.55294389370312]), np.array([0.07206538785026177, -1.066879121199723, 7.411388806354198])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Pratham Orbit, Time : 3 Orbits, J2'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Pratham Orbit, Time : 3 Orbits, J2'

#Pratham - J2 + Drag, 3 Orbits
r_out, v_out = j2.propagate(pratham_orbit['r'], pratham_orbit['v'], 3*pratham_orbit['T'], drag=True)
r_check, v_check = np.array([7087.579092275149, 90.71912418507162, -70.53311043601157]), np.array([0.07204441261305872, -1.066879474760484, 7.41138977384836])
r_error, v_error = ( j2.norm( r_out - r_check ) / j2.norm(r_check) ),( j2.norm( v_out - v_check ) / j2.norm(v_check) ) 
assert(r_error < (error_percentage/100) ), 'Error: Position, Pratham Orbit, Time : 3 Orbits, J2'
assert(v_error < (error_percentage/100) ), 'Error: Velocity, Pratham Orbit, Time : 3 Orbits, J2'

print('Code Successful')