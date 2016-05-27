from EXOSIMS.Observatory.WFIRSTObservatory import WFIRSTObservatory
import astropy.units as u
from astropy.time import Time
import numpy as np
import os, inspect
import scipy.interpolate as interpolate
try:
    import cPickle as pickle
except:
    import pickle

class WFIRSTObservatoryL2(WFIRSTObservatory):
    """ WFIRST Observatory at L2 implementation. 
    Only difference between this and the WFIRSTObservatory implementation
    is the orbit method, and carrying an internal missionStartTime value,
    which should be the same as the one in TimeKeeping (necessary due to 
    orbit method interface, which calls for absolute time).

    Orbit is stored in pickled dictionary on disk (generated by MATLAB
    code adapted from E. Kolemen (2008).  Describes approx. 6 month halo
    which is then patched for the entire mission duration).
    """

    def __init__(self, missionStart=60634., orbit_datapath=None, **specs):
        
        #Run prototype constructor __init__ 
        WFIRSTObservatory.__init__(self,**specs)
        
        #set own missionStart value
        self.missionStart = Time(float(missionStart), format='mjd')
        
        #find and load halo orbit datat        
        #data is al in heliocentric ecliptic coords.  time is 2\pi = 1 sidreal year
        if orbit_datapath is None:
            classpath = os.path.split(inspect.getfile(self.__class__))[0]
            filename = 'L2_halo_orbit_six_month.p'
            orbit_datapath = os.path.join(classpath, filename)
        if not os.path.exists(orbit_datapath):
            raise Exception("Orbit data file not found.")
        halo = pickle.load( open( orbit_datapath, "rb" ) )
        
        #unpack orbit properties
        self.orbit_period = halo['te'].flatten()[0]/(2*np.pi)*u.year
        self.L2_dist = halo['x_lpoint'][0][0]*u.AU
        self.orbit_pos = halo['state'][:,0:3]*u.AU #this is in rotating frame wrt the sun
        self.orbit_pos[:,0] -= self.L2_dist #now with respect to L2, still in rotating frame
        self.orbit_time = halo['t'].flatten()/(2*np.pi)*u.year
        #create interpolant (years & AU units)
        self.orbit_interp = interpolate.interp1d(self.orbit_time.value,\
                self.orbit_pos.value.T,kind='cubic')

    def orbit(self, time):
        """Finds WFIRST L2 Halo orbit position vector
        
        This method finds the WFIRST L2 Halo rbit position 
        vector as 1D numpy array (astropy Quantity with units of km) in the
        heliocentric equatorial frame, stores this vector in self.r_sc,
        and returns True if successful.
        
        Args:
            time (Time):
                current absolute time (astropy Time)
            
        Returns:
            success (bool):
                True if successful, False if not
        
        """
        
        #find time from mission start and interpolated position
        deltime = (time - self.missionStart).to('year')
        cpos = self.orbit_interp(np.mod(deltime,self.orbit_period).value)
        
        #add L2 position to get current ecliptic coord
        th = np.mod(deltime.value,1.)*2*np.pi
        cpos += np.array([np.cos(th),np.sin(th),0])*self.L2_dist.to('AU').value
        
        #finally, rotate into equatorial plane
        obe = self.obe(self.cent(time))
        cpos = (np.dot(self.rot(np.radians(-obe),1),cpos)*u.AU).to('km')
        
        self.r_sc = cpos
        
        return np.all(np.isfinite(self.r_sc))

