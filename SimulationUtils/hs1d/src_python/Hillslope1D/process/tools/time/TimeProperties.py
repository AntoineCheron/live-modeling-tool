import numpy as np
import TimeUnit as TU

class TimeProperties(object):
    """
        Properties management for time (unit and values) used to describe the
        hillslope and its paramters variations
        #######################################################################
        @param
            tmin : minimal time value (usually 0)
            tmax : maximal time value
            Nt : number of time steps from tmin to tmax
            unit : string corresponding to the time unit of the time serie
                choices : year, days, hour, min, sec
            time_custom : -1 or a vector of time values if databased recharge

        @attributes
          t : time vector
          Nt : number of time values (=len(t))
          tmax : maximal time value
          tmin : minimal time value
          unit : time unit
          TU : TimeUnit class
        #######################################################################
    """


    def __init__(self, tmin=0, tmax=35, Nt=35*24*10, unit='days', time_custom=-1):
        self.tmin = tmin
        self.Nt = Nt
        self.unit = unit
        self.TU = TU.TimeUnit(tmax,unit)
        self.tmax = self.TU.time_to_seconds(tmax)
        if isinstance(time_custom,int):
            self.time_properties()
        else:
            self.t = time_custom

    def time_properties(self):
        """
            Creates a time vector in sec, based on tmin, tmax and Nt (after
            conversion from unit to seconds)
        """
        self.t = np.arange(self.tmin, self.tmax, (self.tmax-self.tmin)/(self.Nt-1))
        self.t = np.reshape(self.t,(len(self.t),1))
        return self.t


