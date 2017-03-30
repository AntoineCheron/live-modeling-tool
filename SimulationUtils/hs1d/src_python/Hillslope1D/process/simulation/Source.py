import TimeProperties as TP
import numpy as np

class Source(object):
  """
      Class managing temporal aspect and source terms (recharge)

      @param
        period : period used to compute the recharge chronicle (if not 'databased')
        recharge_type : type of used recharge ('periodical, 'square', 'steady', 'random'
                                               'databased')
        recharge_rate : either a float or a vector of recharge (m/s)
        tmin : minimal time value of the serie
        tmax : maximal time value of the serie
        Nt : number of time values
        unit : time unit ('days', 'hour', 'min', 'sec', 'year')
        time_custom : -1 or a vector containing time values if 'databased' recharge

      @attributes
        recharge : recharge value at the asked time
        period : period used to compute the recharge chronicle (if not 'databased')
        tmax : maximal time value of the chronicle
        recharge_type : type of used recharge
        t_chronicle : a vector containing time values of the chronicle
        recharge_chronicle : all recharge values
        TP : TimeProperties class

  """

  def __init__(self, period, recharge_type, recharge_rate, tmin, tmax, Nt, unit, time_custom):
      self.recharge_type = recharge_type
      self.TP = TP.TimeProperties(tmin, tmax, Nt, unit, time_custom)
      self.tmax = self.TP.TU.time_to_seconds(tmax)
      self.period = period
      self.period = self.source_terms()
      self.set_recharge_chronicle(recharge_rate)

  def source_terms(self):

      if self.recharge_type == 'periodical' or self.recharge_type == 'square':
          if self.period is None:
              self.period = self.TP.TU.time_to_seconds(5)

      elif self.recharge_type == 'steady' or self.recharge_type is None:
          self.period = 'inf'

      elif self.recharge_type == 'random':
            self.period = 0

      elif self.recharge_type == 'databased':
          self.period = -1

      return self.period

  def set_recharge_chronicle(self, recharge_rate):

  #        print('time=',self.TP.t)
      if self.period == 'inf' or self.period is None:
          recharge_rate = (recharge_rate*10**-3)/86400
          self.recharge_rate = recharge_rate

          self.recharge_chronicle = self.recharge_rate*np.ones(self.TP.t)

      elif self.period == 0:
          self.recharge_chronicle = self.set_random_recharge()

      elif self.period == -1:
          self.recharge_chronicle = recharge_rate

      else:
          recharge_rate = (recharge_rate*10**-3)/86400
          self.recharge_rate = recharge_rate
          if self.recharge_type == 'periodical':
              self.recharge_chronicle = self.recharge_rate*(1+np.cos(2*np.pi*(self.TP.t/self.period)))
          elif self.recharge_type == 'square':
              Int_ = np.floor(self.TP.t/self.period)
              Odd_rest = Int_%2
              self.recharge_chronicle = self.recharge_rate*Odd_rest
              rem_stiff = self.TP.t%(2*self.period)
              bool1 = rem_stiff < 60
              self.recharge_chronicle[bool1] = self.recharge_rate*(1-(rem_stiff[bool1])/60)
              rem_stiff = (self.TP.t + self.period)%(2*self.period)
              bool1 = rem_stiff < 60
              self.recharge_chronicle[bool1] = self.recharge_rate*((rem_stiff[bool1])/60)

      return self.recharge_chronicle

  def compute_recharge_rate(self, t):
    if self.period == 'inf':
        self.recharge = self.recharge_rate

    elif self.period == 0:
        self.t_chronicle = self.TP.time_properties()

        self.recharge = np.interp(t, self.t_chronicle, self.recharge_chronicle)

    elif self.period == -1:
        self.t_chronicle = self.TP.t
        self.recharge = np.interp(t, self.t_chronicle, self.recharge_chronicle)

    elif self.recharge_type == 'periodical':
        self.recharge = self.recharge_rate*(1+np.cos(2*np.pi*(t/self.period)))
    elif self.recharge_type == 'square':
        Int_ = np.floor(t/self.period)
        Odd_rest = Int_%2
        recharge = self.recharge_rate*Odd_rest
        rem_stiff = t%(2*self.period)
        bool1 = rem_stiff < 60
        if bool1 is True:
            self.recharge = self.recharge_rate*(1-(rem_stiff)/60)
        rem_stiff = (t + self.period)%(2*self.period)
        bool1 = rem_stiff < 60
        if bool1 is True:
            self.recharge = self.recharge_rate*((rem_stiff)/60)
    return self.recharge
