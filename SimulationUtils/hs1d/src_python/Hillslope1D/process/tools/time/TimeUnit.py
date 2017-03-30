class TimeUnit(object):
    """
        Class managing the conversion between time units dor the calculation
        ######################################################################
        @param
          tmax : the maximal time of the timeserie
          unit : string corresponding to the time unit
            choices : year, days, hour, min, sec

        @attributes
          tmax : the maximal time of the timeserie
          unit : string corresponding to the time unit
        ######################################################################
    """

    def __init__(self, tmax=86400*35, unit='sec'):
        self.tmax = tmax
        self.unit = unit


    def time_to_seconds(self,tmax):
        """
            Conversion between each time unit to seconds
        """
        sec_in_min = 60
        sec_in_hour = sec_in_min*60
        sec_in_day = sec_in_hour*24
        sec_in_year = sec_in_day*365

        if self.unit == '':
            self.unit = ('enter time units - seconds, minutes, days, year')
            if self.unit == 'days':
                tmax = tmax*sec_in_day
            elif self.unit == 'hour':
                tmax = tmax*sec_in_hour
            elif self.unit == 'min':
                tmax = tmax*sec_in_min
            elif self.unit == 'sec':
                tmax = tmax
            elif self.unit == 'year':
                tmax = tmax*sec_in_year
            else:
                print('UNIT UNKNOWN #s\n', self.unit)

        elif self.unit == 'year':
            tmax = self.tmax*sec_in_year
        elif self.unit == 'days':
            tmax = tmax * sec_in_day
        elif self.unit == 'hour':
            tmax = tmax * sec_in_hour
        elif self.unit == 'min':
            tmax = tmax * sec_in_min
        elif self.unit == 'sec':
            tmax = tmax
        else:
            print('UNIT UNKNOWN #s\n', self.unit)

        return tmax

    def time_to_days(self):
        sec_in_day = 86400
        min_in_day = 1440
        hour_in_day = 24
        day_in_year = 365


        if self.unit == '':
            self.unit = ('enter time units - seconds, minutes, days, year')
            if self.unit == 'hour':
                self.tmax = self.tmax/hour_in_day
            elif self.unit == 'min':
                self.tmax = self.tmax/min_in_day
            elif self.unit == 'sec':
                self.tmax = self.tmax/sec_in_day
            elif self.unit == 'year':
                self.tmax = self.tmax*day_in_year
            else:
                print('UNIT UNKNOWN #s\n', self.unit)

        elif self.unit == 'sec':
            self.tmax = self.tmax/self.sec_in_day
        elif self.unit == 'min':
            self.tmax = self.tmax / min_in_day
        elif self.unit == 'hour':
            self.tmax = self.tmax / hour_in_day
        elif self.unit == 'year':
            self.tmax = self.tmax*day_in_year
        else:
            print('UNIT UNKNOWN #s\n', self.unit)

        return self.tmax

    def time_to_years(self):
        days_in_year = 365
        hours_in_year = 24*days_in_year
        min_in_year = 60*hours_in_year
        sec_in_year = 60*min_in_year


        if self.unit == '':
            self.unit = ('enter time units - seconds, minutes, days, year')
            if self.unit == 'days':
                self.tmax = self.tmax*days_in_year
            elif self.unit == 'hour':
                self.tmax = self.tmax*hours_in_year
            elif self.unit == 'min':
                self.tmax = self.tmax*min_in_year
            elif self.unit == 'sec':
                self.tmax = self.tmax
            else:
                print('UNIT UNKNOWN #s\n', self.unit)

        elif self.unit == 'sec':
            self.tmax = self.tmax / sec_in_year
        elif self.unit == 'min':
            self.tmax = self.tmax / min_in_year
        elif self.unit == 'hour':
            self.tmax = self.tmax / hours_in_year
        elif self.unit == 'days':
            self.tmax = self.tmax / days_in_year
        else:
            print('UNIT UNKNOWN #s\n', self.unit)

        return self.tmax

    def time_to_hours(self):
        hours_in_day = 24
        hours_in_year = 365*hours_in_day
        min_in_hour = 60
        sec_in_hour = 60*min_in_hour

        if self.unit == '':
            self.unit = ('enter time units - seconds, minutes, days, year')
            if self.unit == 'year':
                self.tmax = self.tmax*hours_in_year
            elif self.unit == 'days':
                self.tmax = self.tmax*hours_in_day
            elif self.unit == 'hour':
                self.tmax = self.tmax
            elif self.unit == 'min':
                self.tmax = self.tmax/min_in_hour
            elif self.unit == 'sec':
                self.tmax = self.tmax/sec_in_hour
            else:
                print('UNIT UNKNOWN #s\n', self.unit)

        elif self.unit == 'sec':
            self.tmax = self.tmax / sec_in_hour
        elif self.unit == 'min':
            self.tmax = self.tmax / min_in_hour
        elif self.unit == 'hour':
            self.tmax = self.tmax
        elif self.unit == 'days':
            self.tmax = self.tmax * hours_in_day
        elif self.unit == 'year':
            self.tmax = self.tmax * hours_in_year
        else:
            print('UNIT UNKNOWN #s\n', self.unit)

        return self.tmax
