
import numpy as np
import datetime


class Time:
    
    
    def __init__(self, year, month, day):
        
        '''
        Parameters:
            
            year: int
                  year in yyyy format
    
            month: int
                   month in m format (do not include extra left 0, e.g.: January is 1, not 01)
    
            day: int
                 day in d format (do not include extra left 0, e.g.: 1 not 01)        
            
        '''
        
        self.year = year
        self.month = month
        self.day = day


        

    
    def dech2datetime(self, dec_h):
        
        '''
        Transforms times in decimal hour format to datetime.datetime
        yyyy-mm-dd hh:mm:ss.f
        
        Parameters:
            
            dec_h: float or array (N,1)
                   times in dec hour
       
        Returns:
            
            dt: list
                times in yyyy-mm-dd hh:mm:ss.f format
                list elements are datetime.datetime objects
                
        '''
        
        if (type(dec_h) == float) or (type(dec_h) == np.float64):
            
            dec_h = np.asarray([dec_h])
        
        l = len(dec_h)
        
        
        dt = []
        
        for i in range(l):
            
            
            hour = int(dec_h[i])
            
            m = 60*(dec_h - hour)
            minute = int(m[i])
            
            s = 60*(m - minute)
            second = int(s[i])
            
            ms = (1e6)*(s - second)
            microsecond = int(ms[i])
            
            
            dt_i = datetime.datetime(self.year, self.month, self.day,
                                      
                                      hour, minute, second, microsecond)
            
            dt.append(dt_i)
        
        return dt
    
    
    
    
    
    
    def dech2time(dec_h):
        
        '''
        Transforms times in decimal hour format to datetime.datetime hh:mm:ss.f
        
        Parameters:
            
            dec_h: float or array (N,1)
                   times in dec hour
       
        Returns:
            
            dt: list
                times in hh:mm:ss.f format
                list elements are datetime.time objects
                
        '''
        
        if (type(dec_h) == float) or (type(dec_h) == np.float64):
            
            dec_h = np.asarray([dec_h])
        
        l = len(dec_h)
        
        
        dt = []
        
        for i in range(l):
            
            
            hour = int(dec_h[i])
            
            m = 60*(dec_h - hour)
            minute = int(m[i])
            
            s = 60*(m - minute)
            second = int(s[i])
            
            ms = (1e6)*(s - second)
            microsecond = int(ms[i])
            
            
            dt_i = datetime.time(hour, minute, second, microsecond)
            
            dt.append(dt_i)
        
        return dt
    
    
    
    
    


    
    
    
    def strtime2dech(time):
        
        '''
        Transfroms time in format 'hh:mm:ss.f' to decimal hour
        
        Parameters:
            
            time = str or list
                   time value in format 'hh:mm:ss.f'
        
        Returns:
            
            dech: float or array (N,1)
                  time in decimal hour
            
        '''
        
        if type(time) == str: time = [time]
        print(time)
        
        l = len(time)
        
        dech = []
        
        for i in range(l):
            
            t = datetime.datetime.strptime(time[i],'%H:%M:%S.%f')
            
            h = t.hour
            m = t.minute
            s = t.second
            ms = t.microsecond
            
            dech_i = h + m/60 + s/3600 + ms/(2.7778e10)
            
            dech.append(dech_i)
        
        return dech
    
    
    
    

    def strdatetime2dech(datetime):
        
        '''
        Transfroms date-time in format 'yyyy-mm-dd hh:mm:ss.f' to decimal hour
        
        Parameters:
            
            time = str or list
                   time value in format 'yyyy-mm-dd hh:mm:ss.f'
        
        Returns:
            
            dech: float or array (N,1)
                  time in decimal hour
            
        '''
        
        t = datetime.datetime.strptime(datetime,'%H:%M:%S.%f')
        
        h = t.hour
        m = t.minute
        s = t.second
        ms = t.microsecond
        
        dech = h + m/60 + s/3600 + ms/(2.7778e10)
        
        return dech

    
    
    
    def hms2dec(h, m, s):
        
        '''
        Gives decimal time value.
        
        Parameters:
            
            h = float orlist of floats
                hour value/s
            
            m = float orlist of floats
                minute value/s
            
            s = float orlist of floats
                second value/s
        
        Returns:
            
            time in decimal hour format
            
        '''
        
        if (type(h) == float) or (type(h) == int or type(h) == np.float64):
            h = [h]
            m = [m]
            s = [s]
        
        else:
            h = list(h)
            m = list(m)
            s = list(s)
            
        t = [h[i] + m[i]/60 + s[i]/3600 for i in range(len(h))]
        
        
        return t   

        




    def doy2date(year, doy):
    
        '''
        Obtains month and day from day of year.
        
        Parameters:
            year: float
            
            doy: float
        
        Returns:
            month: float
            
            day: float
        '''
        
        
        d = datetime.date(year, 1, 1) + datetime.timedelta(doy - 1)
        
        month = d.month
        day = d.day
        
        return month, day
    
    
    
    
    
    
    # def date2doy(year, month, day):
        
    #     '''
    #     Returns day of year from date
        
    #     Parameters:
            
    #         year: float
            
    #         month: float
            
    #         day: float
        
    #     Returns:
            
    #         doy: foat
    #     '''
        
    #     doy = datetime.datetime(year,month,day).timetuple().tm_yday
        
    #     return doy
        
    
    
    
    
    def date2doy(self):
        
        '''
        Returns day of year from date
        
        
        Returns:
            
            doy: foat
        '''
        
        doy = datetime.datetime(self.year,
                                
                                self.month,
                                
                                self.day).timetuple().tm_yday
        
        return doy
        
        
        
        
 