


def import_MAG_data(year, month, day, resolution='1s'):
    
    '''
    Imports 1-day MAVEN MAG L2 calibrated magnetic field and spacecraft
    position time series in MSO (Sun State) coordinates.
    
    Parameters:
        
        year: int
        
        month: int
               without leading zeros
        
        day: int
             without leading zeros
        
        resolution: str
                    '1s' for low resolution data
                    '32ms' for high resolution data
    
    Returns:
        
        t_dec: N-1 array
               time vector in decimal hour
       
        B_mso: N-3 array
               magnetic field cartesian vector in MSO coordinates
               
        B: N-1 array
           magnetic field magnitude
           
       r_mso: N-3 array
              spacecraft position cartesian vector in MSO coordinates 
          
    '''
    
    
    import numpy as np
    from pathlib import Path
    from func_time import Time
    ftm = Time
    ftm0 = Time(year,month,day)
    
    
    
    doy = ftm0.date2doy()
    
    
    if resolution == '32ms': resolution = ''
    
    
    # complete with leading zeros for file name
    month = str(month).zfill(2)
    day = str(day).zfill(2)
    doy = str(doy).zfill(3)
    
   
    path1 = r'/Users/sofiaburne/Desktop/BS-MPB/maven_mag_calibrated/'
    
    path2 = r'{}/{}/'.format(year,month)
    
    path = path1+path2
   
    file_name = r'mvn_mag_l2_{}{}ss{}_{}{}{}_v01_r01.sts'.format(year,
                                                                 
                                                                 doy,
                                                                 
                                                                 resolution,
                                                                 
                                                                 year,
                                                                 
                                                                 month,
                                                                 
                                                                 day)
    
    
    # change r01 for r02 when necessary
   
    my_file = Path(path+file_name)
   
    if my_file.is_file() == False:
        
        file_name = r'mvn_mag_l2_{}{}ss{}_{}{}{}_v01_r02.sts'.format(year,
                                                                    
                                                                    doy,
                                                                    
                                                                    resolution,
                                                                    
                                                                    year,
                                                                    
                                                                    month,
                                                                    
                                                                    day)
   
    
    
    h,m,s,ms, Bx,By,Bz, x,y,z = np.loadtxt( path+file_name,
                                           
                                           skiprows=155,
                                           
                                           usecols=(2,3,4,5,7,8,9,11,12,13),
                                           
                                           unpack=True)
    
    
    t_dec = np.asarray( ftm.hms2dec(h,m,s+ms*1e-3) )
    
    B_mso = np.asarray([Bx,By,Bz]).T
    B = np.linalg.norm(B_mso, axis=1)
    
    R_Mars = 3389.5
    x, y, z = x/R_Mars, y/R_Mars, z/R_Mars
    r_mso = np.asarray([x,y,z]).T
    
    
    return t_dec, B_mso, B, r_mso
    


    
