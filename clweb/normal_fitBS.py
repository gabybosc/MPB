import numpy as np



def norm_fit(R_cross, L, eps, x0):
    
    '''
    Calculates bowshock external normal vector (using a given bowshock model) at the crossing point.
    
    Parameters:
        
        R_cross: array
                 bowshock crossing point in aberrated MSO coordinates
        
        L: float
           semi-latus rectum of the conic
        
        eps: float
             eccentricity of the conic
        
        x0: float
            focus point of the conic
    
    
    Returns:
        
        n_fit: array
               bowshock external normal vector at the crossing point in MSO aberrated coordinates. 
        
        normal: array
                normal vector n_fit before normalization

    '''


    a = L/(eps**2 - 1)
    b = L/(eps**2 - 1)**(1/2)
    c = x0 + L*eps/(eps**2 - 1)
    
    normal = np.array([ float(2*(R_cross[0]-c)/a**2), float(-2*R_cross[1]/b**2), float(- 2*R_cross[2]/b**2) ])
    n = np.linalg.norm(normal)
    n_fit = normal/n
    
    # define external normal
    if n_fit[0] < 0:
        n_fit = - n_fit
        
    return n_fit, normal












def err_norm_fit(R_cross, err_R_cross, L, eps, x0, param, err_param, normal):
    
    '''
    Returns the normal vector error, for the normal to the bowshock curve at the crossing point.
    
    Parameters:
        
        R_cross: array
                 bowshock crossing point in aberrated MSO coordinates
        
        err_R_cross: array
                     error of R_cross
        
        L: float
           semilactus rectum of the bowshock curve
           if this is the parameter that has been readjusted from the original model,
           then the re-adjusted value should be entered, and param should be 'L'.

        eps: float
             eccentricity of the bowshock curve
             if this is the parameter that has been readjusted from the original model,
             then the re-adjusted value should be entered, and param should be 'eps'.    
 
        x0: float
            focus of the bowshock curve
            if this is the parameter that has been readjusted from the original model,
            then the re-adjusted value should be entered, and param should be 'x0'.    
        
        parm: str
              parameter of the bowshock curve that has been readjusted from orginal model
              it take values: 'L' or 'eps' or 'x0', or 'None' if no adjustments were done.
        
        err_param: float
                   error value of the parameter that has been readjusted (0 if no adjustments were done)
        
        normal: array
                normal vector before normalization
    
    Returns:
        
        err_N: array
               error values of the normal vector in MSO aberrated coordinates
    '''
    
    
    # rename
    Rc = R_cross
    err_Rc = err_R_cross
    
    
    # using that N = [nx,ny,nz]/n
    
    # derivatives of nx, ny, nz and n on x, y, z
    
    dxnx = 2/L**2 * (eps**2 - 1)**2
    dyny = -2/L**2 * (eps - 1)
    dznz = dyny
    
    arg = (eps**2 - 1)**4 * (x0 - Rc[0] + L*eps/(eps**2 - 1))**2 + (eps - 1)**2*Rc[1]**2 + (eps - 1)**2*Rc[2]**2
    
    dxn = 2 * (eps**2 - 1)**4 * (Rc[0] - x0 - L*eps/(eps**2 - 1)) / ( L**2 * np.sqrt(np.float64(arg)) )
    dyn = 2 * (eps - 1)**2 * Rc[1] / ( L**2 * np.sqrt(np.float64(arg)) )
    dzn = 2 * (eps - 1)**2 * Rc[2] / ( L**2 * np.sqrt(np.float64(arg)) )     
    
    
    # derivatives of nx, ny, nz and n on L / eps / x0

    if (param == 'L'):  
    
        dpnx = -2/L**3 * ( (eps**2 - 1) * ( 2*x0 + 2*eps**2 *(Rc[0] - x0) - 2*Rc[0] - eps*L ) ) 
        dpny = 4*Rc[1]/L**3 * (eps - 1)
        dpnz = 4*Rc[2]/L**3 * (eps - 1)
        dpn = 2/L**3 * ( -eps*L*(eps**2 - 1)**3 * (Rc[0] - x0 - eps*L/(eps**2 - 1)) -2 *arg ) / np.sqrt(np.float64(arg))

    if (param == 'eps'):  
    
        dpnx = 2/L**2 * ( L - 3*eps**2*L + 4*eps*(eps**2 - 1) *(Rc[0] - x0) )
        dpny = -2*Rc[1]/L**2
        dpnz = -2*Rc[2]/L**2
        dpn = 2/L**2 * ( L*(eps**4 - 1) *(-eps*L - Rc[0] + eps**2 *(Rc[0] - x0) + x0) + 4*eps*(eps**2 - 1) *(Rc[0] - x0 + eps*(L + eps*(x0 - Rc[0]) ) ) + (Rc[1]**2 + Rc[2]**2)*(eps-1) ) / np.sqrt(np.float64(arg))
    
    if (param == 'x0'):  
    
        dpnx = -2/L**2 * (eps**2 - 1)**2
        dpny = 0
        dpnz = 0
        dpn = -2/L**2 * (eps**2 - 1)**4 * (Rc[0] - x0 - L*eps/(eps**2 - 1)) / np.sqrt(np.float64(arg))
    
    if (param == 'None'):  
    
        dpnx = 0
        dpny = 0
        dpnz = 0
        dpn = 0
    
    
    # errors of nx, ny, nz and n
    err_nx = np.sqrt( np.float64( (dxnx * err_Rc[0])**2 + (dpnx * err_param)**2 ) )
    err_ny = np.sqrt( np.float64( (dyny * err_Rc[1])**2 + (dpny * err_param)**2 ) )
    err_nz = np.sqrt( np.float64( (dznz * err_Rc[2])**2 + (dpnz * err_param)**2 ) )
    err_n = np.sqrt( np.float64( (dxn * err_Rc[0])**2 + (dyn * err_Rc[1])**2 + (dzn * err_Rc[2])**2 + (dpn * err_param)**2 ) )    
    
    nx = normal[0]
    ny = normal[1]
    nz = normal[2]
    n = np.linalg.norm(normal)
    
    # errors of each component of N
    err_Nx = np.sqrt( np.float64( (err_nx/n)**2 + (nx*err_n)**2 ) )
    err_Ny = np.sqrt( np.float64( (err_ny/n)**2 + (ny*err_n)**2 ) )
    err_Nz = np.sqrt( np.float64( (err_nz/n)**2 + (nz*err_n)**2 ) )
        
    err_N = np.array([err_Nx, err_Ny, err_Nz])
    
    return err_N
