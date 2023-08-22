import numpy as np


def cartesian2polar(x,y,z,r0,sym_axis='x'):
    
    '''
    Obtain polar coordinates from cartesian components.
    
    Parameters:
        x,y,z: float
               position of 1 point in conic surface measured from origin
               
        r0: (1,3) array
            focus of conic surface
            
        sym_axis: str
                  axis of symmetry for conic surface (e.g.: 'x')
        
    Returns:
        r:  float
            distance of point in conic surface measured from focus
            
        theta: float
               polar angle measured from symmetry axis
               
        phi: float
             azimuth angle measured in plane perpendicular to symmetry axis
    '''
    
    X = x - r0[0]
    Y = y - r0[1]
    Z = z - r0[2]
    
    r = np.sqrt( X**2 + Y**2 + Z**2 )
    
    if sym_axis == 'x':
        theta = np.arccos(X/r)
        phi = np.arctan2(Z,Y)
        
    elif sym_axis == 'y':
        theta = np.arccos(Y/r)
        phi = np.arctan2(X,Z)
    
    elif sym_axis == 'z':
        theta = np.arccos(Z/r)
        phi = np.arctan2(Y,X)
    
    
    return r, theta, phi




def fit_L(r,theta,eps=1.03):
    
    '''
    Gives L from conic surface fitting a given point.
    
    Parameters:
        r: float
           distance of point in conic surface measured from focus
           in units of Martian radii
           
        theta: float
               polar angle measured from symmetry axis
               in radians
               
        eps: float
             eccentricity of conic surface
    
    Returns:
        L: float
           semi-lactus rectum for conic surface containing given point
           in units of Martian radii
    
    '''
    
    L = r * (1 + eps*np.cos(theta))
    
    return L




def Rsd(L,x0=0.64,eps=1.03):
    
    '''
    Gives standoff distance from Vignes model.
    
    Parameters:
        L: float
           semi-lactus rectum for conic surface containing specific point
           in units of Martian radii
           
    Returns:
        Rsd: float
             standoff distance
             in units of Martian radii
    
    '''
    
    Rsd = x0 + L/(1 + eps)
    
    return Rsd
    
    
    
def Rtd(L,y,z,x0=0.64,eps=1.03):
    
    '''
    Gives terminator standoff distance from Vignes model.
    
    Parameters:
        L: float
           semi-lactus rectum for conic surface containing specific point
           in units of Martian radii
        
        y,z: float
             position of point in conic surface measured from origin
             in units of Martian radii
               
    Returns:
        Rtd: float
             terminator standoff distance
             in units of Martian radii
    '''

    #phi = np.arccos( x0/np.sqrt(y**2 + z**2) ) # phi = theta(x=0)
    phi = np.pi - np.arctan2( np.sqrt(y**2 + z**2), x0 )
    
    Rtd = np.sin(phi) * L/( 1 + eps*np.cos(phi) )
    
    return Rtd 
    
    