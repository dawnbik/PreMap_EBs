Danbi Section

def download_lightcurve(name):
    """Downloads the lightcurve for a given object. Returns the first table in the search result.
    Uses lightkurve package."""
    search = lk.search_lightcurve(name)
    # search.table["dataURL"] = search.table["dataURI"]
    return search[0].download()

# def find_period(table):
#     """Returns the period of a binary system given the lightcurve. Uses LombScargle from astropy.timeseries"""
#     freq, power = LombScargle(table.time.value, table.flux.value).autopower()
#     max_power = np.max(power)
#     freq_at_max_power = freq[np.argmax(power)]
#     period = 1.0 / freq_at_max_power
#     return period

def find_period(name):
    data = ascii.read('rnaasac6e42t1_mrt.txt')
    period = data[name]['BLS-Period']
    return period






Amelia Section
def distance_from_simbad(tic_id):
    '''
    Finds the distance to the system using Simbad and the tic id.
    
    Parameters
    ----------
    tic_id : tic id of the star system in Simbad
    
    Returns
    -------
    Parsecs
        The distance to the star in Parsecs.
        
    '''
    custom_simbad = Simbad()
    
    custom_simbad.add_votable_fields('parallax')
    
    dist_arcsec = (custom_simbad.query_object(tic_id)['PLX_VALUE'][0])/1000
    
    distance = 1/dist_arcsec
    
    return distance*u.parsec

def Lum_from_Tess_Flux(Tess_flux, distance):
    """
    Take Flux values from TESS, and convert into intrinsic Luminosities.
   
    Parameters
    ----------
    Tess_flux : float
                Given in Counts (or electrons) per second, pass the float value
    distance   : units of parsecs
                distance to the object in units of parsecs
               
    Returns
    -------
    Watts
        The luminosity of the object in terms of Watts
       
    """
    T_mag = (-2.5*np.log10(Tess_flux)) + 20.4403
   
    Physical_units= (2416*u.Jansky)*10**(12.5065/(T_mag))
   
    hz=(800*u.nm).to(u.Hz, equivalencies=u.spectral())

    dist_in_meters = distance.to(u.m)
   
    Lum=Physical_units.to(u.W/(u.Hz*u.m**2))*hz*dist_in_meters**2
   
    return Lum

def mass_from_lum(lum):
    '''
    Converts the luminosity of the star to mass.
    
    Parameters
    ----------
    lum: Units of Watts
    
    Returns
    -------
    Solar Mass
        Mass of the star in units of Solar Mass.
    
    '''
    sun_lum = 3.827E26 * u.watt
    
    mass = (lum/(1.4*sun_lum))**(1/3.5)
    
    return mass*u.Msun

def axis_from_masses(mass_star1, mass_star2, period):
    '''
    Finds the semi-major axis of the orbit using the mass of both binary stars and the period of the system.
    
    Parameters
    ----------
    mass_star1: Units of solar mass
    
    mass_star2: Units of solar mass
    
    period: float
    
    Returns
    -------
    Parsecs
        The semi-major axis of the orbit in Parsecs.
        
    Astronomical Units
        The semi_major axis of the orbit in astronomical units.
        
    '''
    
    p_days = period*u.day
    
    constant = ac.G
    
    axis = ((p_days.to(u.second)**2)/(4*(np.pi**2)/(constant*(mass_star1.to(u.kg)+mass_star2.to(u.kg)))))**(1/3)
    
    return  axis.to(u.pc), axis.to(u.au)
