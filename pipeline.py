import numpy as np
from astropy.table import Table, join, MaskedColumn, vstack, Column
from astropy.io import ascii
from numpy import inf
import matplotlib.pyplot as plt
import lightkurve as lk
from astroquery.simbad import Simbad
from astropy import units as u


#Danbi Section
def download_lightcurve(name):
    """Downloads the lightcurve for a given object. Returns the first table in the search result.
    Uses lightkurve package.
    
    Parameters
    ----------
    name : string
           id of the target object
    
    Returns
    -------
    The first (index 0) SearchResult (lightkurve class) object containing data table.
    """
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

def find_period(name, data):
    """Returns the period of the given object in days from the catalog.
    
    Parameters
    ----------
    name : string
            id of the target object
    
    Returns
    -------
    period : float
             Orbital period (in days)"""
    period = data[np.where(data['Obj-ID'] == name)]['BLS-Period']
    return period

def find_fluxes(search, period, range):
    """Takes in a SearchResult object, the corresponding period, and desired fractional width, and returns a pair of fluxes for each object.
    Returns fluxes of the dimmer star, the brighter star, and total flux, in that order. Units will be preserved from the input objects.

    Parameters
    ----------
    search : SearchResult
             contains data for the target object
    period : float
             Orbital period of the target system
    range : float
            fractional value for how narrow one wishes the range of datapoints to include. Must be less than 1.

    Returns
    -------
    flux_A : float32
             flux of the dimmer star in units of electrons per seconds.
    flux_B : float32
             flux of the brighter star in units of electrons per seconds.
    flux_tot : float32
            total flux of the system, in units of electrons per seconds.
    """
    half_folded = search.time.value % (period/2.0)
    folded = search.time.value % period
    # find fractional range, gather datapoints, find median
    width = period * range # some value around 0.05?
    left = (period/4.0) - (width/2.0) # divisor 4 bc 1/2 * 1/2
    right = (period/4.0) + (width/2.0)
    flux_tot = half_folded[left:right]
    # update left:right for full period
    left = (period/2.0) - (width/2.0)
    right = (period/2.0) + (width/2.0)
    flux_B = folded[left:right]
    flux_A = flux_tot - flux_B
    return flux_A, flux_B, flux_tot


#Amelia Section
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


class Properties_of_EBs:
    def __init__(self, tic_id, data_table):
        
        #required inputs
        self.tic_id = tic_id
        self.data_table = data_table

        #parameters that a user could change
        self.range = 0.05 #The fraction of the period to search for total flux

        self.showPlots = True

    def download_lightcurve(self, tic_id):
        """Downloads the lightcurve for a given object. Returns the first table in the search result.
        Uses lightkurve package.
        
        Parameters
        ----------
        tic_id : string
            id of the target object
        
        Returns
        -------
        The first (index 0) SearchResult (lightkurve class) object containing data table from SPOC.
        """
        search = lk.search_lightcurve(tic_id, author='SPOC')

        self.lightcurve_table = search[0].download()
        # search.table["dataURL"] = search.table["dataURI"]
        return self.lightcurve_table
    
    def plot_lc(self, lightcurve_table):
        fig=plt.figure()
        plt.title(self.lightcurve_table.label)
        plt.scatter(self.lightcurve_table['time'].value, self.lightcurve_table['flux'].value, color='b')
        plt.xlabel("Time [J Day]", fontsize=16)
        plt.ylabel("Flux [e/s]", fontsize=16)
        plt.show()
        plt.close(fig)
    
    def find_period(self, tic_id, data_table):
        """Returns the period of the given object in days from the catalog.
        
        Parameters
        ----------
        tic_id : string
                id of the target object
        data_table : string
                table of systems we are searching for
        
        Returns
        -------
        period : float
                Orbital period (in days)"""
        
        self.period = self.data_table[np.where(self.data_table['Obj-ID'] == self.tic_id)]['BLS-Period']
        return self.period

    def find_fluxes(self, lightcurve_table, period):
        """Takes in a SearchResult object, the corresponding period, and desired fractional width, and returns a pair of fluxes for each object.
        Returns fluxes of the dimmer star, the brighter star, and total flux, in that order. Units will be preserved from the input objects.

        Parameters
        ----------
        lightcurve_table : SearchResult
                contains data for the target object
        period : float
                Orbital period of the target system
        range : float
                fractional value for how narrow one wishes the range of datapoints to include. Must be less than 1.

        Returns
        -------
        flux_A : float32
                flux of the dimmer star in units of electrons per seconds.
        flux_B : float32
                flux of the brighter star in units of electrons per seconds.
        flux_tot : float32
                total flux of the system, in units of electrons per seconds.
        """
        half_folded = self.lightcurve_table.time.value % (self.period/2.0)
        folded = self.lightcurve_table.time.value % self.period

        flux=self.lightcurve_table.flux.value
       
        # find fractional range, gather datapoints, find median
        width = self.period * self.range # some value around 0.05?
        left = (self.period/4.0) - (width/2.0) # divisor 4 bc 1/2 * 1/2
        right = (self.period/4.0) + (width/2.0)
        self.flux_tot = np.median(flux[np.where((half_folded>left) & (half_folded < right))])
        
        # update left:right for full period
        left = (self.period/2.0) - (width/2.0)
        right = (self.period/2.0) + (width/2.0)
        self.flux_B = np.median(flux[np.where((folded>left) & (folded < right))])
        self.flux_A = self.flux_tot - self.flux_B
        return self.flux_A, self.flux_B, self.flux_tot
    
    def distance_from_simbad(self, tic_id):
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

        self.distance = distance*u.parsec
        
        return self.distance


    def Lum_from_Tess_Flux(self, Tess_flux, distance):
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

        dist_in_meters = self.distance.to(u.m)
    
        Lum=Physical_units.to(u.W/(u.Hz*u.m**2))*hz*dist_in_meters**2
    
        return Lum

    def mass_from_lum(self, lum):
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

    def axis_from_masses(self, mass_star1, mass_star2, period):
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
        
        p_days = self.period*u.day
        
        constant = ac.G
        
        axis = ((p_days.to(u.second)**2)/(4*(np.pi**2)/(constant*(mass_star1.to(u.kg)+mass_star2.to(u.kg)))))**(1/3)
        
        return  axis.to(u.pc), axis.to(u.au)
    
    def runAll(self):
        self.download_lightcurve(self.tic_id)
        self.find_period(self.tic_id, self.data_table)
        self.distance_from_simbad(self.tic_id)
        self.find_fluxes(self.lightcurve_table, self.period)
        
        self.lum_A = Lum_from_Tess_Flux(self.flux_A, self.distance)
        self.lum_B = Lum_from_Tess_Flux(self.flux_B, self.distance)

        self.M_A = mass_from_lum(self.lum_A) 
        self.M_B = mass_from_lum(self.lum_B)  

        self.pc_separation, self.au_separation = axis_from_masses(self.M_A, self.M_B, self.period)

        if self.showPlots:
            self.plot_lc(self.lightcurve_table)