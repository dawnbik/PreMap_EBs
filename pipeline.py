Danbi Section
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
    








Amelia Section
