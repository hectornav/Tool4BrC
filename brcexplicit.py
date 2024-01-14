
import pandas as pd

from modules import conc2abs as c2a
from modules import constants
from modules import data_retrieval as dr
from modules import utils as ut

#model = dr.get_explct4brcmass('Payerne', remove_negatives=True)
#observation = dr.get_obsabs370('Payerne', remove_negatives=True)
ri_brc_strng = 0.1
ri_brc_blchd = 0.01

#abs_brc = c2a.get_absorption4brc(model, WAVELENGTH=constants.WAVELENGTH, REL_HUM=constants.RELATIVE_HUMIDITY, ri_brc_strng=ri_brc_strng, ri_brc_blchd=ri_brc_blchd, SA=False)
# Load the observed data from the CSV file.
observed_data_clean = pd.read_csv('../absorption/NInventory/obs/absorption/absorption_brc370.csv')

# Get a list of unique stations from the observed data.
unique_stations = observed_data_clean[['station_name', 'category_station']].drop_duplicates()

# Get a list of station names from the unique stations dataframe.
list_of_stations = unique_stations['station_name'].unique()

for stn in list_of_stations:
    ut.plot_line_seasonal_abs_brc(stn, ri_brc_strng, ri_brc_blchd)

