
import json 
import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import table

dict_of_names = {}
with open('modules/conv_stat_names.json') as file:
    dict_of_names = json.load(file)

#print(dict_of_names['stations']['Ispra'])

values_path = 'num_val/'

dfs =[]
for station in os.listdir(values_path):
    df = pd.read_csv(values_path+station)
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)
#rename station names as in the json file
df['station'] = df['station'].replace(dict_of_names['stations'])
#change column name from num_val to num_pts
df = df.rename(columns={'num_val':'pts'})
df = df.rename(columns={'station':'stn'})
#sort by station name
df = df.sort_values(by=['stn'])
#ocultar index
sns.set_style("whitegrid")
fig, ax = plt.subplots(figsize=(8, 3)) 
# Ocultar los ejes
ax.xaxis.set_visible(False)  
ax.yaxis.set_visible(False)  
ax.set_frame_on(False)  

# Crear la tabla y eliminar la celda de color por defecto
tabla = table(ax, df.T, loc='center', cellLoc='center', rowLoc='center')

# Personalizar cada celda de la tabla
for key, cell in tabla.get_celld().items():
    cell.set_linewidth(0.5)
    cell.set_edgecolor('black')
    
    # Ocultar índice
    if key[0] == 0 or key[1] == -1:
        cell.set_visible(False)
        
    if key[1] == 0:  # Encabezado de la tabla
        cell.set_fontsize(12)
        #cell.set_text_props(weight='bold')
        #cell.set_facecolor('lightgrey')
    else:  # Cuerpo de la tabla
        cell.set_fontsize(12)
        cell.set_facecolor('white')

#fit the image 
fig.tight_layout()
plt.savefig('num_points.png', bbox_inches='tight', dpi=300)

# Cerrar el gráfico para liberar memoria
plt.close()