# read in the libraries
import pandas as pd
import numpy as np
import seaborn 
from datetime import datetime
import datetime
from sqlalchemy import create_engine
import sys
sys.path.append(r'C:\work\src')
import settings as s
import matplotlib
import matplotlib.pyplot as plt
import shapely.wkt as shwk
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.io.img_tiles import MapQuestOpenAerial
from cartopy.io.img_tiles import MapQuestOSM
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


#receive credentials via settings file
cs = getattr(s, 'localhost_SurfaceBook')
#cs = getattr(s, 'AWS_localhost')
password = str(cs['password'])
password = password[2:len(password) - 1] #accounts for conversion of encoded password and type cast

# connect to the database
engine = create_engine('postgresql://' + cs['user'] + ':' + str(cs['password'])[2:-1] + '@' + cs['host'] + ':' + cs['port'] + '/' + cs['dbname'])

# get the GRACE data

Gv12 = pd.read_sql('SELECT mascon AS mascon, values_filter1d AS mass, date FROM mascon_solution WHERE version = 12', engine, index_col = 'date')
Gv16 = pd.read_sql('SELECT mascon AS mascon, values_filter1d AS mass, date FROM mascon_solution WHERE version = 16', engine, index_col = 'date')
Gv1601 = pd.read_sql('SELECT mascon AS mascon, values_filter1d AS mass, date FROM mascon_solution WHERE version = 1601', engine, index_col = 'date')

masconAreasv12 = pd.read_sql('SELECT mascon, area_km2 FROM mascon_fit WHERE version = 12', engine)
masconAreasv16 = pd.read_sql('SELECT mascon, area_km2 FROM mascon_fit WHERE version = 16', engine)

Gv12['date'] = Gv12.index
Gv16['date'] = Gv16.index
Gv1601['date'] = Gv1601.index

Gv12_mascons = pd.pivot_table(Gv12, index = 'date', columns='mascon', values='mass')
Gv12_mascons = Gv12_mascons * masconAreasv12.mean()[1]/1e5 

Gv16_mascons = pd.pivot_table(Gv16, index = 'date', columns='mascon', values='mass')
Gv16_mascons = Gv16_mascons * masconAreasv16.mean()[1]/1e5

Gv1601_mascons = pd.pivot_table(Gv1601, index = 'date', columns='mascon', values='mass')
Gv1601_mascons = Gv1601_mascons * masconAreasv16.mean()[1]/1e5

Gv12_sum = Gv12['mass'].groupby(Gv12.index).mean()* masconAreasv12.sum()[1]/1e5 
Gv16_sum = Gv16['mass'].groupby(Gv16.index).mean()* masconAreasv16.sum()[1]/1e5
Gv1601_sum = Gv1601['mass'].groupby(Gv1601.index).mean()* masconAreasv16.sum()[1]/1e5

# pickled data are in meters total per mascon. So no need for area correction. If you divide by 1000 you get km 
# and the sum is over the area of a mascon, so the km^2 is implied. Thus "modeled" comes out in Gt when divided by 1000

modeled = ((pd.read_pickle(r'C:\work\datadrive\SnowModel\CFSR_final.p'))/1000).cumsum() 

# resampling modeled data to the time index of Luthcke's data
modeled_mascons = modeled.reindex(Gv16_sum.index, method = 'ffill')

# pickled data are also using the older mascon numbers so let's renumber here
modeled_mascons.columns = Gv16_mascons.columns

modeled_sum=modeled_mascons.sum(axis = 1)

all = Gv16_mascons.join(Gv1601_mascons, how='inner', rsuffix='_fm').join(modeled_mascons, how='inner', rsuffix='_mod')
all_sum = Gv16_sum.to_frame().join(Gv1601_sum.to_frame(), how='inner', rsuffix='_fm').join(modeled_sum.to_frame(), how='inner')
all_sum['d'] = all_sum['mass_fm']+all_sum[0]
all_sum.columns = ['GRACE','delta','mod','mod+delta']

# now we will plot both a map and the data for each of 3 different regions in Alaska

# ____________________________
# Western Gulf of Alaska : MAP
# ____________________________

masconGeoms = pd.read_sql('SELECT mascon, ST_AsText(geom) AS geom FROM mascon_fit WHERE version = 14 AND mascon IN (1484, 1485, 1480, 1481, 1482, 1483, 1473, 1474, 1475, 1463, 1464, 1465, 1466, 1453, 1454, 1455, 1447)', engine)

masconCentroids = pd.read_sql('SELECT mascon, ST_X(ST_Centroid(geom)) AS longitude, ST_Y(ST_Centroid(geom)) AS latitude FROM mascon_fit WHERE version = 14 AND mascon IN (1484, 1485, 1480, 1481, 1482, 1483, 1473, 1474, 1475, 1463, 1464, 1465, 1466, 1453, 1454, 1455, 1447)', engine)

masconShapely = []
for i in range(len(masconGeoms)):
    masconShapely.append(shwk.loads(masconGeoms.geom[i]))


# Gets the longitude/latitude bounds for a glacier and assigns it to local variables. This is used in the plotting 
# below to define the plot extent
longLatQuery = ("SELECT MIN(g.lowLong) AS lowLong, MAX(g.highLong) AS highLong, MIN(g.lowLat) \
                AS lowLat, MAX(g.highLat) AS highLat FROM (SELECT ST_Xmin(geom) as lowLong, \
                ST_XMax(geom) as highLong, ST_YMin(geom) as lowLat, ST_YMax(geom) as highLat \
                FROM mascon_fit WHERE version = 14 AND mascon IN (1484, 1485, 1480, 1481, 1482, \
                1483, 1473, 1474, 1475, 1463, 1464, 1465, 1466, 1453, 1454, 1455, 1447)) AS g")
                
longLatDataFrame = pd.read_sql_query(longLatQuery, engine)

plt.figure(figsize=(9,5))
geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
extent = longLatDataFrame.values.tolist()[0]
# this is to put a bit of a buffer around the map
correction = 1.0
for i in range(4):
    if i % 2 == 0:
        extent[i] -= correction
    else:
        extent[i] += correction
imagery = MapQuestOSM()
ax = plt.axes(projection=imagery.crs)
ax.set_extent(extent, geodetic)
#ax.coastlines(resolution='100m')
ax.add_image(imagery, 7)
ax.add_feature(cf.ShapelyFeature(masconShapely,geodetic,edgecolor='red', facecolor='none'))

gl = ax.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

for i in range(len(masconCentroids)): 
   longitudeVal = masconCentroids.longitude[i] 
   latitudeVal = masconCentroids.latitude[i]  
   label = str(masconCentroids.mascon[i]) 
   plt.text(longitudeVal, latitudeVal, str(label), horizontalalignment='center',transform=ccrs.Geodetic())

fig = ax.get_figure()
fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\westernMap.png')


# ______________________________
# Western Gulf of Alaska : PLOTS
# ______________________________

matplotlib.style.use('seaborn-poster')
nr = 6
nc = 4
fig, axes = plt.subplots(nrows = nr, ncols = nc, figsize = (16,16))
fig.tight_layout()
subset = {1:1484, 2:1485, 4:1480, 5:1481, 6:1482, 7:1483, 8:1473, 9:1474, 
           10:1475, 12:1463, 13:1464, 14:1465, 15:1466, 16:1453, 
           17:1454, 18:1455, 20:1447}

for i in range(nr*nc):
    plotPos = np.unravel_index(i,(nr,nc))
    try:
        mascon = subset[i]
        GRACE = all[str(mascon)] - all[str(mascon)][0]
        FM = all[str(mascon)+'_fm']-all[str(mascon)+'_fm'][0]
        mod = all[mascon]-all[mascon][0]
        GRACE.plot(ax = axes[plotPos[0],plotPos[1]])
        mod.plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        (mod + FM).plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        axes[plotPos[0],plotPos[1]].set_title(str(mascon))
        axes[plotPos[0],plotPos[1]].get_xaxis().set_visible(False)
        if plotPos == (3,3):
            print('true')
            axes[plotPos[0],plotPos[1]].legend([r'GRACE estimated','SnowModel', 'SnowModel + $\Delta$ forward model'],loc='lower left', bbox_to_anchor=(-0.2, -1.1))
    except:
        axes[plotPos[0],plotPos[1]].axis('off')
 
fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\westernPlot.png')


# ____________________________
# Eastern Gulf of Alaska : MAP
# ____________________________

masconGeoms = pd.read_sql('SELECT mascon, ST_AsText(geom) AS geom FROM mascon_fit WHERE version = 14 AND mascon IN (1476, 1477, 1478, 1467, 1468, 1469, 1470, 1456, 1457, 1458, 1459, 1448, 1449)', engine)


masconShapely = []
for i in range(len(masconGeoms)):
    masconShapely.append(shwk.loads(masconGeoms.geom[i]))

# Gets the longitude/latitude bounds for a glacier and assigns it to local variables. This is used in the plotting 
# below to define the plot extent
longLatQuery = ("SELECT MIN(g.lowLong) AS lowLong, MAX(g.highLong) AS highLong, MIN(g.lowLat) \
                AS lowLat, MAX(g.highLat) AS highLat FROM (SELECT ST_Xmin(geom) as lowLong, \
                ST_XMax(geom) as highLong, ST_YMin(geom) as lowLat, ST_YMax(geom) as highLat \
                FROM mascon_fit WHERE version = 14 AND mascon IN (1476, 1477, 1478, 1467, 1468, \
                1469, 1470, 1456, 1457, 1458, 1459, 1448, 1449)) AS g")
                
longLatDataFrame = pd.read_sql_query(longLatQuery, engine)

plt.figure(figsize=(9,5))
geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
extent = longLatDataFrame.values.tolist()[0]
# this is to put a bit of a buffer around the map
correction = 1.0
for i in range(4):
    if i % 2 == 0:
        extent[i] -= correction
    else:
        extent[i] += correction
imagery = MapQuestOSM()
ax = plt.axes(projection=imagery.crs)
ax.set_extent(extent, geodetic)
#ax.coastlines(resolution='100m')
ax.add_image(imagery, 7)
ax.add_feature(cf.ShapelyFeature(masconShapely,geodetic,edgecolor='red', facecolor='none'))

gl = ax.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

fig = ax.get_figure()
fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\easternMap.png')

# ______________________________
# Eastern Gulf of Alaska : PLOTS
# ______________________________

matplotlib.style.use('seaborn-poster')
nr = 4
nc = 4
fig, axes = plt.subplots(nrows = nr, ncols = nc, figsize = (16,16))
fig.tight_layout()
subset = {0:1476, 1:1477, 2: 1478, 4:1467, 5:1468, 6:1469, 7:1470,
          8:1456, 9:1457, 10:1458, 11:1459, 14:1448,15:1449}

for i in range(nr*nc):
    plotPos = np.unravel_index(i,(nr,nc))
    try:
        mascon = subset[i]
        GRACE = all[str(mascon)] - all[str(mascon)][0]
        FM = all[str(mascon)+'_fm']-all[str(mascon)+'_fm'][0]
        mod = all[mascon]-all[mascon][0]
        GRACE.plot(ax = axes[plotPos[0],plotPos[1]])
        mod.plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        (mod + FM).plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        axes[plotPos[0],plotPos[1]].set_title(str(mascon))
        axes[plotPos[0],plotPos[1]].get_xaxis().set_visible(False)
    except:
        axes[plotPos[0],plotPos[1]].axis('off')

fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\easternPlot.png')


# _________________________________
# Southeastern Gulf of Alaska : MAP
# _________________________________

masconGeoms = pd.read_sql('SELECT mascon, ST_AsText(geom) AS geom FROM mascon_fit WHERE version = 14 AND mascon IN (1460, 1461, 1450, 1451, 1444, 1445, 1439, 1440, 1441, 1436, 1437, 1438, 1434, 1435, 1432)', engine)


masconShapely = []
for i in range(len(masconGeoms)):
    masconShapely.append(shwk.loads(masconGeoms.geom[i]))

# Gets the longitude/latitude bounds for a glacier and assigns it to local variables. This is used in the plotting 
# below to define the plot extent
longLatQuery = ("SELECT MIN(g.lowLong) AS lowLong, MAX(g.highLong) AS highLong, MIN(g.lowLat) \
                AS lowLat, MAX(g.highLat) AS highLat FROM (SELECT ST_Xmin(geom) as lowLong, \
                ST_XMax(geom) as highLong, ST_YMin(geom) as lowLat, ST_YMax(geom) as highLat \
                FROM mascon_fit WHERE version = 14 AND mascon IN (1460, 1461, 1450, 1451, 1444, \
                1445, 1439, 1440, 1441, 1436, 1437, 1438, 1434, 1435, 1432)) AS g")
                
longLatDataFrame = pd.read_sql_query(longLatQuery, engine)

plt.figure(figsize=(9,5))
geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
extent = longLatDataFrame.values.tolist()[0]
# this is to put a bit of a buffer around the map
correction = 1.0
for i in range(4):
    if i % 2 == 0:
        extent[i] -= correction
    else:
        extent[i] += correction
imagery = MapQuestOSM()
ax = plt.axes(projection=imagery.crs)
ax.set_extent(extent, geodetic)
#ax.coastlines(resolution='100m')
ax.add_image(imagery, 7)
ax.add_feature(cf.ShapelyFeature(masconShapely,geodetic,edgecolor='red', facecolor='none'))

gl = ax.gridlines(draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

fig = ax.get_figure()
fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\southeasternMap.png')

# ___________________________________
# Southeastern Gulf of Alaska : PLOTS
# ___________________________________

matplotlib.style.use('seaborn-poster')
nr = 7
nc = 3
fig, axes = plt.subplots(nrows = nr, ncols = nc, figsize = (16,16))
fig.tight_layout()
subset = {0:1460, 1:1461, 3:1450, 4:1451, 6:1444, 7:1445, 9:1439,
          10:1440, 11:1441, 12:1436, 13:1437,14:1438,16:1434,
          17:1435,20:1432}

for i in range(nr*nc):
    plotPos = np.unravel_index(i,(nr,nc))
    try:
        mascon = subset[i]
        GRACE = all[str(mascon)] - all[str(mascon)][0]
        FM = all[str(mascon)+'_fm']-all[str(mascon)+'_fm'][0]
        mod = all[mascon]-all[mascon][0]
        GRACE.plot(ax = axes[plotPos[0],plotPos[1]])
        mod.plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        (mod + FM).plot(ax = axes[plotPos[0],plotPos[1]], legend = False)
        axes[plotPos[0],plotPos[1]].set_title(str(mascon))
        axes[plotPos[0],plotPos[1]].get_xaxis().set_visible(False)
    except:
        axes[plotPos[0],plotPos[1]].axis('off')

fig.savefig(r'C:\Users\Anthony Arendt\Google Drive\forwardModelingGRACE\figures\southeasternPlot.png')