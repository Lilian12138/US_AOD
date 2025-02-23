import pandas as pd
import glob
import geopandas as gpd
import numpy as np
import math
import matplotlib.pyplot as plt
from pykrige.ok import OrdinaryKriging
import xarray as xr
import salem
import warnings
warnings.filterwarnings("ignore")

def boundary(df):
    js_box = df.geometry.total_bounds
    # create grid points the spatial space is equal
    N_lon = abs(math.ceil((int(js_box[0]-1) - int(js_box[2]+1))/0.1))
    N_lat = abs(math.ceil((int(js_box[1]-1) - int(js_box[3])+1)/0.1))
    grid_lon = np.linspace(js_box[0],js_box[2],N_lon, endpoint=False)
    grid_lat = np.linspace(js_box[1],js_box[3],N_lat, endpoint=False)
    print(f'the spatial box is: \n{js_box}')
    print(f'the spatial space is: \n{round(grid_lon[1]-grid_lon[0],2)} x {round(grid_lat[1]-grid_lat[0],2)}')
    return grid_lon, grid_lat

if __name__ == '__main__':

    #### 1. open the air quality csv files and merge to panel table. 
    AOD_list = glob.glob(r'AOD_interpolation_processing\UScounty_AOD_*.csv')
    # merge modis AOD
    AOD_df_list = []
    for path in AOD_list:
        AOD_df = pd.read_csv(path)
        AOD_df.drop(['system:index','Shape_Area', 'Shape_Leng', '.geo'], axis = 1, inplace=True)
        AOD_df = AOD_df.melt(id_vars=['GEOID'], var_name='Date', value_name='AOD')
        AOD_df['Date'] = pd.to_datetime(AOD_df['Date'].str.split('_', expand=True)[1]).dt.strftime('%Y-%m-%d')
        AOD_df_list.append(AOD_df)
    AOD_df = pd.concat(AOD_df_list)
    AOD_df = AOD_df.groupby(by = ['GEOID', 'Date']).mean().reset_index()
    # add the entir time df to check the empty time and GEOID
    AOD_df['Date'] =  AOD_df['Date'].astype('datetime64[ns]')
    GEOID_lst = AOD_df['GEOID'].drop_duplicates().to_list()
    df = pd.DataFrame(GEOID_lst,columns=['GEOID'])
    df = pd.DataFrame(pd.date_range(start=AOD_df['Date'].min(), end=AOD_df['Date'].max()), columns=['Date']).merge(df, how ='cross')
    AOD_df = df.merge(AOD_df, on = ['GEOID','Date'], how ='left')

    #### 2. time series interpolation using nearest valid time by conty id
    time_interp_AOD_df_lst = []
    for geoitem in GEOID_lst:
        df = AOD_df.loc[AOD_df['GEOID'] == geoitem].sort_values(by='Date', ascending=True)
        df['AOD_spline'] = df['AOD'].interpolate(method='spline', order=3, limit=5)
        time_interp_AOD_df_lst.append(df)
    time_interp_AOD_df = pd.concat(time_interp_AOD_df_lst)
    # time_interp_AOD_df.to_csv(r'AOD_interpolation_processing\scratch\time_interp_AOD_df.csv', index=False)

    #### 3. spatial interpolation using kriging
    # open selected us county
    uscounty_gdf = gpd.read_file(r'AOD_interpolation_processing\uscounty_select.shp')
    time_interp_AOD_df = time_interp_AOD_df.merge(uscounty_gdf, left_on='GEOID', right_on='geoid_j', how = 'left')
    time_interp_AOD_df.drop(columns=['GEOID_y','geoid_j'],inplace=True)
    time_interp_AOD_df.rename(columns={'GEOID_x':'GEOID'},inplace=True)
    grid_lon, grid_lat = boundary(uscounty_gdf)
    # get the date list
    Date_list = time_interp_AOD_df['Date'].drop_duplicates().to_list()
    
    for date in Date_list:
        df_daily = time_interp_AOD_df.loc[time_interp_AOD_df["Date"]==date].dropna(axis=0,how='any')
        # Add check for empty DataFrame
        if len(df_daily) == 0:
            print(f"No valid data for date: {date}")
            continue
        know_lon = df_daily["X_lon"].values
        know_lat = df_daily["Y_lat"].values
        know_z = df_daily["AOD"].values
        OK = OrdinaryKriging(know_lon,know_lat, know_z,variogram_model='gaussian',nlags=6)
        z,ss = OK.execute('grid', grid_lon, grid_lat)
        proj= 'EPSG:4326'

        # create dataarray from np.array
        da_z = xr.DataArray(
            data=z,
            dims=["lat","lon"],
            coords=dict(
                lon = (['lon'], grid_lon),
                lat = (['lat'], grid_lat))
        )
        pd_na = uscounty_gdf
        pd_na['Date'] = date
        for idx, row in pd_na.iterrows():
            data_roi = da_z.salem.roi(geometry=row['geometry'],crs=proj)
            data_z_roi = data_roi.mean(skipna=True,keep_attrs=True).values
            pd_na.loc[idx,'Filled_AOD'] = data_z_roi
        print('The number of NA is', len(pd_na.loc[pd_na['Filled_AOD'].isna()]), 'for ', date)
        # saving in scrach folder to decrease the merroy storage using.
        pd_na.loc[:,['GEOID','Filled_AOD','Date']].to_csv(f'AOD_interpolation_processing\scratch\spatial_interp_AOD_df_{date}.csv', index=False)
    #### 4. merge
    spatial_interp_AOD_df_lst = []
    for file in glob.glob(f'AOD_interpolation_processing\scratch\spatial_interp_AOD_df_*.csv'):
        df = pd.read_csv(file)
        spatial_interp_AOD_df_lst.append(df)
    spatial_interp_AOD_df=pd.concat(spatial_interp_AOD_df_lst)
    AOD_df = time_interp_AOD_df.merge(spatial_interp_AOD_df, on = ['GEOID', 'Date'], how ='left')
    AOD_df['Merge_AOD'] = AOD_df['AOD_spline']
    AOD_df.loc[AOD_df['Merge_AOD'].isna(), 'Merge_AOD'] = AOD_df['Filled_AOD']
    AOD_df.to_csv('AOD_interpolation_processing\AOD_merge.csv', index=False)
    print(AOD_df)
