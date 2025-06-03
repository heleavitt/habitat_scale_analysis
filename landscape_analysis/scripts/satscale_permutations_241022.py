# process_permutation.py

import sys
import geopandas as gpd
import pandas as pd
import numpy as np
import gc
from pathlib import Path

# Import helper functions (assign_class_num, site_edge_filter, calculate_metrics, etc.)

def site_edge_filter(df, man, manmin, m, mmin):
    df['site_type'] = df.apply(lambda row: (
        'mangrove' if row['edge_l.mangrove'] >= man else
        'marsh' if row['edge_l.marsh'] >= m else
        'mixed' if (row['edge_l.marsh'] >= mmin and row['edge_l.mangrove'] >= manmin) else
        'reject'
    ), axis=1)
    return df[df['site_type'] != 'reject']

def calculate_metrics(pffw_sites, habitat_poly, buffer_distance, edge_distance):
    print('Calculating metrics for buffer', buffer_distance, 'and edge', edge_distance)
    pffw_sites_buffered = pffw_sites.copy()
    pffw_sites_buffered['geometry'] = pffw_sites.buffer(buffer_distance)
    habitat_poly = habitat_poly.to_crs(pffw_sites_buffered.crs)
    # Transform to the desired CRS

    intersected = gpd.overlay(pffw_sites_buffered, habitat_poly, how='intersection')
    intersected['area'] = intersected.geometry.area
    intersected['land_water_ratio'] = intersected.apply(
        lambda row: row['area'] if row['Class_name'] not in ['Water', 'Mud Flat'] else 0, axis=1
    )
    gc.collect()
    # Instead, calculate the mud area in `intersected`
    intersected['mud'] = intersected.apply(
        lambda row: row['area'] if row['Class_name'] == 'Mud Flat' else 0, axis=1
    )
    
    group_by_site = intersected.groupby('site_date_key').agg({
        'area': 'sum',
        'land_water_ratio': 'sum',
        'fetch_distance': 'mean',
        'mud': 'sum'
    }).reset_index()

    group_by_site['land_water_ratio'] = group_by_site['land_water_ratio'] / group_by_site['area']
    group_by_site['perc'] = group_by_site['area'] / (np.pi * buffer_distance**2)
    # Adjust water_poly to include 'Mudflat'
    # Combine `Water` and `Mudflat` into a single `MultiPolygon`
    water_and_mud_poly_gdf = intersected[intersected['Class_name'].isin(['Water', 'Mud Flat'])][['site_date_key','Class_name', 'geometry']]

    # Convert back to a GeoDataFrame with a single combined geometry
    water_and_mud_poly_gdf['Class_name'] = "Water_and_Mudflat"
    
    # Use the combined water and mud polygon for buffering and intersection
    water_in_buffer = gpd.clip(water_and_mud_poly_gdf, pffw_sites_buffered)
    water_in_buffer['geometry'] = water_in_buffer.buffer(edge_distance)

    # Update land_poly to exclude `Water_and_Mudflat`
    land_poly = intersected[~intersected['Class_name'].isin(['Water', 'Mud Flat'])]

    # Perform the intersection with the buffered combined water and mud polygon
    edge_land_intersection = gpd.clip(land_poly, water_in_buffer)
    gc.collect()

    edge_land_intersection = gpd.clip(edge_land_intersection, pffw_sites_buffered)
        
    edge_land_intersection['area'] = edge_land_intersection.geometry.area

    hab_100 = edge_land_intersection.groupby(['site_date_key', 'Class_name']).agg({
        'area': 'sum'
    }).reset_index()
    hab_100['e.perc'] = hab_100['area'] / (np.pi * buffer_distance**2)
    del edge_land_intersection

    gc.collect()
    hc100 = hab_100.pivot(index='site_date_key', columns='Class_name', values='e.perc').fillna(0)
    hc100['edge_man'] = hc100.get('Mangrove', 0)
    hc100['edge_mar'] = hc100.get('Saltmarsh', 0)

    hc100['edge_l.mangrove'] = hc100['edge_man'] / (hc100['edge_mar'] + hc100['edge_man'])
    hc100['edge_l.marsh'] = hc100['edge_mar'] / (hc100['edge_mar'] + hc100['edge_man'])
    
    hc100 = hc100.reset_index().merge(group_by_site[['site_date_key', 'land_water_ratio','mud', 'fetch_distance']], on='site_date_key')
    hc100 = site_edge_filter(hc100, man=0.75, manmin=0.25, m=0.75, mmin=0.25)
    
    return hc100

def process_permutation(edge_distance, buffer_distance, output_folder):
    # Load preprocessed data
    habitat_poly = gpd.read_file("output/satscale/preprocessed_data/habitat_poly.gpkg")
    pffw_sites = gpd.read_file("output/satscale/preprocessed_data/pffw_sites.gpkg")

    hc100 = calculate_metrics(pffw_sites, habitat_poly, buffer_distance, edge_distance)
    output_path = Path(output_folder) / f'google2022_edge{int(edge_distance)}_buf{int(buffer_distance)}.csv'

    gc.collect()
    # Check if hc100 is not empty
    if not hc100.empty:
        hc100.to_csv(output_path, index=False)
        print(f"Saved CSV for buffer {buffer_distance} and edge {edge_distance} to {output_path}")
    else:
        print(f"No data to save for buffer {buffer_distance} and edge {edge_distance}")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Process permutations of edge and buffer distances.')
    parser.add_argument('edge_distance', type=float, help='Edge distance')
    parser.add_argument('buffer_distance', type=float, help='Buffer distance')

    args = parser.parse_args()

    output_folder = "output/satscale"
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    process_permutation(args.edge_distance, args.buffer_distance, output_folder)
