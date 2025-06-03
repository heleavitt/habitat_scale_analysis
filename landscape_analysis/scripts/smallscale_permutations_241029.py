import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Point, LineString
from shapely.ops import unary_union
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import glob
import os
# Helper functions
def assign_class_num(row):
    if row['Class_name'] == 'Water':
        return 1
    elif row['Class_name'] == 'Mangrove':
        return 3
    elif row['Class_name'] == 'Port and Sand':
        return 4
    else:
        return 0

def deg2rad(degrees):
    return np.deg2rad(degrees)

def rad2deg(radians):
    return np.rad2deg(radians)

def site_edge_filter(df, man, manmin, m, mmin):
    df['site_type'] = df.apply(lambda row: (
        'mangrove' if row['edge_l.mangrove'] >= man else
        'marsh' if row['edge_l.marsh'] >= m else
        'mixed' if (row['edge_l.marsh'] >= mmin and row['edge_l.mangrove'] >= manmin) else
        'reject'
    ), axis=1)
    return df[df['site_type'] != 'reject']

def calculate_fetch(site_point, wind_direction, land_poly, distance_m=10000):
    """Calculate the fetch distance from a site point in a given wind direction."""
    wind_direction_rad = np.deg2rad(wind_direction)

    # Calculate the end point of the wind line
    far_x = site_point.x + distance_m * np.cos(wind_direction_rad)
    far_y = site_point.y + distance_m * np.sin(wind_direction_rad)

    # Create a LineString representing the wind line
    wind_line = LineString([site_point, (far_x, far_y)])

    # Ensure that intersection returns individual geometries properly
    intersection = land_poly.intersection(wind_line)

    if isinstance(intersection, gpd.GeoSeries):
        intersection = intersection.unary_union  # Merge to a single geometry if multiple

    # Handle intersection cases properly
    if intersection.is_empty:
        return np.nan  # No intersection
    elif hasattr(intersection, 'geoms'):
        # Multiple intersections: Find the nearest one
        distances = [site_point.distance(geom) for geom in intersection.geoms]
        return min(distances) if distances else np.nan
    else:
        # Single intersection
        return site_point.distance(intersection)




# Metric calculation and shapefile processing
def calculate_metrics(pffw_sites, habitat_poly, buffer_distance, edge_distance):
    pffw_sites_buffered = pffw_sites.copy()
    pffw_sites_buffered['geometry'] = pffw_sites.buffer(buffer_distance)
    
    intersected = gpd.overlay(pffw_sites_buffered, habitat_poly, how='intersection')
    intersected['area'] = intersected.geometry.area
    intersected['land_water_ratio'] = intersected.apply(
        lambda row: row['area'] if row['Class_name'] != 'Water' else 0, axis=1
    )
    
    group_by_site = intersected.groupby('site_date_key').agg({
        'area': 'sum',
        'land_water_ratio': 'sum',
        'fetch_distance': 'mean'
    }).reset_index()

    group_by_site['land_water_ratio'] = group_by_site['land_water_ratio'] / group_by_site['area']
    group_by_site['perc'] = group_by_site['area'] / (np.pi * buffer_distance**2)
    
    water_poly = habitat_poly[habitat_poly['Class_name'] == 'Water']
    water_in_buffer = gpd.clip(water_poly, pffw_sites_buffered)
    water_in_buffer['geometry'] = water_in_buffer.buffer(edge_distance)
        
    land_poly = habitat_poly[habitat_poly['Class_name'] != 'Water']
    land_in_buffer = gpd.sjoin(land_poly, pffw_sites_buffered, predicate='intersects')
    edge_land_intersection = gpd.clip(land_in_buffer, water_in_buffer)
    edge_land_intersection = gpd.clip(edge_land_intersection, pffw_sites_buffered) # ensures that the edge buffer doesn't come outside of the radius buffer

    edge_land_intersection['area'] = edge_land_intersection.geometry.area

            
    hab_100 = edge_land_intersection.groupby(['site_date_key', 'Class_name']).agg({
        'area': 'sum'
    }).reset_index()
    
    hab_100['e.perc'] = hab_100['area'] / (np.pi * buffer_distance**2)
    
    hc100 = hab_100.pivot(index='site_date_key', columns='Class_name', values='e.perc').fillna(0)
    hc100['edge_man'] = hc100.get('Mangrove', 0)
    hc100['edge_mar'] = hc100.get('Spartina', 0)
    hc100['edge_l.mangrove'] = hc100['edge_man'] / (hc100['edge_mar'] + hc100['edge_man'])
    hc100['edge_l.marsh'] = hc100['edge_mar'] / (hc100['edge_mar'] + hc100['edge_man'])
    
    hc100 = hc100.reset_index().merge(group_by_site[['site_date_key', 'land_water_ratio', 'fetch_distance']], on='site_date_key')
    hc100 = site_edge_filter(hc100, man=0.75, manmin=0.25, m=0.75, mmin=0.25)
    
    return hc100

# Process shapefile function
def process_shapefile(shapefile_path, pffw_sites, buffer_distance, edge_distance, output_folder):
    habitat_poly = gpd.read_file(shapefile_path)
    habitat_poly = habitat_poly[habitat_poly['Class_name'] != 'Border']
    habitat_poly['geometry'] = habitat_poly.geometry.buffer(0)


        ## Ensure that only sites within a safe distance of the habitat polygons are processed
    habitat_poly = habitat_poly[~habitat_poly.is_empty & habitat_poly.is_valid]

    # Perform the union operation
    habitat_union = habitat_poly.unary_union
    # Convert habitat_union back to a GeoDataFrame for a spatial join
    habitat_union_gdf = gpd.GeoDataFrame(geometry=[habitat_union], crs=habitat_poly.crs)

    # Use a spatial join to find intersecting points
    pffw_sites_within = gpd.sjoin(pffw_sites, habitat_union_gdf, predicate='intersects')
    

    # Filter the original pffw_sites to match the indices of sites within the safe distance
    pffw_sites_within = pffw_sites[pffw_sites.index.isin(pffw_sites_within.index)]
    
    # Proceed only if there are sites within the safe distance
    if not pffw_sites_within.empty:
        hc100 = calculate_metrics(pffw_sites_within, habitat_poly, buffer_distance, edge_distance)
        output_path = Path(output_folder) / f'{Path(shapefile_path).stem}_edge{edge_distance}_buf{buffer_distance}.csv'
        hc100.to_csv(output_path, index=False)
        print(f"Saved CSV for {shapefile_path} with buffer {buffer_distance} and edge {edge_distance} to {output_path}")

import time

def merge_csvs_for_iteration(buffer_distance, edge_distance, shapefile_paths, output_folder):
    combined_df = []

    for shapefile_path in shapefile_paths:
        csv_path = Path(output_folder) / f'{Path(shapefile_path).stem}_edge{edge_distance}_buf{buffer_distance}.csv'

        # Wait for the file to exist with a timeout
        wait_time = 0
        max_wait_time = 60  # Maximum time to wait (e.g., 120 seconds)
        while not csv_path.exists() and wait_time < max_wait_time:
            time.sleep(5)  # Wait for 5 seconds before checking again
            wait_time += 5

        if not csv_path.exists():
            print(f"File {csv_path} not found after waiting.")
            continue

        df = pd.read_csv(csv_path)
        combined_df.append(df)

    if combined_df:
        merged_df = pd.concat(combined_df, ignore_index=True)
        final_output_path = Path(output_folder) / f'combined_edge{edge_distance}_buf{buffer_distance}.csv'
        merged_df.to_csv(final_output_path, index=False)
        print(f"Saved merged CSV for buffer {buffer_distance} and edge {edge_distance} to {final_output_path}")
    else:
        print(f"No files found for merging with buffer {buffer_distance} and edge {edge_distance}.")

# Main function to process all shapefiles
# Main function to process all shapefiles
def process_multiple_shapefiles(pffw_sites, shapefile_paths, buffer_distances, shapefile_merge_paths, edge_distances, output_folder):
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    with ProcessPoolExecutor(max_workers=50) as executor:
        futures = []
        for shapefile_path in shapefile_paths:
            for buffer_distance in buffer_distances:
                for edge_distance in edge_distances:
                    futures.append(
                        executor.submit(
                            process_shapefile, shapefile_path, pffw_sites, buffer_distance, edge_distance, output_folder
                        )
                    )

        # Process each future as it completes
        for future in as_completed(futures):
            try:
                future.result()  # This will raise an exception if the processing failed
                print("Processing completed successfully.")
            except Exception as e:
                print(f"An error occurred: {e}")

    # Merge results after all parallel tasks are complete
    for buffer_distance in buffer_distances:
        for edge_distance in edge_distances:
            merge_csvs_for_iteration(buffer_distance, edge_distance, shapefile_merge_paths, output_folder)

sample_trip = 2209
# Load your marsh polygon shapefile (equivalent to marshpoly.sf in R)
marshpoly_sf = gpd.read_file("input/google2022.shp")
marshpoly_sf = marshpoly_sf.to_crs(epsg = 32615)
# Load the wind data
wind_data = pd.read_csv("input/CO-OPS_8761724_met.csv")
# Load the CSV file into a GeoDataFrame
pffw_sites = pd.read_csv("input/drop_field.csv")
pffw_sites = pffw_sites[pffw_sites['sample_trip'] == sample_trip]
# Create 'land only' dataframe by excluding water
terhab = marshpoly_sf[~marshpoly_sf['Class_name'].isin(['Water', 'Mud Flat'])][['Class_name', 'geometry']]



# Isolate and consolidate the water polygons
water_poly = marshpoly_sf[marshpoly_sf['Class_name'] == 'Water']
water_poly = unary_union(water_poly.geometry)  # Consolidate into a single multipolygon
water_poly = gpd.GeoDataFrame(geometry=[water_poly], crs=marshpoly_sf.crs)  # Convert back to GeoDataFrame
water_poly['Class_name'] = "Water"  # Re-add the 'Class_name' column

# Bind the water polygon and terrestrial habitats back into one layer
habitat_poly = pd.concat([water_poly, terhab], ignore_index=True)

# Ensure it's cast as MULTIPOLYGON if needed
habitat_poly['geometry'] = habitat_poly['geometry'].apply(
    lambda geom: geom if geom.geom_type == 'MultiPolygon' else MultiPolygon([geom])
)
habitat_poly = habitat_poly.explode().reset_index(drop=True)  # This ensures proper MULTIPOLYGON casting

# Assign class numbers
habitat_poly['Class_num'] = habitat_poly.apply(assign_class_num, axis=1)




pffw_sites = gpd.GeoDataFrame(pffw_sites, geometry=gpd.points_from_xy(pffw_sites.longitude, pffw_sites.latitude), crs="EPSG:4326")
pffw_sites = pffw_sites.to_crs(epsg=32615)  # Transform to the desired CRS


# Convert wind direction from degrees to radians
wind_data['Wind Dir (rad)'] = wind_data['Wind Dir (deg)'].apply(deg2rad)

# Calculate the vector components of wind direction
wind_data['u'] = wind_data['Wind Speed (kn)'] * np.cos(wind_data['Wind Dir (rad)'])
wind_data['v'] = wind_data['Wind Speed (kn)'] * np.sin(wind_data['Wind Dir (rad)'])

# Calculate the mean of the vector components
mean_u = wind_data['u'].mean(skipna=True)
mean_v = wind_data['v'].mean(skipna=True)

# Calculate the prevailing wind direction in radians
mean_wind_dir_rad = np.arctan2(mean_v, mean_u)

# Convert the prevailing wind direction back to degrees and adjust
mean_wind_dir_deg = rad2deg(mean_wind_dir_rad) + 180
if mean_wind_dir_deg < 0:
    mean_wind_dir_deg += 360

# Apply the fetch calculation to each site point
fetch_distances = pffw_sites.geometry.apply(
    lambda geom: calculate_fetch(geom, mean_wind_dir_deg, terhab)
)

# Store the results in the GeoDataFrame
pffw_sites['fetch_distance'] = fetch_distances
# Changed files paths a bit to fix the interupption without redoing everything 
buffer_distances = [20, 30, 50, 70, 100, 120, 150]
edge_distances = [1, 3, 5]
shapefile_paths = glob.glob(os.path.join('input', '**', '*.shp'), recursive=True)


shapefile_merge_paths = shapefile_paths

# Process all shapefiles and save the output
process_multiple_shapefiles(
    pffw_sites=pffw_sites,
    shapefile_paths=shapefile_paths,
    buffer_distances=buffer_distances,
    shapefile_merge_paths=shapefile_merge_paths,
    edge_distances=edge_distances,
    output_folder="output/smallscale")
