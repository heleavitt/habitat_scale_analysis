
import geopandas as gpd
import pandas as pd
import numpy as np
import gc

from shapely import LineString, MultiPolygon, unary_union

# Your helper functions (assign_class_num, deg2rad, rad2deg, etc.)

def deg2rad(degrees):
    return np.deg2rad(degrees)

def rad2deg(radians):
    return np.rad2deg(radians)

def assign_class_num(row):
    if row['Class_name'] == 'Water':
        return 1
    elif row['Class_name'] == 'Saltmarsh':
        return 2
    elif row['Class_name'] == 'Mangrove':
        return 3
    elif row['Class_name'] == 'Manmade':
        return 4
    elif row['Class_name'] == 'Mud Flat':
        return 5
    else:
        return 0

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
    
    gc.collect()

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
        
# Load your marsh polygon shapefile (equivalent to marshpoly.sf in R)
marshpoly_sf = gpd.read_file("input/google2022.shp")
marshpoly_sf =  marshpoly_sf.to_crs(epsg=32615)

# Create 'land only' dataframe by excluding water
terhab = marshpoly_sf[~marshpoly_sf['Class_name'].isin(['Water', 'Mud Flat'])][['Class_name', 'geometry']]

# Isolate and consolidate the water polygons
water_poly = marshpoly_sf[marshpoly_sf['Class_name'].isin(['Water', 'Mud Flat'])][['Class_name', 'geometry']]
water_poly = unary_union(water_poly.geometry)  # Consolidate into a single multipolygon
water_poly = gpd.GeoDataFrame(geometry=[water_poly], crs=marshpoly_sf.crs)  # Convert back to GeoDataFrame
water_poly['Class_name'] = "Water"  # Re-add the 'Class_name' column

# Bind the water polygon and terrestrial habitats back into one layer
habitat_poly = marshpoly_sf
del marshpoly_sf

# Ensure it's cast as MULTIPOLYGON if needed
habitat_poly['geometry'] = habitat_poly['geometry'].apply(
    lambda geom: geom if geom.geom_type == 'MultiPolygon' else MultiPolygon([geom])
)
habitat_poly = habitat_poly.explode(index_parts=False).reset_index(drop=True)

# Assign class numbers
habitat_poly['Class_num'] = habitat_poly.apply(assign_class_num, axis=1)

# Save the preprocessed habitat data
habitat_poly.to_file("output/satscale/preprocessed_data/habitat_poly.gpkg", driver="GPKG")
water_poly.to_file("output/satscale/preprocessed_data/water_poly.gpkg", driver="GPKG")
terhab.to_file("output/satscale/preprocessed_data/terhab.gpkg", driver="GPKG")

# Load the CSV file into a GeoDataFrame
pffw_sites = pd.read_csv("raw_input/drop_field.csv")
pffw_sites = pffw_sites[pffw_sites['sample_trip'] == 2209]

pffw_sites = gpd.GeoDataFrame(pffw_sites, geometry=gpd.points_from_xy(pffw_sites.longitude, pffw_sites.latitude), crs="EPSG:4326")
pffw_sites = pffw_sites.to_crs(epsg=32615)  # Transform to the desired CRS

# Load the wind data
wind_data = pd.read_csv("input/CO-OPS_8761724_met.csv")

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

fetch_distances = pffw_sites.geometry.apply(
    lambda geom: calculate_fetch(geom, mean_wind_dir_deg, terhab)
)
del wind_data
del terhab
pffw_sites['fetch_distance'] = fetch_distances

# Save the preprocessed site data
pffw_sites.to_file("output/satscale/preprocessed_data/pffw_sites.gpkg", driver="GPKG")
