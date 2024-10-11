
import time
import math 
import pandas as pd
import folium
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import wkt
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Point, Polygon
import geopandas as gpd
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap  # Ensure this line is included


def read_and_project_transport_data(filepath, line_numbers=None, crs='EPSG:32633', transport_type=None):
    """
    Reads public transport data from a shapefile, optionally filters by line numbers and transport type,
    and reprojects it to a specified CRS.

    Parameters:
        filepath (str): Path to the shapefile.
        crs (str): Coordinate Reference System to project the data to. Default is 'EPSG:32633'.
        line_numbers (int, list, tuple, or None): Line numbers to filter the data by. Default is None (no filtering).
        transport_type (str, list, or None): Transport type(s) to filter the data by. Default is None (no filtering).

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the projected transport data.
    """
    # Read shapefile
    gdf = gpd.read_file(filepath)
    
    # Apply filtering based on line_numbers if provided
    if line_numbers is not None:
        if isinstance(line_numbers, (list, tuple)):
            gdf = gdf[gdf["Lijn_Numbe"].isin(line_numbers)]
        elif isinstance(line_numbers, str): #int):
            gdf = gdf[gdf["Lijn_Numbe"] == line_numbers]
        else:
            raise ValueError("line_numbers must be an int, list, or tuple")
    
    #Apply filtering based on transport_type if provided
    if (transport_type is not None) and (line_numbers is None):
    #if transport_type is not None:
        if isinstance(transport_type, list):
            gdf = gdf[gdf["type"].isin(transport_type)]
        elif isinstance(transport_type, str):
            gdf = gdf[gdf["type"] == transport_type]
        else:
            raise ValueError("transport_type must be a str or list of str")

    # Reproject to the specified CRS
    gdf_projected = gdf.to_crs(crs)
    
    return gdf_projected



def calculate_buffer(gdf, buffer_distance, crs='EPSG:3857'):
    """
    Calculates a buffer around each geometry in the GeoDataFrame.

    Parameters:
        gdf (gpd.GeoDataFrame): GeoDataFrame with the geometries to buffer.
        buffer_distance (float): Buffer distance in meters. Default is 100 meters.
        crs (str): Coordinate Reference System to use for buffering. Default is 'EPSG:3857' (meters).

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the buffered geometries.
    """
    # Reproject to a CRS that uses meters for distance calculations
    gdf_meters = gdf.to_crs(crs)
    
    # Create a buffer around each geometry
    gdf_meters['geometry'] = gdf_meters['geometry'].buffer(buffer_distance)
    
    return gdf_meters


def load_and_check_cbs_data(filepath, lines):
    """
    Loads CBS population data from a shapefile, performs data checks, and ensures CRS consistency with transport data.

    Parameters:
        filepath (str): Path to the CBS shapefile.
        trams_gdf (gpd.GeoDataFrame): GeoDataFrame containing transport data to match CRS with.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the loaded CBS data with matching CRS to the transport data.
        pd.Series: Series containing counts of NaN values for each column.
        pd.Series: Series containing counts of zero values for each column.
    """
    # Load CBS data
    cbs_gdf = gpd.read_file(filepath)
    
    # Ensure CRS consistency with transport data
    if cbs_gdf.crs != lines.crs:
        cbs_gdf = cbs_gdf.to_crs(lines.crs)
    
    # Count NaN values in each column
    nan_counts = cbs_gdf.isna().sum()
    
    # Count 0 values in each column
    zero_counts = (cbs_gdf == 0).sum()
     
    
    return cbs_gdf


def spatial_join_and_remove_duplicates(cbs_gdf, buffered_lines, how='inner', predicate='intersects'):
    """
    Performs a spatial join to find intersecting geometries between CBS and transport data, and removes duplicates.

    Parameters:
        cbs_gdf (gpd.GeoDataFrame): GeoDataFrame containing CBS data.
        trams_gdf (gpd.GeoDataFrame): GeoDataFrame containing transport data.
        how (str): Type of join. Default is 'inner'.
        predicate (str): Spatial predicate to use for the join. Default is 'intersects'.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the result of the spatial join with duplicates removed.
    """

    # Ensure both GeoDataFrames have the same CRS
    if cbs_gdf.crs != buffered_lines.crs:
        buffered_lines = buffered_lines.to_crs(cbs_gdf.crs)
    
    # Perform spatial join
    joined_gdf = gpd.sjoin(cbs_gdf, buffered_lines, how=how, predicate=predicate)
    
    # Remove duplicates
    
    #joined_gdf = joined_gdf.drop_duplicates()

    #joined_gdf = joined_gdf[joined_gdf.duplicated(subset='crs28992', keep=False)]

    joined_gdf = joined_gdf.drop_duplicates(subset='crs28992', keep='first')  # Keep the first occurrence and remove others


    
    return joined_gdf


def generate_summary_statistics(cbs_gdf):
    """
    Calculates summary statistics from the CBS GeoDataFrame, excluding specified columns,
    and formats the results into a summary DataFrame.

    Parameters:
        cbs_gdf (gpd.GeoDataFrame): GeoDataFrame containing CBS data.

    Returns:
        pd.DataFrame: A DataFrame containing the summary statistics.
    """
    # Calculate the sum for numeric columns excluding 'G_woz_woni' and 'geometry'
    sum_values = cbs_gdf.drop(columns=['geometry', 'G_woz_woni']).sum()
    
    # Calculate the average for 'G_woz_woni'
    average_woz_woni = cbs_gdf['G_woz_woni'].mean()
    
    # Round the sums and average to zero decimal places
    rounded_sum_values = sum_values.round(0)
    rounded_average_woz_woni = round(average_woz_woni, 0)
    
    # Create the summary row DataFrame
    summary_row_ams = pd.DataFrame({
        'Area': ['Amsterdam'],
        'A_inhab': [rounded_sum_values.get('A_inhab', 0)],
        'A_0_15': [rounded_sum_values.get('A_0_15', 0)],
        'A_15_25': [rounded_sum_values.get('A_15_25', 0)],
        'A_25_45': [rounded_sum_values.get('A_25_45', 0)],
        'A_45_65': [rounded_sum_values.get('A_45_65', 0)],
        'A_65+': [rounded_sum_values.get('A_65+', 0)],
        'A_woning': [rounded_sum_values.get('A_woning', 0)],
        'G_woz_woni': [rounded_average_woz_woni],
        'A_nederlan': [rounded_sum_values.get('A_nederlan', 0)],
        'A_west_mig': [rounded_sum_values.get('A_west_mig', 0)],
        'A_n_west_m': [rounded_sum_values.get('A_n_west_m', 0)],
    })

    return summary_row_ams


def generate_summary_statistics_for_sensed_data(joined_gdf):
    """
    Generates summary statistics for a GeoDataFrame by calculating sums for numeric columns,
    calculating the average for a specific column, and formatting the results into a summary DataFrame.
    
    The function excludes the columns 'crs28992', 'G_woz_woni', and 'geometry' from summation,
    and calculates the average for the 'G_woz_woni' column.
    
    Parameters:
        gdf (gpd.GeoDataFrame): GeoDataFrame containing the data.

    Returns:
        pd.DataFrame: A DataFrame containing the summary statistics.
    """

    # Calculate the sum for numeric columns excluding specified columns
    sum_values = joined_gdf.drop(columns=['geometry', 'G_woz_woni']).sum()
    
    # Calculate the average for the specified column
    average_value = joined_gdf['G_woz_woni'].mean()
    
    # Round the sums and average to zero decimal places
    rounded_sum_values = sum_values.round(0)
    rounded_average_value = round(average_value, 0)
    
    # Create the summary row DataFrame
    summary_row = pd.DataFrame({
        'Area': ['Sensed Area'],  # Customize this label if needed
        'A_inhab': [rounded_sum_values.get('A_inhab', 0)],
        'A_0_15': [rounded_sum_values.get('A_0_15', 0)],
        'A_15_25': [rounded_sum_values.get('A_15_25', 0)],
        'A_25_45': [rounded_sum_values.get('A_25_45', 0)],
        'A_45_65': [rounded_sum_values.get('A_45_65', 0)],
        'A_65+': [rounded_sum_values.get('A_65+', 0)],
        'A_woning': [rounded_sum_values.get('A_woning', 0)],
        'G_woz_woni': [rounded_average_value],
        'A_nederlan': [rounded_sum_values.get('A_nederlan', 0)],
        'A_west_mig': [rounded_sum_values.get('A_west_mig', 0)],
        'A_n_west_m': [rounded_sum_values.get('A_n_west_m', 0)],
    })

    return summary_row

def concatenate_dataframes(dfs, ignore_index=True):
    """
    Concatenates a list of DataFrames into a single DataFrame.

    Parameters:
        dfs (list of pd.DataFrame): List of DataFrames to concatenate.
        ignore_index (bool): Whether to ignore the index and reset it in the resulting DataFrame.

    Returns:
        pd.DataFrame: A single DataFrame resulting from the concatenation of the input DataFrames.
    """
    # Concatenate the DataFrames
    merged_df = pd.concat(dfs, ignore_index=ignore_index)
    
    return merged_df


def calculate_and_compare_sums(cbs_gdf, sensed_gdf):
    """
    Calculates and compares the sums of columns for two GeoDataFrames, and computes the percentage
    of the sensed values relative to the city's total values.

    Parameters:
        cbs_gdf (gpd.GeoDataFrame): GeoDataFrame containing the full city's data.
        sensed_gdf (gpd.GeoDataFrame): GeoDataFrame containing the sensed data.

    Returns:
        pd.DataFrame: A DataFrame containing the sums for both datasets and the percentage of sensed values.
    """
    # Calculate the sum of all columns except 'G_woz_woni', 'geometry', 'index_right', and 'Lijn_Numbe' for the city-wide data
    cbs_sums = cbs_gdf.drop(columns=['crs28992', 'G_woz_woni', 'geometry']).sum()
    
    # Extract values from the Series
    values_t = cbs_sums.values

    # Calculate the sum of all columns except 'G_woz_woni', 'geometry', 'index_right', and 'Lijn_Numbe' for the sensed data
    sensed_sums = sensed_gdf.drop(columns=['crs28992', 'G_woz_woni', 'geometry', 'index_right', 'Lijn_Numbe','type']).sum()
    
    # Extract values from the Series
    values_s = sensed_sums.values

    # Extract keys (index) from the Series
    keys = sensed_sums.index

    # Create a new DataFrame
    data = {
        'Sociodemo': keys,
        'Sums_sensed': values_s,
        'Sums_total': values_t
    }

    sums = pd.DataFrame(data)
    
    # Calculate the percentage of sensed values relative to the city's total values
    sums['Sensed_%'] = ((sums['Sums_sensed'] / sums['Sums_total']) * 100).round(2)
    # Calculate exclusion
    sums['Excluded!'] = (sums['Sums_total'] - sums['Sums_sensed']).round(0)


    # Convert sums to integer values for cleaner display
    sums['Sums_sensed'] = sums['Sums_sensed'].astype(int)
    sums['Sums_total'] = sums['Sums_total'].astype(int)
    sums['Sensed_%'] = sums['Sensed_%'].astype(int)
    sums['Excluded!'] = sums['Excluded!'].astype(int)
    
    return sums


def normalize_statistics(merged_df):
    """
    Normalizes the columns in the merged_statistics DataFrame based on 'A_inhab'.
    Drops 'A_woning' and 'A_inhab', and rounds the results to two decimal places.
    
    Parameters:
        merged_statistics (pd.DataFrame): The DataFrame containing the statistics to normalize.

    Returns:
        pd.DataFrame: A DataFrame with normalized statistics.
    """
    # Create a copy of the DataFrame to avoid modifying the original
    average_stats = merged_df.copy()

    # Drop the 'A_woning' column
    average_stats = average_stats.drop(columns=['A_woning'])

    # Identify columns to normalize
    columns_to_normalize = average_stats.columns.difference(['Area', 'G_woz_woni', 'A_inhab'])

    # Normalize the identified columns
    average_stats[columns_to_normalize] = average_stats[columns_to_normalize].div(average_stats['A_inhab'], axis=0)

    # Drop the 'A_inhab' column
    average_stats = average_stats.drop(columns=['A_inhab'])

    # Round to two decimal places
    average_stats = average_stats.round(2)

    return average_stats


def plot_sums_and_percentages(df, buffer_distance):
    """
    Plots a bar chart comparing 'Sums_sensed' and 'Sums_total' from the given DataFrame,
    and optionally plots the percentage of sensed values.

    Parameters:
        df (pd.DataFrame): DataFrame containing 'Sums_sensed', 'Sums_total', and 'Sociodemo' columns.
        buffer_distance (float): Buffer distance used for the plot title.
    """
    # Create a figure and axis
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    # Define bar width and positions
    bar_width = 0.4
    index = range(len(df))
    
    # Bar chart for Sums_sensed and Sums_total
    bar1 = ax1.bar(index, df['Sums_sensed'], bar_width, label='Pop Sensed', color='#85b66f')
    bar2 = ax1.bar([i + bar_width for i in index], df['Sums_total'], bar_width, label='Pop Total', color='#ffa3c4')
    
    # Remove axis lines and ticks but keep labels
    for spine in ax1.spines.values():
        spine.set_visible(False)
    ax1.tick_params(axis='both', which='both', length=0)
    ax1.xaxis.set_tick_params(width=0)
    ax1.yaxis.set_tick_params(width=0)

    # Set labels and title
    ax1.set_xlabel('Sociodemographics (Inclusion)', fontweight='bold')
    #ax1.set_ylabel('Population', fontweight='bold')
    #ax1.set_title(f'public transport sensing {buffer_distance}m buffer', fontweight='bold', fontsize=14)
    ax1.set_xticks([i + bar_width / 2 for i in index])
    ax1.set_xticklabels(df['Sociodemo'], rotation=45, ha='right')
    
    # Add legends
    ax1.legend(loc='upper right')
    
    # Adjust layout and show plot
    plt.tight_layout()
    
    plt.show()

    # Save plot as PNG
    plt.savefig('sums_percentages_plot.png', bbox_inches='tight')
    plt.close(fig)




def plot_transport_and_population(lines, cbs_gdf, sensed_gdf, ams_gdf, buffer_distance, transport_type=None, lijn=None):
    """
    Plots public transport data, CBS population data, and sensed population with a buffer on a map.
    
    Parameters:
        lines (gpd.GeoDataFrame): GeoDataFrame containing the buffered transport data.
        cbs_gdf (gpd.GeoDataFrame): GeoDataFrame containing the CBS data.
        sensed_gdf (gpd.GeoDataFrame): GeoDataFrame containing the sensed data.
        ams_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary of the area to plot.
        buffer_distance (float): Distance used for buffering the transport data.
    """

 
    # Filter by transport type if provided
    if transport_type is not None:
        if isinstance(transport_type, list):
            lines = lines[lines["type"].isin(transport_type)]
        elif isinstance(transport_type, str):
            lines = lines[lines["type"] == transport_type]
        else:
            raise ValueError("transport_type must be a str or list of str")
    
    # Filter by line numbers if provided
    if lijn is not None:
        if isinstance(lijn, list):
            lines = lines[lines["Lijn_Numbe"].isin(lijn)]
        elif isinstance(lijn, int):
            lines = lines[lines["Lijn_Numbe"] == lijn]
        else:
            raise ValueError("lijn must be an int or list of int")

    # Change projection
    lines = lines.to_crs(ams_gdf.crs)
    cbs_gdf = cbs_gdf.to_crs(ams_gdf.crs)
    sensed_gdf = sensed_gdf.to_crs(ams_gdf.crs)
    
    # Set up the plot with a dark background
    plt.style.use('default')

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 10))
    
    
    # Plot the boundary
    ams_gdf.boundary.plot(ax=ax, linewidth=0.5, edgecolor='black')

    
    # Plot the CBS population
    cbs_gdf.plot(ax=ax, markersize=5, color='#e15989', edgecolor='white', linewidth=0.35)
    
    # Plot the sensed population
    sensed_gdf.plot(ax=ax, markersize=5, color='#85b66f', edgecolor='white', linewidth=0.35)

    # Plot the buffered trams on top
    lines.plot(ax=ax,color='black')  # Adjust linewidth as needed
    
    # Add title and labels
    ax.set_title(f'Public Transport Sensing {buffer_distance}m Buffer', fontweight='bold', fontsize=12)
    
    # Remove X and Y axes
    ax.set_axis_off()
    
    # Create custom legend handles
    line_handle = mlines.Line2D([], [], color='black', linestyle='-', linewidth=2, label='Lines')
    red_dot = mlines.Line2D([], [], color='#e15989', marker='o', markersize=5, linestyle='None', label='CBS Population')
    green_dot = mlines.Line2D([], [], color='#85b66f', marker='o', markersize=5, linestyle='None', label='Sensed Population')
    
    # Add legend
    ax.legend(handles=[line_handle, red_dot, green_dot], loc='upper right')

    # Show plot

    plt.show()
    
    # Save plot as PNG
    plt.savefig('transport_population_plot.png', bbox_inches='tight')
    plt.close(fig)


    
def plot_comparison(average_stats):
    """
    Plots a comparison of metrics for different areas using a scatter plot.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing the normalized statistics with columns ['Area', 'Metric', 'Value'].
        
    Returns:
        None
    """

    df = average_stats
    # Drop the column that should not be included in the plot
    df = df.drop(columns=['G_woz_woni'], errors='ignore')

    # Pivot the DataFrame for easier plotting
    df_melted = df.melt(id_vars='Area', var_name='Metric', value_name='Value')

    # Plotting
    plt.figure(figsize=(14, 6))
    ax = plt.gca()  # Get current axis

    # Define colors for the points
    colors = {
        'Amsterdam': '#ffa3c4',  # Pink color
        'Sensed Area': '#85b66f'  # Green color
    }

    # Plot scatter points for each area
    for area in df['Area'].unique():
        df_area = df_melted[df_melted['Area'] == area]
        plt.scatter(df_area['Metric'], df_area['Value'], color=colors.get(area, '#000000'), s=150, label=area, marker='o')

    # Customizing the plot
    plt.xlabel('Sociodemographics', fontweight='bold')
    plt.ylabel('Percentage %', fontweight='bold')
    plt.xticks(rotation=45)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    # Remove axis lines but keep tick labels
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')

    # Remove grid lines
    ax.grid(False)

    # Show plot
    #plt.show()

    plt.gcf()
    #plt.show()

    # Save plot as PNG
    plt.savefig('comparison.png', bbox_inches='tight')
    plt.close()

def master_function(transport_filepath, cbs_filepath, buffer_distance, line_number=None, transport_type=None, crs='EPSG:32633'):
    """
    Master function to handle the entire process from reading data to generating comparison statistics.

    Parameters:
        transport_filepath (str): Path to the transport data shapefile.
        cbs_filepath (str): Path to the CBS data shapefile.
        buffer_distance (float): Buffer distance for the transport data.
        line_number (int or None): Optional line number to filter transport data.
        crs (str): Coordinate Reference System for reprojecting data.

    Returns:
        buffered_trams_gdf (gpd.GeoDataFrame): GeoDataFrame with buffered transport data.
        sensed_cbs (gpd.GeoDataFrame): GeoDataFrame resulting from spatial join.
        summary_df (pd.DataFrame): DataFrame with concatenated summary statistics.
        comparison_df (pd.DataFrame): DataFrame comparing sums and percentages.
        average_stats (pd.DataFrame): DataFrame with normalized statistics.
    """
    
    # Step 1: Read and project transport data
    buffered_trams_gdf = read_and_project_transport_data(
        filepath=transport_filepath, 
        line_numbers=line_number, 
        crs=crs, 
        transport_type=transport_type
    )
    
    # Step 2: Calculate buffer
    buffered_trams_gdf = calculate_buffer(buffered_trams_gdf, buffer_distance, crs)
    
    # Step 3: Load and check CBS data
    cbs_gdf = load_and_check_cbs_data(cbs_filepath, buffered_trams_gdf)
    
    # Step 4: Perform spatial join
    sensed_cbs = spatial_join_and_remove_duplicates(cbs_gdf, buffered_trams_gdf)
    
    # Step 5: Generate summary statistics
    summary_cbs = generate_summary_statistics(cbs_gdf)
    summary_sensed = generate_summary_statistics_for_sensed_data(sensed_cbs)
    
    # Step 6: Concatenate summary statistics
    summary_df = concatenate_dataframes([summary_cbs, summary_sensed])
    
    # Step 7: Calculate and compare sums
    comparison_df = calculate_and_compare_sums(cbs_gdf, sensed_cbs)

    # Step 8: Normalize statistics
    average_stats = normalize_statistics(summary_df)
    
    return buffered_trams_gdf, sensed_cbs, summary_df, comparison_df, average_stats

def visualization_master_function(lines, cbs, sensed_cbs, ams_gdf, buffer_distance, comparison_df, average_stats, transport_type):
    """
    Master function for visualizing data comparisons.

    Parameters:
        comparison_df (pd.DataFrame): DataFrame containing 'Sums_sensed', 'Sums_total', and 'Sociodemo' columns.
        buffer_distance (float): Buffer distance used for the plot title.

    Returns:
        None
    """
    # Create plots
    fig1 = plot_sums_and_percentages(comparison_df, buffer_distance)
    fig2 = plot_transport_and_population(lines, cbs, sensed_cbs, ams_gdf, buffer_distance, transport_type)
    fig3 = plot_comparison(average_stats)

    return fig1, fig2, fig3






def process_frequencies(file_path):
    """
    Processes the spatial sensing data, calculates weighted sums, and returns the modified GeoDataFrame 
    with weighted columns (including geometry) and a separate DataFrame with ratios of weighted sums 
    to normal sums.
    
    Parameters:
    file_path (str): The path to the shapefile containing sensing data.
    
    Returns:
    tuple: A tuple containing:
        - gpd.GeoDataFrame: The modified GeoDataFrame with weighted columns, including geometry.
        - pd.DataFrame: A DataFrame containing the ratios of weighted sums to normal sums, rounded to 1 decimal point.
    """
    # Load the shapefile
    freq_cbs = gpd.read_file(file_path)
    
    # Define the columns to analyze (including 'geometry')
    columns_to_weight = ['A_inhab', 'A_0_15', 'A_15_25', 'A_25_45', 'A_45_65', 'A_65+', 'A_woning', 
                         'A_nederlan', 'A_west_mig', 'A_n_west_m']

    # Create a new GeoDataFrame for weighted values, preserving crs28992, points_cou, and geometry
    weighted_freq_cbs = freq_cbs[['crs28992', 'count', 'geometry']].copy()

    # Calculate weighted values for the specified columns
    for col in columns_to_weight:
        weighted_freq_cbs[col] = freq_cbs[col] * freq_cbs['count']

    # Drop the 'G_woz_woni' column if it exists
    weighted_freq_cbs = weighted_freq_cbs.drop(columns=['G_woz_woni'], errors='ignore')

    # Define the columns to sum
    columns_to_sum = columns_to_weight
    
    # Initialize dictionaries for weighted sums, normal sums, and ratios
    weighted_sums = {}
    normal_sums = {}
    ratios = {}

    # Loop over the columns to calculate the weighted sum and total sum
    for col in columns_to_sum:
        weighted_sums[col] = (freq_cbs[col] * freq_cbs['count']).sum()
        normal_sums[col] = freq_cbs[col].sum()
        ratios[col] = weighted_sums[col] / normal_sums[col]

    # Convert the ratios dictionary into a DataFrame
    ratios_df = pd.DataFrame([ratios])

    # Round the ratios to one decimal point
    ratios_df = ratios_df.round(1)
    
    # Ensure weighted_freq_cbs remains a GeoDataFrame
    weighted_freq_cbs = gpd.GeoDataFrame(weighted_freq_cbs, geometry='geometry')

    # Replace NaN with 0 for new weighted columns and round to 0 decimal places
    weighted_freq_cbs[columns_to_weight] = weighted_freq_cbs[columns_to_weight].fillna(0).round(0).astype(int)

    # Change names
    weighted_freq_cbs.rename(columns=lambda x: x.replace('A_', 'Weight_') if x.startswith('A_') else x, inplace=True)

    # Sorting
    weighted_freq_cbs = weighted_freq_cbs.sort_values(by='Weight_inhab', ascending=False)

    # Fill NA values
    weighted_freq_cbs = weighted_freq_cbs.fillna(0)


    # Reorder the columns
    columns_ordered = [col for col in weighted_freq_cbs.columns if col not in ['geometry', 'count']] + ['count', 'geometry']
    weighted_freq_cbs = weighted_freq_cbs[columns_ordered]

    # Return the weighted_freq_cbs GeoDataFrame and ratios_df DataFrame separately
    return weighted_freq_cbs, ratios_df



def plot_counts(weighted_freq_cbs, ams_gdf):
    """
    Plots the counts from the weighted_freq_cbs GeoDataFrame with the boundary of the area (ams_gdf) on a map.
    
    Parameters:
        weighted_freq_cbs (gpd.GeoDataFrame): GeoDataFrame containing the weighted data including counts.
        ams_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary of the area to plot.
    """
    # Filter out rows where count is zero
    weighted_freq_cbs = weighted_freq_cbs[weighted_freq_cbs['count'] > 0]
    
    # Change projection
    weighted_freq_cbs = weighted_freq_cbs.to_crs(ams_gdf.crs)
    
    # Set up the plot with a dark background
    plt.style.use('default')

    # Define custom colormap from #e15989 to #85b66f
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ['#e15989', '#85b66f'])

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 10))
    
    # Plot the boundary
    ams_gdf.boundary.plot(ax=ax, linewidth=0.5, edgecolor='black')

    # Create a classifier for equal counts (quantiles)
    classifier = mapclassify.Quantiles(weighted_freq_cbs['count'], k=8)

    # Plot the counts from weighted_freq_cbs
    weighted_freq_cbs.plot(column='count', ax=ax, markersize=5, cmap=custom_cmap, 
                           legend=True, scheme='quantiles', classification_kwds={'k': 8},
                           edgecolor='white', linewidth=0.35)
    
    # Add title and labels
    ax.set_title('Frequency Map - Amount of Measurements', fontweight='bold', fontsize=12)
    
    # Remove X and Y axes
    ax.set_axis_off()
    
    # Show plot
    plt.show()


def plot_weighted_column(weighted_freq_cbs, ams_gdf, column_to_plot):
    """
    Plots the selected column from the weighted_freq_cbs GeoDataFrame with the boundary of the area (ams_gdf) on a map.
    
    Parameters:
        weighted_freq_cbs (gpd.GeoDataFrame): GeoDataFrame containing the weighted data including counts.
        ams_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary of the area to plot.
        column_to_plot (str): The column name from weighted_freq_cbs to visualize.
    """
    # Ensure the column exists in the GeoDataFrame
    if column_to_plot not in weighted_freq_cbs.columns:
        raise ValueError(f"Column '{column_to_plot}' not found in the provided GeoDataFrame.")
    
    # Filter out rows where the selected column is zero
    weighted_freq_cbs = weighted_freq_cbs[weighted_freq_cbs[column_to_plot] > 0]

    # Change projection to match the boundary
    weighted_freq_cbs = weighted_freq_cbs.to_crs(ams_gdf.crs)
    
    # Set up the plot with a dark background
    plt.style.use('default')

    # Define custom colormap from #e15989 to #85b66f
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", ['#e15989', '#85b66f'])

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 10))
    
    # Plot the boundary
    ams_gdf.boundary.plot(ax=ax, linewidth=0.5, edgecolor='black')

    # Create a classifier for equal counts (quantiles)
    classifier = mapclassify.Quantiles(weighted_freq_cbs[column_to_plot], k=8)

    # Plot the selected column using the custom colormap and 8 quantile classes (equal count)
    weighted_freq_cbs.plot(column=column_to_plot, ax=ax, markersize=5, cmap=custom_cmap, 
                           legend=True, scheme='quantiles', classification_kwds={'k': 8},
                           edgecolor='white', linewidth=0.35)

    # Add title and labels
    ax.set_title(f'Frequencies x Units ({column_to_plot})', fontweight='bold', fontsize=12)
    
    # Remove X and Y axes
    ax.set_axis_off()
    
    # Show plot
    plt.show()


def plot_ratios_comparison(ratios_df):
    """
    Plots a comparison of demographic ratios using a scatter plot.
    
    Parameters:
        ratios_df (pd.DataFrame): DataFrame containing the demographic ratio columns.

    Returns:
        None
    """
    # Prepare the DataFrame for plotting
    ratios_df_melted = ratios_df.melt(var_name='Metric', value_name='Value')

    # Create the plot
    plt.figure(figsize=(14, 6))
    ax = plt.gca()  # Get current axis

    # Define colors for the points
    colors = ['#85b66f' if metric != 'A_inhab' else '#ffa3c4' for metric in ratios_df.columns]

    # Scatter plot for demographic ratios
    ax.scatter(ratios_df.columns, ratios_df.iloc[0], color=colors, s=150)

    # Draw a horizontal pink line at the height of A_inhab
    ax.axhline(y=ratios_df['A_inhab'].iloc[0], color='#ffa3c4', linestyle='--', linewidth=2, label='A_inhab Line')

    # Customizing the plot
    plt.xlabel('Sociodemographics', fontweight='bold')
    plt.ylabel('Measurements', fontweight='bold')
    plt.xticks(rotation=45)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Remove axis lines but keep tick labels
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none')

    # Add title
    ax.set_title('Frequencies per Person/Unit', fontweight='bold', fontsize=12)

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()


# ### Master Function Visualistions

# In[20]:


def visualize_frequencies(weighted_freq_cbs, ams_gdf, column_to_plot, ratios_df):
    """
    Master function to visualize frequencies by calling respective plotting functions.

    Parameters:
        weighted_freq_cbs (gpd.GeoDataFrame): GeoDataFrame containing the weighted data including counts.
        ams_gdf (gpd.GeoDataFrame): GeoDataFrame containing the boundary of the area to plot.
        column_to_plot (str): The column name from weighted_freq_cbs to visualize.
        ratios_df (pd.DataFrame): DataFrame containing the demographic ratio columns. Defaults to None.
    """

    # Plot all three visualizations
    plot_counts(weighted_freq_cbs, ams_gdf)
    plot_weighted_column(weighted_freq_cbs, ams_gdf, column_to_plot)
    plot_ratios_comparison(ratios_df)

