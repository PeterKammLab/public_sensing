import streamlit as st
from functions import *  
from functions import visualize_frequencies

import time
import math 
import pandas as pd
import folium
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import wkt
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Point, Polygon
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap
import mapclassify

# Load initial data
ams_gdf = gpd.read_file("gemeente_T.shp")
transport_filepath = 'public_amsterdam.shp'
cbs_filepath = 'cbs_amsterdam_2021_clean.shp'

# Streamlit interface
st.title("Public Transport Sensing")

# Use tabs instead of sidebar
tab1, tab2, tab3 = st.tabs(["**Spatial Analysis**", "**Frequencies**", "**Weights & Index**"])

with tab1:
    # User input for transport type
    transport_type = st.selectbox(
        "Select Transport Type",
        (None, 'Bus', 'Tram', 'Night Bus')
    )

    # User input for line type 
    line_number = st.selectbox(
        "Select Your Favourite Line",
        (None, '1', '2', '3', '4', '5', '7', '12', '13', '14', '17', '19', 
         '24', '26', '18', '21', '22', '34', '35', '36', '37', '38', 
         '41', '43', '44', '47', '49', '61', '63', '65', '68', '231', 
         '232', '233', 'N81', 'N82', 'N83', 'N86', 'N87', 'N88', 
         'N89', 'N91', 'N93', '461', '463', '464', '48', '66', 
         '267', '369', '62', 'N85', '15', 'N84', '40', '245', 
         '246', '247', '76')
    )

    # User input for buffer distance
    buffer_distance = st.slider(
        "Select Buffer Distance (meters)",
        min_value=0,
        max_value=1000,
        value=100
    )

    # Button to perform analysis
    if st.button("Run Analysis"):
        # Call the analysis function
        buffered_trams_gdf, sensed_cbs, merged_df, comparison_df, average_stats = master_function(
            transport_filepath, cbs_filepath, buffer_distance, line_number, transport_type
        )
        lines = read_and_project_transport_data("public_amsterdam.shp")
        buffered_lines = calculate_buffer(lines, buffer_distance)
        cbs = load_and_check_cbs_data("cbs_amsterdam_2021_clean.shp", buffered_lines)

        # Display results
        plot_transport_and_population(lines, cbs, sensed_cbs, ams_gdf, buffer_distance, transport_type)
        st.image('transport_population_plot.png', use_column_width=True)
        
        plot_sums_and_percentages(comparison_df, buffer_distance)
        st.image('sums_percentages_plot.png', use_column_width=True)

        plot_comparison(average_stats)
        st.image('comparison.png', use_column_width=True)

        st.markdown("#### Who do we sense for?", unsafe_allow_html=True)
        st.write(comparison_df)
        st.write(merged_df)

with tab2:
    # Add an h6 Markdown title
    # Add the text with specified bold formatting
     # Display initial map
    
    st.markdown("based on real-time GTFS data, **GVB trams** * **buses** for **13th March 2024** in Amsterdam")

    st.subheader("")
    
    # Define the fixed file path
    file_path = "freq_cbs_1304_fullday.shp"

    # Process frequencies using the fixed file path
    weighted_freq_cbs, ratios_df, freq_cbs = process_frequencies(file_path)

    
    plot_counts(weighted_freq_cbs, ams_gdf)
    st.image('counts_cbs.png', use_column_width=True)

    # Add a description under the map
    st.markdown(
        "<p style='color: grey; font-size: 12.5px;'>This plot shows amount of measurements for each cell on a given day</p>",
        unsafe_allow_html=True
    )

    #<br>
    
    # Display initial map and ratios DataFrame
    st.subheader("")
    plot_ratios_comparison(ratios_df)
    st.image('ratios.png', use_column_width=True)

     # Add a description under the graph
    st.markdown(
    "<p style='color: grey; font-size: 12.5px;'>This graph shows the number of measurements for all cells for each individual or unit, <br> categorized by a specific sociodemographic attribute on a given day.</p>",
    unsafe_allow_html=True
    )
    st.subheader("")
    
    st.markdown("<h6 style='color: black;'>Number of Measurements per Person/Unit (work in progress)</h5>", unsafe_allow_html=True)


    label_mapping = {
        'Inhabitants': 'inhab_index',
        'Age 0-15': '0_15_index',
        'Age 15-25': '15_25_index',
        'Age 25-45': '25_45_index',
        'Age 45-65': '45_65_index',
        'Age 65+': '65+_index',
        'Housing Units': 'woning_index',
        'Dutch': 'nederlan_index',
        'West. Migration': 'west_mig_index',
        'Non-West. Migration': 'n_west_m_index'
    }

    # Create options for the selectbox with user-friendly labels
    options = list(label_mapping.keys())

    

    # Add a unique key to the selectbox
    column_to_plot_label = st.selectbox("Select Attribute", options=options, key="unique_key_for_selectbox")


# Add a description under the map 2
    st.markdown(
    "<p style='color: grey; font-size: 12.5px;'>This map shows the number of measurements per cell for each individual or unit, <br> categorized by a specific sociodemographic attribute on a given day.</p>",
    unsafe_allow_html=True
    )
    
   # Normalize the weights after loading the frequencies
    normalized_gdf = normalize_weights_and_merge(weighted_freq_cbs, freq_cbs)



# Assuming you have your Streamlit tab setup before this

# Assuming you have your Streamlit tab setup before this
with tab3:

    st.markdown("based on real-time GTFS data, **GVB trams** * **buses** for **13th March 2024** in Amsterdam")

    st.subheader("")
    # Define the mapping of user-friendly names to actual column names for index analysis
    index_label_mapping = {
        'Inhabitants': 'inhab_index',
        'Age 0-15': '0_15_index',
        'Age 15-25': '15_25_index',
        'Age 25-45': '25_45_index',
        'Age 45-65': '45_65_index',
        'Age 65+': '65+_index',
        'Housing Units': 'woning_index',
        'Dutch': 'nederlan_index',
        'West. Migration': 'west_mig_index',
        'Non-West. Migration': 'n_west_m_index'
    }

    # Define the mapping of user-friendly names to actual column names for weighted analysis
    weighted_label_mapping = {
        'Inhabitants': 'Weight_inhab',
        'Age 0-15': 'Weight_0_15',
        'Age 15-25': 'Weight_15_25',
        'Age 25-45': 'Weight_25_45',
        'Age 45-65': 'Weight_45_65',
        'Age 65+': 'Weight_65+',
        'Housing Units': 'Weight_woning',
        'Dutch': 'Weight_nederlan',
        'West. Migration': 'Weight_west_mig',
        'Non-West. Migration': 'Weight_n_west_m'
    }

    # Create options for the selectbox with user-friendly labels
    options = list(index_label_mapping.keys())

    #st.subheader("Weights & Index for Sensing")

    # Select column to plot using the user-friendly names
    column_to_plot_label = st.selectbox("Select Attribute", options=options)

    # Get the actual column name for index analysis
    column_to_plot_index = index_label_mapping[column_to_plot_label]
    
    # Get the actual column name for weighted analysis
    column_to_plot_weighted = weighted_label_mapping[column_to_plot_label]

    # Button to perform weighted frequency analysis and display results
    if st.button("Run Weighted Analysis"):
        plot_weighted_column(weighted_freq_cbs, ams_gdf, column_to_plot_weighted)  # Use the normalized GeoDataFrame
        st.image('weights_cbs.png', use_column_width=True)
        # Add a description under the map 2
        st.markdown(
        "<p style='color: grey; font-size: 12.5px;'>Weighted value per cell, calculated as the product of the number of measurements and the number of persons/units <br> categorized by a specific sociodemographic attribute on a given day.</p>",
        unsafe_allow_html=True
        )

    # Button to perform index analysis and display results
    if st.button("Run Index Analysis"):
        plot_weighted_column(normalized_gdf, ams_gdf, column_to_plot_index)  # Use weighted_freq_cbs as input
        st.image('weights_cbs.png', use_column_width=True)
        # Add a description under the map 2
        st.markdown(
        "<p style='color: grey; font-size: 12.5px;'>**'Sensing Index'** representing sensing potential based on the number of persons/units and amount of measurements per cell</p>",
        unsafe_allow_html=True
        )
