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
tab1, tab2 = st.tabs(["Analysis", "Frequencies"])

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
    # Define the fixed file path
    file_path = "freq_cbs_1304_fullday.shp"

    # Process frequencies using the fixed file path
    weighted_freq_cbs, ratios_df = process_frequencies(file_path)

    # Display initial map and ratios DataFrame
    st.subheader("Ratios Comparison")
    plot_ratios_comparison(ratios_df)
    st.image('ratios.png', use_column_width=True)
    
    # Display initial map
    st.subheader("Counts Visualization")
    plot_counts(weighted_freq_cbs, ams_gdf)
    st.image('counts_cbs.png', use_column_width=True)

    # User inputs for frequency visualization
    options = ['Weight_inhab'] + [col for col in weighted_freq_cbs.columns[3:-1] if col != 'count'] 
    column_to_plot = st.selectbox("Select Column to Plot", options=options)
    
    # Button to perform frequency analysis and display results
    if st.button("Run Weighted Analysis"):
        st.subheader("Frequency Visualization")
        plot_weighted_column(weighted_freq_cbs, ams_gdf, column_to_plot)
        st.image('weights_cbs.png', use_column_width=True)
