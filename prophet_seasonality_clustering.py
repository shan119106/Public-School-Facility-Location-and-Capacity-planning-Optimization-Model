"""
Prophet Forecasting with Seasonality-Based Clustering

This script performs time series forecasting using Prophet and clusters weeks
based on the seasonality components of ATV and bookings rather than raw values.
"""

import pandas as pd
import matplotlib.pyplot as plt
from prophet import Prophet
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    from prophet.plot import plot_plotly, plot_components_plotly
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


def prepare_prophet_data(pandas_data, value_column):
    """
    Prepare data for Prophet model.
    
    Args:
        pandas_data: DataFrame with year, week, and value columns
        value_column: Name of the column to forecast
        
    Returns:
        DataFrame with 'ds' (datetime) and prepared data
    """
    pandas_data['week'] = pandas_data['week'].astype(int)
    pandas_data['year'] = pandas_data['year'].astype(int)
    pandas_data['ds'] = pd.to_datetime(
        pandas_data['year'].astype(str) + pandas_data['week'].astype(str) + '1', 
        format='%Y%W%w'
    )
    return pandas_data


def fit_prophet_model(prophet_df, value_column):
    """
    Fit a Prophet model to the data.
    
    Args:
        prophet_df: DataFrame with 'ds' and value columns
        value_column: Name of the value column
        
    Returns:
        Fitted Prophet model and forecast DataFrame
    """
    # Prepare Prophet DataFrame
    prophet_input = prophet_df[['ds', value_column]].rename(columns={value_column: 'y'})
    
    # Fit Prophet model
    model = Prophet(weekly_seasonality=True, yearly_seasonality=True)
    model.fit(prophet_input)
    
    # Make future dataframe and forecast
    future = model.make_future_dataframe(periods=52, freq='W')  # 1 year ahead
    forecast = model.predict(future)
    
    return model, forecast


def extract_seasonality_components(forecast, pandas_data):
    """
    Extract seasonality components from Prophet forecast for the historical period.
    
    Args:
        forecast: Prophet forecast DataFrame
        pandas_data: Original data with 'ds' column
        
    Returns:
        DataFrame with seasonality components aligned with original data
    """
    # Get only the historical dates (matching original data)
    historical_forecast = forecast[forecast['ds'].isin(pandas_data['ds'])].copy()
    
    # Extract seasonality components
    # Prophet stores weekly_seasonality and yearly_seasonality separately
    seasonality_df = historical_forecast[['ds', 'weekly', 'yearly']].copy()
    
    # Combine weekly and yearly seasonality into a total seasonality measure
    seasonality_df['total_seasonality'] = (
        seasonality_df['weekly'] + seasonality_df['yearly']
    )
    
    return seasonality_df


def plot_forecast(model, forecast, title, y_label):
    """
    Plot Prophet forecast results.
    
    Args:
        model: Fitted Prophet model
        forecast: Forecast DataFrame
        title: Plot title
        y_label: Y-axis label
    """
    if PLOTLY_AVAILABLE:
        # Use interactive plotly plots if available
        # Main forecast plot
        fig1 = plot_plotly(model, forecast)
        fig1.update_layout(
            yaxis_title=y_label,
            title=title
        )
        fig1.show()
        
        # Components plot (trend and seasonality breakdown)
        fig2 = plot_components_plotly(model, forecast)
        fig2.update_yaxes(title_text=y_label, row=1, col=1)
        fig2.update_layout(title=f"{y_label} Components (Trend and Seasonality)")
        fig2.show()
    else:
        # Use matplotlib plots as fallback
        from prophet.plot import plot, plot_components
        
        # Main forecast plot
        fig1 = model.plot(forecast)
        fig1.suptitle(title)
        plt.ylabel(y_label)
        plt.show()
        
        # Components plot
        fig2 = model.plot_components(forecast)
        fig2.suptitle(f"{y_label} Components (Trend and Seasonality)")
        plt.show()


def perform_seasonality_clustering(atv_seasonality, bookings_seasonality, n_clusters=3):
    """
    Perform clustering based on seasonality components.
    
    Args:
        atv_seasonality: DataFrame with ATV seasonality components
        bookings_seasonality: DataFrame with bookings seasonality components
        n_clusters: Number of clusters
        
    Returns:
        DataFrame with cluster assignments and visualization
    """
    # Merge seasonality data
    merged_seasonality = pd.merge(
        atv_seasonality[['ds', 'total_seasonality']].rename(
            columns={'total_seasonality': 'ATV_seasonality'}
        ),
        bookings_seasonality[['ds', 'total_seasonality']].rename(
            columns={'total_seasonality': 'bookings_seasonality'}
        ),
        on='ds'
    )
    
    # Scale the seasonality features
    scaler = StandardScaler()
    merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']] = (
        scaler.fit_transform(
            merged_seasonality[['ATV_seasonality', 'bookings_seasonality']]
        )
    )
    
    # Perform K-means clustering on seasonality
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    merged_seasonality['cluster'] = kmeans.fit_predict(
        merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']]
    )
    
    # Visualize clusters
    plt.figure(figsize=(10, 6))
    for i in range(n_clusters):
        cluster_data = merged_seasonality[merged_seasonality['cluster'] == i]
        plt.scatter(
            cluster_data['ATV_seasonality'], 
            cluster_data['bookings_seasonality'], 
            label=f'Cluster {i}',
            alpha=0.6
        )
    
    plt.xlabel('ATV Seasonality')
    plt.ylabel('Bookings Seasonality')
    plt.title('Weekly Clusters based on ATV and Bookings Seasonality')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # Print cluster centers
    print("Cluster Centers (scaled):")
    print(kmeans.cluster_centers_)
    
    # Print cluster statistics
    print("\nCluster Statistics:")
    for i in range(n_clusters):
        cluster_data = merged_seasonality[merged_seasonality['cluster'] == i]
        print(f"\nCluster {i}:")
        print(f"  Size: {len(cluster_data)}")
        print(f"  ATV Seasonality - Mean: {cluster_data['ATV_seasonality'].mean():.4f}, "
              f"Std: {cluster_data['ATV_seasonality'].std():.4f}")
        print(f"  Bookings Seasonality - Mean: {cluster_data['bookings_seasonality'].mean():.4f}, "
              f"Std: {cluster_data['bookings_seasonality'].std():.4f}")
    
    return merged_seasonality, kmeans


def main_workflow(atv_query_result, bookings_query_result):
    """
    Main workflow for Prophet forecasting and seasonality-based clustering.
    
    Args:
        atv_query_result: DataFrame with columns ['year', 'week', 'ATV']
        bookings_query_result: DataFrame with columns ['year', 'week', 'bookings']
    """
    print("=" * 60)
    print("Prophet Forecasting with Seasonality-Based Clustering")
    print("=" * 60)
    
    # ========== ATV Forecasting ==========
    print("\n1. Processing ATV data...")
    pandas_data_atv = prepare_prophet_data(atv_query_result.copy(), 'ATV')
    model_atv, forecast_atv = fit_prophet_model(pandas_data_atv, 'ATV')
    print("   ATV model fitted successfully.")
    
    # Plot ATV forecast
    plot_forecast(model_atv, forecast_atv, "ATV Forecast with Prophet", "ATV")
    
    # Extract ATV seasonality
    atv_seasonality = extract_seasonality_components(forecast_atv, pandas_data_atv)
    print(f"   Extracted ATV seasonality for {len(atv_seasonality)} data points.")
    
    # ========== Bookings Forecasting ==========
    print("\n2. Processing bookings data...")
    pandas_data_bookings = prepare_prophet_data(bookings_query_result.copy(), 'bookings')
    model_bookings, forecast_bookings = fit_prophet_model(pandas_data_bookings, 'bookings')
    print("   Bookings model fitted successfully.")
    
    # Plot bookings forecast
    plot_forecast(model_bookings, forecast_bookings, "Bookings Forecast with Prophet", "Bookings")
    
    # Extract bookings seasonality
    bookings_seasonality = extract_seasonality_components(forecast_bookings, pandas_data_bookings)
    print(f"   Extracted bookings seasonality for {len(bookings_seasonality)} data points.")
    
    # ========== Seasonality-Based Clustering ==========
    print("\n3. Performing seasonality-based clustering...")
    merged_seasonality, kmeans = perform_seasonality_clustering(
        atv_seasonality, 
        bookings_seasonality, 
        n_clusters=3
    )
    
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    
    return {
        'atv_model': model_atv,
        'bookings_model': model_bookings,
        'atv_forecast': forecast_atv,
        'bookings_forecast': forecast_bookings,
        'atv_seasonality': atv_seasonality,
        'bookings_seasonality': bookings_seasonality,
        'merged_seasonality': merged_seasonality,
        'kmeans': kmeans
    }


# Example usage with Spark SQL queries:
"""
# Step 1: Run Spark queries to get ATV data
query_atv = '''
SELECT 
  year(pnr.flt_dptr_date_d) AS year,
  weekofyear(pnr.flt_dptr_date_d) AS week,
  SUM(pnr.pax_ct_i * pnr.rbk_od_avg_fare) / SUM(pnr.pax_ct_i) AS ATV
FROM 
  rm_workspace.mcla_mxbus_seasonality_pnr_data_p2 pnr
GROUP BY 
  year(pnr.flt_dptr_date_d), weekofyear(pnr.flt_dptr_date_d)
ORDER BY 
  year, week
'''
spark_data_atv = spark.sql(query_atv)
pandas_data_atv = spark_data_atv.toPandas()

# Step 2: Run Spark queries to get bookings data
query_bookings = '''
SELECT 
  year(pnr.flt_dptr_date_d) AS year,
  weekofyear(pnr.flt_dptr_date_d) AS week,
  SUM(pnr.pax_ct_i) AS bookings
FROM 
  rm_workspace.mcla_mxbus_seasonality_pnr_data_p2 pnr
GROUP BY 
  year(pnr.flt_dptr_date_d), weekofyear(pnr.flt_dptr_date_d)
ORDER BY 
  year, week
'''
spark_data_bookings = spark.sql(query_bookings)
pandas_data_bookings = spark_data_bookings.toPandas()

# Step 3: Run the main workflow
results = main_workflow(pandas_data_atv, pandas_data_bookings)

# Access results:
# - results['merged_seasonality']: DataFrame with seasonality values and cluster assignments
# - results['kmeans']: Fitted KMeans model
# - results['atv_forecast']: ATV forecast with all components
# - results['bookings_forecast']: Bookings forecast with all components
"""
