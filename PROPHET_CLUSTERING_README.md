# Prophet Forecasting with Seasonality-Based Clustering

This module provides functionality to perform time series forecasting using Facebook Prophet and cluster data points based on their **seasonality components** rather than raw values.

## Problem Context

When analyzing time series data like Average Transaction Value (ATV) and bookings, clustering based on raw values may not capture the underlying seasonal patterns. This implementation extracts the seasonality components from Prophet models and uses them for clustering, which better identifies periods with similar seasonal behavior.

## Key Features

- **Prophet Forecasting**: Fits Prophet models to ATV and bookings data with weekly and yearly seasonality
- **Seasonality Extraction**: Extracts weekly and yearly seasonality components from the forecast
- **Seasonality-Based Clustering**: Performs K-means clustering on the seasonality components instead of raw values
- **Visualization**: Provides plots for forecasts, components, and cluster assignments

## Files

- `prophet_seasonality_clustering.py`: Main module with all functionality
- `example_seasonality_clustering.py`: Example script with sample data generation

## Installation

Install required dependencies:

```bash
pip install pandas numpy matplotlib prophet scikit-learn
```

## Usage

### With Spark SQL Data

```python
from prophet_seasonality_clustering import main_workflow

# Query ATV data from Spark
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
pandas_data_atv = spark.sql(query_atv).toPandas()

# Query bookings data from Spark
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
pandas_data_bookings = spark.sql(query_bookings).toPandas()

# Run the analysis
results = main_workflow(pandas_data_atv, pandas_data_bookings)
```

### With Sample Data

```python
# Run the example script
python example_seasonality_clustering.py
```

## How It Works

### 1. Data Preparation
- Converts year/week to datetime format
- Prepares data in Prophet's expected format

### 2. Prophet Modeling
- Fits separate Prophet models for ATV and bookings
- Enables weekly and yearly seasonality components
- Generates forecasts for 1 year ahead

### 3. Seasonality Extraction
- Extracts `weekly` and `yearly` seasonality components from forecast
- Combines them into a `total_seasonality` measure
- Aligns seasonality with historical data points

### 4. Clustering
- Standardizes seasonality features using StandardScaler
- Applies K-means clustering (default: 3 clusters)
- Clusters based on seasonality patterns, not raw values

### 5. Visualization
- Forecast plots with confidence intervals
- Component breakdowns (trend, weekly, yearly seasonality)
- Scatter plot of clusters in seasonality space

## Key Difference from Original Approach

### Original Approach (Raw Values)
```python
# Clusters based on raw ATV and bookings values
merged_data[['ATV_scaled', 'bookings_scaled']] = scaler.fit_transform(
    merged_data[['ATV', 'bookings']]
)
merged_data['cluster'] = kmeans.fit_predict(
    merged_data[['ATV_scaled', 'bookings_scaled']]
)
```

### New Approach (Seasonality Components)
```python
# Clusters based on seasonality components from Prophet
merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']] = scaler.fit_transform(
    merged_seasonality[['ATV_seasonality', 'bookings_seasonality']]
)
merged_seasonality['cluster'] = kmeans.fit_predict(
    merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']]
)
```

## Output

The `main_workflow` function returns a dictionary with:

- `atv_model`: Fitted Prophet model for ATV
- `bookings_model`: Fitted Prophet model for bookings
- `atv_forecast`: Full forecast DataFrame for ATV
- `bookings_forecast`: Full forecast DataFrame for bookings
- `atv_seasonality`: Seasonality components for ATV
- `bookings_seasonality`: Seasonality components for bookings
- `merged_seasonality`: Combined seasonality data with cluster assignments
- `kmeans`: Fitted K-means clustering model

## Example Output

```
Cluster Statistics:

Cluster 0:
  Size: 35
  ATV Seasonality - Mean: 5.2341, Std: 2.1234
  Bookings Seasonality - Mean: 45.6789, Std: 15.3456

Cluster 1:
  Size: 42
  ATV Seasonality - Mean: -3.4567, Std: 1.8765
  Bookings Seasonality - Mean: -32.4567, Std: 12.6789

Cluster 2:
  Size: 27
  ATV Seasonality - Mean: 0.1234, Std: 1.2345
  Bookings Seasonality - Mean: 2.3456, Std: 8.9012
```

## Benefits of Seasonality-Based Clustering

1. **Pattern Recognition**: Identifies weeks with similar seasonal behavior regardless of absolute values
2. **Trend Independence**: Removes the influence of long-term trends
3. **Better Segmentation**: Groups periods by seasonal characteristics (e.g., holiday patterns, seasonal demand)
4. **Actionable Insights**: Helps identify when to apply similar strategies based on seasonal patterns

## License

This code is part of the Public School Facility Location and Capacity Planning Optimization Model project.
