# Quick Start Guide: Prophet Seasonality Clustering with Spark SQL

This guide shows how to use the Prophet seasonality clustering with your Spark SQL data.

## Installation

```bash
pip install -r requirements_prophet.txt
```

## Usage with Spark SQL

### Step 1: Query Your Data

```python
# Query ATV data
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

# Query bookings data
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
```

### Step 2: Run Seasonality Clustering

```python
from prophet_seasonality_clustering import main_workflow

# Run the complete analysis
results = main_workflow(pandas_data_atv, pandas_data_bookings)

# Access the results
merged_seasonality = results['merged_seasonality']
print(merged_seasonality[['ds', 'ATV_seasonality', 'bookings_seasonality', 'cluster']])
```

### Step 3: Use the Cluster Assignments

```python
# Get the seasonality data with clusters
cluster_df = results['merged_seasonality'][['ds', 'cluster', 'ATV_seasonality', 'bookings_seasonality']]

# You can now use these clusters for various purposes:
# 1. Identify weeks with similar seasonal behavior
# 2. Apply different strategies for different seasonal clusters
# 3. Forecast future weeks based on which cluster they belong to

# Example: Find all weeks in cluster 0
cluster_0_weeks = cluster_df[cluster_df['cluster'] == 0]
print(f"Cluster 0 has {len(cluster_0_weeks)} weeks")
print(cluster_0_weeks)
```

## What's the Difference from the Original?

### Original Approach (Clustering on Raw Values)
```python
# This clusters based on absolute ATV and bookings values
merged_data[['ATV_scaled', 'bookings_scaled']] = scaler.fit_transform(
    merged_data[['ATV', 'bookings']]
)
merged_data['cluster'] = kmeans.fit_predict(
    merged_data[['ATV_scaled', 'bookings_scaled']]
)
```

**Problem**: This approach is heavily influenced by trends. Early weeks with low ATV/bookings get clustered together, late weeks with high values get clustered together, even if they have similar seasonal patterns.

### New Approach (Clustering on Seasonality)
```python
# This clusters based on seasonal components extracted by Prophet
merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']] = scaler.fit_transform(
    merged_seasonality[['ATV_seasonality', 'bookings_seasonality']]
)
merged_seasonality['cluster'] = kmeans.fit_predict(
    merged_seasonality[['ATV_seasonality_scaled', 'bookings_seasonality_scaled']]
)
```

**Benefit**: This approach removes the trend and focuses purely on seasonal patterns. Weeks with similar seasonal behavior get clustered together regardless of when they occur or their absolute values.

## Examples

### Example 1: Testing with Sample Data
```bash
# Run the example with synthetic data
python example_seasonality_clustering.py
```

### Example 2: Comparing Approaches
```bash
# See the difference between raw and seasonality clustering
python comparison_clustering.py
```

### Example 3: Running Tests
```bash
# Verify everything works correctly
python test_seasonality_clustering.py
```

## Understanding the Results

The `main_workflow` function returns a dictionary with:

- `atv_model`: Prophet model fitted to ATV data
- `bookings_model`: Prophet model fitted to bookings data
- `atv_forecast`: Full forecast DataFrame for ATV (including trend, seasonality)
- `bookings_forecast`: Full forecast DataFrame for bookings
- `atv_seasonality`: Extracted seasonality components for ATV
- `bookings_seasonality`: Extracted seasonality components for bookings
- `merged_seasonality`: Combined data with cluster assignments
- `kmeans`: Fitted KMeans clustering model

### Key Columns in `merged_seasonality`

- `ds`: Date (datetime)
- `ATV_seasonality`: Total seasonality for ATV (weekly + yearly)
- `bookings_seasonality`: Total seasonality for bookings (weekly + yearly)
- `ATV_seasonality_scaled`: Standardized ATV seasonality
- `bookings_seasonality_scaled`: Standardized bookings seasonality
- `cluster`: Cluster assignment (0, 1, 2, ...)

## Visualization

The workflow automatically generates plots:
1. Forecast plots with confidence intervals
2. Component breakdowns (trend, weekly, yearly seasonality)
3. Cluster scatter plots in seasonality space

## Advanced Usage

### Customize Number of Clusters
```python
from prophet_seasonality_clustering import perform_seasonality_clustering

# Use 5 clusters instead of 3
merged_seasonality, kmeans = perform_seasonality_clustering(
    atv_seasonality, 
    bookings_seasonality, 
    n_clusters=5
)
```

### Access Individual Components
```python
# Get just the Prophet models and forecasts
from prophet_seasonality_clustering import prepare_prophet_data, fit_prophet_model

pandas_data_atv = prepare_prophet_data(your_data, 'ATV')
model_atv, forecast_atv = fit_prophet_model(pandas_data_atv, 'ATV')

# Now you can work with the Prophet model directly
# For example, make longer-term forecasts
future = model_atv.make_future_dataframe(periods=104, freq='W')  # 2 years
forecast = model_atv.predict(future)
```

## Troubleshooting

### Plotly Not Available
The code works with or without plotly. If plotly is not installed, it falls back to matplotlib for plotting.

### Memory Issues
If you have a very large dataset, consider:
1. Reducing the forecast horizon (change `periods=52` in `fit_prophet_model`)
2. Processing data in batches
3. Using only the necessary columns

### Week Number Issues
The code uses ISO week numbers (Monday as start of week). If your data uses a different convention, adjust the date conversion in `prepare_prophet_data`.

## Further Reading

For more details, see:
- `PROPHET_CLUSTERING_README.md` - Comprehensive documentation
- `prophet_seasonality_clustering.py` - Source code with detailed docstrings
