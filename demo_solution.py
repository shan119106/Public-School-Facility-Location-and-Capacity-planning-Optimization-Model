"""
Simple demonstration showing the exact solution to the problem statement.

This script shows how to cluster on seasonality of ATV and seasonality of bookings
from Prophet models, NOT on raw ATV and bookings values.
"""

import pandas as pd
import numpy as np
from prophet_seasonality_clustering import main_workflow

# Generate simple sample data
np.random.seed(42)
n_weeks = 52
dates = pd.date_range(start='2023-01-02', periods=n_weeks, freq='W-MON')

# Create ATV data with seasonal pattern
atv_data = pd.DataFrame({
    'year': dates.year,
    'week': dates.isocalendar().week,
    'ATV': 100 + 10 * np.sin(2 * np.pi * np.arange(n_weeks) / 52) + np.random.normal(0, 2, n_weeks)
})

# Create bookings data with seasonal pattern
bookings_data = pd.DataFrame({
    'year': dates.year,
    'week': dates.isocalendar().week,
    'bookings': (1000 + 100 * np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4) + 
                np.random.normal(0, 20, n_weeks)).astype(int)
})

print("=" * 80)
print("DEMONSTRATION: Clustering on Prophet Seasonality Components")
print("=" * 80)
print()
print("Problem Statement:")
print("  'clustering is happening for the ATV and bookings but I want it for")
print("   the seasonality of ATV and seasonality of bookings from the prophet model'")
print()
print("Solution:")
print("  1. Fit Prophet models to ATV and bookings data")
print("  2. Extract seasonality components from Prophet forecasts")
print("  3. Cluster on seasonality components (NOT raw values)")
print()
print("=" * 80)

# Run the main workflow
results = main_workflow(atv_data, bookings_data)

# Show the key difference
print()
print("=" * 80)
print("KEY OUTPUTS")
print("=" * 80)
print()
print("Seasonality data with cluster assignments:")
print("-" * 80)
output_df = results['merged_seasonality'][
    ['ds', 'ATV_seasonality', 'bookings_seasonality', 'cluster']
].head(10)
print(output_df.to_string(index=False))
print()
print("Cluster summary:")
print("-" * 80)
for cluster in sorted(results['merged_seasonality']['cluster'].unique()):
    cluster_data = results['merged_seasonality'][
        results['merged_seasonality']['cluster'] == cluster
    ]
    print(f"Cluster {cluster}: {len(cluster_data)} weeks")
    print(f"  ATV seasonality mean: {cluster_data['ATV_seasonality'].mean():.4f}")
    print(f"  Bookings seasonality mean: {cluster_data['bookings_seasonality'].mean():.4f}")
    print()

print("=" * 80)
print("SUCCESS: Clustering is now based on seasonality components!")
print("=" * 80)
