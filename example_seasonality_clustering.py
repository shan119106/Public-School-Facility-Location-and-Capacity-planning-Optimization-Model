"""
Example script demonstrating Prophet forecasting with seasonality-based clustering.

This script creates sample data to demonstrate the functionality without requiring
a Spark connection or database access.
"""

import pandas as pd
import numpy as np
from prophet_seasonality_clustering import main_workflow


def generate_sample_data(n_weeks=104):
    """
    Generate sample ATV and bookings data with seasonal patterns.
    
    Args:
        n_weeks: Number of weeks to generate
        
    Returns:
        Tuple of (atv_data, bookings_data) DataFrames
    """
    np.random.seed(42)
    
    # Generate date range
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    
    # Extract year and week
    years = dates.year
    weeks = dates.isocalendar().week
    
    # Generate ATV with yearly and weekly seasonality
    # Base trend: increasing over time
    trend = np.linspace(100, 120, n_weeks)
    
    # Yearly seasonality (peak in summer, low in winter)
    yearly_cycle = 10 * np.sin(2 * np.pi * np.arange(n_weeks) / 52)
    
    # Weekly seasonality (within year patterns)
    weekly_cycle = 5 * np.sin(2 * np.pi * np.arange(n_weeks) / 4)
    
    # Add noise
    noise_atv = np.random.normal(0, 3, n_weeks)
    
    # Combine components
    atv = trend + yearly_cycle + weekly_cycle + noise_atv
    
    # Generate bookings with different seasonal pattern
    # Base trend: slightly decreasing over time
    bookings_trend = np.linspace(1000, 950, n_weeks)
    
    # Yearly seasonality (different phase than ATV)
    bookings_yearly = 100 * np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4)
    
    # Weekly seasonality
    bookings_weekly = 50 * np.sin(2 * np.pi * np.arange(n_weeks) / 4 + np.pi/6)
    
    # Add noise
    noise_bookings = np.random.normal(0, 30, n_weeks)
    
    # Combine components
    bookings = bookings_trend + bookings_yearly + bookings_weekly + noise_bookings
    
    # Ensure bookings are positive integers
    bookings = np.maximum(bookings, 100).astype(int)
    
    # Create DataFrames
    atv_data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'ATV': atv
    })
    
    bookings_data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'bookings': bookings
    })
    
    return atv_data, bookings_data


def main():
    """
    Main function to demonstrate the seasonality-based clustering workflow.
    """
    print("Generating sample data...")
    atv_data, bookings_data = generate_sample_data(n_weeks=104)  # 2 years of data
    
    print(f"Generated {len(atv_data)} weeks of data")
    print("\nSample ATV data:")
    print(atv_data.head())
    print("\nSample bookings data:")
    print(bookings_data.head())
    
    print("\n" + "=" * 60)
    print("Starting Prophet Forecasting and Seasonality Clustering")
    print("=" * 60)
    
    # Run the main workflow
    results = main_workflow(atv_data, bookings_data)
    
    print("\n" + "=" * 60)
    print("Results Summary")
    print("=" * 60)
    
    # Display sample of merged seasonality data with clusters
    print("\nSample of seasonality data with cluster assignments:")
    print(results['merged_seasonality'][
        ['ds', 'ATV_seasonality', 'bookings_seasonality', 'cluster']
    ].head(10))
    
    # Cluster distribution
    print("\nCluster distribution:")
    print(results['merged_seasonality']['cluster'].value_counts().sort_index())
    
    print("\nDone! Check the generated plots for visual analysis.")


if __name__ == "__main__":
    main()
