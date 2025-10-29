"""
Test script to verify the Prophet seasonality clustering functionality.
This runs in non-interactive mode for testing purposes.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt


def test_prepare_prophet_data():
    """Test data preparation function."""
    from prophet_seasonality_clustering import prepare_prophet_data
    
    # Create sample data
    data = pd.DataFrame({
        'year': [2022, 2022, 2022],
        'week': [1, 2, 3],
        'value': [100, 105, 110]
    })
    
    result = prepare_prophet_data(data, 'value')
    
    assert 'ds' in result.columns, "ds column should be created"
    assert result['ds'].dtype == 'datetime64[ns]', "ds should be datetime"
    print("✓ test_prepare_prophet_data passed")


def test_fit_prophet_model():
    """Test Prophet model fitting."""
    from prophet_seasonality_clustering import fit_prophet_model, prepare_prophet_data
    
    # Generate sample data with trend
    np.random.seed(42)
    n_weeks = 52
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    years = dates.year
    weeks = dates.isocalendar().week
    values = np.linspace(100, 120, n_weeks) + np.random.normal(0, 2, n_weeks)
    
    data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'value': values
    })
    
    prepared_data = prepare_prophet_data(data, 'value')
    model, forecast = fit_prophet_model(prepared_data, 'value')
    
    assert model is not None, "Model should be fitted"
    assert forecast is not None, "Forecast should be generated"
    assert 'yhat' in forecast.columns, "Forecast should contain yhat"
    assert 'weekly' in forecast.columns, "Forecast should contain weekly seasonality"
    assert 'yearly' in forecast.columns, "Forecast should contain yearly seasonality"
    print("✓ test_fit_prophet_model passed")


def test_extract_seasonality_components():
    """Test seasonality extraction."""
    from prophet_seasonality_clustering import (
        extract_seasonality_components, 
        fit_prophet_model, 
        prepare_prophet_data
    )
    
    # Generate sample data
    np.random.seed(42)
    n_weeks = 52
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    years = dates.year
    weeks = dates.isocalendar().week
    values = np.linspace(100, 120, n_weeks) + np.random.normal(0, 2, n_weeks)
    
    data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'value': values
    })
    
    prepared_data = prepare_prophet_data(data, 'value')
    model, forecast = fit_prophet_model(prepared_data, 'value')
    seasonality = extract_seasonality_components(forecast, prepared_data)
    
    assert 'total_seasonality' in seasonality.columns, "Should have total_seasonality"
    assert 'weekly' in seasonality.columns, "Should have weekly component"
    assert 'yearly' in seasonality.columns, "Should have yearly component"
    assert len(seasonality) == len(prepared_data), "Should match original data length"
    print("✓ test_extract_seasonality_components passed")


def test_perform_seasonality_clustering():
    """Test clustering on seasonality components."""
    from prophet_seasonality_clustering import perform_seasonality_clustering
    
    # Create sample seasonality data
    np.random.seed(42)
    n_weeks = 52
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    
    atv_seasonality = pd.DataFrame({
        'ds': dates,
        'weekly': np.sin(2 * np.pi * np.arange(n_weeks) / 52),
        'yearly': np.cos(2 * np.pi * np.arange(n_weeks) / 52),
        'total_seasonality': np.sin(2 * np.pi * np.arange(n_weeks) / 52) + 
                            np.cos(2 * np.pi * np.arange(n_weeks) / 52)
    })
    
    bookings_seasonality = pd.DataFrame({
        'ds': dates,
        'weekly': np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4),
        'yearly': np.cos(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4),
        'total_seasonality': np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4) + 
                            np.cos(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4)
    })
    
    merged_seasonality, kmeans = perform_seasonality_clustering(
        atv_seasonality, 
        bookings_seasonality, 
        n_clusters=3
    )
    
    assert 'cluster' in merged_seasonality.columns, "Should have cluster assignments"
    assert 'ATV_seasonality' in merged_seasonality.columns, "Should have ATV seasonality"
    assert 'bookings_seasonality' in merged_seasonality.columns, "Should have bookings seasonality"
    assert len(merged_seasonality['cluster'].unique()) <= 3, "Should have at most 3 clusters"
    assert kmeans is not None, "KMeans model should be returned"
    print("✓ test_perform_seasonality_clustering passed")


def test_full_workflow():
    """Test the complete workflow with synthetic data."""
    from prophet_seasonality_clustering import main_workflow
    
    # Generate sample data
    np.random.seed(42)
    n_weeks = 52
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    years = dates.year
    weeks = dates.isocalendar().week
    
    # ATV data with seasonality
    trend = np.linspace(100, 120, n_weeks)
    yearly_cycle = 10 * np.sin(2 * np.pi * np.arange(n_weeks) / 52)
    atv = trend + yearly_cycle + np.random.normal(0, 3, n_weeks)
    
    atv_data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'ATV': atv
    })
    
    # Bookings data with different seasonality
    bookings_trend = np.linspace(1000, 950, n_weeks)
    bookings_yearly = 100 * np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4)
    bookings = bookings_trend + bookings_yearly + np.random.normal(0, 30, n_weeks)
    bookings = np.maximum(bookings, 100).astype(int)
    
    bookings_data = pd.DataFrame({
        'year': years,
        'week': weeks,
        'bookings': bookings
    })
    
    # Run workflow
    results = main_workflow(atv_data, bookings_data)
    
    assert 'atv_model' in results, "Should return ATV model"
    assert 'bookings_model' in results, "Should return bookings model"
    assert 'merged_seasonality' in results, "Should return merged seasonality"
    assert 'kmeans' in results, "Should return KMeans model"
    assert 'cluster' in results['merged_seasonality'].columns, "Should have clusters"
    
    print("✓ test_full_workflow passed")
    print("\nFull workflow results summary:")
    print(f"  - Number of data points: {len(results['merged_seasonality'])}")
    print(f"  - Clusters found: {results['merged_seasonality']['cluster'].nunique()}")
    print(f"  - Cluster distribution:")
    for i in range(results['merged_seasonality']['cluster'].nunique()):
        count = (results['merged_seasonality']['cluster'] == i).sum()
        print(f"    Cluster {i}: {count} weeks")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing Prophet Seasonality Clustering Module")
    print("=" * 60)
    print()
    
    tests = [
        test_prepare_prophet_data,
        test_fit_prophet_model,
        test_extract_seasonality_components,
        test_perform_seasonality_clustering,
        test_full_workflow
    ]
    
    for test in tests:
        print(f"Running {test.__name__}...")
        try:
            test()
        except Exception as e:
            print(f"✗ {test.__name__} failed: {e}")
            import traceback
            traceback.print_exc()
            return False
        print()
    
    print("=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
