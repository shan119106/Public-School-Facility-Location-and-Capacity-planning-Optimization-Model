"""
Comparison script: Raw Values vs. Seasonality-Based Clustering

This script demonstrates the difference between clustering on raw ATV/bookings
values versus clustering on seasonality components extracted from Prophet models.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from prophet_seasonality_clustering import (
    prepare_prophet_data,
    fit_prophet_model,
    extract_seasonality_components
)


def generate_sample_data_with_trend(n_weeks=104):
    """
    Generate sample data with strong trend and seasonal patterns.
    This makes the difference between raw and seasonality clustering clear.
    """
    np.random.seed(42)
    dates = pd.date_range(start='2022-01-03', periods=n_weeks, freq='W-MON')
    years = dates.year
    weeks = dates.isocalendar().week
    
    # ATV with strong upward trend
    trend = np.linspace(100, 200, n_weeks)  # Strong upward trend
    yearly_cycle = 15 * np.sin(2 * np.pi * np.arange(n_weeks) / 52)
    weekly_cycle = 5 * np.sin(2 * np.pi * np.arange(n_weeks) / 4)
    atv = trend + yearly_cycle + weekly_cycle + np.random.normal(0, 3, n_weeks)
    
    # Bookings with downward trend but similar seasonality
    bookings_trend = np.linspace(1200, 800, n_weeks)  # Downward trend
    bookings_yearly = 120 * np.sin(2 * np.pi * np.arange(n_weeks) / 52 + np.pi/4)
    bookings_weekly = 50 * np.sin(2 * np.pi * np.arange(n_weeks) / 4 + np.pi/6)
    bookings = bookings_trend + bookings_yearly + bookings_weekly + np.random.normal(0, 30, n_weeks)
    bookings = np.maximum(bookings, 100).astype(int)
    
    return pd.DataFrame({
        'year': years,
        'week': weeks,
        'ds': dates,
        'ATV': atv,
        'bookings': bookings
    })


def cluster_raw_values(data, n_clusters=3):
    """
    Cluster based on raw ATV and bookings values (original approach).
    """
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(data[['ATV', 'bookings']])
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    data['cluster_raw'] = kmeans.fit_predict(scaled_features)
    
    return data, kmeans


def cluster_seasonality(data, n_clusters=3):
    """
    Cluster based on seasonality components from Prophet (new approach).
    """
    # Fit Prophet models
    atv_prepared = prepare_prophet_data(data.copy(), 'ATV')
    model_atv, forecast_atv = fit_prophet_model(atv_prepared, 'ATV')
    atv_seasonality = extract_seasonality_components(forecast_atv, atv_prepared)
    
    bookings_prepared = prepare_prophet_data(data.copy(), 'bookings')
    model_bookings, forecast_bookings = fit_prophet_model(bookings_prepared, 'bookings')
    bookings_seasonality = extract_seasonality_components(forecast_bookings, bookings_prepared)
    
    # Merge seasonality
    merged = pd.merge(
        atv_seasonality[['ds', 'total_seasonality']].rename(
            columns={'total_seasonality': 'ATV_seasonality'}
        ),
        bookings_seasonality[['ds', 'total_seasonality']].rename(
            columns={'total_seasonality': 'bookings_seasonality'}
        ),
        on='ds'
    )
    
    # Cluster on seasonality
    scaler = StandardScaler()
    scaled_seasonality = scaler.fit_transform(
        merged[['ATV_seasonality', 'bookings_seasonality']]
    )
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    merged['cluster_seasonality'] = kmeans.fit_predict(scaled_seasonality)
    
    # Merge back to original data
    data = data.merge(
        merged[['ds', 'ATV_seasonality', 'bookings_seasonality', 'cluster_seasonality']], 
        on='ds'
    )
    
    return data, kmeans


def visualize_comparison(data, output_file='cluster_comparison.png'):
    """
    Create side-by-side comparison of raw vs seasonality clustering.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Raw values clustering
    ax1 = axes[0, 0]
    for cluster in data['cluster_raw'].unique():
        cluster_data = data[data['cluster_raw'] == cluster]
        ax1.scatter(
            cluster_data['ATV'], 
            cluster_data['bookings'],
            label=f'Cluster {cluster}',
            alpha=0.6,
            s=50
        )
    ax1.set_xlabel('ATV (Raw)', fontsize=12)
    ax1.set_ylabel('Bookings (Raw)', fontsize=12)
    ax1.set_title('Clustering on Raw Values\n(Influenced by Trend)', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Time series of raw values with clusters
    ax2 = axes[0, 1]
    for cluster in data['cluster_raw'].unique():
        cluster_data = data[data['cluster_raw'] == cluster]
        ax2.scatter(
            cluster_data['ds'], 
            cluster_data['ATV'],
            label=f'Cluster {cluster}',
            alpha=0.6,
            s=30
        )
    ax2.set_xlabel('Date', fontsize=12)
    ax2.set_ylabel('ATV', fontsize=12)
    ax2.set_title('Raw Values Over Time\n(Clusters follow trend)', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='x', rotation=45)
    
    # Seasonality clustering
    ax3 = axes[1, 0]
    for cluster in data['cluster_seasonality'].unique():
        cluster_data = data[data['cluster_seasonality'] == cluster]
        ax3.scatter(
            cluster_data['ATV_seasonality'], 
            cluster_data['bookings_seasonality'],
            label=f'Cluster {cluster}',
            alpha=0.6,
            s=50
        )
    ax3.set_xlabel('ATV Seasonality', fontsize=12)
    ax3.set_ylabel('Bookings Seasonality', fontsize=12)
    ax3.set_title('Clustering on Seasonality Components\n(Trend Removed)', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Time series of seasonality with clusters
    ax4 = axes[1, 1]
    for cluster in data['cluster_seasonality'].unique():
        cluster_data = data[data['cluster_seasonality'] == cluster]
        ax4.scatter(
            cluster_data['ds'], 
            cluster_data['ATV_seasonality'],
            label=f'Cluster {cluster}',
            alpha=0.6,
            s=30
        )
    ax4.set_xlabel('Date', fontsize=12)
    ax4.set_ylabel('ATV Seasonality', fontsize=12)
    ax4.set_title('Seasonality Over Time\n(Clusters capture patterns)', fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Comparison plot saved to {output_file}")
    return fig


def print_cluster_analysis(data):
    """
    Print detailed comparison of clustering results.
    """
    print("\n" + "=" * 80)
    print("CLUSTERING COMPARISON ANALYSIS")
    print("=" * 80)
    
    print("\n### RAW VALUES CLUSTERING ###")
    print("-" * 80)
    for cluster in sorted(data['cluster_raw'].unique()):
        cluster_data = data[data['cluster_raw'] == cluster]
        print(f"\nRaw Cluster {cluster}:")
        print(f"  Size: {len(cluster_data)} weeks")
        print(f"  ATV Range: {cluster_data['ATV'].min():.2f} - {cluster_data['ATV'].max():.2f}")
        print(f"  ATV Mean: {cluster_data['ATV'].mean():.2f} ± {cluster_data['ATV'].std():.2f}")
        print(f"  Bookings Range: {cluster_data['bookings'].min():.0f} - {cluster_data['bookings'].max():.0f}")
        print(f"  Bookings Mean: {cluster_data['bookings'].mean():.0f} ± {cluster_data['bookings'].std():.0f}")
        print(f"  Date Range: {cluster_data['ds'].min().date()} to {cluster_data['ds'].max().date()}")
    
    print("\n### SEASONALITY-BASED CLUSTERING ###")
    print("-" * 80)
    for cluster in sorted(data['cluster_seasonality'].unique()):
        cluster_data = data[data['cluster_seasonality'] == cluster]
        print(f"\nSeasonality Cluster {cluster}:")
        print(f"  Size: {len(cluster_data)} weeks")
        print(f"  ATV Seasonality Mean: {cluster_data['ATV_seasonality'].mean():.4f} ± {cluster_data['ATV_seasonality'].std():.4f}")
        print(f"  Bookings Seasonality Mean: {cluster_data['bookings_seasonality'].mean():.4f} ± {cluster_data['bookings_seasonality'].std():.4f}")
        print(f"  Raw ATV Mean: {cluster_data['ATV'].mean():.2f} ± {cluster_data['ATV'].std():.2f}")
        print(f"  Raw Bookings Mean: {cluster_data['bookings'].mean():.0f} ± {cluster_data['bookings'].std():.0f}")
        
        # Show date distribution to demonstrate seasonal patterns
        months = cluster_data['ds'].dt.month.value_counts().sort_index()
        print(f"  Month Distribution: {dict(months)}")
    
    print("\n" + "=" * 80)
    print("KEY INSIGHT:")
    print("=" * 80)
    print("Raw clustering groups weeks by absolute values (early/mid/late periods).")
    print("Seasonality clustering groups weeks by similar seasonal patterns")
    print("(e.g., holiday periods, summer/winter patterns) regardless of when they occur.")
    print("=" * 80)


def main():
    """
    Main function to run the comparison.
    """
    print("=" * 80)
    print("COMPARING RAW VALUES vs SEASONALITY-BASED CLUSTERING")
    print("=" * 80)
    
    # Generate data
    print("\nGenerating sample data with strong trends...")
    data = generate_sample_data_with_trend(n_weeks=104)
    print(f"Generated {len(data)} weeks of data (2 years)")
    
    # Cluster on raw values
    print("\nPerforming clustering on raw values...")
    data, kmeans_raw = cluster_raw_values(data.copy(), n_clusters=3)
    print("Raw values clustering complete.")
    
    # Cluster on seasonality
    print("\nPerforming clustering on seasonality components...")
    data, kmeans_seasonality = cluster_seasonality(data, n_clusters=3)
    print("Seasonality clustering complete.")
    
    # Print analysis
    print_cluster_analysis(data)
    
    # Visualize
    print("\nGenerating comparison visualization...")
    fig = visualize_comparison(data)
    
    print("\n" + "=" * 80)
    print("COMPARISON COMPLETE!")
    print("=" * 80)
    print("\nCheck cluster_comparison.png for visual comparison.")


if __name__ == "__main__":
    main()
