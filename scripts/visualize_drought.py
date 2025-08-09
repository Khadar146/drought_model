#!/usr/bin/env python3
"""
Somaliland Drought Analysis Visualization
Author: Khadar - Enhanced drought monitoring for climate research
Purpose: Visualize SPI time series and drought events for scientific analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import seaborn as sns

# Set scientific plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_spi_data(filepath, var_name='SPI_12'):
    """Load SPI data from NetCDF file"""
    with nc.Dataset(filepath, 'r') as ds:
        spi = ds.variables[var_name][:]
        time = ds.variables['time'][:]
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
        time_units = ds.variables['time'].units
        
    return spi, time, lat, lon, time_units

def create_time_axis(time_vals, time_units):
    """Convert NetCDF time to datetime objects"""
    if 'since' in time_units:
        base_date_str = time_units.split('since ')[-1]
        base_date = datetime.strptime(base_date_str, '%Y-%m-%d')
        
        if 'seconds' in time_units:
            dates = [base_date + timedelta(seconds=int(t)) for t in time_vals]
        elif 'days' in time_units:
            dates = [base_date + timedelta(days=int(t)) for t in time_vals]
        else:
            # Fallback: assume monthly data starting from 1980-01
            dates = [datetime(1980, 1, 1) + timedelta(days=30*i) for i in range(len(time_vals))]
    else:
        dates = [datetime(1980, 1, 1) + timedelta(days=30*i) for i in range(len(time_vals))]
    
    return dates

def compute_regional_average(spi_data):
    """Compute spatial average, handling NaN values"""
    return np.nanmean(spi_data, axis=(0, 1))

def identify_drought_events(spi_ts, threshold=-1.0):
    """Identify drought events from SPI time series"""
    events = []
    in_drought = False
    start_idx = None
    
    for i, spi in enumerate(spi_ts):
        if not np.isnan(spi):
            if spi <= threshold and not in_drought:
                # Start of drought
                in_drought = True
                start_idx = i
            elif spi > threshold and in_drought:
                # End of drought
                in_drought = False
                events.append((start_idx, i-1))
        elif in_drought and np.isnan(spi):
            # End drought if we hit NaN
            in_drought = False
            events.append((start_idx, i-1))
    
    # Handle drought extending to end of series
    if in_drought:
        events.append((start_idx, len(spi_ts)-1))
    
    return events

def plot_drought_analysis():
    """Create comprehensive drought analysis plots"""
    
    # Load SPI data for different accumulation periods
    spi_1, time_vals, lat, lon, time_units = load_spi_data('output/spi/spi_01.nc', 'SPI_01')
    spi_3, _, _, _, _ = load_spi_data('output/spi/spi_03.nc', 'SPI_03')
    spi_12, _, _, _, _ = load_spi_data('output/spi/spi_12.nc', 'SPI_12')
    
    # Create time axis
    dates = create_time_axis(time_vals, time_units)
    years = [d.year + d.month/12.0 for d in dates]
    
    # Compute regional averages
    spi_1_avg = compute_regional_average(spi_1)
    spi_3_avg = compute_regional_average(spi_3)
    spi_12_avg = compute_regional_average(spi_12)
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(4, 1, figsize=(15, 16))
    fig.suptitle('Somaliland Drought Analysis (1980-2024)\nHistorical SPI Analysis for Climate Research', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Multi-period SPI comparison
    ax1 = axes[0]
    ax1.plot(years, spi_1_avg, label='SPI-1 (1-month)', alpha=0.7, linewidth=1)
    ax1.plot(years, spi_3_avg, label='SPI-3 (3-month)', alpha=0.8, linewidth=1.5)
    ax1.plot(years, spi_12_avg, label='SPI-12 (12-month)', alpha=0.9, linewidth=2)
    
    # Add drought threshold lines
    ax1.axhline(y=-1.0, color='orange', linestyle='--', alpha=0.7, label='Moderate Drought')
    ax1.axhline(y=-1.5, color='red', linestyle='--', alpha=0.7, label='Severe Drought')
    ax1.axhline(y=-2.0, color='darkred', linestyle='--', alpha=0.7, label='Extreme Drought')
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    ax1.set_ylabel('SPI Value')
    ax1.set_title('Regional Average SPI - Multiple Accumulation Periods')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-4, 3)
    
    # Plot 2: Focus on SPI-12 with drought events highlighted
    ax2 = axes[1]
    ax2.plot(years, spi_12_avg, color='navy', linewidth=2, label='SPI-12')
    
    # Identify and highlight drought events
    drought_events = identify_drought_events(spi_12_avg, threshold=-1.0)
    for start, end in drought_events:
        ax2.axvspan(years[start], years[end], alpha=0.3, color='red', label='Drought Event' if start == drought_events[0][0] else "")
        
        # Highlight extreme droughts
        if np.min(spi_12_avg[start:end+1]) <= -2.0:
            ax2.axvspan(years[start], years[end], alpha=0.2, color='darkred')
    
    ax2.axhline(y=-1.0, color='orange', linestyle='--', alpha=0.7)
    ax2.axhline(y=-2.0, color='darkred', linestyle='--', alpha=0.7)
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    ax2.set_ylabel('SPI-12 Value')
    ax2.set_title('Annual Drought Index (SPI-12) with Drought Events Highlighted')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-4, 3)
    
    # Plot 3: Drought event statistics
    ax3 = axes[2]
    
    # Calculate decade-wise drought frequency
    decades = np.arange(1980, 2030, 10)
    drought_freq = []
    for decade in decades[:-1]:
        decade_events = [e for e in drought_events 
                        if decade <= years[e[0]] < decade + 10]
        drought_freq.append(len(decade_events))
    
    bars = ax3.bar([f"{d}s" for d in decades[:-1]], drought_freq, 
                   color=['skyblue', 'lightcoral', 'lightgreen', 'orange', 'purple'][:len(drought_freq)])
    ax3.set_ylabel('Number of Drought Events')
    ax3.set_title('Drought Event Frequency by Decade')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, freq in zip(bars, drought_freq):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                str(freq), ha='center', va='bottom', fontweight='bold')
    
    # Plot 4: Running trend analysis
    ax4 = axes[3]
    
    # Calculate 10-year running mean
    window = 120  # 10 years * 12 months
    running_mean = pd.Series(spi_12_avg).rolling(window=window, center=True).mean()
    
    ax4.plot(years, spi_12_avg, color='lightblue', alpha=0.5, linewidth=1, label='Monthly SPI-12')
    ax4.plot(years, running_mean, color='darkblue', linewidth=3, label='10-year Running Mean')
    
    # Add linear trend line
    valid_idx = ~np.isnan(spi_12_avg)
    if np.sum(valid_idx) > 10:
        trend_coef = np.polyfit(np.array(years)[valid_idx], spi_12_avg[valid_idx], 1)
        trend_line = np.polyval(trend_coef, years)
        ax4.plot(years, trend_line, color='red', linestyle='--', linewidth=2, 
                label=f'Linear Trend ({trend_coef[0]:.4f}/year)')
    
    ax4.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax4.set_xlabel('Year')
    ax4.set_ylabel('SPI-12 Value')
    ax4.set_title('Long-term Drought Trends in Somaliland')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('output/spi/somaliland_drought_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('output/spi/somaliland_drought_analysis.pdf', bbox_inches='tight')
    
    print("✓ Drought analysis plots saved:")
    print("  - output/spi/somaliland_drought_analysis.png")
    print("  - output/spi/somaliland_drought_analysis.pdf")
    
    # Print key statistics
    print("\n=== KEY DROUGHT STATISTICS ===")
    print(f"Analysis period: {years[0]:.1f} - {years[-1]:.1f}")
    print(f"Total drought events (SPI-12 <= -1.0): {len(drought_events)}")
    
    extreme_events = [e for e in drought_events if np.min(spi_12_avg[e[0]:e[1]+1]) <= -2.0]
    print(f"Extreme drought events (SPI-12 <= -2.0): {len(extreme_events)}")
    
    if drought_events:
        durations = [(e[1] - e[0] + 1) for e in drought_events]
        print(f"Average drought duration: {np.mean(durations):.1f} months")
        print(f"Longest drought: {np.max(durations)} months")
    
    print(f"Drought frequency: {len(drought_events) / (years[-1] - years[0]) * 10:.1f} events per decade")
    
    if valid_idx.sum() > 0:
        print(f"Mean SPI-12 (1980-2024): {np.nanmean(spi_12_avg):.3f}")
        print(f"SPI-12 standard deviation: {np.nanstd(spi_12_avg):.3f}")
    
    plt.show()

if __name__ == "__main__":
    print("Generating Somaliland Drought Analysis Visualizations...")
    print("=" * 50)
    plot_drought_analysis()
    print("\nVisualization complete! Open the PNG/PDF files to view detailed analysis.")
