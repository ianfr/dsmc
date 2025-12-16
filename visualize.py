#!/usr/bin/env python3
"""
Interactive 3D visualization of DSMC particle data using PyVista.

Usage: 
    python visualize.py <output_folder> [timestep]         # Single timestep
    python visualize.py <output_folder> --loop             # Loop all timesteps
    python visualize.py <output_folder> --loop --fps 10    # Loop at 10 FPS
    
Examples:
    python visualize.py DSMC_OUT 50           # Visualize timestep 50
    python visualize.py build/DSMC_OUT        # Visualize timestep 0
    python visualize.py DSMC_OUT --loop       # Animate all timesteps
    python visualize.py DSMC_OUT --loop --fps 10        Animate at 10 FPS
    python visualize.py DSMC_OUT --vmin 0 --vmax 1000   Set color scale from 0 to 1000 m/s
    python visualize.py DSMC_OUT 5 --clip-axis z --clip-start 0.4       # Clip Z axis from 40%% to 100%%
    python visualize.py DSMC_OUT 5 --clip-axis x --clip-start 0.25 --clip-end 0.75  # Clip X axis middle 50%%
"""

import sys
import argparse
import glob
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import pyvista as pv


def load_particle_data(filepath: str) -> tuple[np.ndarray, np.ndarray]:
    """Load particle positions and velocities from CSV file using pandas."""
    df = pd.read_csv(filepath)
    df.columns = df.columns.str.strip()  # Remove whitespace from column names
    df = df.fillna(0.0)  # Replace NaN with zeros
    
    positions = df[['x', 'y', 'z']].values
    speeds = df['velmag'].values
    
    return positions, speeds


def get_available_timesteps(data_dir: str, include_initial: bool = False) -> list[tuple[int, str]]:
    """Get list of available timesteps from output directory.
    
    Args:
        data_dir: Path to directory containing CSV files
        include_initial: If True, include afterCreate file as timestep -1
    
    Returns:
        List of (timestep_number, filepath) tuples sorted by timestep.
    """
    files = glob.glob(f"{data_dir}/*-particle.csv")
    timesteps = []
    
    for f in files:
        basename = os.path.basename(f)
        if basename.startswith('afterCreate'):
            if include_initial:
                timesteps.append((-1, f))
            continue
        try:
            ts = int(basename.split('-')[0])
            timesteps.append((ts, f))
        except ValueError:
            continue
    
    return sorted(timesteps)


def create_plotter():
    """Create and configure a PyVista plotter with consistent styling."""
    plotter = pv.Plotter()
    plotter.set_background('#2b2b2b')  # Dark gray background
    return plotter


def clip_particles(positions: np.ndarray, speeds: np.ndarray, domain_bounds: tuple, 
                   clip_axis: str = None, clip_start: float = 0.0, clip_end: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Clip particles to a box along specified axis.
    
    Args:
        positions: Nx3 array of particle positions
        speeds: N array of particle speeds
        domain_bounds: Tuple of (x_min, x_max, y_min, y_max, z_min, z_max)
        clip_axis: Axis to clip along ('x', 'y', or 'z'), None for no clipping
        clip_start: Start position as proportion of domain length [0-1]
        clip_end: End position as proportion of domain length [0-1]
    
    Returns:
        Tuple of (filtered_positions, filtered_speeds)
    """
    if clip_axis is None:
        return positions, speeds
    
    # Map axis to index
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    if clip_axis.lower() not in axis_map:
        print(f"Warning: Invalid clip axis '{clip_axis}', ignoring clipping")
        return positions, speeds
    
    axis_idx = axis_map[clip_axis.lower()]
    x_min, x_max, y_min, y_max, z_min, z_max = domain_bounds
    
    # Get domain bounds for the specified axis
    if axis_idx == 0:
        domain_min, domain_max = x_min, x_max
    elif axis_idx == 1:
        domain_min, domain_max = y_min, y_max
    else:
        domain_min, domain_max = z_min, z_max
    
    # Convert proportions to absolute positions
    domain_length = domain_max - domain_min
    abs_start = domain_min + clip_start * domain_length
    abs_end = domain_min + clip_end * domain_length
    
    # Filter particles
    axis_positions = positions[:, axis_idx]
    mask = (axis_positions >= abs_start) & (axis_positions <= abs_end)
    
    return positions[mask], speeds[mask]


def get_domain_bounds(data_dir: str) -> tuple[float, float, float, float, float, float]:
    """Get domain bounds from afterCreate-particle.csv file.
    
    Returns:
        Tuple of (x_min, x_max, y_min, y_max, z_min, z_max)
    """
    # Look for afterCreate file
    aftercreate_files = glob.glob(f"{data_dir}/afterCreate-particle.csv")
    
    if not aftercreate_files:
        print("Warning: afterCreate-particle.csv not found, using default bounds")
        return (0, 1, 0, 1, 0, 1)
    
    # Load the initial particle positions
    positions, _ = load_particle_data(aftercreate_files[0])
    
    # Compute bounds
    x_min, y_min, z_min = positions.min(axis=0)
    x_max, y_max, z_max = positions.max(axis=0)
    
    return (x_min, x_max, y_min, y_max, z_min, z_max)


def visualize_timestep(timestep: int, data_dir: str, plotter=None, domain_bounds=None, vmin=None, vmax=None,
                       clip_axis=None, clip_start=0.0, clip_end=1.0):
    """Create interactive 3D visualization of particles at given timestep.
    
    Args:
        timestep: Timestep number to visualize (-1 for initial state)
        data_dir: Directory containing particle data
        plotter: Existing PyVista plotter to reuse (optional)
        domain_bounds: Tuple of (x_min, x_max, y_min, y_max, z_min, z_max)
        vmin: Minimum velocity for color scale (optional)
        vmax: Maximum velocity for color scale (optional)
        clip_axis: Axis to clip along ('x', 'y', or 'z'), None for no clipping
        clip_start: Start position as proportion of domain length [0-1]
        clip_end: End position as proportion of domain length [0-1]
    """
    
    # Build filename
    if timestep == -1:
        # Use afterCreate file
        aftercreate_files = glob.glob(f"{data_dir}/afterCreate-particle.csv")
        if not aftercreate_files:
            print("Error: afterCreate-particle.csv not found")
            return None
        filename = aftercreate_files[0]
    else:
        filename = f"{data_dir}/{timestep:010d}-particle.csv"
    
    print(f"Loading: {filename}")
    
    try:
        positions, speeds = load_particle_data(filename)
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        available = get_available_timesteps(data_dir, include_initial=True)
        if available:
            print(f"Available timesteps: {available[0][0]} - {available[-1][0]}")
        return None
    
    # Get domain bounds if not provided
    if domain_bounds is None:
        domain_bounds = get_domain_bounds(data_dir)
    
    # Apply clipping if specified
    if clip_axis is not None:
        positions, speeds = clip_particles(positions, speeds, domain_bounds, clip_axis, clip_start, clip_end)
    
    n_particles = len(positions)
    print(f"  Particles: {n_particles:,}")
    print(f"  Speed range: {speeds.min():.2e} - {speeds.max():.2e}")
    
    # Create point cloud
    cloud = pv.PolyData(positions)
    cloud['speed'] = speeds
    
    # Create or use existing plotter
    if plotter is None:
        plotter = create_plotter()
        show_after = True
    else:
        plotter.clear()
        show_after = False
    
    # Add mesh with turbo colormap (works well on dark backgrounds)
    mesh_kwargs = {
        'scalars': 'speed',
        'cmap': 'turbo',
        'point_size': 3,
        'render_points_as_spheres': True,
        'scalar_bar_args': {
            'title': 'Speed (m/s)',
            'title_font_size': 16,
            'label_font_size': 14,
            'color': 'white',
            'font_family': 'arial'
        }
    }
    
    # Apply custom color limits if provided
    if vmin is not None or vmax is not None:
        clim = [vmin if vmin is not None else speeds.min(),
                vmax if vmax is not None else speeds.max()]
        mesh_kwargs['clim'] = clim
    
    plotter.add_mesh(cloud, **mesh_kwargs)
    
    # Add domain bounding box
    x_min, x_max, y_min, y_max, z_min, z_max = domain_bounds
    box = pv.Box(bounds=(x_min, x_max, y_min, y_max, z_min, z_max))
    plotter.add_mesh(box, style='wireframe', color='white', line_width=2, opacity=1.0)
    
    # Add axes
    plotter.add_axes(
        line_width=4,
        color='white',
        x_color='red',
        y_color='green',
        z_color='blue',
        xlabel='X',
        ylabel='Y',
        zlabel='Z'
    )
    
    # Add title
    title = f"DSMC Particles - Initial State" if timestep == -1 else f"DSMC Particles - Timestep {timestep}"
    plotter.add_title(
        title,
        font_size=14,
        color='white',
        font='arial'
    )
    
    # Add text annotation with statistics
    stats_text = f"Particles: {n_particles:,}\nSpeed: {speeds.min():.2e} - {speeds.max():.2e}"
    plotter.add_text(
        stats_text,
        position='lower_left',
        font_size=10,
        color='white',
        font='arial'
    )
    
    if show_after:
        print("\nControls:")
        print("  Left-click + drag: Rotate")
        print("  Middle-click + drag: Pan")
        print("  Scroll: Zoom")
        print("  'r': Reset camera")
        print("  'q': Quit")
        plotter.show()
    
    return plotter


def animate_timesteps(data_dir: str, fps: int = 5, vmin=None, vmax=None,
                      clip_axis=None, clip_start=0.0, clip_end=1.0):
    """Animate all available timesteps in a loop.
    
    Args:
        data_dir: Directory containing particle data
        fps: Frames per second for animation
        vmin: Minimum velocity for color scale (optional)
        vmax: Maximum velocity for color scale (optional)
        clip_axis: Axis to clip along ('x', 'y', or 'z'), None for no clipping
        clip_start: Start position as proportion of domain length [0-1]
        clip_end: End position as proportion of domain length [0-1]
    """
    
    timesteps = get_available_timesteps(data_dir, include_initial=True)
    
    if not timesteps:
        print(f"Error: No particle data found in {data_dir}")
        return
    
    print(f"Found {len(timesteps)} timesteps")
    print(f"Animating at {fps} FPS")
    print("\nControls:")
    print("  Left-click + drag: Rotate")
    print("  Middle-click + drag: Pan")
    print("  Scroll: Zoom")
    print("  'r': Reset camera")
    print("  'q': Quit")
    print()
    
    # Get domain bounds once from afterCreate file
    domain_bounds = get_domain_bounds(data_dir)
    
    plotter = create_plotter()
    plotter.show(interactive_update=True, auto_close=False)
    
    frame_delay = 1.0 / fps
    
    try:
        while True:
            for ts, filepath in timesteps:
                start_time = time.time()
                
                # Load data
                positions, speeds = load_particle_data(filepath)
                
                # Apply clipping if specified
                if clip_axis is not None:
                    positions, speeds = clip_particles(positions, speeds, domain_bounds, clip_axis, clip_start, clip_end)
                
                # Clear and update
                plotter.clear()
                
                # Create point cloud
                cloud = pv.PolyData(positions)
                cloud['speed'] = speeds
                
                # Add mesh
                mesh_kwargs = {
                    'scalars': 'speed',
                    'cmap': 'turbo',
                    'point_size': 3,
                    'render_points_as_spheres': True,
                    'scalar_bar_args': {
                        'title': 'Speed (m/s)',
                        'title_font_size': 16,
                        'label_font_size': 14,
                        'color': 'white',
                        'font_family': 'arial'
                    }
                }
                
                # Apply custom color limits if provided
                if vmin is not None or vmax is not None:
                    clim = [vmin if vmin is not None else speeds.min(),
                            vmax if vmax is not None else speeds.max()]
                    mesh_kwargs['clim'] = clim
                
                plotter.add_mesh(cloud, **mesh_kwargs)
                
                # Add domain bounding box
                x_min, x_max, y_min, y_max, z_min, z_max = domain_bounds
                box = pv.Box(bounds=(x_min, x_max, y_min, y_max, z_min, z_max))
                plotter.add_mesh(box, style='wireframe', color='white', line_width=2, opacity=1.0)
                
                # Add axes
                plotter.add_axes(
                    line_width=4,
                    color='white',
                    x_color='red',
                    y_color='green',
                    z_color='blue',
                    xlabel='X',
                    ylabel='Y',
                    zlabel='Z'
                )
                
                # Add title
                title = f"DSMC Particles - Initial State" if ts == -1 else f"DSMC Particles - Timestep {ts}"
                plotter.add_title(
                    title,
                    font_size=14,
                    color='white',
                    font='arial'
                )
                
                # Add stats
                n_particles = len(positions)
                stats_text = f"Particles: {n_particles:,}\nSpeed: {speeds.min():.2e} - {speeds.max():.2e}"
                plotter.add_text(
                    stats_text,
                    position='lower_left',
                    font_size=10,
                    color='white',
                    font='arial'
                )
                
                # Update display
                plotter.update()
                
                # Sleep to maintain FPS
                elapsed = time.time() - start_time
                sleep_time = max(0, frame_delay - elapsed)
                if sleep_time > 0:
                    time.sleep(sleep_time)
                
                # Check if window was closed
                if not plotter.render_window:
                    return
                    
    except KeyboardInterrupt:
        print("\nAnimation stopped by user")
    finally:
        plotter.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize DSMC particle simulation data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python visualize.py DSMC_OUT 50                                   # Visualize timestep 50
    python visualize.py build/DSMC_OUT 0                                # Visualize timestep 0
    python visualize.py DSMC_OUT --loop                               # Animate all timesteps
    python visualize.py DSMC_OUT --loop --fps 10                      # Animate at 10 FPS
    python visualize.py DSMC_OUT -1 --vmin 0 --vmax 1000                 # Set color scale from 0 to 1000 m/s
    python visualize.py DSMC_OUT 5 --clip-axis z --clip-start 0.4       # Clip Z axis from 40%% to 100%%
    python visualize.py DSMC_OUT 5 --clip-axis x --clip-start 0.25 --clip-end 0.75  # Clip X axis middle 50%%
"""
    )
    
    parser.add_argument(
        'data_dir',
        help='Path to directory containing particle CSV files'
    )
    
    parser.add_argument(
        'timestep',
        type=int,
        nargs='?',
        default=0,
        help='Timestep to visualize (default: 0, ignored if --loop is used)'
    )
    
    parser.add_argument(
        '--loop',
        action='store_true',
        help='Animate all timesteps in a loop'
    )
    
    parser.add_argument(
        '--fps',
        type=int,
        default=5,
        help='Frames per second for animation (default: 5)'
    )
    
    parser.add_argument(
        '--vmin',
        type=float,
        default=None,
        help='Minimum velocity for color scale (default: auto from data)'
    )
    
    parser.add_argument(
        '--vmax',
        type=float,
        default=None,
        help='Maximum velocity for color scale (default: auto from data)'
    )
    
    parser.add_argument(
        '--clip-axis',
        type=str,
        choices=['x', 'y', 'z', 'X', 'Y', 'Z'],
        default=None,
        help='Axis to clip along (x, y, or z)'
    )
    
    parser.add_argument(
        '--clip-start',
        type=float,
        default=0.0,
        help='Start position for clipping as proportion of domain [0-1] (default: 0.0)'
    )
    
    parser.add_argument(
        '--clip-end',
        type=float,
        default=1.0,
        help='End position for clipping as proportion of domain [0-1] (default: 1.0)'
    )
    
    args = parser.parse_args()
    
    # Validate data directory
    if not os.path.isdir(args.data_dir):
        print(f"Error: Directory not found: {args.data_dir}")
        sys.exit(1)
    
    # Run visualization
    if args.loop:
        animate_timesteps(args.data_dir, args.fps, vmin=args.vmin, vmax=args.vmax,
                         clip_axis=args.clip_axis, clip_start=args.clip_start, clip_end=args.clip_end)
    else:
        visualize_timestep(args.timestep, args.data_dir, vmin=args.vmin, vmax=args.vmax,
                          clip_axis=args.clip_axis, clip_start=args.clip_start, clip_end=args.clip_end)


if __name__ == "__main__":
    main()
