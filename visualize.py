#!/usr/bin/env python3
"""
Interactive 3D visualization of DSMC particle data using PyVista.
Usage: python visualize.py [timestep]
       python visualize.py              # defaults to timestep 0
       python visualize.py 10           # visualize timestep 10
"""

import sys
import glob
import numpy as np
import pyvista as pv


def load_particle_data(filepath: str) -> tuple[np.ndarray, np.ndarray]:
    """Load particle positions and velocities from CSV file."""
    data = np.loadtxt(filepath, delimiter=',', skiprows=1)
    data = np.nan_to_num(data, nan=0.0)  # Replace NaN with zeros
    positions = data[:, :3]  # x, y, z
    speeds = data[:, 3]      # velmag
    return positions, speeds


def get_available_timesteps(data_dir: str = "DSMC_OUT") -> list[int]:
    """Get list of available timesteps from output directory."""
    files = glob.glob(f"{data_dir}/*-particle.csv")
    timesteps = []
    for f in files:
        basename = f.split('/')[-1]
        if basename.startswith('afterCreate'):
            continue
        try:
            ts = int(basename.split('-')[0])
            timesteps.append(ts)
        except ValueError:
            continue
    return sorted(timesteps)


def visualize_timestep(timestep: int, data_dir: str = "DSMC_OUT"):
    """Create interactive 3D visualization of particles at given timestep."""
    
    # Build filename
    filename = f"{data_dir}/{timestep:010d}-particle.csv"
    
    print(f"Loading particle data from: {filename}")
    
    try:
        positions, speeds = load_particle_data(filename)
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        available = get_available_timesteps(data_dir)
        if available:
            print(f"Available timesteps: {available[0]} - {available[-1]}")
        return
    
    n_particles = len(positions)
    print(f"Loaded {n_particles:,} particles")
    print(f"Speed range: {speeds.min():.2e} - {speeds.max():.2e}")
    print(f"Position bounds:")
    print(f"  X: [{positions[:,0].min():.2e}, {positions[:,0].max():.2e}]")
    print(f"  Y: [{positions[:,1].min():.2e}, {positions[:,1].max():.2e}]")
    print(f"  Z: [{positions[:,2].min():.2e}, {positions[:,2].max():.2e}]")
    
    # Create point cloud
    cloud = pv.PolyData(positions)
    cloud['speed'] = speeds
    
    # Create plotter
    plotter = pv.Plotter()
    plotter.add_mesh(
        cloud,
        scalars='speed',
        cmap='plasma',
        point_size=2,
        render_points_as_spheres=False,
        scalar_bar_args={'title': 'Speed (m/s)'}
    )
    
    # Add bounding box
    bounds = cloud.bounds
    plotter.add_bounding_box(color='white', line_width=1)
    
    # Set up camera and display
    plotter.add_axes()
    plotter.set_background('black')
    plotter.add_title(f"DSMC Particles - Timestep {timestep}", font_size=12)
    
    print("\nControls:")
    print("  Left-click + drag: Rotate")
    print("  Middle-click + drag: Pan")
    print("  Scroll: Zoom")
    print("  'r': Reset camera")
    print("  'q': Quit")
    
    plotter.show()


def main():
    # Parse command line arguments
    if len(sys.argv) > 1:
        try:
            timestep = int(sys.argv[1])
        except ValueError:
            print(f"Error: Invalid timestep '{sys.argv[1]}'. Must be an integer.")
            sys.exit(1)
    else:
        timestep = 0
    
    visualize_timestep(timestep)


if __name__ == "__main__":
    main()
