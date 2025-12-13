# DSMC (Direct Simulation Monte Carlo)

This project uses Metal **GPU acceleration** on Apple Silicon for parallelized a DSMC code for dilute gas. Metal acceleration structures are used for embedding arbitrary 3D objects in the computational domain. 

## Example

Heat diffusion through Argon gas at STP when the bottom boundary is hot:
* There are 10^5 simulated particles, each representing ~270 million Argon atoms
* The side length of the cubic domain is 100 microns
* The computational grid is 2x2
* Reflection off of the domain boundaries is diffuse

Initial snapshot  |  . |  .| .
:-------------------------:|:-------------------------:|:-----------:|---------------:
![image](https://user-images.githubusercontent.com/49919175/221903145-bac8047f-21c9-45de-88ff-ddcd78b39015.png)  |  ![image](https://user-images.githubusercontent.com/49919175/221903182-cca674c3-abd5-4006-af05-1fca7cdef571.png) | ![image](https://user-images.githubusercontent.com/49919175/221903233-a3799aa3-87a3-49f3-87b7-924b0ec8faf9.png) | ![image](https://user-images.githubusercontent.com/49919175/221903292-aac2b597-c7b0-4e08-b09b-366759c60e35.png)


### Features
- **GPU-accelerated particle updates**: Position integration runs on Metal compute shaders
- **GPU-accelerated collisions**: Rejection sampling collision detection parallelized per cell
- **Spatial hashing on GPU**: Particles sorted into cells for efficient neighbor lookups
- **Ray-traced mesh collisions**: Embed arbitrary 3D objects using Metal Acceleration Structures
- **Diffuse reflection**: Particles reflect off domain boundaries and mesh surfaces

### Building (CMake & VS Code)

```bash
cd build
cmake ..
make
./dsmc
```

### Building (XCode)

```bash
cd build-xcode
cmake -G Xcode ..
```

Note that if using XCode, 'config.json' has to be in the build-xcode folder and not the top-level repo folder, unlike with CMake.

Now you can open build-xcode/dsmc.xcodeproj like a normal project, using scheme "dsmc" for building.

For best performance, select Product -> Build For -> Profiling which will generate a Release folder inside build-xcode.

### Configuration

Edit `config.json` to configure the simulation.

#### Mesh Options
- `type`: "none", "sphere", or "box" for built-in primitives
- `size`: Size of the primitive (radius for sphere, side length for box)
- `file`: Path to an OBJ or STL file for custom geometry

## Architecture

### GPU Pipeline (Metal)

1. **Collision Detection Phase**
   - Compute cell indices for each particle (spatial hashing)
   - Count particles per cell (histogram)
   - Compute prefix sum for cell offsets
   - Reorder particles by cell for coalesced memory access
   - Calculate inter-particle collisions using rejection sampling

2. **Position Update Phase**
   - Integrate particle positions using velocity

3. **Mesh Intersection Phase** (if mesh loaded)
   - Cast rays from particles in velocity direction
   - Use Metal Acceleration Structure for fast triangle intersection
   - Apply diffuse reflection on hit

4. **Boundary Enforcement Phase**
   - Reflect particles off domain boundaries

### File Structure

```
├── Metal/
│   └── Shaders.metal      # GPU compute kernels
├── MetalCompute.h/.mm     # Metal API wrapper
├── Grid.h/.mm             # GPU-accelerated grid
├── Mesh.h/.cpp            # OBJ/STL mesh loader
├── main.cpp               # Entry point
└── config.json            # Simulation configuration
```

## Overview

This project uses Direct Simulation Monte Carlo (DSMC) methods for modeling (dilute) gas flow.

See [this page](https://ianfriedri.notion.site/DSMC-Readme-fe1d611d1d7d49e09fda3d6477bea094) for mathematical background and more information.

## Performance

On Apple Silicon (M1/M2/M3), the Metal version achieves significant speedups:
- 100K particles: ~10-50x faster than CPU
- 1M particles: ~50-100x faster than CPU

The acceleration structure-based mesh intersection adds minimal overhead while enabling complex geometry interactions.