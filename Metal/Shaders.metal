//
// Shaders.metal
// DSMC Metal Compute Shaders
//

#include <metal_stdlib>
#include <metal_raytracing>

using namespace metal;
using namespace raytracing;

// =============================================================================
// Data Structures
// =============================================================================

struct Particle {
    packed_float3 pos;
    float  pad0;  // padding for alignment
    packed_float3 vel;
    float  pad1;
};

struct SimulationParams {
    float delta_t;       // timestep
    float v_max;         // maximum velocity
    float d;             // particle diameter
    float cell_length;   // cell side length
    float domain_min;    // domain minimum (0)
    float domain_max;    // domain maximum (dim * cell_length)
    float v_mult;        // velocity multiplier for reflections
    float eps;           // small epsilon for boundary enforcement
    uint  num_particles;
    uint  num_cells_per_dim;
    uint  f_n;           // number of physical particles per simulated particle
    float V_c;           // cell volume
};

struct CollisionParams {
    uint cell_start;     // start index in sorted particle buffer
    uint cell_end;       // end index in sorted particle buffer
    uint m_cand;         // number of collision candidates
    uint seed;           // random seed for this cell
};

// =============================================================================
// Random Number Generation (PCG-based)
// =============================================================================

// PCG random number generator
uint pcg_hash(uint input) {
    uint state = input * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

// Generate random float in [0, 1)
float rand_float(thread uint& seed) {
    seed = pcg_hash(seed);
    return float(seed) / float(0xFFFFFFFFu);
}

// Generate random float in [min, max)
float rand_range(thread uint& seed, float min_val, float max_val) {
    return min_val + rand_float(seed) * (max_val - min_val);
}

// Generate random unit vector on sphere
float3 random_unit_sphere(thread uint& seed) {
    float theta = 2.0f * M_PI_F * rand_float(seed);
    float u = 1.0f - 2.0f * rand_float(seed);  // uniform in [-1, 1]
    u = clamp(u, -1.0f, 1.0f);  // ensure valid range for acos
    float phi = acos(u);
    return float3(
        sin(phi) * cos(theta),
        sin(phi) * sin(theta),
        cos(phi)
    );
}

// Generate random unit vector in hemisphere (pointing in +normal direction)
float3 random_hemisphere(thread uint& seed, float3 normal) {
    float3 v = random_unit_sphere(seed);
    if (dot(v, normal) < 0.0f) {
        v = -v;
    }
    return v;
}

// =============================================================================
// Phase 1: Basic GPU Particle System
// =============================================================================

// Kernel 1: Update particle positions (replaces Cell::updatePositions)
kernel void updatePositions(
    device Particle* particles [[buffer(0)]],
    constant SimulationParams& params [[buffer(1)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    particles[gid].pos += particles[gid].vel * params.delta_t;
}

// Kernel 2: Enforce domain boundaries (replaces Grid::enforceDomain)
kernel void enforceDomain(
    device Particle* particles [[buffer(0)]],
    constant SimulationParams& params [[buffer(1)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    device Particle& p = particles[gid];
    uint seed = gid * 1099087573u + params.num_particles;
    
    float3 pos = p.pos;
    float3 vel = p.vel;
    float speed = length(vel);
    
    // Handle zero-speed case - give particle a small random velocity
    if (speed < 0.0001f) {
        speed = params.v_mult;
        vel = speed * random_unit_sphere(seed);
    }
    
    // X boundaries
    if (pos.x >= params.domain_max) {
        pos.x = params.domain_max - params.eps;
        vel = speed * random_hemisphere(seed, float3(-1, 0, 0));
    }
    if (pos.x <= params.domain_min) {
        pos.x = params.domain_min + params.eps;
        vel = speed * random_hemisphere(seed, float3(1, 0, 0));
    }
    
    // Y boundaries (with energy modification at min Y - hot wall)
    if (pos.y >= params.domain_max) {
        pos.y = params.domain_max - params.eps;
        vel = speed * random_hemisphere(seed, float3(0, -1, 0));
    }
    if (pos.y <= params.domain_min) {
        pos.y = params.domain_min + params.eps;
        vel = speed * random_hemisphere(seed, float3(0, 1, 0));
        vel.y *= 1.05f; // gain 5% energy (hot wall)
    }
    
    // Z boundaries
    if (pos.z >= params.domain_max) {
        pos.z = params.domain_max - params.eps;
        vel = speed * random_hemisphere(seed, float3(0, 0, -1));
    }
    if (pos.z <= params.domain_min) {
        pos.z = params.domain_min + params.eps;
        vel = speed * random_hemisphere(seed, float3(0, 0, 1));
    }
    
    p.pos = pos;
    p.vel = vel;
}

// =============================================================================
// Phase 2: Spatial Hashing & Collision Detection
// =============================================================================

// Compute cell index from position
inline uint3 getCellIndex3D(float3 pos, float cell_length, uint num_cells) {
    uint3 idx;
    idx.x = clamp(uint(pos.x / cell_length), 0u, num_cells - 1u);
    idx.y = clamp(uint(pos.y / cell_length), 0u, num_cells - 1u);
    idx.z = clamp(uint(pos.z / cell_length), 0u, num_cells - 1u);
    return idx;
}

inline uint getCellIndex(float3 pos, float cell_length, uint num_cells) {
    uint3 idx = getCellIndex3D(pos, cell_length, num_cells);
    return idx.x * num_cells * num_cells + idx.y * num_cells + idx.z;
}

// Kernel 3: Compute cell indices for each particle (spatial hashing)
kernel void computeCellIndices(
    device Particle* particles [[buffer(0)]],
    device uint* cellIndices [[buffer(1)]],
    device uint* particleIndices [[buffer(2)]],
    constant SimulationParams& params [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    uint cellIdx = getCellIndex(particles[gid].pos, params.cell_length, params.num_cells_per_dim);
    cellIndices[gid] = cellIdx;
    particleIndices[gid] = gid;
}

// Kernel 4: Count particles per cell (histogram)
kernel void countParticlesPerCell(
    device const uint* cellIndices [[buffer(0)]],
    device atomic_uint* cellCounts [[buffer(1)]],
    constant SimulationParams& params [[buffer(2)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    uint cellIdx = cellIndices[gid];
    atomic_fetch_add_explicit(&cellCounts[cellIdx], 1u, memory_order_relaxed);
}

// Kernel 5: Compute cell start/end offsets (prefix sum done on CPU or separate kernel)
kernel void computeCellOffsets(
    device const uint* cellCounts [[buffer(0)]],
    device uint* cellStarts [[buffer(1)]],
    device uint* cellEnds [[buffer(2)]],
    constant uint& numCells [[buffer(3)]],
    constant uint* prefixSum [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= numCells) return;
    
    cellStarts[gid] = (gid == 0) ? 0 : prefixSum[gid - 1];
    cellEnds[gid] = prefixSum[gid];
}

// Kernel 6: Reorder particles by cell (for coalesced access)
kernel void reorderParticles(
    device const Particle* particlesIn [[buffer(0)]],
    device Particle* particlesOut [[buffer(1)]],
    device const uint* particleIndices [[buffer(2)]],
    constant SimulationParams& params [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    uint srcIdx = particleIndices[gid];
    if (srcIdx < params.num_particles) {
        particlesOut[gid] = particlesIn[srcIdx];
    }
}

// Kernel 7: Inter-particle collisions using rejection sampling (per-cell)
kernel void calculateCollisions(
    device Particle* particles [[buffer(0)]],
    device const uint* cellStarts [[buffer(1)]],
    device const uint* cellEnds [[buffer(2)]],
    constant SimulationParams& params [[buffer(3)]],
    constant uint& frameNumber [[buffer(4)]],
    uint gid [[thread_position_in_grid]]  // one thread per cell
) {
    uint numCells = params.num_cells_per_dim * params.num_cells_per_dim * params.num_cells_per_dim;
    if (gid >= numCells) return;
    
    uint start = cellStarts[gid];
    uint end = cellEnds[gid];
    uint n_c = end - start;
    
    if (n_c < 2) return;  // need at least 2 particles for collision
    
    // Initialize random seed for this cell
    uint seed = gid * 1099087573u + frameNumber * 747796405u;
    
    // Calculate number of collision candidates
    float pi = M_PI_F;
    uint m_cand = uint((float(n_c) * float(n_c - 1) * float(params.f_n) * pi *
                        params.d * params.d * params.v_max * params.delta_t) / (2.0f * params.V_c));
    
    // Limit to reasonable number of collision candidates
    m_cand = min(m_cand, n_c * 10u);
    
    for (uint i = 0; i < m_cand; i++) {
        // Pick 2 random particles in this cell
        uint idx1 = start + uint(rand_float(seed) * float(n_c));
        uint idx2 = start + uint(rand_float(seed) * float(n_c));
        
        // Make sure they're different
        if (idx1 == idx2) continue;
        
        // Clamp to valid range
        idx1 = min(idx1, end - 1);
        idx2 = min(idx2, end - 1);
        
        // Pick random proportion of maximum speed
        float v_tmp = rand_float(seed) * params.v_max;
        
        // Relative speed
        float3 vel1 = particles[idx1].vel;
        float3 vel2 = particles[idx2].vel;
        float v_r = length(vel1 - vel2);
        
        // Skip if relative velocity is too small (avoid numerical issues)
        if (v_r < 0.0001f) continue;
        
        if (v_r > v_tmp) {
            // Collision accepted - calculate new velocities
            float3 v_star = random_unit_sphere(seed);
            
            // Center of mass velocity
            float3 v_cm = 0.5f * (vel1 + vel2);
            
            // Update velocities (momentum conserving)
            particles[idx1].vel = v_cm + 0.5f * v_star * v_r;
            particles[idx2].vel = v_cm - 0.5f * v_star * v_r;
        }
    }
}

// =============================================================================
// Phase 3: Ray-Object Intersection with Acceleration Structures
// =============================================================================

// Ray structure for intersection testing
struct Ray {
    float3 origin;
    float3 direction;
    float min_distance;
    float max_distance;
};

// Kernel 8: Test particle trajectories against mesh using acceleration structure
kernel void intersectMesh(
    device Particle* particles [[buffer(0)]],
    constant SimulationParams& params [[buffer(1)]],
    instance_acceleration_structure accelStruct [[buffer(2)]],
    constant uint& frameNumber [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    device Particle& p = particles[gid];
    float speed = length(p.vel);
    if (speed < 0.0001f) return;
    
    // Create ray from particle position in velocity direction
    float3 velDir = p.vel / speed;  // safe since we checked speed above
    
    ray r;
    r.origin = p.pos;
    r.direction = velDir;
    r.min_distance = 0.0001f;
    r.max_distance = speed * params.delta_t * 2.0f;  // look ahead
    
    // Create intersector
    intersector<triangle_data, instancing> intersector;
    intersector.accept_any_intersection(false);  // we want closest hit
    
    // Perform intersection
    intersection_result<triangle_data, instancing> result = intersector.intersect(r, accelStruct);
    
    if (result.type != intersection_type::none) {
        // Hit! Perform diffuse reflection
        float hit_distance = result.distance;
        
        // Get triangle normal (from primitive data if available, otherwise compute)
        // For now, use the geometric normal
        uint primIdx = result.primitive_id;
        
        // Move particle to just before hit point
        float3 hit_point = r.origin + r.direction * (hit_distance - params.eps);
        p.pos = hit_point;
        
        // Get surface normal (simplified - assumes we have it in instance data)
        // In a full implementation, you'd interpolate vertex normals
        float3 normal = result.triangle_front_facing ?
                        -normalize(r.direction) : normalize(r.direction);
        
        // Initialize random seed
        uint seed = gid * 1099087573u + frameNumber * 747796405u + primIdx;
        
        // Diffuse reflection: random direction in hemisphere around normal
        float speed = length(p.vel);
        p.vel = speed * random_hemisphere(seed, -r.direction);
    }
}

// Alternative kernel using primitive acceleration structure (no instancing)
kernel void intersectMeshPrimitive(
    device Particle* particles [[buffer(0)]],
    constant SimulationParams& params [[buffer(1)]],
    primitive_acceleration_structure accelStruct [[buffer(2)]],
    device const packed_float3* triangleNormals [[buffer(3)]],
    constant uint& frameNumber [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    device Particle& p = particles[gid];
    float speed = length(p.vel);
    if (speed < 0.0001f) return;
    
    // Create ray from particle position in velocity direction
    ray r;
    r.origin = p.pos;
    r.direction = normalize(p.vel);
    r.min_distance = 0.0001f;
    r.max_distance = speed * params.delta_t * 2.0f;
    
    // Create intersector
    intersector<triangle_data> inter;
    inter.accept_any_intersection(false);
    
    // Perform intersection
    intersection_result<triangle_data> result = inter.intersect(r, accelStruct);
    
    if (result.type != intersection_type::none) {
        float hit_distance = result.distance;
        uint primIdx = result.primitive_id;
        
        // Move particle to just before hit point
        float3 hit_point = r.origin + r.direction * (hit_distance * 0.99f);
        p.pos = hit_point;
        
        // Get pre-computed triangle normal
        float3 normal = float3(triangleNormals[primIdx]);
        if (!result.triangle_front_facing) {
            normal = -normal;
        }
        
        // Initialize random seed
        uint seed = gid * 1099087573u + frameNumber * 747796405u + primIdx;
        
        // Diffuse reflection
        p.vel = speed * random_hemisphere(seed, normal);
    }
}

// =============================================================================
// Utility Kernels
// =============================================================================

// Initialize particles with random positions and velocities
kernel void initializeParticles(
    device Particle* particles [[buffer(0)]],
    constant SimulationParams& params [[buffer(1)]],
    constant uint& seed_offset [[buffer(2)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    uint seed = gid * 1099087573u + seed_offset;
    
    // Random position in domain
    particles[gid].pos = float3(
        rand_range(seed, params.domain_min, params.domain_max),
        rand_range(seed, params.domain_min, params.domain_max),
        rand_range(seed, params.domain_min, params.domain_max)
    );
    
    // Random velocity direction with magnitude v_mult
    particles[gid].vel = params.v_mult * random_unit_sphere(seed);
    particles[gid].pad0 = 0;
    particles[gid].pad1 = 0;
}

// Copy particle data back for output (convert to double precision on CPU)
kernel void copyParticlePositions(
    device const Particle* particles [[buffer(0)]],
    device float4* positions [[buffer(1)]],
    constant SimulationParams& params [[buffer(2)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= params.num_particles) return;
    
    positions[gid] = float4(particles[gid].pos, length(particles[gid].vel));
}

// Clear buffer to zero
kernel void clearBuffer(
    device atomic_uint* buffer [[buffer(0)]],
    constant uint& count [[buffer(1)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= count) return;
    atomic_store_explicit(&buffer[gid], 0u, memory_order_relaxed);
}
