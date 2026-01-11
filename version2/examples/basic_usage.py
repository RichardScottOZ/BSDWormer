"""
Basic usage example for BSDWormer.

This example demonstrates how to:
1. Load a geophysical raster dataset
2. Process it at multiple depth levels
3. Export results as images and VTK files
"""

import numpy as np
from pathlib import Path

# Import BSDWormer
from bsdwormer import Wormer, load_config
from bsdwormer.visualization import write_vtk_worm_levels


def main():
    """Run a basic worming workflow."""
    
    # Configuration (optional)
    # config = load_config('config.yaml')
    # wormer = Wormer(config=config)
    
    # Or use defaults
    wormer = Wormer()
    
    # Load input raster
    input_file = "path/to/your/magnetic_data.tif"
    
    print(f"Loading data from {input_file}...")
    wormer.import_raster(input_file)
    
    print(f"Data shape: {wormer.base_grid.shape}")
    print(f"Cell size: dx={wormer.dx}, dy={wormer.dy}")
    
    # Build padded raster for FFT processing
    # Padded shape should be larger than input, preferably power of 2
    rows, cols = wormer.base_grid.shape
    padded_rows = 2 ** int(np.ceil(np.log2(rows * 1.5)))
    padded_cols = 2 ** int(np.ceil(np.log2(cols * 1.5)))
    
    print(f"Building padded raster: ({padded_rows}, {padded_cols})...")
    wormer.build_padded_raster(
        padded_shape=(padded_rows, padded_cols),
        rolloff_size=100,
        pad_type='hann'
    )
    
    # Process at multiple depth/height levels
    # Units depend on your input data (typically meters)
    depths = [50, 100, 200, 500, 1000]
    
    print("\nProcessing worm levels...")
    for i, dz in enumerate(depths, 1):
        print(f"  Level {i}/{len(depths)}: dz = {dz}")
        
        # Compute worms as points (sub-pixel precision)
        wormer.worm_level_as_points(dz)
        
        # Build graph structure and segments
        wormer.build_worm_segs(
            dz=dz,
            clipped=True,
            log_vals=True,
            from_image=False
        )
        
        # Build VTK representation
        wormer.build_level_for_vtk(
            dz,
            invert_z=True,
            delta_z_in_units=dz
        )
        
        # Optionally export each level as an image
        worm_image = wormer.worm_level_as_image(dz)
        output_image = f"worms_level_{dz}m.tif"
        wormer.export_raster(worm_image, output_image)
        print(f"    Exported: {output_image}")
    
    # Export all levels to VTK for 3D visualization
    output_vtk = "worms_multilevel"
    print(f"\nExporting VTK files: {output_vtk}_level_*.vtk")
    write_vtk_worm_levels(
        output_vtk,
        wormer.all_points,
        wormer.all_lines,
        wormer.all_vals
    )
    
    print("\nProcessing complete!")
    print(f"Output files:")
    print(f"  - Raster images: worms_level_*m.tif")
    print(f"  - VTK files: {output_vtk}_level_*.vtk")
    print("\nVTK files can be visualized in ParaView or similar software.")


if __name__ == "__main__":
    main()
