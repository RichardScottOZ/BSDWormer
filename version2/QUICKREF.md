# BSDWormer Quick Reference

## Installation

```bash
cd version2
pip install -e .
```

## Basic Usage

### Simple Processing

```python
from bsdwormer import Wormer

# Create wormer
wormer = Wormer()

# Load data
wormer.import_raster('magnetic_data.tif')

# Build padded grid
wormer.build_padded_raster((2048, 2048))

# Process at depth
wormer.worm_level_as_image(dz=100)

# Export
wormer.export_raster(wormer.worm_image, 'output.tif')
```

### With Configuration

```python
from bsdwormer import Wormer, load_config

config = load_config('config.yaml')
wormer = Wormer(config=config)
wormer.import_raster('data.tif')
```

### Multiple Levels

```python
depths = [50, 100, 200, 500]
for dz in depths:
    wormer.worm_level_as_points(dz)
    wormer.build_worm_segs(dz=dz)
    wormer.build_level_for_vtk(dz)

# Export VTK
from bsdwormer.visualization import write_vtk_worm_levels
write_vtk_worm_levels('output', wormer.all_points, 
                      wormer.all_lines, wormer.all_vals)
```

## Key Classes

### Wormer
Main class for worm detection.

**Methods:**
- `import_raster(filename)` - Load raster data
- `build_padded_raster(shape, rolloff=100, pad_type='hann')` - Prepare for FFT
- `worm_level(dz)` - Process at depth dz
- `worm_level_as_image(dz)` - Get worm image
- `worm_level_as_points(dz)` - Get sub-pixel points
- `build_worm_segs(dz=dz)` - Build graph structure
- `build_level_for_vtk(dz)` - Prepare VTK output
- `export_raster(array, filename)` - Save raster

### FourierDomainGrid
Manages spatial/Fourier domain data.

**Methods:**
- `set_spatial_grid(grid)` - Set spatial data
- `set_hat_grid(grid)` - Set Fourier data
- `simple_fft(grid)` - Forward FFT
- `simple_ifft(grid)` - Inverse FFT
- `build_wavenumbers(grid)` - Compute wavenumbers

### FourierDomainOps
Fourier domain operations.

**Methods:**
- `build_mod_k()` - Modulus of wavenumber
- `build_upward_continuation_op(dz)` - Upward continuation
- `build_dx_op()` - X derivative operator
- `build_dy_op()` - Y derivative operator
- `build_grad_vector(fdg)` - Gradient vector
- `canny_edge_detect(fdg)` - Edge detection

## Configuration File

```yaml
padding:
  rolloff_size: 100
  pad_type: 'hann'

processing:
  nodata_value: -100
  log_vals: true
  clipped: true

output:
  format: 'GTiff'
  compression: 'LZW'
```

## Command Line (Future)

```bash
# Process single file
bsdwormer process input.tif -o output.tif --dz 100

# Multiple levels
bsdwormer process input.tif --levels 50,100,200,500

# With config
bsdwormer process input.tif --config config.yaml
```

## Common Patterns

### Batch Processing

```python
from pathlib import Path
from bsdwormer import Wormer

input_dir = Path('data/')
output_dir = Path('results/')
output_dir.mkdir(exist_ok=True)

for input_file in input_dir.glob('*.tif'):
    wormer = Wormer()
    wormer.import_raster(str(input_file))
    wormer.build_padded_raster((2048, 2048))
    wormer.worm_level_as_image(100)
    
    output_file = output_dir / f"{input_file.stem}_worms.tif"
    wormer.export_raster(wormer.worm_image, str(output_file))
```

### Custom Parameters

```python
wormer = Wormer()
wormer.import_raster('data.tif')

# Custom padding
wormer.build_padded_raster(
    padded_shape=(4096, 4096),
    rolloff_size=200,
    pad_type='hamming'
)

# Custom worm detection
wormer.build_worm_segs(
    nodata_in_worm_image=-999,
    clipped=True,
    log_vals=False,
    dz=100
)
```

## Troubleshooting

### GDAL Not Found
```bash
# Ubuntu/Debian
sudo apt-get install python3-gdal

# macOS
brew install gdal
```

### Import Errors
```bash
# Reinstall in development mode
pip install -e .
```

### Memory Issues
- Use smaller padded sizes
- Process one level at a time
- Clear intermediate results

### Poor Results
- Check input data quality
- Adjust rolloff_size
- Try different pad_type
- Experiment with depth levels

## Tips

1. **Padded Size**: Use power of 2 for faster FFT
2. **Depth Levels**: Start with coarse, refine based on features
3. **Log Values**: Usually better for visualization
4. **VTK Output**: Use ParaView for 3D visualization
5. **Configuration**: Create config file for reproducibility

## References

- Full documentation: `docs/`
- Examples: `examples/`
- Tests: `tests/`
- Original code: `../src/`
