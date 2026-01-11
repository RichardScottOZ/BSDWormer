"""
Test suite for FourierDomainGrid.

These tests validate the basic functionality of the FourierDomainGrid class.
"""

import numpy as np
import pytest
from bsdwormer.core.fourier_domain_grid import FourierDomainGrid


class TestFourierDomainGridInitialization:
    """Tests for FourierDomainGrid initialization."""
    
    def test_default_initialization(self):
        """Test that FourierDomainGrid initializes with correct defaults."""
        fdg = FourierDomainGrid()
        assert fdg.dx == 1.0
        assert fdg.dy == 1.0
        assert fdg.spatial_grid is None
        assert fdg.hat_grid is None
    
    def test_custom_spacing(self):
        """Test initialization with custom grid spacing."""
        fdg = FourierDomainGrid(dx=2.5, dy=3.0)
        assert fdg.dx == 2.5
        assert fdg.dy == 3.0


class TestGridOperations:
    """Tests for grid setting and manipulation."""
    
    def test_set_spatial_grid(self):
        """Test setting spatial grid."""
        fdg = FourierDomainGrid()
        grid = np.zeros((64, 64), dtype=complex)
        fdg.set_spatial_grid(grid)
        assert np.array_equal(fdg.spatial_grid, grid)
    
    def test_set_hat_grid(self):
        """Test setting Fourier grid."""
        fdg = FourierDomainGrid()
        grid = np.zeros((64, 64), dtype=complex)
        fdg.set_hat_grid(grid)
        assert np.array_equal(fdg.hat_grid, grid)


class TestWavenumbers:
    """Tests for wavenumber computation."""
    
    def test_build_wavenumbers(self):
        """Test wavenumber array construction."""
        fdg = FourierDomainGrid(dx=1.0, dy=1.0)
        grid = np.zeros((512, 512), dtype=complex)
        fdg.build_wavenumbers(grid)
        
        assert fdg.grid_shape == (512, 512)
        assert fdg.grid_x_len == 512
        assert fdg.grid_y_len == 512
        assert len(fdg.kx) == 512
        assert len(fdg.ky) == 512
        assert fdg.kx[0] == 0.0
        assert fdg.ky[0] == 0.0
    
    def test_wavenumber_symmetry(self):
        """Test that wavenumbers are symmetric around zero."""
        fdg = FourierDomainGrid(dx=1.0, dy=1.0)
        grid = np.zeros((128, 128), dtype=complex)
        fdg.build_wavenumbers(grid)
        
        # Check Nyquist frequency
        assert np.isclose(np.max(fdg.kx), 0.5 - 1.0/128)
        assert np.isclose(np.min(fdg.kx), -0.5)


class TestFFTOperations:
    """Tests for FFT and inverse FFT."""
    
    def test_fft_delta_function(self):
        """Test FFT of delta function at origin."""
        fdg = FourierDomainGrid()
        grid = np.zeros((512, 512), dtype=complex)
        grid[0, 0] = 1.0 + 0.0j
        
        hat = fdg.simple_fft(grid)
        # FFT of delta at origin should be constant (1+0j) everywhere
        assert np.allclose(hat, 1.0 + 0.0j)
    
    def test_ifft_delta_function(self):
        """Test inverse FFT of delta function at origin."""
        fdg = FourierDomainGrid()
        hat = np.zeros((512, 512), dtype=complex)
        hat[0, 0] = 1.0 + 0.0j
        
        spatial = fdg.simple_ifft(hat)
        # IFFT of delta should be constant everywhere (scaled by 1/N^2)
        expected = (1.0 + 0.0j) / (512.0 * 512.0)
        assert np.allclose(spatial, expected)
    
    def test_fft_roundtrip(self):
        """Test that FFT followed by IFFT recovers original."""
        fdg = FourierDomainGrid()
        original = np.random.rand(64, 64) + 1j * np.random.rand(64, 64)
        
        hat = fdg.simple_fft(original)
        recovered = fdg.simple_ifft(hat)
        
        assert np.allclose(original, recovered)
    
    def test_fft_real_input(self):
        """Test FFT with real-valued input."""
        fdg = FourierDomainGrid()
        real_grid = np.random.rand(64, 64)
        
        hat = fdg.simple_fft(real_grid)
        recovered = fdg.simple_ifft(hat)
        
        # Should recover real part
        assert np.allclose(real_grid, recovered.real)
        # Imaginary part should be near zero
        assert np.allclose(recovered.imag, 0, atol=1e-10)


class TestRepr:
    """Tests for string representation."""
    
    def test_repr_empty(self):
        """Test repr with no grids set."""
        fdg = FourierDomainGrid(dx=2.0, dy=3.0)
        repr_str = repr(fdg)
        assert "dx=2.0" in repr_str
        assert "dy=3.0" in repr_str
        assert "None" in repr_str
    
    def test_repr_with_grids(self):
        """Test repr with grids set."""
        fdg = FourierDomainGrid()
        grid = np.zeros((64, 128))
        fdg.set_spatial_grid(grid)
        
        repr_str = repr(fdg)
        assert "(64, 128)" in repr_str


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
