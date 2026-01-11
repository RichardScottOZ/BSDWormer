"""
Fourier Domain Grid operations for BSDWormer.

This module provides classes for working with grids in both spatial
and Fourier (frequency) domains.
"""

from typing import Optional, Tuple
import numpy as np
from numpy.typing import NDArray


class FourierDomainGrid:
    """Manages spatial and Fourier domain representations of gridded data.
    
    This class handles transformations between spatial and wavenumber domains,
    and maintains the necessary metadata for proper scaling and interpretation.
    
    Attributes:
        spatial_grid: Data in spatial domain
        hat_grid: Data in Fourier (wavenumber) domain
        dx: Grid spacing in x direction
        dy: Grid spacing in y direction
        kx: Wavenumber array in x direction
        ky: Wavenumber array in y direction
        grid_shape: Shape of the grid (rows, columns)
    
    Example:
        >>> import numpy as np
        >>> fdg = FourierDomainGrid(dx=1.0, dy=1.0)
        >>> grid = np.zeros((512, 512), dtype=complex)
        >>> fdg.set_spatial_grid(grid)
        >>> fdg.set_hat_grid(fdg.simple_fft(fdg.spatial_grid))
        >>> fdg.build_wavenumbers(grid)
    """
    
    def __init__(self, dx: float = 1.0, dy: float = 1.0):
        """Initialize FourierDomainGrid.
        
        Args:
            dx: Grid spacing in x direction (default: 1.0)
            dy: Grid spacing in y direction (default: 1.0)
        """
        self.spatial_grid: Optional[NDArray] = None
        self.hat_grid: Optional[NDArray] = None
        self.dx = dx
        self.dy = dy
        self.kx: Optional[NDArray] = None
        self.ky: Optional[NDArray] = None
        self.grid_shape: Optional[Tuple[int, int]] = None
        self.grid_x_len: Optional[int] = None
        self.grid_y_len: Optional[int] = None
        
    def build_wavenumbers(self, grid: NDArray) -> None:
        """Build wavenumber arrays based on grid size.
        
        Computes the discrete Fourier transform sample frequencies
        for both x and y directions.
        
        Args:
            grid: 2D array to determine grid dimensions
        """
        self.grid_shape = grid.shape
        self.grid_y_len = self.grid_shape[0]  # rows
        self.grid_x_len = self.grid_shape[1]  # columns
        
        # fftfreq returns the DFT sample frequencies
        self.kx = np.fft.fftfreq(self.grid_x_len, d=self.dx)
        self.ky = np.fft.fftfreq(self.grid_y_len, d=self.dy)

    def set_spatial_grid(self, grid: NDArray) -> None:
        """Set the spatial domain grid.
        
        Args:
            grid: 2D array representing data in spatial domain
        """
        self.spatial_grid = grid
    
    def set_hat_grid(self, grid: NDArray) -> None:
        """Set the Fourier (wavenumber) domain grid.
        
        Args:
            grid: 2D array representing data in Fourier domain
        """
        self.hat_grid = grid
    
    def simple_fft(self, spatial_grid: NDArray) -> NDArray:
        """Perform 2D Fast Fourier Transform.
        
        Args:
            spatial_grid: Complex array in spatial domain
            
        Returns:
            Complex array in Fourier domain
        """
        return np.fft.fft2(spatial_grid)

    def simple_ifft(self, hat_grid: NDArray) -> NDArray:
        """Perform 2D Inverse Fast Fourier Transform.
        
        Args:
            hat_grid: Complex array in Fourier domain
            
        Returns:
            Complex array in spatial domain
        """
        return np.fft.ifft2(hat_grid)
    
    def __repr__(self) -> str:
        """String representation of FourierDomainGrid."""
        spatial_shape = self.spatial_grid.shape if self.spatial_grid is not None else None
        hat_shape = self.hat_grid.shape if self.hat_grid is not None else None
        return (
            f"FourierDomainGrid(dx={self.dx}, dy={self.dy}, "
            f"spatial_shape={spatial_shape}, hat_shape={hat_shape})"
        )
