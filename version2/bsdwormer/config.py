"""Configuration management for BSDWormer."""

from typing import Any, Dict, Optional
import yaml
import json
from pathlib import Path


class Config:
    """Configuration class for BSDWormer.
    
    Attributes:
        padding: Padding configuration
        processing: Processing parameters
        output: Output settings
    """
    
    def __init__(self, config_dict: Optional[Dict[str, Any]] = None):
        """Initialize configuration.
        
        Args:
            config_dict: Dictionary of configuration parameters
        """
        self.padding = {
            'rolloff_size': 100,
            'pad_type': 'hann'
        }
        self.processing = {
            'nodata_value': -100,
            'log_vals': True,
            'clipped': True
        }
        self.output = {
            'format': 'GTiff',
            'compression': 'LZW'
        }
        
        if config_dict:
            self.update(config_dict)
    
    def update(self, config_dict: Dict[str, Any]) -> None:
        """Update configuration from dictionary.
        
        Args:
            config_dict: Dictionary of configuration parameters
        """
        if 'padding' in config_dict:
            self.padding.update(config_dict['padding'])
        if 'processing' in config_dict:
            self.processing.update(config_dict['processing'])
        if 'output' in config_dict:
            self.output.update(config_dict['output'])
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary.
        
        Returns:
            Dictionary of configuration parameters
        """
        return {
            'padding': self.padding,
            'processing': self.processing,
            'output': self.output
        }


def load_config(config_path: str) -> Config:
    """Load configuration from file.
    
    Args:
        config_path: Path to configuration file (YAML or JSON)
        
    Returns:
        Config object
        
    Raises:
        ValueError: If file format is not supported
        FileNotFoundError: If config file doesn't exist
    """
    path = Path(config_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(path, 'r') as f:
        if path.suffix in ['.yaml', '.yml']:
            config_dict = yaml.safe_load(f)
        elif path.suffix == '.json':
            config_dict = json.load(f)
        else:
            raise ValueError(
                f"Unsupported config format: {path.suffix}. "
                "Use .yaml, .yml, or .json"
            )
    
    return Config(config_dict)


def save_config(config: Config, config_path: str) -> None:
    """Save configuration to file.
    
    Args:
        config: Config object to save
        config_path: Path to save configuration file (YAML or JSON)
        
    Raises:
        ValueError: If file format is not supported
    """
    path = Path(config_path)
    config_dict = config.to_dict()
    
    with open(path, 'w') as f:
        if path.suffix in ['.yaml', '.yml']:
            yaml.dump(config_dict, f, default_flow_style=False)
        elif path.suffix == '.json':
            json.dump(config_dict, f, indent=2)
        else:
            raise ValueError(
                f"Unsupported config format: {path.suffix}. "
                "Use .yaml, .yml, or .json"
            )
