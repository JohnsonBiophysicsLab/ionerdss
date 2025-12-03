"""
Plotting configuration and styles.
"""

import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any

class PlotStyle:
    """Manages plot styling."""
    
    @staticmethod
    def set_style(style: str = "whitegrid", context: str = "paper"):
        """
        Sets the global plotting style.
        
        Args:
            style: Seaborn style ('whitegrid', 'darkgrid', etc.)
            context: Seaborn context ('paper', 'notebook', 'talk', 'poster')
        """
        sns.set_style(style)
        sns.set_context(context)
        plt.rcParams['figure.figsize'] = (10, 6)
        plt.rcParams['font.size'] = 12

    @staticmethod
    def get_default_kwargs() -> Dict[str, Any]:
        return {
            "linewidth": 2,
            "markersize": 8,
        }

