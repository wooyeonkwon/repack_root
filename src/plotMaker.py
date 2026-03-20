"""Utilities for plotting dual-axis graphs.

This module adds explicit y-axis range controls for GraphDualAxis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt

Range = Optional[Tuple[float, float]]


@dataclass
class Series:
    x: Sequence[float]
    y: Sequence[float]
    label: str = ""
    color: Optional[str] = None


class PlotMaker:
    """Small plotting helper.

    `GraphDualAxis` supports `y1_range` and `y2_range` so callers can
    explicitly control the left/right y-axis ranges.
    """

    def GraphDualAxis(
        self,
        x: Sequence[float],
        y1: Sequence[float],
        y2: Sequence[float],
        *,
        y1_label: str = "y1",
        y2_label: str = "y2",
        x_label: str = "x",
        title: str = "",
        y1_range: Range = None,
        y2_range: Range = None,
        y1_color: str = "tab:blue",
        y2_color: str = "tab:red",
        y1_style: str = "-",
        y2_style: str = "-",
    ):
        """Draw a dual-axis line chart.

        Args:
            x, y1, y2: Data arrays with matching lengths.
            y1_range: Optional (min, max) for left y-axis.
            y2_range: Optional (min, max) for right y-axis.
        """
        if not (len(x) == len(y1) == len(y2)):
            raise ValueError("x, y1, y2 must have the same length")

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(x, y1, linestyle=y1_style, color=y1_color, label=y1_label)
        ax2.plot(x, y2, linestyle=y2_style, color=y2_color, label=y2_label)

        ax1.set_xlabel(x_label)
        ax1.set_ylabel(y1_label, color=y1_color)
        ax2.set_ylabel(y2_label, color=y2_color)
        if title:
            ax1.set_title(title)

        if y1_range is not None:
            self._validate_range(y1_range, "y1_range")
            ax1.set_ylim(*y1_range)
        if y2_range is not None:
            self._validate_range(y2_range, "y2_range")
            ax2.set_ylim(*y2_range)

        fig.tight_layout()
        return fig, (ax1, ax2)

    @staticmethod
    def _validate_range(axis_range: Tuple[float, float], name: str) -> None:
        if len(axis_range) != 2:
            raise ValueError(f"{name} must be a tuple/list with exactly 2 values")
        lo, hi = axis_range
        if lo >= hi:
            raise ValueError(f"{name} min must be smaller than max")
