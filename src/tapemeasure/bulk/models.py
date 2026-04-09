from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BulkinessConfig:
    mode: str
    input_csv: Path
    output_dir: Path
    radii: list[float] = field(default_factory=lambda: [10.0])
    k_neighbors: int = 10
    subtract_baseline: bool = False
    baseline_spacing: float = 3.8
    histogram_bins: int = 40
    make_plots: bool = False
