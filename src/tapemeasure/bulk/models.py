from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BulkinessConfig:
    mode: str
    input_path: Path
    output_dir: Path
    radius: float = 10.0
    max_neighbors: int = 15
    make_plots: bool = False
