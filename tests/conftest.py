from __future__ import annotations

import sys
from pathlib import Path

# Make tools/ importable when pytest is run from the repo root.
sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
