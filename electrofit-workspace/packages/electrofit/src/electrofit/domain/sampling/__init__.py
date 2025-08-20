"""Sampling domain logic (frame selection, conformer extraction helpers).

Public API kept deliberately small for orchestrators:

select_frame_indices(traj, n, method, seed) -> list[int]
prepare_conformer_directory(...)

Implementation details (RMSD farthest point etc.) live here to keep pipeline
step modules thin and side-effect free on import.
"""
from .frames import select_frame_indices  # noqa: F401
# Backwards compatibility alias for legacy tests importing workflows.step4_extract_conforms._select_indices
_select_indices = select_frame_indices  # noqa: F401
from .conformer_io import prepare_conformer_directory  # noqa: F401
