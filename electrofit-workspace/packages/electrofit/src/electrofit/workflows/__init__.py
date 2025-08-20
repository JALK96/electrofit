"""Removed legacy workflows package.

All previous entry points have migrated to ``electrofit.pipeline.steps``.
This package is intentionally inert and raises an ImportError to surface
necessary migration. See documentation for updated usage examples.
"""

raise ImportError(
	"electrofit.workflows has been removed; use electrofit.pipeline.steps.<stepN> modules (e.g. electrofit.pipeline.steps.step4)."
)


