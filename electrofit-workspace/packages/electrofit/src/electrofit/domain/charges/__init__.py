"""Domain logic for computing and applying atomic partial charges.

Current status: placeholder package created during migration.
Upcoming steps:
- Factor pure functions out of core.process_conform (geometry prep, gaussian run, resp fitting, mol2 update).
- Introduce service-level orchestration separate from CLI / scratch management.

Guidelines:
- No direct subprocess invocations here long-term (delegate to adapters layer)
- Keep I/O side-effects explicit via function parameters.
"""

__all__ = []
