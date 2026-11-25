#!/usr/bin/env python3
"""The primer pair an rCRUX release was built with, from the files it commits.

rcrux-py records it in the marker's build.log (``--forward-primer X`` in the
command it echoes, or ``Primers: fwd=X, rev=Y``) and, for preset markers, in
src/rcrux/config.py (``PrimerConfig(name="<marker>", forward="X", reverse="Y")``).
Where both have it they agree; some releases have neither.

Standard library only: the release launcher runs it with the system python3.

Usage:
    python3 rcrux_primers.py MARKER BUILD_LOG CONFIG_PY   # prints "FORWARD REVERSE"
"""

import re
import sys


def primers(marker: str, build_log: str, config_py: str) -> tuple[str, str] | None:
    """(forward, reverse), or None when neither file records them."""
    found = []
    for flag, key in (("forward", "fwd"), ("reverse", "rev")):
        match = re.search(rf"--{flag}-primer\s+([A-Za-z]+)", build_log) or re.search(
            rf"\b{key}=([A-Za-z]+)", build_log
        )
        found.append(match.group(1) if match else None)
    if None not in found:
        return found[0], found[1]
    preset = re.search(
        rf'name="{re.escape(marker)}",\s*forward="([A-Za-z]+)",\s*reverse="([A-Za-z]+)"',
        config_py,
    )
    return (preset.group(1), preset.group(2)) if preset else None


def main(argv: list[str]) -> int:
    marker, build_log, config_py = argv
    with open(build_log) as log, open(config_py) as config:
        pair = primers(marker, log.read(), config.read())
    if pair is None:
        return 1
    print(*pair)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
