"""Optional checkout provenance for simulation logs."""

from pathlib import Path
import subprocess


def get_git_info():
    repository = Path(__file__).resolve().parent.parent
    unavailable = ('git hash: unavailable', 'git branch: unavailable')
    if not (repository / '.git').exists():
        return unavailable
    try:
        values = [subprocess.check_output(
            ['git', '-C', str(repository), 'rev-parse', *args],
            text=True, stderr=subprocess.DEVNULL).strip()
            for args in [('HEAD',), ('--abbrev-ref', 'HEAD')]]
    except (OSError, subprocess.CalledProcessError):
        return unavailable
    return f'git hash: {values[0]}', f'git branch: {values[1]}'
