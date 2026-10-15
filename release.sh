#!/bin/bash
set -euo pipefail; IFS=$'\n\t'

NAME="PyECLOUD"
VER=$( python -c 'exec(open("PyECLOUD/_version.py").read()); print(__version__)' )

echo "========================================================================"
echo "Building $NAME v$VER source archive"
echo "========================================================================"

# Build before tagging, so a failed build does not leave a pushed tag behind.
# Only upload this release's source archive, even if dist/ contains old builds.
DIST_DIR=$(mktemp -d)
trap 'rm -r -- "$DIST_DIR"' EXIT
python -m build --sdist --outdir "$DIST_DIR"

echo "========================================================================"
echo "Tagging $NAME v$VER"
echo "========================================================================"

git tag "v$VER"
git push origin "v$VER"

echo "========================================================================"
echo "Releasing $NAME v$VER on PyPI"
echo "========================================================================"

twine upload "$DIST_DIR"/*.tar.gz
