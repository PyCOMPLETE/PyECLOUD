#!/bin/bash
set -euo pipefail; IFS=$'\n\t'

NAME=$( python setup.py --name )
VER=$( python setup.py --version )

echo "========================================================================"
echo "Tagging $NAME v$VER"
echo "========================================================================"

git tag "v$VER"
git push origin "v$VER"

echo "========================================================================"
echo "Releasing $NAME v$VER on PyPI"
echo "========================================================================"

# Only upload this release's source archive, even if dist/ contains old builds.
DIST_DIR=$(mktemp -d)
trap 'rm -r -- "$DIST_DIR"' EXIT
python setup.py sdist --dist-dir "$DIST_DIR"
twine upload "$DIST_DIR"/*.tar.gz
