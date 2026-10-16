#!/bin/bash
# Tag a release. Pushing the tag to the main repository triggers the
# "Publish to PyPI" workflow (.github/workflows/publish.yml), which builds and
# uploads the source distribution.
set -euo pipefail; IFS=$'\n\t'

NAME="PyECLOUD"
REMOTE="${REMOTE:-origin}"
VER=$( python -c 'exec(open("PyECLOUD/_version.py").read()); print(__version__)' )

# The publish workflow only runs for vX.Y.Z tags; refuse anything else.
if [[ ! "$VER" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Version '$VER' is not of the form X.Y.Z; the publish workflow would not run." >&2
    exit 1
fi

echo "========================================================================"
echo "Checking that $NAME v$VER builds"
echo "========================================================================"

# Build and check locally first, so a broken release is caught before tagging.
DIST_DIR=$(mktemp -d)
trap 'rm -r -- "$DIST_DIR"' EXIT
python -m build --sdist --outdir "$DIST_DIR"
python -m twine check --strict "$DIST_DIR"/*.tar.gz

echo "========================================================================"
echo "Tagging $NAME v$VER and pushing to $REMOTE"
echo "========================================================================"

git tag "v$VER"
git push "$REMOTE" "v$VER"

echo "========================================================================"
echo "Pushed v$VER. Follow the \"Publish to PyPI\" workflow in the Actions tab."
echo "========================================================================"
