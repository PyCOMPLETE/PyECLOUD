#!/usr/bin/env python3
"""Build, check and publish a source-only release, then push its Git tag."""

import argparse
import ast
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import tomllib


ROOT = Path(__file__).resolve().parent


def run(*command):
    print("+ " + shlex.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, check=True)


def git_output(*args):
    return subprocess.check_output(
        ["git", *args], cwd=ROOT, text=True
    ).strip()


def release_metadata():
    with (ROOT / "pyproject.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    version_attr = metadata["tool"]["setuptools"]["dynamic"]["version"]["attr"]
    module, attribute = version_attr.rsplit(".", 1)
    version_file = ROOT / (module.replace(".", "/") + ".py")
    for node in ast.parse(version_file.read_text()).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == attribute
            for target in node.targets
        ):
            return metadata["project"]["name"], ast.literal_eval(node.value)
    raise ValueError(f"Cannot find {attribute} in {version_file}")


def check_release_checkout(tag):
    if git_output("status", "--porcelain"):
        raise ValueError("Commit or stash checkout changes before publishing.")
    if git_output("tag", "--list", tag):
        raise ValueError(f"Tag {tag} already exists locally.")
    remote = subprocess.run(
        ["git", "ls-remote", "--exit-code", "--tags", "origin", f"refs/tags/{tag}"],
        cwd=ROOT, capture_output=True, text=True,
    )
    if remote.returncode == 0:
        raise ValueError(f"Tag {tag} already exists on origin.")
    if remote.returncode != 2:  # Git returns 2 when no matching ref exists.
        raise ValueError(f"Cannot check origin tags: {remote.stderr.strip()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-only", action="store_true",
                        help="build and check the source archive; do not upload or tag")
    parser.add_argument("--no-isolation", action="store_true",
                        help="use build dependencies from the active environment")
    args = parser.parse_args()

    name, version = release_metadata()
    tag = f"v{version}"
    if not args.build_only:
        check_release_checkout(tag)
    print(f"Preparing {name} {version} (source distribution only)", flush=True)

    # A fresh directory prevents old wheels or other releases from being uploaded.
    with tempfile.TemporaryDirectory(prefix=f"{name}-release-") as temporary:
        command = [sys.executable, "-m", "build", "--sdist", "--outdir", temporary]
        if args.no_isolation:
            command.append("--no-isolation")
        run(*command, str(ROOT))
        artifacts = list(Path(temporary).iterdir())
        if len(artifacts) != 1 or not artifacts[0].name.endswith(".tar.gz"):
            raise ValueError(f"Expected exactly one source archive, got {artifacts}")
        run(sys.executable, "-m", "twine", "check", "--strict", artifacts[0])
        destination = ROOT / "dist" / artifacts[0].name
        destination.parent.mkdir(exist_ok=True)
        shutil.copy2(artifacts[0], destination)

    print(f"Checked source archive: {destination}", flush=True)
    if args.build_only:
        return

    run(sys.executable, "-m", "twine", "upload", "--non-interactive", destination)
    # Tag only after a successful upload, so failed uploads leave no release tag.
    run("git", "tag", "-a", tag, "-m", f"Release {name} {version}")
    run("git", "push", "origin", f"refs/tags/{tag}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as exc:
        raise SystemExit(str(exc))
