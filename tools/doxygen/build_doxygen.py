#!/usr/bin/env python3
"""Build the complete CCPP Physics GitHub Pages site locally."""

import argparse
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


NCAR_PHYSICS_URL = "https://github.com/NCAR/ccpp-physics.git"
UFS_PHYSICS_URL = "https://github.com/ufs-community/ccpp-physics.git"
GENERATE_METADATA_HTML = Path(__file__).with_name("generate_metadata_html.py")


def run(*command, cwd=None):
    """Run a command while showing the reproducible command line."""
    print("+ " + shlex.join(map(str, command)), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def require_tools():
    missing = [
        command
        for command in ("git", "doxygen", "dot", "perl", "bibtex")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(
            "Missing build commands: "
            + ", ".join(missing)
            + ". See tools/doxygen/README.md for installation requirements."
        )


def require_paths(ref_name, root, paths):
    """Report an incompatible repository ref before starting Doxygen."""
    missing = [str(path) for path in paths if not (root / path).is_file()]
    if missing:
        raise SystemExit(
            f"{ref_name} does not contain files required by this build stage: "
            + ", ".join(missing)
        )


def reset_directory(directory):
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)


def move_contents(source, destination):
    destination.mkdir(parents=True, exist_ok=True)
    contents = list(source.iterdir()) if source.is_dir() else []
    if not contents:
        raise RuntimeError(f"No generated files found in {source}")
    for path in contents:
        shutil.move(str(path), destination / path.name)


def copy_metadata(metadata, html):
    tables = list(metadata.glob("*.html"))
    if not tables:
        raise RuntimeError(f"No metadata tables found in {metadata}")
    html.mkdir(parents=True, exist_ok=True)
    for table in tables:
        shutil.copy2(table, html / table.name)


def checkout_physics(physics, remote_url, ref):
    """Check out a ccpp-physics branch, tag, or commit with its submodules."""
    if not (physics / ".git").exists():
        if physics.exists():
            shutil.rmtree(physics)
        run("git", "init", "--quiet", physics)
    run("git", "fetch", "--depth", "1", remote_url, ref, cwd=physics)
    run("git", "checkout", "--force", "--detach", "FETCH_HEAD", cwd=physics)
    run("git", "clean", "-ffdx", "--quiet", cwd=physics)
    run("git", "submodule", "sync", "--recursive", cwd=physics)
    run(
        "git",
        "submodule",
        "update",
        "--init",
        "--recursive",
        "--force",
        "--depth",
        "1",
        cwd=physics,
    )


def build_section(physics, site, ref_name, doxyfile, section):
    """Build one documentation section from the checked-out Physics ref."""
    docs = physics / "physics/docs"
    require_paths(ref_name, docs, (doxyfile,))
    metadata = docs / "metadata_html"
    run(
        sys.executable,
        GENERATE_METADATA_HTML,
        "--physics",
        physics,
        "--output",
        metadata,
    )
    doc_html = docs / "doc/html"
    reset_directory(doc_html)
    copy_metadata(metadata, doc_html)
    run("doxygen", doxyfile, cwd=docs)
    move_contents(doc_html, site / section)


def build_landing_page(physics, site, ref_name):
    """Build the site landing page from the checked-out Physics ref."""
    docs = physics / "physics/docs"
    require_paths(ref_name, docs, ("mainpage_doxyfile",))
    landing = docs / "html/main"
    reset_directory(landing)
    run("doxygen", "mainpage_doxyfile", cwd=docs)
    move_contents(landing, site)


def build(workspace, site, ncar_physics_ref, v7_physics_ref, ufs_physics_ref):
    """Build every documentation section, one Physics ref at a time."""
    require_tools()
    physics = workspace / "physics_repo"
    if site.exists():
        shutil.rmtree(site)
    site.mkdir(parents=True)

    ncar = f"NCAR/ccpp-physics ref {ncar_physics_ref!r}"
    checkout_physics(physics, NCAR_PHYSICS_URL, ncar_physics_ref)
    for doxyfile, section in (
        ("ccpphsd_doxyfile", "HSD"),
        ("ccppgw_doxyfile", "GWv1"),
        ("ccpplandda_doxyfile", "LandDA"),
    ):
        build_section(physics, site, ncar, doxyfile, section)
    build_landing_page(physics, site, ncar)

    checkout_physics(physics, NCAR_PHYSICS_URL, v7_physics_ref)
    build_section(
        physics,
        site,
        f"NCAR/ccpp-physics ref {v7_physics_ref!r}",
        "ccpp_doxyfile",
        "V7",
    )

    checkout_physics(physics, UFS_PHYSICS_URL, ufs_physics_ref)
    build_section(
        physics,
        site,
        f"ufs-community/ccpp-physics ref {ufs_physics_ref!r}",
        "ccppsrw_doxyfile",
        "SRWv3",
    )

    index = site / "index.html"
    if not index.is_file():
        raise RuntimeError(f"The build did not produce {index}")
    print(f"Site built in {site}")


def main():
    default_workspace = Path(__file__).resolve().parent / ".build/doxygen"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ncar-physics-ref",
        default="main",
        help="NCAR/ccpp-physics ref for HSD, GWv1, LandDA, and the landing "
        "page (default: main)",
    )
    parser.add_argument(
        "--v7-physics-ref",
        default="v7.0.1",
        help="NCAR/ccpp-physics ref for the V7 docs (default: v7.0.1)",
    )
    parser.add_argument(
        "--ufs-physics-ref",
        default="release/srw-v3",
        help="ufs-community/ccpp-physics ref for the SRWv3 docs "
        "(default: release/srw-v3)",
    )
    parser.add_argument("--workspace", type=Path, default=default_workspace)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    workspace = arguments.workspace.resolve()
    site = arguments.output.resolve() if arguments.output else workspace / "site"
    build(
        workspace,
        site,
        arguments.ncar_physics_ref,
        arguments.v7_physics_ref,
        arguments.ufs_physics_ref,
    )


if __name__ == "__main__":
    main()
