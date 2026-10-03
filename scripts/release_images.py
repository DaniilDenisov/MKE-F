"""Build and optionally push traceable MKE-F Docker release images.

The release version comes from the newest CHANGELOG.md heading. Both images get
OCI version and revision labels and an additional git-<short SHA> tag.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHANGELOG_VERSION_RE = re.compile(r"(?m)^## \[([^]]+)]")
REVISION_LABEL = "org.opencontainers.image.revision"
VERSION_LABEL = "org.opencontainers.image.version"

IMAGES = (
    ("solver", "denisovds/mkef-solver", "docker/solver.Dockerfile", "runtime"),
    ("web", "denisovds/mkef-web", "docker/web.Dockerfile", None),
)


def run(command: list[str], *, dry_run: bool = False) -> None:
    print(" ".join(command))
    if not dry_run:
        subprocess.run(command, cwd=ROOT, check=True)


def output(command: list[str]) -> str:
    result = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def release_version() -> str:
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    match = CHANGELOG_VERSION_RE.search(changelog)
    if not match:
        raise SystemExit("CHANGELOG.md has no '## [version]' release heading")
    return match.group(1)


def require_clean_worktree() -> None:
    status = output(["git", "status", "--porcelain=v1"])
    if status:
        raise SystemExit(
            "Working tree is not clean. Commit or stash all changes before "
            "building a release."
        )


def require_release_compose_version(version: str) -> None:
    compose = (ROOT / "compose.release.yaml").read_text(encoding="utf-8")
    missing = [
        image
        for _, image, _, _ in IMAGES
        if f"image: {image}:{version}" not in compose
    ]
    if missing:
        joined = ", ".join(missing)
        raise SystemExit(
            f"compose.release.yaml does not reference version {version} for: {joined}"
        )


def require_pushed_head(revision: str) -> None:
    try:
        upstream = output(
            ["git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{upstream}"]
        )
        upstream_revision = output(["git", "rev-parse", upstream])
    except subprocess.CalledProcessError as error:
        raise SystemExit("Current branch has no readable upstream branch") from error
    if upstream_revision != revision:
        raise SystemExit(
            f"HEAD {revision} is not the commit at {upstream} ({upstream_revision}). "
            "Push the source commit before pushing its images."
        )


def ensure_git_tag(tag: str, revision: str, *, dry_run: bool, disabled: bool) -> None:
    if disabled:
        return
    existing = output(["git", "tag", "--list", tag])
    if existing:
        tagged_revision = output(["git", "rev-list", "-n", "1", tag])
        if tagged_revision != revision:
            raise SystemExit(
                f"Git tag {tag} points to {tagged_revision}, not HEAD {revision}"
            )
        print(f"Git tag {tag} already points to HEAD")
        return
    run(
        ["git", "tag", "-a", tag, "-m", f"Release {tag.removeprefix('v')}"],
        dry_run=dry_run,
    )


def inspect_label(image: str, label: str) -> str:
    template = f'{{{{ index .Config.Labels "{label}" }}}}'
    return output(["docker", "image", "inspect", image, "--format", template])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build traceable solver and web Docker images for the current release."
    )
    parser.add_argument("--push", action="store_true", help="Push all Docker tags")
    parser.add_argument(
        "--no-latest", action="store_true", help="Do not build or push the latest tag"
    )
    parser.add_argument(
        "--no-git-tag", action="store_true", help="Do not create the local v<version> tag"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print Docker and tag commands only"
    )
    args = parser.parse_args()

    require_clean_worktree()
    version = release_version()
    require_release_compose_version(version)
    revision = output(["git", "rev-parse", "--verify", "HEAD"])
    short_revision = revision[:12]

    if args.push:
        require_pushed_head(revision)

    git_tag = f"v{version}"
    ensure_git_tag(
        git_tag,
        revision,
        dry_run=args.dry_run,
        disabled=args.no_git_tag,
    )

    print(f"Release {version} from Git commit {revision}")
    built_images: list[tuple[str, list[str]]] = []

    for name, repository, dockerfile, target in IMAGES:
        tags = [f"{repository}:{version}", f"{repository}:git-{short_revision}"]
        if not args.no_latest:
            tags.append(f"{repository}:latest")

        command = [
            "docker",
            "build",
            "--pull",
            "--build-arg",
            f"MKEF_VERSION={version}",
            "--build-arg",
            f"VCS_REF={revision}",
            "-f",
            dockerfile,
        ]
        if target:
            command.extend(["--target", target])
        for tag in tags:
            command.extend(["-t", tag])
        command.append(".")

        print(f"Building {name} image")
        run(command, dry_run=args.dry_run)
        built_images.append((tags[0], tags))

        if not args.dry_run:
            actual_revision = inspect_label(tags[0], REVISION_LABEL)
            actual_version = inspect_label(tags[0], VERSION_LABEL)
            if actual_revision != revision or actual_version != version:
                raise SystemExit(
                    f"Image metadata mismatch for {tags[0]}: "
                    f"version={actual_version!r}, revision={actual_revision!r}"
                )

    if args.push:
        for _, tags in built_images:
            for tag in tags:
                run(["docker", "push", tag], dry_run=args.dry_run)

    print("Verified Docker metadata:" if not args.dry_run else "Planned metadata:")
    for primary_tag, _ in built_images:
        print(f"  {primary_tag}: version={version}, revision={revision}")
    if not args.no_git_tag:
        print(f"Push the Git tag with: git push origin {git_tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
