#!/usr/bin/env python3
"""Build a deterministic GitHub Wiki tree from canonical docs sources."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPOSITORY_ROOT / "docs" / "wiki.yml"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "build" / "wiki"
SLUG_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]*$")
MARKDOWN_LINK_PATTERN = re.compile(r"(!?)\[([^\]]+)\]\(([^)]+)\)")
MANAGED_FILE_INVENTORY = ".postgwas-wiki-files"
PACKAGED_DEFAULTS_DIRECTORY = REPOSITORY_ROOT / "src" / "postgwas" / "config" / "defaults"
CONFIGURATION_DEFAULTS_TOKEN = "<!-- GENERATED: PACKAGED CONFIGURATION DEFAULTS -->"
WIKI_IMAGE_SUFFIXES = frozenset(
    {".gif", ".jpeg", ".jpg", ".png", ".svg", ".webp"}
)
GENERATED_NOTICE = (
    "<!-- Generated from {source}. Edit the source file, not the Wiki copy. -->\n\n"
)


class WikiBuildError(ValueError):
    """Raised when the Wiki manifest or a canonical page is invalid."""


@dataclass(frozen=True)
class WikiPage:
    title: str
    slug: str
    source: Path
    source_relative: Path
    section: str | None = None


@dataclass(frozen=True)
class WikiSite:
    title: str
    footer_source: Path
    footer_source_relative: Path
    pages: tuple[WikiPage, ...]


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise WikiBuildError(f"{label} must be a mapping")
    return value


def _required_text(mapping: dict[str, Any], key: str, label: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise WikiBuildError(f"{label}.{key} must be non-empty text")
    return value.strip()


def _repository_file(value: str, label: str) -> tuple[Path, Path]:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise WikiBuildError(f"{label} must be a repository-relative path")
    resolved = (REPOSITORY_ROOT / relative).resolve()
    expected_root = (REPOSITORY_ROOT / "docs").resolve()
    if not resolved.is_relative_to(expected_root):
        raise WikiBuildError(f"{label} must remain under {expected_root}")
    if not resolved.is_file():
        raise WikiBuildError(f"{label} does not exist: {relative}")
    return resolved, relative


def load_manifest(path: Path = DEFAULT_MANIFEST) -> WikiSite:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise WikiBuildError(f"Cannot read Wiki manifest {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise WikiBuildError(f"Invalid Wiki manifest {path}: {exc}") from exc

    root = _mapping(document, "manifest")
    unknown_root = set(root) - {"wiki"}
    if unknown_root:
        raise WikiBuildError(f"Unknown manifest keys: {', '.join(sorted(unknown_root))}")
    wiki = _mapping(root.get("wiki"), "wiki")
    unknown_wiki = set(wiki) - {"title", "footer_source", "pages"}
    if unknown_wiki:
        raise WikiBuildError(f"Unknown wiki keys: {', '.join(sorted(unknown_wiki))}")

    title = _required_text(wiki, "title", "wiki")
    footer_value = _required_text(wiki, "footer_source", "wiki")
    footer_source, footer_relative = _repository_file(footer_value, "wiki.footer_source")
    raw_pages = wiki.get("pages")
    if not isinstance(raw_pages, list) or not raw_pages:
        raise WikiBuildError("wiki.pages must be a non-empty list")

    pages = []
    slugs: set[str] = set()
    sources: set[Path] = set()
    for index, raw_page in enumerate(raw_pages):
        label = f"wiki.pages[{index}]"
        page = _mapping(raw_page, label)
        unknown_page = set(page) - {"title", "slug", "source", "section"}
        if unknown_page:
            raise WikiBuildError(
                f"Unknown {label} keys: {', '.join(sorted(unknown_page))}"
            )
        page_title = _required_text(page, "title", label)
        slug = _required_text(page, "slug", label)
        if not SLUG_PATTERN.fullmatch(slug):
            raise WikiBuildError(
                f"{label}.slug must contain only letters, numbers, and hyphens: {slug}"
            )
        slug_key = slug.casefold()
        if slug_key in slugs:
            raise WikiBuildError(f"Duplicate Wiki slug: {slug}")
        slugs.add(slug_key)

        source_value = _required_text(page, "source", label)
        source, source_relative = _repository_file(source_value, f"{label}.source")
        if source.suffix.lower() != ".md":
            raise WikiBuildError(f"{label}.source must be a Markdown file")
        if source in sources:
            raise WikiBuildError(f"Canonical source is published more than once: {source_relative}")
        sources.add(source)

        section = page.get("section")
        if section is not None and (not isinstance(section, str) or not section.strip()):
            raise WikiBuildError(f"{label}.section must be non-empty text when provided")
        pages.append(
            WikiPage(
                title=page_title,
                slug=slug,
                source=source,
                source_relative=source_relative,
                section=section.strip() if isinstance(section, str) else None,
            )
        )

    if pages[0].slug != "Home":
        raise WikiBuildError("The first Wiki page must use the slug 'Home'")
    return WikiSite(title, footer_source, footer_relative, tuple(pages))


def _split_link_target(target: str) -> tuple[str, str]:
    path, marker, anchor = target.partition("#")
    return path, f"#{anchor}" if marker else ""


def _yaml_mapping(path: Path) -> dict[str, Any]:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise WikiBuildError(f"Cannot read canonical configuration {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise WikiBuildError(f"Invalid canonical configuration {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise WikiBuildError(f"Canonical configuration must be a mapping: {path}")
    return document


def render_packaged_configuration_defaults() -> str:
    """Render the same packaged-default hierarchy assembled by config.loader."""
    application = _yaml_mapping(PACKAGED_DEFAULTS_DIRECTORY / "application.yaml")
    application["pipeline"] = _yaml_mapping(
        PACKAGED_DEFAULTS_DIRECTORY / "pipeline.yaml"
    )
    application["resources"] = _yaml_mapping(
        PACKAGED_DEFAULTS_DIRECTORY / "resources.yaml"
    )

    modules: dict[str, Any] = {}
    module_directory = PACKAGED_DEFAULTS_DIRECTORY / "modules"
    for path in sorted(module_directory.glob("*.yaml")):
        document = _yaml_mapping(path)
        module = document.get("module", document)
        if not isinstance(module, dict):
            raise WikiBuildError(
                f"Packaged module configuration must be a mapping: {path}"
            )
        modules[path.stem] = module
    if not modules:
        raise WikiBuildError(f"No packaged module defaults found in {module_directory}")
    application["modules"] = modules

    rendered = yaml.safe_dump(
        application,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
    ).rstrip()
    return (
        "Current packaged settings:\n\n"
        "```yaml\n"
        f"{rendered}\n"
        "```"
    )


def render_page(
    page: WikiPage,
    source_to_slug: dict[Path, str],
    assets: dict[str, Path],
) -> str:
    text = page.source.read_text(encoding="utf-8")
    if not text.startswith("# "):
        raise WikiBuildError(
            f"Canonical Wiki page must start with an H1: {page.source_relative}"
        )
    if text.count(CONFIGURATION_DEFAULTS_TOKEN) > 1:
        raise WikiBuildError(
            f"Configuration-defaults token appears more than once: {page.source_relative}"
        )
    if CONFIGURATION_DEFAULTS_TOKEN in text:
        text = text.replace(
            CONFIGURATION_DEFAULTS_TOKEN,
            render_packaged_configuration_defaults(),
        )

    def replace_link(match: re.Match[str]) -> str:
        image_marker, label, target = match.groups()
        if target.startswith(("#", "http://", "https://", "mailto:")):
            return match.group(0)
        raw_path, anchor = _split_link_target(target)
        if not raw_path:
            return match.group(0)
        local_target = (page.source.parent / raw_path).resolve()
        if image_marker:
            docs_root = (REPOSITORY_ROOT / "docs").resolve()
            if (
                not local_target.is_relative_to(docs_root)
                or not local_target.is_file()
            ):
                raise WikiBuildError(
                    f"Broken local image in {page.source_relative}: {target}"
                )
            if local_target.suffix.lower() not in WIKI_IMAGE_SUFFIXES:
                raise WikiBuildError(
                    f"Unsupported Wiki image type in {page.source_relative}: {target}"
                )
            published_name = local_target.name
            previous_source = assets.get(published_name)
            if previous_source is not None and previous_source != local_target:
                raise WikiBuildError(
                    f"Wiki image filename collision: {previous_source} and {local_target}"
                )
            assets[published_name] = local_target
            return f"![{label}]({published_name}{anchor})"
        if local_target.suffix.lower() != ".md":
            if local_target.exists():
                raise WikiBuildError(
                    f"Wiki non-image asset publication is not implemented: "
                    f"{page.source_relative} -> {target}"
                )
            raise WikiBuildError(
                f"Broken local link in {page.source_relative}: {target}"
            )
        slug = source_to_slug.get(local_target)
        if slug is None:
            raise WikiBuildError(
                f"Local Markdown link targets an unpublished page: "
                f"{page.source_relative} -> {target}"
            )
        return f"[{label}]({slug}{anchor})"

    rendered = MARKDOWN_LINK_PATTERN.sub(replace_link, text).rstrip() + "\n"
    return GENERATED_NOTICE.format(source=page.source_relative.as_posix()) + rendered


def render_sidebar(site: WikiSite) -> str:
    sections: OrderedDict[str, list[WikiPage]] = OrderedDict()
    unsectioned = []
    for page in site.pages:
        if page.section is None:
            unsectioned.append(page)
        else:
            sections.setdefault(page.section, []).append(page)

    lines = [f"# {site.title}", ""]
    lines.extend(f"- [{page.title}]({page.slug})" for page in unsectioned)
    for section, pages in sections.items():
        if lines[-1]:
            lines.append("")
        lines.append(f"## {section}")
        lines.append("")
        lines.extend(f"- [{page.title}]({page.slug})" for page in pages)
    return "\n".join(lines).rstrip() + "\n"


def render_wiki(site: WikiSite) -> dict[str, str | bytes]:
    source_to_slug = {page.source.resolve(): page.slug for page in site.pages}
    assets: dict[str, Path] = {}
    rendered = {
        f"{page.slug}.md": render_page(page, source_to_slug, assets)
        for page in site.pages
    }
    rendered["_Sidebar.md"] = render_sidebar(site)
    footer = site.footer_source.read_text(encoding="utf-8").rstrip() + "\n"
    rendered["_Footer.md"] = GENERATED_NOTICE.format(
        source=site.footer_source_relative.as_posix()
    ) + footer
    for filename, source in assets.items():
        if filename in rendered or filename == MANAGED_FILE_INVENTORY:
            raise WikiBuildError(
                f"Wiki image conflicts with a generated file: {filename}"
            )
        rendered[filename] = source.read_bytes()
    return rendered


def write_wiki(rendered: dict[str, str | bytes], output: Path) -> None:
    output = output.resolve()
    if output.exists() and not output.is_dir():
        raise WikiBuildError(f"Wiki output path is not a directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    inventory = output / MANAGED_FILE_INVENTORY
    previous_files: set[str] = set()
    if inventory.exists():
        for filename in inventory.read_text(encoding="utf-8").splitlines():
            if not filename or Path(filename).name != filename:
                raise WikiBuildError(f"Invalid managed-file inventory entry: {filename!r}")
            previous_files.add(filename)
    elif any(output.iterdir()):
        raise WikiBuildError(
            f"Wiki output directory is not empty and is not managed by this tool: {output}"
        )

    for filename in rendered:
        if (output / filename).exists() and filename not in previous_files:
            raise WikiBuildError(
                f"Refusing to replace unmanaged output file: {output / filename}"
            )

    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        for filename, content in rendered.items():
            target = temporary / filename
            if isinstance(content, bytes):
                target.write_bytes(content)
            else:
                target.write_text(content, encoding="utf-8")
        (temporary / MANAGED_FILE_INVENTORY).write_text(
            "\n".join(sorted(rendered)) + "\n", encoding="utf-8"
        )
        for filename in previous_files - set(rendered):
            stale_file = output / filename
            if stale_file.is_file():
                stale_file.unlink()
        for filename in rendered:
            (temporary / filename).replace(output / filename)
        (temporary / MANAGED_FILE_INVENTORY).replace(inventory)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build or validate the generated PostGWAS GitHub Wiki tree."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Wiki manifest (default: docs/wiki.yml).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Generated Wiki directory (default: build/wiki).",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate all sources and links without writing generated files.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        site = load_manifest(args.manifest)
        rendered = render_wiki(site)
        if not args.check:
            write_wiki(rendered, args.output)
    except WikiBuildError as exc:
        raise SystemExit(f"Wiki build failed: {exc}") from exc
    action = "validated" if args.check else f"built at {args.output}"
    print(f"Wiki {action}: {len(site.pages)} pages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
