# Documentation maintenance

The PostGWAS user guide is maintained as canonical Markdown in the main
repository. The root [`README.md`](../../README.md) introduces installation,
available analyses and the first workflow. The [user-guide home](home.md)
contains the complete manifest-ordered documentation index, so all pages can
be browsed directly without GitHub Wiki access.

## Source contract

[`docs/wiki.yml`](../wiki.yml) is the documentation manifest. Its page order is
the order used by the user-guide index and optional generated navigation, and each
entry defines:

- the user-facing page title;
- a unique, stable Wiki slug;
- the canonical Markdown source under `docs/`;
- the optional sidebar section.

`docs/wiki/home.md` and `docs/wiki/footer.md` are source files. The special Wiki
filenames `Home.md`, `_Sidebar.md`, and `_Footer.md` exist only in generated
output. Existing canonical pages elsewhere under `docs/` can be published by
adding them to the manifest; they should not be copied into this directory.

## Validate sources

Run validation without writing generated files:

```console
python tools/docs/build_wiki.py --check
```

Validation rejects malformed manifest fields, unsafe or missing source paths,
duplicate slugs, broken local links, links to unpublished Markdown pages, and
local assets that do not yet have a defined publication rule.
It also checks the README's local destinations and the four stable navigation
sections for every published analysis guide, including guides outside `wiki/`.

```console
python tools/docs/validate_wiki_cli.py
```

The CLI checker inspects live help for accepted options/choices and usage
arguments also labelled `Required:`, and reuses the method-aware planner/registry for pipeline
resource requirements. Inline command names are navigation, not complete
recipes. Fenced CLI-only examples must supply the declared requirements;
examples using a placeholder `--run-config` require manual configuration review.
This check does not open references, validate every conditional input/schema,
execute analyses or establish that their outputs are correct.

Pull requests and pushes to the default branch run the same contract through
`.github/workflows/wiki-docs.yml`.

## Build an optional Wiki-compatible tree

Generate the publishable files locally:

```console
python tools/docs/build_wiki.py
```

Output is written to the ignored `build/wiki/` directory. The builder records
the files it owns in `.postgwas-wiki-files`; a later build removes only stale
files listed in that inventory and preserves unrelated files. The inventory is
part of the generated publication tree and must be retained in the Wiki Git
repository so automated rebuilds can distinguish managed files from unrelated
content. GitHub does not render it as a Wiki page.

An alternative output directory can be selected explicitly:

```console
python tools/docs/build_wiki.py --output /path/to/generated-wiki
```

The builder refuses to use an existing non-empty directory unless that
directory contains its managed-file inventory. This prevents accidental
replacement of a Wiki clone or another user-owned directory.

## Publication

GitHub Wiki publication is not configured for this private repository. Users
browse the canonical pages through the complete documentation index in the
user-guide home. The workflow validates the Markdown, generated navigation, links,
and documented commands, but does not attempt to push a `.wiki.git` repository.

If Wiki access is enabled in the future, the optional generated tree can be
reviewed and published as a separate change. The canonical Markdown and the
user-guide index must remain the source of truth.

## Add a page

1. Create or update the canonical Markdown file under `docs/`.
2. Start the page with one level-one heading.
3. Add the page to `docs/wiki.yml` with a stable, unique slug.
4. Use relative links between canonical Markdown sources.
5. Add its canonical source link to the complete index in `docs/wiki/home.md`
   in manifest order. Add or update its README catalogue entry when introducing
   a public analysis command.
6. Run the documentation validation tools listed above.
7. Inspect the source page and optional generated navigation.

For module documentation, begin with
[`docs/templates/module-page.md`](../templates/module-page.md). Configuration and
pipeline sections marked as generated must not be filled by copying defaults or
dependencies manually.

## Maintain the harmonisation policy reference

The detailed source/step contract and generated policy tables live in
[`docs/wiki/harmonisation/policies.md`](harmonisation/policies.md), not the
landing README. The generated block is derived from canonical YAML and schema
constraints; do not edit its rows by hand.

```console
python tools/docs/update_harmonisation_policies.py --check
```

After an intentional change to policy documentation or schema, run the same
command without `--check`, review the generated diff, and run the documentation
checks. The local tests retain complete sample-sheet field coverage and the seven-stage,
29-step processing-order contract in the detailed guide while checking the
README's introductory navigation separately.

## Local regression suite

`tests/` is deliberately excluded from Git and preserved only in maintainer
working copies. A fresh clone does not contain it. Maintainers with those files
can run `python -m pytest -q`; do not describe that as a command available from
every public checkout. Published CI uses the documentation tools above,
installer shell-syntax checks and the complete Linux installation job. Removing
the tracked tests reduces published regression coverage; these remaining checks
are not a replacement for the full suite or reference-backed analysis tests.
