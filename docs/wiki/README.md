# Wiki maintenance

The public GitHub Wiki is a generated view of canonical Markdown maintained in
the main PostGWAS repository. Do not maintain an independent copy of the same
documentation in the Wiki repository.

## Source contract

[`docs/wiki.yml`](../wiki.yml) is the publication manifest. Its page order is
the sidebar order, and each entry defines:

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

Pull requests and pushes to the default branch run the same contract through
`.github/workflows/wiki-docs.yml`.

## Build the Wiki tree

Generate the publishable files locally:

```console
python tools/docs/build_wiki.py
```

Output is written to the ignored `build/wiki/` directory. The builder records
the Markdown files it owns in `.postgwas-wiki-files`; a later build removes only
stale files listed in that inventory and preserves unrelated files. Publish only
the generated Markdown files, not the inventory.

An alternative output directory can be selected explicitly:

```console
python tools/docs/build_wiki.py --output /path/to/generated-wiki
```

The builder refuses to use an existing non-empty directory unless that
directory contains its managed-file inventory. This prevents accidental
replacement of a Wiki clone or another user-owned directory.

## Publish

Wiki publication is deliberately separate from validation. After the repository
Wiki is enabled and its initial Home page exists:

1. Build the Wiki tree.
2. Clone the repository's `.wiki.git` repository into a separate directory.
3. Copy the generated Markdown files into that clone.
4. Review the Wiki repository diff.
5. Commit and push the Wiki repository.

Automated publication should run only after changes reach the default branch.
Use a repository-scoped write credential selected for the Wiki; do not add a
broad personal token to the workflow.

## Add a page

1. Create or update the canonical Markdown file under `docs/`.
2. Start the page with one level-one heading.
3. Add the page to `docs/wiki.yml` with a stable, unique slug.
4. Use relative links between canonical Markdown sources.
5. Run the Wiki validation tests.
6. Inspect the generated page and navigation before publication.

For module documentation, begin with
[`docs/templates/module-page.md`](../templates/module-page.md). Configuration and
pipeline sections marked as generated must not be filled by copying defaults or
dependencies manually.
