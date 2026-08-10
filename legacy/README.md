# Legacy material

This directory preserves historical, duplicate, experimental, or incomplete
PostGWAS files that are intentionally excluded from the installable
`postgwas` package.

Files here are retained for comparison and recovery only. They are not public
APIs, pipeline targets, supported command-line entry points, or test fixtures.
Names include the original module or purpose so similarly named historical
files remain distinguishable.

## Contents

- `backups/`: superseded implementations kept for reference.
- `cli/`: duplicate module/pipeline CLI implementations.
- `config/`: deprecated configuration interfaces and examples.
- `containers/`: incomplete or superseded container definitions.
- `dependency_specs/`: module-local dependency snapshots replaced by project metadata.
- `examples/`: historical commands containing environment-specific paths.
- `experimental/`: providers and resources not integrated into production services.
- `prototypes/`: incomplete or internally unreferenced utilities.
- `test_data/`: large historical inputs and generated outputs.
- `vendor/`: architecture-specific external binaries not distributed with Python.

Production code must never import from this directory.
