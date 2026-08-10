# Harmonisation test data

`manifest_v2.csv` uses the canonical version-2 manifest names. Its ADHD input
uses the exact absolute fixture path requested for this local test workspace.
The manifest loader also supports relative paths for portable user manifests.

For quantitative traits, `control_count` or `control_count_column` represents
total N and all case-count fields must be empty. For case-control traits, both
control and case counts are required.

The compressed ADHD fixture is retained for header validation and scientific
regression tests. It is test input, not generated output.
