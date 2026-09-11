# CALDERA runtime validation

CALDERA preflight uses the same shared R runtime resolver as SuSiE. It discovers
the configured Rscript's actual `.Library` and `.libPaths()`, puts that native
library first in the child process's `R_LIBS` and `R_LIBS_USER`, and retains other
discovered libraries afterwards. The parent process and user settings are not
modified, and no packages are installed automatically.

The required namespaces `data.table` and `dplyr` are native software-contract
invariants from the pinned upstream `caldera()` implementation, not adjustable
model settings. Missing or incompatible required packages fail preflight before
the model runs. Both preflight and the checked CALDERA child use `--vanilla`, so
site/user profiles or environment files cannot change the validated library
selection during startup. The validated environment is passed to the R child;
nonzero native exits and missing outputs remain fatal. The canonical log records
R version, required package versions and the library search order, not the full
inherited environment. Pipeline resource reuse carries the validated runtime
alongside unchanged reference identities. No gene scores, credible sets, model
parameters, references or significance thresholds are altered by this correction.

Sources:

- [R library search-path documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html)
- [R startup and `--vanilla` documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html)
- [Pinned CALDERA implementation](https://github.com/kheilbron/caldera/blob/81a8a0308741ae986660f711bbf6a7abbd3bca19/z_caldera.R)
