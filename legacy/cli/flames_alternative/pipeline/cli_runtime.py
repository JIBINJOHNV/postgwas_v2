# postgwas/flames/clis/pipeline/cli_runtime.py


def add_runtime_arguments(parser):
    runtime = parser.add_argument_group("5. Runtime Options")

    runtime.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous partial results."
    )

    runtime.add_argument(
        "--dry_run",
        action="store_true",
        help="Validate configuration without executing."
    )

    return parser
