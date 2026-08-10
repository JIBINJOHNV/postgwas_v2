from postgwas.config.models.common import StrictModel


class PipelineConfig(StrictModel):
    modules: list[str]
    fail_fast: bool
    stop_after: str | None = None
