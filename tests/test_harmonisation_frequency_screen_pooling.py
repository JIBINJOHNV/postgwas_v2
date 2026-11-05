"""Selected-frequency screens pool GWAS-row counts rather than shard fractions."""

import polars as pl
import pytest

from postgwas.modules.harmonisation.shared import statistics


def test_unequal_shards_use_pooled_counts_instead_of_mean_percentages():
    screens = [
        statistics.frequency_maf_screen(pl.DataFrame({"AF": [.1] * 94 + [.9] * 6}), "AF", .95),
        statistics.frequency_maf_screen(pl.DataFrame({"AF": [.1]}), "AF", .95),
    ]

    pooled = statistics.pool_frequency_maf_screens(iter(screens), .95)

    assert pooled["rows"] == pooled["usable"] == 101
    assert pooled["unusable"] == 0
    assert pooled["at_or_below_0_5"] == 95
    assert pooled["low_fraction_of_usable"] == pytest.approx(95 / 101)
    assert pooled["maf_like"] is False
    assert pooled["column_present"] is True


@pytest.mark.parametrize("low_count, expected", [(190, False), (191, True)])
def test_pooled_screen_uses_strict_cutoff_with_half_counted_as_low(low_count, expected):
    frame = pl.DataFrame({"AF": [.5] * low_count + [1.] * (200 - low_count)})
    screens = [statistics.frequency_maf_screen(frame.slice(start, 100), "AF", .1)
               for start in (0, 100)]

    pooled = statistics.pool_frequency_maf_screens(screens, .95)

    assert pooled["at_or_below_0_5"] == low_count
    assert pooled["usable"] == 200
    assert pooled["decision_cutoff"] == .95
    assert pooled["maf_like"] is expected


def test_pooled_screen_keeps_endpoints_and_excludes_invalid_or_missing_values():
    frame = pl.DataFrame({"AF": ["0", "0.5", "1", None, "nonnumeric", "NaN", "inf", "-0.1", "1.1"]})
    screens = [statistics.frequency_maf_screen(frame.slice(0, 4), "AF", .65),
               statistics.frequency_maf_screen(frame.slice(4), "AF", .65)]

    pooled = statistics.pool_frequency_maf_screens(screens, .65)

    assert pooled == {
        "column_present": True, "rows": 9, "usable": 3, "unusable": 6,
        "at_or_below_0_5": 2, "low_fraction_of_usable": 2 / 3,
        "decision_cutoff": .65, "maf_like": True,
    }


def test_empty_pool_does_not_claim_a_present_or_usable_frequency_column():
    pooled = statistics.pool_frequency_maf_screens(iter(()), .95)

    assert pooled == {
        "column_present": False, "rows": 0, "usable": 0, "unusable": 0,
        "at_or_below_0_5": 0, "low_fraction_of_usable": 0.,
        "decision_cutoff": .95, "maf_like": False,
    }


def test_missing_shard_column_is_reported_without_losing_other_shard_counts():
    screens = [statistics.frequency_maf_screen(pl.DataFrame({"OTHER": [1, 2]}), "AF", .95),
               statistics.frequency_maf_screen(pl.DataFrame({"AF": [1.]}), "AF", .95)]

    pooled = statistics.pool_frequency_maf_screens(screens, .95)

    assert pooled["column_present"] is False
    assert pooled["rows"] == 3
    assert pooled["usable"] == 1
    assert pooled["unusable"] == 2
    assert pooled["low_fraction_of_usable"] == 0.
    assert pooled["maf_like"] is False


def test_zero_usable_values_remain_explicit_for_the_dataset_stop_decision():
    screens = [statistics.frequency_maf_screen(pl.DataFrame({"AF": [None, "nan", "bad"]}), "AF", .95)]

    pooled = statistics.pool_frequency_maf_screens(screens, .95)

    assert pooled["column_present"] is True
    assert pooled["usable"] == 0
    assert pooled["unusable"] == 3
    assert pooled["maf_like"] is False
