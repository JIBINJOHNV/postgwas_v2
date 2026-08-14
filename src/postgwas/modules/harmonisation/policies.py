"""
policies.py — the single source of truth for every configurable value in the
harmonisation pipeline.

Governing principle (plan v3, Part 6):

  * every value the pipeline uses is a policy,
  * every policy default is loaded from the canonical harmonisation YAML,
  * every policy carries a plain-English explanation used by canonical config
    export, while run logs record only the resolved value a step consulted.

Where the plan recommends a different value (the rows marked with a warning
triangle), the option is implemented but the CURRENT value stays the default and
the recommendation is spelled out in the help text.

Public API
----------
    POLICIES                       frozen registry: dotted key -> Policy
    Policy                         one registry entry (frozen dataclass)
    PolicyError                    raised with ALL problems at once
    load_policies(mapping=None)    -> Policies
    default_policies()             -> Policies  (same as load_policies({}))
    describe(keys, policies=None)  -> [(key, value, help), ...]
    format_value(value)            -> str, the rendering used in the log

    Policies.get('eaf.out_of_range')      dotted-key access
    Policies.eaf.out_of_range             attribute access
    Policies.help('eaf.out_of_range')     plain-English text
    Policies.describe([...])              values plus config-export explanations
    Policies.resolved_dict()              flat dict for the run manifest

This module imports nothing from postgwas.* and nothing outside the standard
library, so it is importable and testable on its own.

Python 3.8 compatible.
"""

import copy
import dataclasses
import difflib
import os
import re
import types
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

__all__ = [
    "POLICIES",
    "Policy",
    "PolicyError",
    "PolicyProblem",
    "Policies",
    "load_policies",
    "default_policies",
    "describe",
    "format_value",
    "policy_keys",
    "group_keys",
]


# =============================================================================
# 1. Validators
# =============================================================================

_TRUE_WORDS = ("true", "yes", "on", "y", "t", "1")
_FALSE_WORDS = ("false", "no", "off", "n", "f", "0")


class Validator(object):
    """Base validator.  `coerce` normalises, `check` returns an error string or
    None, `describe` explains the constraint in plain English."""

    kind = "any"

    def coerce(self, value):
        return value

    def check(self, value):
        return None

    def describe(self):
        return "any value"

    def __repr__(self):
        return "<%s %s>" % (type(self).__name__, self.describe())


class AnyValue(Validator):
    kind = "any"


class Enum(Validator):
    """One of a fixed set of members.  Matching is case-insensitive for strings
    so that `CLIP` in a YAML file resolves to `clip`."""

    kind = "enum"

    def __init__(self, members, allow_none=False):
        self.members = tuple(members)
        self.allow_none = bool(allow_none)

    def coerce(self, value):
        if value is None:
            return None
        if isinstance(value, str):
            stripped = value.strip()
            if self.allow_none and stripped.lower() in ("none", "null", "disabled", ""):
                return None
            for member in self.members:
                if isinstance(member, str) and member.lower() == stripped.lower():
                    return member
            return stripped
        return value

    def check(self, value):
        if value is None:
            if self.allow_none:
                return None
            return "must be one of %s (got nothing)" % self._members_text()
        if value in self.members:
            return None
        suggestion = ""
        if isinstance(value, str):
            close = difflib.get_close_matches(
                value.lower(), [str(m).lower() for m in self.members], n=1, cutoff=0.6
            )
            if close:
                suggestion = "  Did you mean '%s'?" % close[0]
        return "must be one of %s (got %r)%s" % (self._members_text(), value, suggestion)

    def _members_text(self):
        return ", ".join("'%s'" % m for m in self.members)

    def describe(self):
        text = "one of " + self._members_text()
        if self.allow_none:
            text += ", or null"
        return text


class Flag(Validator):
    """A true/false switch."""

    kind = "bool"

    def coerce(self, value):
        if isinstance(value, str):
            word = value.strip().lower()
            if word in _TRUE_WORDS:
                return True
            if word in _FALSE_WORDS:
                return False
        return value

    def check(self, value):
        if isinstance(value, bool):
            return None
        return "must be true or false (got %r)" % (value,)

    def describe(self):
        return "true or false"


class Num(Validator):
    """A number, optionally an integer, optionally bounded."""

    def __init__(
        self,
        minimum=None,
        maximum=None,
        integer=False,
        allow_none=False,
        exclusive_minimum=False,
    ):
        self.minimum = minimum
        self.maximum = maximum
        self.integer = bool(integer)
        self.allow_none = bool(allow_none)
        self.exclusive_minimum = bool(exclusive_minimum)
        self.kind = "int" if integer else "float"

    def coerce(self, value):
        if isinstance(value, bool):
            return value
        if value is None:
            return None
        if isinstance(value, str):
            word = value.strip()
            if self.allow_none and word.lower() in ("none", "null", "disabled", ""):
                return None
            try:
                if self.integer:
                    return int(word, 10)
                return float(word)
            except ValueError:
                return value
        if self.integer and isinstance(value, float) and float(value).is_integer():
            return int(value)
        if (not self.integer) and isinstance(value, int):
            return float(value)
        return value

    def check(self, value):
        if value is None:
            if self.allow_none:
                return None
            return "must be a number (got nothing)"
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return "must be %s (got %r)" % (
                "a whole number" if self.integer else "a number",
                value,
            )
        if self.integer and not float(value).is_integer():
            return "must be a whole number (got %r)" % (value,)
        if self.minimum is not None:
            if self.exclusive_minimum and value <= self.minimum:
                return "must be greater than %s (got %r)" % (self.minimum, value)
            if (not self.exclusive_minimum) and value < self.minimum:
                return "must be at least %s (got %r)" % (self.minimum, value)
        if self.maximum is not None and value > self.maximum:
            return "must be at most %s (got %r)" % (self.maximum, value)
        return None

    def describe(self):
        what = "whole number" if self.integer else "number"
        if self.minimum is not None and self.maximum is not None:
            span = "%s between %s and %s" % (what, self.minimum, self.maximum)
        elif self.minimum is not None:
            span = "%s of at least %s" % (what, self.minimum)
        elif self.maximum is not None:
            span = "%s of at most %s" % (what, self.maximum)
        else:
            span = "any %s" % what
        if self.allow_none:
            span += ", or null"
        return span


class Text(Validator):
    """A string, optionally constrained by a regular expression."""

    kind = "str"

    def __init__(self, pattern=None, allow_none=False, description=None):
        self.pattern = re.compile(pattern) if pattern else None
        self.pattern_text = pattern
        self.allow_none = bool(allow_none)
        self.description = description

    def coerce(self, value):
        if value is None:
            return None
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return str(value)
        if isinstance(value, str):
            word = value.strip()
            if self.allow_none and word.lower() in ("none", "null", ""):
                return None
            return word
        return value

    def check(self, value):
        if value is None:
            if self.allow_none:
                return None
            return "must be text (got nothing)"
        if not isinstance(value, str):
            return "must be text (got %r)" % (value,)
        if self.pattern is not None and not self.pattern.match(value):
            return "must look like %s (got %r)" % (
                self.description or self.pattern_text,
                value,
            )
        return None

    def describe(self):
        if self.description:
            return self.description
        if self.pattern_text:
            return "text matching %s" % self.pattern_text
        return "any text"


class ListOf(Validator):
    """A list.  A comma-separated string is accepted and split."""

    kind = "list"

    def __init__(self, item, allow_empty=True, unique=True, upper=False,
                 nested_alternatives=False):
        self.item = item
        self.allow_empty = bool(allow_empty)
        self.unique = bool(unique)
        self.upper = bool(upper)
        # When true a nested list is allowed and means "at least one of these",
        # e.g. [chr, pos, [beta, zscore]] = chr and pos and (beta or zscore).
        self.nested_alternatives = bool(nested_alternatives)

    def coerce(self, value):
        if value is None:
            return None
        if isinstance(value, str):
            parts = [p.strip() for p in value.split(",")]
            value = [p for p in parts if p != ""]
        if isinstance(value, tuple):
            value = list(value)
        if isinstance(value, list):
            out = []
            for item in value:
                if self.nested_alternatives and isinstance(item, (list, tuple)):
                    inner = []
                    for member in item:
                        member = self.item.coerce(member)
                        if self.upper and isinstance(member, str):
                            member = member.upper()
                        inner.append(member)
                    out.append(inner)
                    continue
                item = self.item.coerce(item)
                if self.upper and isinstance(item, str):
                    item = item.upper()
                out.append(item)
            return out
        return value

    def check(self, value):
        if not isinstance(value, list):
            return "must be a list (got %r)" % (value,)
        if not value and not self.allow_empty:
            return "must not be empty"
        problems = []
        flat = []
        for item in value:
            if self.nested_alternatives and isinstance(item, (list, tuple)):
                if not item:
                    problems.append("an alternative group must not be empty")
                    continue
                for member in item:
                    message = self.item.check(member)
                    if message:
                        problems.append("item %r %s" % (member, message))
                    flat.append(member)
                continue
            message = self.item.check(item)
            if message:
                problems.append("item %r %s" % (item, message))
            flat.append(item)
        if self.unique and len(set(map(_hashable, flat))) != len(flat):
            problems.append("contains repeated entries")
        if problems:
            return "; ".join(problems)
        return None

    def describe(self):
        text = "a list where each %s" % self.item.describe()
        if self.nested_alternatives:
            text += "; a nested list means at least one of its members"
        return text


class ChromosomeSet(ListOf):
    """A list of chromosome labels.  Accepts the `1-22,X,Y,MT` shorthand."""

    def __init__(self, members):
        ListOf.__init__(self, Enum(members), allow_empty=False, unique=True, upper=True)

    def coerce(self, value):
        if isinstance(value, str):
            value = [p.strip() for p in value.split(",") if p.strip() != ""]
        if isinstance(value, (list, tuple)):
            expanded = []
            for item in value:
                text = str(item).strip().upper()
                match = re.match(r"^(\d+)\s*-\s*(\d+)$", text)
                if match:
                    low, high = int(match.group(1)), int(match.group(2))
                    if low <= high:
                        expanded.extend(str(n) for n in range(low, high + 1))
                        continue
                expanded.append(text)
            return expanded
        return value


class StringMap(Validator):
    """A mapping of text to text.  Accepts `23->X,24->Y` and `23:X,24:Y`."""

    kind = "dict"

    def __init__(self, value_members=None):
        self.value_members = tuple(value_members) if value_members else None

    def coerce(self, value):
        if isinstance(value, str):
            pairs = {}
            for chunk in value.replace("→", "->").split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                for arrow in ("->", "=>", ":", "="):
                    if arrow in chunk:
                        left, right = chunk.split(arrow, 1)
                        pairs[left.strip().upper()] = right.strip().upper()
                        break
                else:
                    pairs[chunk] = chunk
            return pairs
        if isinstance(value, dict):
            return dict(
                (str(k).strip().upper(), str(v).strip().upper())
                for k, v in value.items()
            )
        return value

    def check(self, value):
        if not isinstance(value, dict):
            return "must be a mapping such as {23: X, 24: Y} (got %r)" % (value,)
        problems = []
        for key, item in value.items():
            if not isinstance(key, str) or not isinstance(item, str):
                problems.append("entry %r -> %r must be text on both sides" % (key, item))
            elif self.value_members and item not in self.value_members:
                problems.append(
                    "entry %r -> %r: the new name must be one of %s"
                    % (key, item, ", ".join("'%s'" % m for m in self.value_members))
                )
        if problems:
            return "; ".join(problems)
        return None

    def describe(self):
        return "a mapping of old name to new name"


class NumericPair(Validator):
    """Two numbers, low then high.  `disabled` / null turns the check off."""

    kind = "list"

    def __init__(self, minimum=None, maximum=None, allow_none=True):
        self.minimum = minimum
        self.maximum = maximum
        self.allow_none = bool(allow_none)

    def coerce(self, value):
        if value is None:
            return None
        if isinstance(value, str):
            word = value.strip().lower()
            if word in ("disabled", "off", "none", "null", ""):
                return None
            parts = [p.strip() for p in value.replace("-", ",").split(",") if p.strip()]
            value = parts
        if isinstance(value, (list, tuple)):
            out = []
            for item in value:
                try:
                    out.append(float(item))
                except (TypeError, ValueError):
                    out.append(item)
            return out
        return value

    def check(self, value):
        if value is None:
            if self.allow_none:
                return None
            return "must be two numbers"
        if not isinstance(value, list) or len(value) != 2:
            return "must be two numbers, low then high, or 'disabled' (got %r)" % (value,)
        for item in value:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                return "must be two numbers, low then high (got %r)" % (value,)
        if value[0] >= value[1]:
            return "the first number must be smaller than the second (got %r)" % (value,)
        if self.minimum is not None and value[0] < self.minimum:
            return "the low end must be at least %s (got %r)" % (self.minimum, value)
        if self.maximum is not None and value[1] > self.maximum:
            return "the high end must be at most %s (got %r)" % (self.maximum, value)
        return None

    def describe(self):
        return "two numbers, low then high, or 'disabled'"


def _hashable(value):
    if isinstance(value, list):
        return tuple(value)
    if isinstance(value, dict):
        return tuple(sorted(value.items()))
    return value


# =============================================================================
# 2. Registry entry
# =============================================================================


@dataclasses.dataclass(frozen=True)
class Policy(object):
    """One entry in the registry.

    key            dotted key, e.g. 'eaf.out_of_range'
    default        the value used when the run does not provide an override
    type           'enum' | 'bool' | 'int' | 'float' | 'str' | 'list' | 'dict'
    validator      enum members or numeric range, see the Validator classes
    help           plain English; the logger prints this verbatim
    origin         the current implementation that consumes the value
    recommendation the plan's recommended value when it differs from the default
                   (always also stated inside `help`)
    """

    key: str
    default: Any
    type: str
    validator: Validator
    help: str
    origin: str = ""
    recommendation: Optional[str] = None

    @property
    def group(self):
        return self.key.split(".", 1)[0]

    @property
    def name(self):
        return self.key.split(".", 1)[1] if "." in self.key else self.key

    def default_copy(self):
        return copy.deepcopy(self.default)


_REGISTRY = []  # type: List[Policy]


# -----------------------------------------------------------------------------
# 6.1  Input reading
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# The registry lives with all other packaged configuration in
# config/defaults/modules/harmonisation.yaml, NOT here.
#
# It is data: each entry has a default, a type, a validator specification, the
# plain-English help text, and its consuming operation. Keeping it
# out of Python means a user can read it, diff it between releases and copy it
# as a starting config without opening any source.  This module only loads and
# validates it.
# -----------------------------------------------------------------------------

DEFAULTS_RELATIVE_PATH = os.path.join(
    "..", "..", "config", "defaults", "modules", "harmonisation.yaml"
)

#: Keys the registry file marks as the ones users most often change.
COMMON_KEYS = ()  # type: Tuple[str, ...]

#: Groups in the order the pipeline reaches them, each with the stage it
#: applies at.  Drives this file and every generated template.
GROUP_ORDER = ()   # type: Tuple[Tuple[str, str], ...]

#: Within a group, the order the owning step actually reads the settings.
KEY_ORDER = {}     # type: Dict[str, List[str]]

#: Every known field: whether it is required at read and at export, and which
#: step supplies it otherwise.  Explains why a field is absent from
#: columns.mandatory.  Loaded from the registry's field_lifecycle block.
FIELD_LIFECYCLE = {}  # type: Dict[str, Dict[str, Any]]

#: What a dataset must supply, expressed in config column names.  Loaded from
#: the ``required_inputs:`` block of the registry file; see validator.py.
REQUIRED_INPUTS = {}  # type: Dict[str, Dict[str, Any]]


class RegistryError(RuntimeError):
    """The shipped registry file is missing or malformed."""


def _normalise_help(text):
    """Tidy whitespace inside each line but keep the line structure."""
    lines = []
    for line in str(text).split("\n"):
        indent = line[: len(line) - len(line.lstrip())]
        lines.append((indent + " ".join(line.split())).rstrip())
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines)


def _validator_from_spec(spec, key):
    """Rebuild a Validator from its YAML description."""
    if not isinstance(spec, dict):
        raise RegistryError("%s: validator must be a mapping, got %r" % (key, spec))
    kind = spec.get("kind")
    if kind == "enum":
        return Enum(tuple(spec.get("members") or ()), bool(spec.get("allow_none", False)))
    if kind == "flag":
        return Flag()
    if kind == "num":
        return Num(
            minimum=spec.get("minimum"),
            maximum=spec.get("maximum"),
            integer=bool(spec.get("integer", False)),
            allow_none=bool(spec.get("allow_none", False)),
            exclusive_minimum=bool(spec.get("exclusive_minimum", False)),
        )
    if kind == "text":
        return Text(
            pattern=spec.get("pattern"),
            allow_none=bool(spec.get("allow_none", False)),
            description=spec.get("description"),
        )
    if kind == "chromosome_set":
        return ChromosomeSet(tuple(spec.get("members") or ()))
    if kind == "list":
        return ListOf(
            _validator_from_spec(spec.get("item") or {"kind": "any"}, key),
            allow_empty=bool(spec.get("allow_empty", True)),
            unique=bool(spec.get("unique", True)),
            upper=bool(spec.get("upper", False)),
            nested_alternatives=bool(spec.get("nested_alternatives", False)),
        )
    if kind == "numeric_pair":
        return NumericPair(
            minimum=spec.get("minimum"),
            maximum=spec.get("maximum"),
            allow_none=bool(spec.get("allow_none", True)),
        )
    if kind == "string_map":
        members = spec.get("value_members")
        return StringMap(tuple(members) if members else None)
    if kind == "any":
        return AnyValue()
    raise RegistryError("%s: unknown validator kind %r" % (key, kind))


def registry_path():
    """Return the one packaged harmonisation policy registry."""
    return os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), DEFAULTS_RELATIVE_PATH)
    )


def load_registry(path=None):
    """Read the canonical harmonisation YAML and build the policy registry."""
    location = path or registry_path()
    try:
        import yaml
    except ImportError:                                     # pragma: no cover
        raise RegistryError(
            "PyYAML is required to read the policy registry (%s)." % location
        )
    if not os.path.exists(location):
        raise RegistryError(
            "The policy registry is missing.\n"
            "  Expected at: %s\n"
            "  This file ships in config/defaults/modules and holds every default. If "
            "postgwas was installed as a package, it was probably left out of "
            "package_data. Reinstall PostGWAS with its packaged YAML files."
            % location
        )
    with open(location) as handle:
        document = yaml.safe_load(handle) or {}
    block = document.get("policies")
    if not isinstance(block, dict):
        raise RegistryError(
            "%s has no 'policies:' mapping at the top level." % location
        )
    built = []
    for group in sorted(block):
        entries = block[group]
        if not isinstance(entries, dict):
            raise RegistryError("%s: group %r must be a mapping." % (location, group))
        for leaf in sorted(entries):
            entry = entries[leaf]
            key = "%s.%s" % (group, leaf)
            if not isinstance(entry, dict):
                raise RegistryError("%s: %s must be a mapping." % (location, key))
            for required in ("default", "type", "validator", "help"):
                if required not in entry:
                    raise RegistryError(
                        "%s: %s is missing '%s'." % (location, key, required)
                    )
            built.append(
                Policy(
                    key=key,
                    default=entry["default"],
                    type=entry["type"],
                    validator=_validator_from_spec(entry["validator"], key),
                    help=_normalise_help(entry["help"]),
                    origin=entry.get("origin", entry.get("today_at")) or "",
                    recommendation=entry.get("recommendation"),
                )
            )
    if not built:
        raise RegistryError("%s defines no policies." % location)
    known = set(policy.key for policy in built)
    global COMMON_KEYS, REQUIRED_INPUTS
    COMMON_KEYS = tuple(
        key for key in (document.get("common") or []) if key in known
    )
    global GROUP_ORDER
    KEY_ORDER.clear()
    order = document.get("group_order") or []
    pairs = []
    numbered = []
    for entry in order:
        if isinstance(entry, dict) and entry.get("group"):
            numbered.append((entry.get("order", len(numbered) + 1),
                             str(entry["group"]), str(entry.get("stage", ""))))
            KEY_ORDER[str(entry["group"])] = [str(k) for k in (entry.get("keys") or [])]
        elif isinstance(entry, str):
            numbered.append((len(numbered) + 1, entry, ""))
    numbered.sort(key=lambda item: item[0])
    pairs = [(name, stage) for _n, name, stage in numbered]
    seen_groups = set(key.split(".")[0] for key in known)
    named = set(name for name, _stage in pairs)
    missing = sorted(seen_groups - named)
    if missing:
        raise RegistryError(
            "%s: group_order does not mention %s. Every group must be placed so "
            "the generated files stay in pipeline order." % (location, ", ".join(missing))
        )
    GROUP_ORDER = tuple(p for p in pairs if p[0] in seen_groups)

    global FIELD_LIFECYCLE
    lifecycle = document.get("field_lifecycle") or {}
    if not isinstance(lifecycle, dict):
        raise RegistryError("%s: 'field_lifecycle' must be a mapping." % location)
    FIELD_LIFECYCLE = lifecycle

    groups = document.get("required_inputs") or {}
    if not isinstance(groups, dict):
        raise RegistryError("%s: 'required_inputs' must be a mapping." % location)
    for name, spec in groups.items():
        if not isinstance(spec, dict) or not isinstance(spec.get("any_of"), list):
            raise RegistryError(
                "%s: required_inputs.%s needs an 'any_of' list of alternatives."
                % (location, name)
            )
        for alternative in spec["any_of"]:
            if not isinstance(alternative, list) or not alternative:
                raise RegistryError(
                    "%s: required_inputs.%s: each alternative must be a non-empty "
                    "list of config keys." % (location, name)
                )
    REQUIRED_INPUTS = groups
    return built


_REGISTRY.extend(load_registry())



# The frozen registry.
POLICIES = types.MappingProxyType(
    dict((policy.key, policy) for policy in _REGISTRY)
)  # type: Dict[str, Policy]

_GROUP_ORDER = []  # type: List[str]
for _policy in _REGISTRY:
    if _policy.group not in _GROUP_ORDER:
        _GROUP_ORDER.append(_policy.group)
GROUPS = tuple(_GROUP_ORDER)

if len(POLICIES) != len(_REGISTRY):
    raise RuntimeError("policies.py: duplicate policy key in the registry")


def policy_keys():
    # type: () -> List[str]
    """Every policy key, in registry order."""
    return [policy.key for policy in _REGISTRY]


def group_keys(group):
    # type: (str) -> List[str]
    """Every policy key in one group, in registry order."""
    return [policy.key for policy in _REGISTRY if policy.group == group]


# =============================================================================
# 3. Errors
# =============================================================================


@dataclasses.dataclass(frozen=True)
class PolicyProblem(object):
    key: str
    value: Any
    message: str

    def line(self):
        return "%-44s %s" % ("%s = %s" % (self.key, format_value(self.value)), self.message)


class PolicyError(Exception):
    """Raised when one or more policy values are wrong.

    Every problem found is reported at once; the loader never stops at the
    first.  `problems` holds the individual `PolicyProblem` records.
    """

    def __init__(self, problems):
        # type: (Sequence[PolicyProblem]) -> None
        self.problems = list(problems)
        Exception.__init__(self, self._render())

    def _render(self):
        count = len(self.problems)
        head = "%d problem%s with the policy settings:" % (count, "" if count == 1 else "s")
        lines = [head, ""]
        for problem in self.problems:
            lines.append("  " + problem.line())
        lines.append("")
        lines.append("Nothing was run. Fix the settings above and try again.")
        return "\n".join(lines)

    def __str__(self):
        return self._render()


# =============================================================================
# 4. Value formatting (shared with the logger's settings block)
# =============================================================================


def format_value(value):
    # type: (Any) -> str
    """Render a policy value the way the log and the error messages show it."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if value != value:  # NaN
            return "nan"
        if value == 0.0:
            return "0.0"
        if abs(value) < 1e-4 or abs(value) >= 1e6:
            return repr(value)
        text = ("%.10f" % value).rstrip("0")
        if text.endswith("."):
            text += "0"
        return text
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(format_value(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(
            "%s -> %s" % (key, format_value(item)) for key, item in value.items()
        ) + "}"
    return str(value)


# =============================================================================
# 5. The Policies object
# =============================================================================

_MISSING = object()


class _Group(object):
    """Attribute access for one group, e.g. `policies.eaf.out_of_range`."""

    __slots__ = ("_name", "_values")

    def __init__(self, name, values):
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_values", values)

    def __getattr__(self, item):
        if item.startswith("__"):
            raise AttributeError(item)
        values = object.__getattribute__(self, "_values")
        if item in values:
            return values[item]
        name = object.__getattribute__(self, "_name")
        raise AttributeError(_unknown_key_message("%s.%s" % (name, item)))

    def __setattr__(self, item, value):
        raise AttributeError(
            "policy values are read-only; build a new Policies object instead"
        )

    def __getitem__(self, item):
        try:
            return self.__getattr__(item)
        except AttributeError as exc:
            raise KeyError(str(exc))

    def __contains__(self, item):
        return item in object.__getattribute__(self, "_values")

    def __iter__(self):
        return iter(object.__getattribute__(self, "_values"))

    def keys(self):
        return list(object.__getattribute__(self, "_values").keys())

    def items(self):
        return list(object.__getattribute__(self, "_values").items())

    def as_dict(self):
        return copy.deepcopy(dict(object.__getattribute__(self, "_values")))

    def __repr__(self):
        name = object.__getattribute__(self, "_name")
        values = object.__getattribute__(self, "_values")
        return "<policies.%s %d settings>" % (name, len(values))


def _unknown_key_message(key):
    close = difflib.get_close_matches(key, list(POLICIES.keys()), n=3, cutoff=0.5)
    text = "there is no policy called '%s'" % key
    if close:
        text += ". Did you mean %s?" % " or ".join("'%s'" % c for c in close)
    return text


def _rebuild_policies(values, required_inputs):
    return Policies(values, required_inputs=required_inputs)


class Policies(object):
    """Resolved policy values.

    Access a value by attribute, `p.eaf.out_of_range`, or by dotted key,
    `p.get('eaf.out_of_range')`.  `p.help('eaf.out_of_range')` returns the
    plain-English explanation the logger prints.  The object is read-only and
    can be pickled, so it can be handed to worker processes as data.
    """

    def __init__(self, values, required_inputs=None):
        # type: (Dict[str, Any], Optional[Dict[str, Any]]) -> None
        flat = dict(values)
        groups = {}  # type: Dict[str, Dict[str, Any]]
        for key, value in flat.items():
            group, _, name = key.partition(".")
            groups.setdefault(group, {})[name or group] = value
        object.__setattr__(self, "_values", flat)
        object.__setattr__(
            self, "_groups", dict((g, _Group(g, v)) for g, v in groups.items())
        )
        # What a dataset must supply.  Travels with the policies so a worker or
        # a validator always sees the same rules the run was configured with.
        object.__setattr__(
            self, "required_inputs",
            copy.deepcopy(required_inputs if required_inputs is not None
                          else REQUIRED_INPUTS),
        )

    def with_required_inputs(self, override):
        """A copy whose requirement groups are replaced by `override`."""
        merged = copy.deepcopy(dict(self.required_inputs))
        for name, spec in (override or {}).items():
            if isinstance(spec, dict) and isinstance(merged.get(name), dict):
                inner = dict(merged[name]); inner.update(spec); merged[name] = inner
            else:
                merged[name] = spec
        return Policies(dict(self._values), required_inputs=merged)

    # -- construction ------------------------------------------------------
    def __reduce__(self):
        return (_rebuild_policies, (dict(self._values), dict(self.required_inputs)))

    # -- access ------------------------------------------------------------
    def __getattr__(self, item):
        if item.startswith("__"):
            raise AttributeError(item)
        groups = object.__getattribute__(self, "_groups")
        if item in groups:
            return groups[item]
        values = object.__getattribute__(self, "_values")
        if item in values:
            return values[item]
        raise AttributeError(_unknown_key_message(item))

    def __setattr__(self, item, value):
        raise AttributeError(
            "policy values are read-only; use load_policies(...) to build a new object"
        )

    def get(self, key, default=_MISSING):
        # type: (str, Any) -> Any
        """Value for a dotted key, e.g. get('eaf.out_of_range')."""
        values = object.__getattribute__(self, "_values")
        if key in values:
            return values[key]
        if default is not _MISSING:
            return default
        raise KeyError(_unknown_key_message(key))

    def __getitem__(self, key):
        return self.get(key)

    def __contains__(self, key):
        return key in object.__getattribute__(self, "_values")

    def __iter__(self):
        # type: () -> Iterator[str]
        return iter(self.keys())

    def __len__(self):
        return len(object.__getattribute__(self, "_values"))

    def keys(self):
        # type: () -> List[str]
        """Every policy key, in registry order."""
        values = object.__getattribute__(self, "_values")
        ordered = [key for key in policy_keys() if key in values]
        extra = [key for key in values if key not in POLICIES]
        return ordered + sorted(extra)

    def items(self):
        # type: () -> List[Tuple[str, Any]]
        return [(key, self.get(key)) for key in self.keys()]

    def groups(self):
        # type: () -> List[str]
        return [g for g in GROUPS if g in object.__getattribute__(self, "_groups")]

    # -- explanation -------------------------------------------------------
    def help(self, key):
        # type: (str) -> str
        """The plain-English explanation for a dotted key."""
        policy = POLICIES.get(key)
        if policy is None:
            raise KeyError(_unknown_key_message(key))
        return policy.help

    def policy(self, key):
        # type: (str) -> Policy
        """The registry entry for a dotted key."""
        policy = POLICIES.get(key)
        if policy is None:
            raise KeyError(_unknown_key_message(key))
        return policy

    def describe(self, keys=None):
        # type: (Optional[Iterable[str]]) -> List[Tuple[str, Any, str]]
        """(key, value, help) for each key, for the logger's settings block.

        An unknown key is reported in place rather than raising, so a wiring
        mistake shows up in the log instead of killing the run.
        """
        if keys is None:
            keys = self.keys()
        out = []
        for key in keys:
            if key in POLICIES:
                out.append((key, self.get(key, POLICIES[key].default), POLICIES[key].help))
            else:
                out.append((key, self.get(key, None), "(" + _unknown_key_message(key) + ")"))
        return out

    # -- output ------------------------------------------------------------
    def resolved_dict(self):
        # type: () -> Dict[str, Any]
        """Every key and its final value, flat and dotted, for the run manifest."""
        values = object.__getattribute__(self, "_values")
        return dict((key, copy.deepcopy(values[key])) for key in self.keys())

    def nested_dict(self):
        # type: () -> Dict[str, Any]
        """The same values arranged by group, for writing back out as YAML."""
        out = {}  # type: Dict[str, Any]
        for key, value in self.items():
            group, _, name = key.partition(".")
            if name:
                out.setdefault(group, {})[name] = copy.deepcopy(value)
            else:
                out[group] = copy.deepcopy(value)
        return out

    def changed_from_default(self):
        # type: () -> Dict[str, Tuple[Any, Any]]
        """Keys whose value differs from the default, as key -> (default, value)."""
        out = {}
        for key, value in self.items():
            policy = POLICIES.get(key)
            if policy is None:
                continue
            if value != policy.default:
                out[key] = (copy.deepcopy(policy.default), copy.deepcopy(value))
        return out

    def is_default(self, key):
        # type: (str) -> bool
        policy = POLICIES.get(key)
        if policy is None:
            raise KeyError(_unknown_key_message(key))
        return self.get(key) == policy.default

    def with_overrides(self, mapping):
        # type: (Optional[Dict[str, Any]]) -> "Policies"
        """A new Policies with `mapping` merged over these values."""
        return load_policies(mapping, base=self.resolved_dict())

    def __repr__(self):
        changed = len(self.changed_from_default())
        return "<Policies %d settings, %d changed from default>" % (len(self), changed)


# =============================================================================
# 6. Loading
# =============================================================================


def _flatten(mapping, prefix, flat, unknown):
    # type: (Dict[Any, Any], str, Dict[str, Any], List[Tuple[str, Any]]) -> None
    """Walk a possibly-nested mapping into flat dotted keys.

    Dotted keys, nested groups and a mixture of the two are all accepted.  A
    branch stops as soon as the joined key names a real policy, so a policy
    whose value is itself a mapping (chromosome.rename_map) is not walked into.
    """
    for raw_key, value in mapping.items():
        key = str(raw_key).strip()
        full = (prefix + "." + key) if prefix else key
        if full in POLICIES:
            flat[full] = value
        elif isinstance(value, dict) and value:
            _flatten(value, full, flat, unknown)
        elif isinstance(value, dict):
            # An empty group block is simply "no overrides here".
            continue
        else:
            unknown.append((full, value))


def load_policies(mapping=None, base=None):
    # type: (Optional[Dict[str, Any]], Optional[Dict[str, Any]]) -> Policies
    """Merge a user mapping over the defaults and validate everything.

    `mapping` may be nested ({'eaf': {'out_of_range': 'null'}}), dotted
    ({'eaf.out_of_range': 'null'}), or a mixture.  Every problem found is
    collected and reported together in one PolicyError - the loader never stops
    at the first bad value.

    `base` lets an already-resolved dict stand in for the defaults, which is how
    Policies.with_overrides works.
    """
    values = {}  # type: Dict[str, Any]
    for policy in _REGISTRY:
        values[policy.key] = policy.default_copy()
    if base:
        for key, value in base.items():
            if key in values:
                values[key] = copy.deepcopy(value)

    problems = []  # type: List[PolicyProblem]

    if mapping is None:
        mapping = {}
    if not isinstance(mapping, dict):
        raise PolicyError(
            [
                PolicyProblem(
                    "(whole policy block)",
                    mapping,
                    "must be a mapping of settings, such as {'eaf': {'out_of_range': 'null'}}",
                )
            ]
        )

    flat = {}  # type: Dict[str, Any]
    unknown = []  # type: List[Tuple[str, Any]]
    _flatten(mapping, "", flat, unknown)

    for key, value in unknown:
        problems.append(PolicyProblem(key, value, _unknown_key_message(key)))

    for key in sorted(flat, key=lambda k: policy_keys().index(k)):
        policy = POLICIES[key]
        raw = flat[key]
        try:
            coerced = policy.validator.coerce(raw)
        except Exception as exc:  # a bad shape, e.g. a list where text was meant
            problems.append(
                PolicyProblem(key, raw, "could not be read as %s (%s)"
                              % (policy.validator.describe(), exc))
            )
            continue
        message = policy.validator.check(coerced)
        if message:
            problems.append(PolicyProblem(key, raw, message))
            continue
        values[key] = coerced

    # Cross-checks that only make sense once every value is in place.
    problems.extend(_cross_check(values))

    if problems:
        raise PolicyError(problems)
    return Policies(values)


def _cross_check(values):
    # type: (Dict[str, Any]) -> List[PolicyProblem]
    """Relationships between settings that a single validator cannot see."""
    problems = []

    def pair(low_key, high_key, what):
        low, high = values.get(low_key), values.get(high_key)
        if low is None or high is None:
            return
        if low >= high:
            problems.append(
                PolicyProblem(
                    low_key, low,
                    "must be smaller than %s (%s), otherwise the %s band is empty"
                    % (high_key, format_value(high), what),
                )
            )

    pair("pvalue.clip_low", "pvalue.clip_high", "p-value")
    pair("info.clip_min", "info.clip_max", "imputation quality")
    pair("input.delimiter_min_columns", "input.delimiter_max_columns", "column count")

    mlogp_max = values.get("pvalue.mlogp_max")
    if mlogp_max is not None and mlogp_max <= values.get("pvalue.mlogp_min", 0.0):
        problems.append(
            PolicyProblem(
                "pvalue.mlogp_max", mlogp_max,
                "must be larger than pvalue.mlogp_min (%s), or null to turn the cap off"
                % format_value(values.get("pvalue.mlogp_min")),
            )
        )

    tolerance = values.get("eaf.clip_tolerance")
    if tolerance is not None and values.get("eaf.out_of_range") == "fail" and tolerance > 1.0:
        # Not an error, just impossible to reach; leave it alone.
        pass

    warn = values.get("validation.warn_reject_fraction")
    hard = values.get("validation.max_reject_fraction")
    if warn is not None and hard is not None and warn > hard:
        problems.append(
            PolicyProblem(
                "validation.warn_reject_fraction", warn,
                "must not be above validation.max_reject_fraction (%s), or the run fails "
                "before the warning can ever be printed" % format_value(hard),
            )
        )

    duplicate_order = values.get("duplicates.selection_order") or []
    if "input_order" not in duplicate_order:
        problems.append(
            PolicyProblem(
                "duplicates.selection_order", duplicate_order,
                "must contain input_order so scientifically tied duplicate rows "
                "always produce the same retained row",
            )
        )
    elif duplicate_order[-1] != "input_order":
        problems.append(
            PolicyProblem(
                "duplicates.selection_order", duplicate_order,
                "must place input_order last because it is the deterministic final "
                "tie-breaker, not a quality measure",
            )
        )

    return problems


def default_policies():
    # type: () -> Policies
    """Return every policy at its canonical YAML default."""
    return load_policies({})


def describe(keys, policies=None):
    # type: (Iterable[str], Optional[Policies]) -> List[Tuple[str, Any, str]]
    """(key, value, help) for each key, for the logger's settings block.

    With no Policies object the defaults are described.
    """
    if policies is None:
        policies = default_policies()
    return policies.describe(keys)


# =============================================================================
# Discovery helpers
# =============================================================================
#
# Inside a container this file is baked into the image, so a user cannot read it
# to find out what is configurable.  These render the registry as something they
# can see and edit, and the module is runnable directly:
#
#     python -m harmonisation.policies              # readable table
#     python -m harmonisation.policies --yaml       # commented YAML template
#     python -m harmonisation.policies --yaml eaf   # one group only
#
# The YAML it emits is a valid `policies:` block: save it, edit it, and pass it
# with --defaults.


def _ordered_groups(keys):
    """Group names in pipeline order, falling back to alphabetical."""
    present = []
    for name, _stage in GROUP_ORDER:
        if any(k.split(".")[0] == name for k in keys):
            present.append(name)
    for k in sorted(keys):
        head = k.split(".")[0]
        if head not in present:
            present.append(head)
    return present


def _ordered_keys(group, keys):
    """Settings in the order the owning step reads them, else alphabetical."""
    wanted = KEY_ORDER.get(group) or []
    present = set(keys)
    out = [group + "." + leaf for leaf in wanted if group + "." + leaf in present]
    out.extend(sorted(k for k in keys if k not in out))
    return out


def _stage_of(group):
    for name, stage in GROUP_ORDER:
        if name == group:
            return stage
    return ""


def as_table(group=None, policies=None):
    # type: (Optional[str], Any) -> str
    """Every setting with its default and its plain-English meaning."""
    pol = policies if policies is not None else default_policies()
    keys = [k for k in pol.keys() if group is None or k.split(".")[0] == group]
    if not keys:
        known = ", ".join(sorted({k.split(".")[0] for k in pol.keys()}))
        return "No policy group %r. Known groups: %s" % (group, known)
    width = min(46, max(len(k) for k in keys))
    out = []
    for position, head in enumerate(_ordered_groups(keys), 1):
        stage = _stage_of(head)
        out.append("")
        out.append("%02d  [%s]%s" % (position, head, ("   %s" % stage) if stage else ""))
        group_keys = [k for k in keys if k.split(".")[0] == head]
        for key in _ordered_keys(head, group_keys):
            out.append("  %-*s = %s" % (width, key, format_value(pol.get(key))))
            for line in _wrap(str(pol.help(key) or ""), 74):
                out.append(("      %s" % line).rstrip())
    return "\n".join(out).lstrip("\n")


def _short_gloss(pol, key, width=58):
    """A few words describing a setting, for a trailing comment."""
    validator = pol.policy(key).validator if hasattr(pol, "policy") else None
    members = getattr(validator, "members", None)
    if members:
        choices = " | ".join(str(m) for m in members)
        if len(choices) <= width:
            return choices
    text = str(pol.help(key) or "").strip()
    sentence = text.split(". ")[0].rstrip(".")
    if len(sentence) > width:
        sentence = sentence[: width - 1].rsplit(" ", 1)[0] + "\u2026"
    return sentence


def as_yaml_template(group=None, changed_only=False, policies=None,
                     style="minimal", common_only=False,
                     include_header=True, include_required_inputs=True):
    # type: (Optional[str], bool, Any, str, bool) -> str
    """A ``policies:`` block ready to save and pass with --parameter-config.

    ``style='minimal'`` one line per setting with a short trailing comment.
    ``style='full'``    the whole help text as comment lines above each setting.
    ``style='values'``  key-value pairs only, with no comments.
    ``common_only``     just the settings the registry marks as commonly changed.
    """
    if style == "simple":
        style = "minimal"
    if style not in ("full", "minimal", "values"):
        raise ValueError("style must be one of: full, minimal, values")
    pol = policies if policies is not None else default_policies()
    keys = [k for k in pol.keys() if group is None or k.split(".")[0] == group]
    if common_only:
        keys = [k for k in keys if k in COMMON_KEYS]
    if changed_only:
        keys = [k for k in keys if not pol.is_default(k)]
    if not keys:
        return "# no settings matched\n"

    if not include_header or style == "values":
        head = []
    elif common_only:
        head = [
            "# postgwas - the settings most runs need.",
            "# Every value shown is the shipped default; delete any line you are not",
            "# changing.  For the complete set:  python -m harmonisation.policies --yaml",
        ]
    else:
        head = [
            "# postgwas policy settings.",
            "# Every value shown is the shipped default; delete any line you are not",
            "# changing.  Run with --full for the complete explanation of each one.",
        ]
    if head:
        head.append("# Save this output as YAML and pass it with --run-config.")
    out = head + ["policies:"]

    width = max(len(k.split(".", 1)[1]) for k in keys) + 1
    for position, head_group in enumerate(_ordered_groups(keys), 1):
        stage = _stage_of(head_group)
        if style != "values":
            out.append("")
            out.append("  # %02d  %s" % (position, stage) if stage else "  # %02d" % position)
        out.append("  %s:" % head_group)
        for key in _ordered_keys(head_group, [k for k in keys if k.split(".")[0] == head_group]):
            leaf = key.split(".", 1)[1]
            if style == "full":
                for line in _wrap(str(pol.help(key) or ""), 68):
                    out.append("    # %s" % line)
                out.append("    %s: %s" % (leaf, _yaml_scalar(pol.get(key))))
            elif style == "minimal":
                entry = "    %-*s %s" % (width, leaf + ":", _yaml_scalar(pol.get(key)))
                out.append("%-46s # %s" % (entry, _short_gloss(pol, key)))
            else:
                out.append("    %s: %s" % (leaf, _yaml_scalar(pol.get(key))))
    required = _required_inputs_yaml(pol, style) if include_required_inputs else ""
    return required + "\n".join(out) + "\n"


def _required_inputs_yaml(pol, style="simple"):
    """The --config requirement table, as an editable YAML block."""
    groups = getattr(pol, "required_inputs", None) or REQUIRED_INPUTS
    if not groups:
        return ""
    out = [
        "",
        "# ===========================================================================",
        "# What every row of the --config CSV must supply, in the CONFIG COLUMN NAMES",
        "# used in that file.  A group is satisfied when every key of ANY ONE",
        "# alternative is set (not blank, not NA).  A config that fails this is",
        "# rejected before the summary statistics file is opened.",
        "#",
        "# Editing this changes what is accepted.  Set `required: false` to stop",
        "# demanding a group, or add an alternative to accept another source.",
        "# ===========================================================================",
        "required_inputs:",
    ]
    for name in groups:
        spec = groups[name] or {}
        alternatives = spec.get("any_of") or []
        readable = "  or  ".join(" + ".join(a) for a in alternatives)
        out.append("")
        out.append("  # %s%s" % (spec.get("label", name),
                                 ("   ->   " + readable) if readable else ""))
        if style == "full" and spec.get("why"):
            for line in _wrap(" ".join(str(spec["why"]).split()), 68):
                out.append("  # %s" % line)
        out.append("  %s:" % name)
        out.append("    required: %s" % _yaml_scalar(bool(spec.get("required", True))))
        out.append("    any_of:")
        for alternative in alternatives:
            out.append("      - [%s]" % ", ".join(alternative))
    return "\n".join(out) + "\n"


def _wrap(text, width):
    # type: (str, int) -> List[str]
    """Wrap to `width`, keeping any explicit line breaks the author put in."""
    raw = str(text)
    if "\n" in raw:
        lines = []
        for chunk in raw.split("\n"):
            if not chunk.strip():
                lines.append("")
                continue
            # keep the author's indentation so a list of values stays a list
            indent = chunk[: len(chunk) - len(chunk.lstrip())]
            wrapped = _wrap_one(chunk, max(20, width - len(indent)))
            lines.append(indent + wrapped[0])
            for extra in wrapped[1:]:
                lines.append(indent + "  " + extra)
        return lines
    return _wrap_one(raw, width)


def _wrap_one(text, width):
    # type: (str, int) -> List[str]
    words = str(text).split()
    if not words:
        return []
    lines, line = [], words[0]
    for word in words[1:]:
        if len(line) + 1 + len(word) <= width:
            line += " " + word
        else:
            lines.append(line)
            line = word
    lines.append(line)
    return lines


#: Bare words YAML 1.1 reads as something other than a string.  A policy value
#: like ``off`` MUST be quoted or the file will not load back as itself.
_YAML_RESERVED = frozenset("""
y yes n no true false on off null none ~
""".split())


def _yaml_scalar(value):
    # type: (Any) -> str
    """Render a policy value so that reading the file back gives the same value."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return "[%s]" % ", ".join(_yaml_scalar(v) for v in value)
    if isinstance(value, dict):
        return "{%s}" % ", ".join(
            "%s: %s" % (k, _yaml_scalar(v)) for k, v in value.items()
        )
    if isinstance(value, str):
        if value == "":
            return '""'
        if value.lower() in _YAML_RESERVED:
            return '"%s"' % value
        # Anything that would come back as a number must be quoted too.
        try:
            float(value)
            return '"%s"' % value
        except ValueError:
            pass
        if value.strip() != value or any(
            c in value for c in ":#{}[],&*?|<>=!%@`\"'"
        ):
            return '"%s"' % value
        return value
    return repr(value)
