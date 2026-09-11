"""Atomic writers for small structured PostGWAS reports."""

from __future__ import annotations

import csv
from html import escape
from html.parser import HTMLParser
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import yaml


_INTERACTIVE_TABLE_MARKER = "postgwas-table-interactivity"
_INTERACTIVE_TABLE_STYLE = """<style id="postgwas-table-interactivity-style">
.postgwas-column-filters th{position:static!important;padding:.35rem .45rem!important;background:#f8fafc!important}
.postgwas-column-filter{box-sizing:border-box;width:100%;min-width:7rem;padding:.38rem .48rem;border:1px solid #b8c4d6;border-radius:.35rem;background:#fff;color:#172033;font:inherit;font-size:.78rem;text-transform:none;letter-spacing:normal}
.postgwas-sort-button{display:inline-flex!important;align-items:center;gap:.35rem;border:0;background:transparent;color:inherit;font:inherit;font-weight:inherit;text-align:left;cursor:pointer}
.postgwas-sort-button::after{content:"↕";opacity:.45;font-size:.82em}
.postgwas-sort-button[data-direction="ascending"]::after{content:"▲";opacity:1}
.postgwas-sort-button[data-direction="descending"]::after{content:"▼";opacity:1}
@media print{.postgwas-column-filters{display:none!important}}
</style>"""
_INTERACTIVE_TABLE_SCRIPT = r"""<script id="postgwas-table-interactivity">
(function() {
  'use strict';
  const dynamicSelector = '.paginated-results[data-source],.data-browser[data-source],.postgwas-data-browser[data-source]';

  function displayed(value, nullValue) {
    return value === null || value === undefined || value === ''
      ? (nullValue || 'NA')
      : String(value);
  }

  function sortable(value) {
    const text = displayed(value).trim();
    const numeric = Number(text.replaceAll(',', '').replace(/%$/, ''));
    return text !== '' && Number.isFinite(numeric)
      ? {missing: false, numeric: true, value: numeric}
      : {missing: /^(?:NA|N\/A|NAN|NONE|NULL|NOT AVAILABLE|-|\.)$/i.test(text), numeric: false, value: text.toLocaleLowerCase()};
  }

  function compare(left, right, ascending) {
    const a = sortable(left), b = sortable(right);
    if (a.missing !== b.missing) return a.missing ? 1 : -1;
    let order;
    if (a.numeric && b.numeric) order = a.value - b.value;
    else order = String(a.value).localeCompare(String(b.value), undefined, {numeric: true});
    return ascending ? order : -order;
  }

  function stableSort(rows, column, ascending) {
    return rows.map(function(row, index) { return {row: row, index: index}; })
      .sort(function(left, right) {
        return compare(left.row[column], right.row[column], ascending) || left.index - right.index;
      }).map(function(item) { return item.row; });
  }

  function ensureHeader(table, columns) {
    let head = table.tHead;
    if (!head) head = table.createTHead();
    let row = head.rows[0];
    if (!row) row = head.insertRow();
    if (!row.cells.length) {
      columns.forEach(function(column) {
        const cell = document.createElement('th');
        cell.scope = 'col'; cell.textContent = column; row.appendChild(cell);
      });
    }
    return row;
  }

  function makeHeaders(headerRow, columns, sort) {
    Array.from(headerRow.cells).forEach(function(cell, index) {
      const label = columns[index] || cell.textContent.trim() || `Column ${index + 1}`;
      cell.replaceChildren();
      const button = document.createElement('button');
      button.type = 'button'; button.className = 'postgwas-sort-button';
      button.textContent = label; button.setAttribute('aria-label', `Sort by ${label}`);
      button.addEventListener('click', function() { sort(index, button); });
      cell.appendChild(button);
    });
  }

  function makeFilters(table, columns, apply) {
    const existing = table.tHead.querySelector('.postgwas-column-filters');
    if (existing) existing.remove();
    const row = table.tHead.insertRow(-1); row.className = 'postgwas-column-filters';
    columns.forEach(function(column, index) {
      const cell = document.createElement('th'); cell.scope = 'col';
      const input = document.createElement('input');
      input.type = 'search'; input.className = 'postgwas-column-filter';
      input.placeholder = 'Filter…'; input.autocomplete = 'off';
      input.setAttribute('aria-label', `Filter ${column || `column ${index + 1}`}`);
      input.addEventListener('input', apply); cell.appendChild(input); row.appendChild(cell);
    });
    return Array.from(row.querySelectorAll('input'));
  }

  function rowMatches(row, filters, globalTerm, nullValue) {
    if (globalTerm && !row.some(function(value) {
      return displayed(value, nullValue).toLocaleLowerCase().includes(globalTerm);
    })) return false;
    return filters.every(function(input, index) {
      const term = input.value.trim().toLocaleLowerCase();
      return !term || displayed(row[index], nullValue).toLocaleLowerCase().includes(term);
    });
  }

  function cloneControl(scope, selectors) {
    const control = scope.querySelector(selectors);
    if (!control) return null;
    const clone = control.cloneNode(true); control.replaceWith(clone); return clone;
  }

  function initializeDynamic(container) {
    const table = container.querySelector('table');
    const payload = document.getElementById(container.dataset.source);
    if (!table || !payload || table.dataset.postgwasInteractive) return;
    const parsed = JSON.parse(payload.textContent);
    const sourceRows = Array.isArray(parsed) ? parsed : parsed.rows;
    if (!Array.isArray(sourceRows)) return;
    let columns = Array.isArray(parsed.columns) ? parsed.columns.map(String) : [];
    const initialHeader = table.tHead && table.tHead.rows[0];
    if (!columns.length && initialHeader) {
      columns = Array.from(initialHeader.cells).map(function(cell, index) {
        return cell.textContent.trim() || `Column ${index + 1}`;
      });
    }
    if (!columns.length && sourceRows.length) {
      columns = sourceRows[0].map(function(_value, index) { return `Column ${index + 1}`; });
    }
    if (!columns.length) return;
    let header = ensureHeader(table, columns);
    const cleanHeader = header.cloneNode(true);
    header.replaceWith(cleanHeader); header = cleanHeader;
    const body = table.tBodies[0] || table.createTBody();
    const globalSearch = cloneControl(container, '.result-search,.table-search,#search');
    const previous = cloneControl(container, '.previous-page,.previous,#prev');
    const next = cloneControl(container, '.next-page,.next,#next');
    const status = container.querySelector('.page-status,.page-state,#page');
    const pageSize = Math.max(1, Number(container.dataset.pageSize) || sourceRows.length || 1);
    const nullValue = container.dataset.nullValue || 'NA';
    let filtered = sourceRows.slice(), page = 0, sortColumn = null, ascending = true;
    let filters = [];

    function apply() {
      const globalTerm = globalSearch ? globalSearch.value.trim().toLocaleLowerCase() : '';
      filtered = sourceRows.filter(function(row) { return rowMatches(row, filters, globalTerm, nullValue); });
      if (sortColumn !== null) filtered = stableSort(filtered, sortColumn, ascending);
      page = 0; render();
    }
    function sort(column, button) {
      ascending = sortColumn === column ? !ascending : true; sortColumn = column;
      table.querySelectorAll('.postgwas-sort-button').forEach(function(item) { item.removeAttribute('data-direction'); });
      button.dataset.direction = ascending ? 'ascending' : 'descending';
      filtered = stableSort(filtered, column, ascending); page = 0; render();
    }
    function render() {
      const pages = Math.max(1, Math.ceil(filtered.length / pageSize));
      page = Math.max(0, Math.min(page, pages - 1));
      const start = page * pageSize, end = Math.min(start + pageSize, filtered.length);
      body.replaceChildren();
      filtered.slice(start, end).forEach(function(row) {
        const tr = document.createElement('tr');
        columns.forEach(function(_column, index) {
          const td = document.createElement('td'); td.textContent = displayed(row[index], nullValue); tr.appendChild(td);
        }); body.appendChild(tr);
      });
      if (status) status.textContent = filtered.length
        ? `${start + 1}–${end} of ${filtered.length.toLocaleString()} · page ${page + 1} of ${pages}`
        : '0 matching rows';
      if (previous) previous.disabled = page === 0;
      if (next) next.disabled = page + 1 >= pages;
    }
    makeHeaders(header, columns, sort);
    filters = makeFilters(table, columns, apply);
    if (globalSearch) globalSearch.addEventListener('input', apply);
    if (previous) previous.addEventListener('click', function() { page -= 1; render(); });
    if (next) next.addEventListener('click', function() { page += 1; render(); });
    table.dataset.postgwasInteractive = 'dynamic'; render();
  }

  function initializeStatic(table) {
    if (table.dataset.postgwasInteractive) return;
    const bodies = Array.from(table.tBodies);
    const rows = bodies.flatMap(function(body) { return Array.from(body.rows); });
    let header = table.tHead && table.tHead.rows[0];
    const columnCount = header && header.cells.length
      ? header.cells.length
      : rows.length
        ? Math.max.apply(null, rows.map(function(row) { return row.cells.length; }))
        : 0;
    if (!columnCount) return;
    let columns;
    if (header && header.cells.length) {
      columns = Array.from(header.cells).map(function(cell, index) {
        return cell.textContent.trim() || `Column ${index + 1}`;
      });
    } else {
      columns = table.classList.contains('key-value') && columnCount === 2
        ? ['Field', 'Value']
        : Array.from({length: columnCount}, function(_value, index) { return `Column ${index + 1}`; });
      header = ensureHeader(table, columns);
    }
    let sortColumn = null, ascending = true, filters = [];
    function apply() {
      rows.forEach(function(row) {
        const values = Array.from({length: columnCount}, function(_value, index) {
          return row.cells[index] ? row.cells[index].textContent : '';
        });
        row.hidden = !rowMatches(values, filters, '');
      });
    }
    function sort(column, button) {
      ascending = sortColumn === column ? !ascending : true; sortColumn = column;
      table.querySelectorAll('.postgwas-sort-button').forEach(function(item) { item.removeAttribute('data-direction'); });
      button.dataset.direction = ascending ? 'ascending' : 'descending';
      const ordered = stableSort(rows.map(function(row) {
        return [row.cells[column] ? row.cells[column].textContent : '', row];
      }), 0, ascending).map(function(item) { return item[1]; });
      ordered.forEach(function(row) { row.parentNode.appendChild(row); }); apply();
    }
    makeHeaders(header, columns, sort); filters = makeFilters(table, columns, apply);
    table.dataset.postgwasInteractive = 'static';
  }

  document.querySelectorAll(dynamicSelector).forEach(initializeDynamic);
  document.querySelectorAll('table').forEach(initializeStatic);
})();
</script>"""


def _with_table_interactivity(document: str) -> str:
    """Add one offline sort/filter controller to every published HTML table."""
    if "<table" not in document.lower() or _INTERACTIVE_TABLE_MARKER in document:
        return document
    lower = document.lower()
    head = lower.rfind("</head>")
    if head >= 0:
        document = document[:head] + _INTERACTIVE_TABLE_STYLE + document[head:]
    else:
        document = _INTERACTIVE_TABLE_STYLE + document
    lower = document.lower()
    body = lower.rfind("</body>")
    if body >= 0:
        return document[:body] + _INTERACTIVE_TABLE_SCRIPT + document[body:]
    return document + _INTERACTIVE_TABLE_SCRIPT


class _ReportContentParser(HTMLParser):
    """Copy static report content; omit scripts, controls and source-page IDs."""

    # HTML structure is a protocol invariant, not a scientific column mapping.
    tags = {"section", "article", "div", "span", "p", "h1", "h2", "h3", "h4",
            "table", "thead", "tbody", "tfoot", "tr", "th", "td", "strong",
            "em", "b", "i", "code", "pre", "ul", "ol", "li", "dl", "dt", "dd",
            "details", "summary", "br", "hr"}
    void_tags = {"input", "img", "br", "hr", "meta", "link"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.in_main = False
        self.found_main = False
        self.skipped = 0
        self.parts = []

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if tag == "main":
            self.in_main = self.found_main = True
            return
        if not self.in_main:
            return
        if self.skipped:
            if tag not in self.void_tags:
                self.skipped += 1
            return
        if tag in {"script", "style", "nav", "form", "button", "label", "noscript"} or (
            set(attrs.get("class", "").split()) & {"result-tools", "table-tools"}
        ):
            self.skipped = 1
            return
        if tag in self.tags:
            safe_attrs = "".join(
                ' %s="%s"' % (name, escape(value or "", quote=True))
                for name, value in attributes if name in {"class", "colspan", "rowspan", "scope", "open"}
            )
            self.parts.append("<%s%s>" % (tag, safe_attrs))

    def handle_endtag(self, tag):
        if tag in self.void_tags:
            return
        if tag == "main":
            self.in_main = False
        elif self.in_main and self.skipped:
            self.skipped -= 1
        elif self.in_main and tag in self.tags and tag not in self.void_tags:
            self.parts.append("</%s>" % tag)

    def handle_data(self, data):
        if self.in_main and not self.skipped:
            self.parts.append(escape(data))


def read_static_report_content(path: str | Path) -> str:
    """Copy published table text and explanations without executing source HTML.

    Links become plain text; the containing report links the original document.
    """
    parser = _ReportContentParser()
    parser.feed(Path(path).read_text(encoding="utf-8"))
    parser.close()
    return "".join(parser.parts) if parser.found_main else ""


def collect_html_reports(results: Mapping[str, Any] | None) -> list[dict[str, str]]:
    """Collect explicit module-result HTML paths without searching directories.

    Nested target results may refer to the same report (for example, multiple
    formatter outputs). Keep one entry per module and resolved report path.
    File availability is a presentation concern, not proof of stage success.
    """
    reports = []
    seen = set()

    def visit(module, value):
        if isinstance(value, Mapping):
            raw = value.get("html_report")
            if raw is not None and raw != "":
                if not isinstance(raw, (str, Path)):
                    raise ValueError("Recorded html_report must be a file path")
                path = str(Path(raw).expanduser().resolve())
                identity = (module, path)
                if identity not in seen:
                    reports.append({"module": module, "path": path})
                    seen.add(identity)
            for child in value.values():
                if isinstance(child, (Mapping, list, tuple)):
                    visit(module, child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit(module, child)

    for module, value in (results or {}).items():
        visit(str(module), value)
    return reports


def _atomic_text_path(destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    return tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=destination.parent,
        prefix=".%s." % destination.name,
        delete=False,
    )


def write_yaml_report(value: Mapping[str, Any], path: str | Path) -> Path:
    """Write one ordered YAML document without exposing a partial result."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            yaml.safe_dump(
                dict(value), handle, sort_keys=False, allow_unicode=True,
            )
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


def write_delimited_report(
    records: Sequence[Mapping[str, Any]],
    path: str | Path,
    *,
    fieldnames: Sequence[str],
    delimiter: str,
    null_value: str,
) -> Path:
    """Write homogeneous records atomically using one configured delimiter."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fieldnames),
                delimiter=delimiter,
                extrasaction="ignore",
            )
            writer.writeheader()
            for record in records:
                writer.writerow({
                    name: null_value if record.get(name) is None else record.get(name)
                    for name in fieldnames
                })
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


def write_html_report(document: str, path: str | Path) -> Path:
    """Write one interactive, complete UTF-8 HTML document atomically."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            handle.write(_with_table_interactivity(document))
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


def write_text_lines_report(lines: Sequence[Any], path: str | Path) -> Path:
    """Write one value per line atomically, without adding a table header."""
    destination = Path(path)
    temporary_path = None
    try:
        with _atomic_text_path(destination) as handle:
            temporary_path = Path(handle.name)
            for value in lines:
                handle.write("%s\n" % value)
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


__all__ = [
    "read_static_report_content",
    "collect_html_reports",
    "write_delimited_report",
    "write_html_report",
    "write_text_lines_report",
    "write_yaml_report",
]
