/* Run: node tests/test_fine_mapping_report_interactions.cjs REPORT.html
 * Execute the emitted report script against a small DOM fixture. This tests
 * event/state behavior, not browser layout or native disclosure rendering.
 */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');

function element(textContent = '', selectors = {}) {
  const listeners = {};
  return {
    textContent, hidden: false, open: false, value: '',
    querySelectorAll(selector) { return selectors[selector] || []; },
    querySelector(selector) { return this.querySelectorAll(selector)[0] || null; },
    addEventListener(event, callback) { listeners[event] = callback; },
    fire(event) { assert.ok(listeners[event], event); listeners[event](); }
  };
}

const firstVariants = [element('shared_variant PIP 0.7'), element('rare_<variant> PIP 0.26')];
const secondVariants = [element('shared_variant PIP 0.96')];
const firstStats = element('set coverage 0.95');
const secondStats = element('set coverage 0.95');
function card(title, variants, stats) {
  return element('', {
    '.cs-summary': [element(title)], '.provenance': [element('recorded locus')],
    '.set-details': [stats], '.variant-detail': variants
  });
}
const cards = [
  card('chr1 L1 initial warning_sample_size', firstVariants, firstStats),
  card('chr2 L1 additional', secondVariants, secondStats)
];
const search = element();
const status = element();
const expand = element();
const collapse = element();
const controls = element();
controls.hidden = true;
const details = [...cards, ...firstVariants, ...secondVariants, firstStats, secondStats];
const browser = element('', {
  '.set-search': [search], '.set-search-state': [status], '.cs-card': cards,
  '.expand-sets': [expand], '.collapse-sets': [collapse],
  '.table-tools': [controls], details
});
const document = element('', {'.credible-set-browser': [browser]});
const html = fs.readFileSync(process.argv[2], 'utf8');
const script = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(match => match[1]).join('\n');
assert.ok(script.includes("'.credible-set-browser'"), 'combined browser script present');
vm.runInNewContext(script, {document});
assert.equal(controls.hidden, false);
assert.equal(status.textContent, '2 of 2 sets · 3 of 3 variant entries');

// A hit inside a collapsed set opens it; nonmatching members are only hidden.
cards[1].open = true;
secondVariants[0].open = true;
function find(term) { search.value = term; search.fire('input'); }
find('  RARE_<VARIANT>  ');
assert.equal(cards[0].open, true);
assert.equal(cards[1].hidden, true);
assert.equal(firstVariants[0].hidden, true);
assert.equal(firstVariants[1].hidden, false);
assert.equal(status.textContent, '1 of 2 sets · 1 of 3 variant entries');

// Shared variants remain represented in every containing set.
find('shared_variant');
assert.ok(cards.every(card => card.open && !card.hidden));
assert.equal(status.textContent, '2 of 2 sets · 2 of 3 variant entries');
find('chr2');
assert.equal(cards[0].hidden, true);
assert.equal(secondVariants[0].hidden, false);
find('warning_sample_size');
assert.ok(firstVariants.every(variant => !variant.hidden));
find('does not exist');
assert.ok(cards.every(card => card.hidden));
assert.equal(status.textContent, 'No matching sets or variants.');
expand.fire('click');
assert.ok(cards.every(card => !card.open));

// Clearing restores the pre-search state, including expanded SNP details.
find('');
assert.equal(cards[0].open, false);
assert.equal(cards[1].open, true);
assert.equal(secondVariants[0].open, true);
assert.ok(details.every(node => !node.hidden));
expand.fire('click');
assert.ok(cards.every(card => card.open));
collapse.fire('click');
assert.ok(details.every(node => !node.open));

// Empty reports contain no browser and must also execute without an error.
vm.runInNewContext(script, {document: element()});
console.log('PASS: collapsed search, case/whitespace, shared variants, locus/warning matches, no matches, clear/restore, expand/collapse, empty report.');
