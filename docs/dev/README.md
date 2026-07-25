# Development notes

Working documents from the v0.7.0 correctness audit. These are a historical record, not
user-facing documentation, and they are **not** part of the built documentation — they are kept
because the changelog references their bug and plan IDs.

| File | What it is |
| --- | --- |
| [`bugs.md`](bugs.md) | The findings report. Bug IDs `B1`–`B15` (defects) and `C1`–`C4` (consistency issues) are cited throughout `CHANGELOG.md`. |
| [`plan.md`](plan.md) | The remediation plan derived from `bugs.md`, with plan IDs `P1`–`P17`. |
| [`memory.md`](memory.md) | The execution log: which plan item was addressed how, and what was deliberately deferred. |

For the current state of the package, read [`CHANGELOG.md`](../../CHANGELOG.md) instead — in
particular its *Known follow-ups* section, which carries forward the items these files left open.
