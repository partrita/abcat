# ADR-001: Improve CLI feedback

## Decision

The CLI will validate input before analysis. Batch commands will continue after
an item fails and will save the failure in the result. The `analyze` command
will provide a short table by default and a detailed table with `--verbose`.

## Reasons

- Users need a clear cause for input errors.
- One bad FASTA record must not remove valid results.
- A short default table is easier to read.

## Scope

- Validate the 20 standard amino acids.
- Validate the JSON/CSV format and numbering scheme.
- Show batch progress and a success/failure summary.
- Support `--sequence-id` for one result in console output.

## Trade-offs

Batch output can contain failed records. Consumers must check `success` or
`Status` before using analysis fields.
