# Frozen postprocessor results

These ten files were copied byte-for-byte from existing Octave exports before
the recovery-operator refactoring (2026-11-04). They cover JSON versions 1–7,
MPCs, releases, linear loads, catalog sections and spatial static/modal/transient
and partial data. `manifest.json` records the original SHA-256 and format/type.
`.gitattributes` disables line-ending conversion for these JSON files so their
checksums stay identical on Windows and Linux checkouts.

Do not regenerate these files with the current exporter during tests. Their
purpose is to verify the explicit rejection of unsupported saved files, before
building a scene. Add further legacy samples with explicit
provenance when needed; changing a fixture or checksum requires review.

`scripts/test_recovery_operators.cjs` checks the manifest and rejects every file.
`scripts/test_recovery_browser.cjs` also loads all frozen files through
the real file picker and checks the instruction to recalculate the original TXT.
Fresh v8 results and direct solver references live in the isolated output
directory of each Docker pipeline run, not here.
