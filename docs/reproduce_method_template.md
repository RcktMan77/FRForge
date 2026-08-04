# How to reproduce method `{{METHOD}}`

**Snapshot path:** `{{SNAPSHOT_DIR}}`  
**Git ref:** `{{GIT_REF}}`  
**Frozen at:** `{{CREATED_AT}}`  
**Julia:** `{{JULIA_VERSION}}` · **Package:** `{{PACKAGE_VERSION}}`  
**Manifest sha256:** `{{MANIFEST_SHA256}}`

## When this snapshot was created

Only freeze after a method is **short-listed** or has survived a **robustness** look — not after every invent run.

## Steps

1. Check out `{{GIT_REF}}` (or a revision that includes the method sources).
2. Instantiate the project: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
3. **Cheap verify** (default; no invent re-run):

   ```bash
   frforge snapshot verify {{SNAPSHOT_DIR}}
   ```

4. **Full re-run** (explicit; regenerates primary metrics within FP noise):

   ```bash
   frforge snapshot verify {{SNAPSHOT_DIR}} --rerun
   # or:
   frforge invent --method {{METHOD}} --baseline {{BASELINE}}
   ```

5. **Tables** from frozen JSON:

   ```bash
   frforge snapshot tables {{SNAPSHOT_DIR}} --out tables.md --csv tables.csv
   ```

## Primary frozen metrics

See `SNAPSHOT.json` → `primary_metrics` and `reports/compare.json`.
