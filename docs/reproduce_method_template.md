# How to reproduce method `{{METHOD}}`

**Snapshot path:** `{{SNAPSHOT_DIR}}`  
**Git ref:** `{{GIT_REF}}`  
**Frozen at:** `{{CREATED_AT}}`  
**Julia:** `{{JULIA_VERSION}}` · **Package:** `{{PACKAGE_VERSION}}`  
**Manifest sha256:** `{{MANIFEST_SHA256}}`

## When this snapshot was created

Only freeze after a method is **short-listed on coarse invent**, has passed **fine-mesh confirm**, and preferably a **robustness** look — not after every invent run.

**Paper-facing freezes should use:**

```bash
frforge snapshot freeze … --require-confirm
```

That hard-fails unless a `confirmed` log entry (or confirm compare JSON) exists. Default freeze only **warns** if confirm is missing.

## Evaluation path (two-tier)

1. Coarse invent (composite-score history): `frforge invent --method {{METHOD}} --baseline {{BASELINE}}`
2. Fine-mesh confirmation: `frforge confirm --method {{METHOD}} --baseline {{BASELINE}}`
3. Optional robustness matrix: `frforge robustness --method {{METHOD}}`
4. Freeze with `--require-confirm` for publication claims

Confirm re-runs 2D Riemann cfg 6, reduced Double-Mach, and (by default) isentropic vortex order on denser meshes. It does **not** rewrite invent composite scores.

**Serial residual only for official confirm.** A threaded confirm (`--threads N` with \(N>1\)) is informational only and **does not** satisfy fine-mesh confirmation for `publication_grade` or `snapshot freeze --require-confirm`.

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
   # then re-confirm if claiming publication grade:
   frforge confirm --method {{METHOD}} --baseline {{BASELINE}}
   ```

5. **Tables** from frozen JSON:

   ```bash
   frforge snapshot tables {{SNAPSHOT_DIR}} --out tables.md --csv tables.csv
   ```

## Primary frozen metrics

See `SNAPSHOT.json` → `primary_metrics` and `reports/compare.json`.
