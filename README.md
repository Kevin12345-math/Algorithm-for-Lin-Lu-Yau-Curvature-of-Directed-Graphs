# Ricci Curvature Calculator for Cayley Graphs

Computes **Lin-Lu-Yau Ricci curvature** κ(e, g) on Cayley graphs of finitely generated groups,
based on: *Algorithm for Calculating Ricci Curvature* by Kevin Fung (May 2026).

---

## Architecture

The script (`ricci_curvature_all.py`) is organized into six sections:

### 1. Group Classes
Each group class encodes a specific group's element representation and multiplication rules.

| Class | Group | Generators |
|---|---|---|
| `Q4m` | Generalized quaternion group $Q_{4m}$ | $a, a^{-1}, b, b^{-1}$ |
| `DihedralGroup` | Dihedral group $D_n$ | $a, a^{-1}, b$ |
| `CyclicGroupZn` | Cyclic group $\mathbb{Z}_n$ with step $k$ | $\pm 1, \pm k$ |
| `SymmetricGroup` | Symmetric group $S_n$ | adjacent transpositions $\sigma_i$ |
| `AlternatingGroup` | Alternating group $A_n$ | Carmichael 3-cycles $V_i, V_i^{-1}$ |
| `GeneralLinearGroupZ` | $\mathrm{GL}_n(\mathbb{Z})$ | elementary matrices (finite only for $n=1$) |

All classes implement a common interface: `identity()`, `multiply(elem, gen)`, `neighborhood(elem, S)`, `elements()`, `elem_name(elem)`.

### 2. Core Algorithms (group-agnostic)

- **`all_distances(group, S)`** — BFS all-pairs shortest-path distances on the Cayley graph.
- **`min_matching_cost(sources, targets, dist)`** — Minimum-weight perfect matching via brute-force permutations.
- **`ricci_curvature(group, dist, x, y, S)`** — Computes κ(x, y) using the two-case Lin-Lu-Yau formula.

### 3. Computation & Output

- **`compute_results(...)`** — Iterates over a parameter range, computes curvatures for all generators and all generating sets $S$. Optionally auto-extends the parameter range until values stabilize.
- **`print_text_table / print_latex_table`** — Console / single LaTeX table renderers.
- **`to_fraction / to_latex`** — Convert float curvature values to readable fraction strings.

### 4. Group Registry & S-Configs

- **`GROUP_REGISTRY`** — Central registry mapping group names to their class, parameter name/values, and extension limits.
- **`auto_generate_S_configs(group_cls)`** — Automatically enumerates every non-empty subset of a group's generator set as a candidate $S$.

### 5. Combined Table Output

- **`build_combined_tables`** — Groups results by $|S|$ and collects all generator symbols.
- **`write_combined_latex`** — Writes multi-column combined LaTeX tables (one table per $|S|$ size).
- **`write_combined_markdown`** — Writes the same data as Markdown tables.

### 6. Main Entry Point

Orchestrates the full pipeline: compute → console summary → write LaTeX → write Markdown.

---

## Usage

```bash
python ricci_curvature_all.py
```

No command-line arguments are needed. To customize which groups or parameter ranges are computed,
edit **`GROUP_REGISTRY`** near the bottom of the file:

```python
GROUP_REGISTRY = {
    'D_n': {
        'class': DihedralGroup,
        'param_name': 'n',
        'param_values': [3, 4, 5],   # starting values
        'max_param': 10,             # auto-extend up to this
    },
    # add or remove entries here
}
```

To add a new group, implement the six-method interface, then add an entry to `GROUP_REGISTRY`.

---

## Output

After running, three outputs are produced:

### Console
A summary table printed to stdout for each group, showing how many S-configs and generators were processed:

```
######################################################################
#  D_n  —  6 tables
######################################################################
  |S|=1: 3 configs × 3 generators
  |S|=2: 3 configs × 3 generators
  |S|=3: 1 configs × 3 generators
```

### `ricci_tables.tex`
A complete, compilable LaTeX document containing combined multi-column tables.
Each table covers one group and one $|S|$ size (singletons / pairs / triples / …),
with columns for each parameter value and rows for each generator $g$, reporting $\kappa(e, g)$ as an exact fraction.
Columns where the curvature stabilizes are marked $n \geq k$.

### `ricci_tables.md`
The same data rendered as GitHub-flavored Markdown tables with LaTeX math,
grouped by group name and $|S|$.

---

## Dependencies

Standard library only — no third-party packages required.

```
Python ≥ 3.8
itertools, collections, fractions, math
```
