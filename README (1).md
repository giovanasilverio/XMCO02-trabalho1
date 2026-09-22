# README — Minimum s-t Cut DIMACS Instances

## File Format

Each file is a plain-text ASCII file. Lines are read top-to-bottom. Leading whitespace is ignored. Every line begins with a **single-character descriptor** that determines how the rest of the line is parsed.

### Line Types

| Descriptor | Meaning | Format |
|------------|---------|--------|
| `c` | Comment | `c <free text>` |
| `p` | Problem line (exactly one, before any `a`/`n` lines) | `p min <n> <m>` |
| `n` | Node designation (source / sink) | `n <node> s` or `n <node> t` |
| `a` | Arc (directed edge with capacity) | `a <from> <to> <capacity>` |

### Field Descriptions

#### `c` — Comment
Free-form text. Ignored by parsers. Used here to describe the instance structure and node-numbering conventions.

#### `p` — Problem Line
```
p min/mcf <n> <m>
```
- `min` is a literal token identifying the problem as min-cut.
- `mcf` is a literal token identifying the problem as max-flow multicommodity.
- `<n>` = total number of vertices.
- `<m>` = total number of directed arcs.

The `p` line **must appear before any `a` lines** so the solver can pre-allocate memory.

#### `n` — Node Designation
```
n <node> s     ← source
n <node> s 1   ← source for commodity 1
n <node> t     ← sink/target
n <node> t 1   ← sink/targe for commodity 1
```
- `<node>` is the vertex index.
- Exactly one source and one sink per instance (minimum cut).
- For multicommodity max flow, we have one source/sink per commodity/product.

#### `a` — Arc
```
a <from> <to> <capacity>
```
- `<from>`, `<to>` are vertex indices.
- `<capacity>` is a non-negative integer.
- Consider every arc as **undirected** from `<from>` to `<to>`, i.e. it goes both ways.

---

## Vertex Numbering Convention

**All instances use 1-based node numbering**, as required by the DIMACS convention:

- Node `1` is always the **source** `s`.
- The highest-numbered node is always the **sink** `t`.
- All intermediate nodes are labeled `2 … n-1`.

This differs from the 0-based numbering sometimes used in algorithm textbooks. If your code uses 0-based indexing internally, subtract 1 from every node ID when reading.

---

## Capacity Semantics

- Capacities are **integers** in the range `1 … 40` in these instances.
- Capacity `0` arcs are not present (they would be no-ops and are omitted).
- A minimum cut is reported as a partition `(S, T)` with `s ∈ S`, `t ∈ T`, and
  `value = Σ capacity(u,v)` over all `u ∈ S, v ∈ T`.

---

## File Naming

```
instance<N>.min
```

`.min` is the conventional extension for DIMACS min-cut files.
