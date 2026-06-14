"""
Unified Ricci Curvature Calculator for Cayley Graphs of Finitely Generated Groups.

Based on: "Algorithm for Calculating Ricci Curvature" by Kevin Fung, May 2026.

Supported groups:
  - Q_{4m}  (generalized quaternion group)
  - D_n     (dihedral group)

To add a new group:
  1. Write a group class implementing:
     - __init__(self, param)   : constructor with group parameter
     - identity(self) -> elem  : return identity element
     - generator_a(self) -> elem
     - generator_b(self) -> elem
     - multiply(self, elem, gen) -> elem : 'a','A','b','B'
     - neighborhood(self, elem, S) -> list
     - elem_name(self, elem) -> str
     - self.n   : the modulus / order of generator a
     - self.order : total group order
  2. Add entry to GROUP_REGISTRY
  3. Add S_CONFIGS entries

Usage: python ricci_curvature_all.py
"""

from itertools import permutations
from collections import deque
from fractions import Fraction
import math

# ============================================================
# 1. Group Classes
# ============================================================

class Q4m:
    """Generalized quaternion group Q_{4m} = <a, b | a^{2m}=e, b^2=a^m, b^{-1}ab=a^{-1}>.

    Elements: (k, t) where k ∈ {0,...,2m-1}, t ∈ {0,1}.
      (k, 0) = a^k,  (k, 1) = a^k b
    """
    # (generator_symbol, LaTeX_name)
    GENERATORS = [('a', 'a'), ('A', 'a^{-1}'), ('b', 'b'), ('B', 'b^{-1}')]
    def __init__(self, m):
        if m < 2:
            raise ValueError("m must be >= 2")
        self.m = m
        self.n = 2 * m
        self.order = 4 * m

    def identity(self):
        return (0, 0)

    def generator_a(self):
        return (1, 0)

    def generator_b(self):
        return (0, 1)

    def multiply(self, elem, gen):
        k, t = elem
        n, m = self.n, self.m
        if gen == 'a':
            return ((k + 1) % n, 0) if t == 0 else ((k - 1) % n, 1)
        elif gen == 'A':
            return ((k - 1) % n, 0) if t == 0 else ((k + 1) % n, 1)
        elif gen == 'b':
            return (k, 1) if t == 0 else ((k + m) % n, 0)
        elif gen == 'B':
            return ((k + m) % n, 1) if t == 0 else (k, 0)
        else:
            raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        k, t = elem
        if t == 0:
            if k == 0:    return "e"
            elif k == 1: return "a"
            elif k == self.m: return f"a^{self.m}"
            else:        return f"a^{k}"
        else:
            if k == 0:    return "b"
            elif k == 1: return "ab"
            else:        return f"a^{k}b"

    def elements(self):
        """Iterate over all group elements."""
        for k in range(self.n):
            for t in (0, 1):
                yield (k, t)


class DihedralGroup:
    """Dihedral group D_n = <a, b | a^n=e, b^2=e, ba=a^{-1}b>.

    Elements: (k, t) where k ∈ {0,...,n-1}, t ∈ {0,1}.
      (k, 0) = a^k (rotation),  (k, 1) = a^k b (reflection)
    """
    # b is self-inverse in D_n, so only 3 distinct generators
    GENERATORS = [('a', 'a'), ('A', 'a^{-1}'), ('b', 'b')]
    def __init__(self, n):
        if n < 3:
            raise ValueError("n must be >= 3 for D_n")
        self.n = n
        self.order = 2 * n

    def identity(self):
        return (0, 0)

    def generator_a(self):
        return (1, 0)

    def generator_b(self):
        return (0, 1)

    def multiply(self, elem, gen):
        k, t = elem
        n = self.n
        if gen == 'a':
            return ((k + 1) % n, 0) if t == 0 else ((k - 1) % n, 1)
        elif gen == 'A':
            return ((k - 1) % n, 0) if t == 0 else ((k + 1) % n, 1)
        elif gen == 'b' or gen == 'B':  # b is self-inverse in D_n
            return (k, 1) if t == 0 else (k, 0)
        else:
            raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        k, t = elem
        if t == 0:
            if k == 0: return "e"
            elif k == 1: return "a"
            else: return f"a^{k}"
        else:
            if k == 0: return "b"
            elif k == 1: return "ab"
            else: return f"a^{k}b"

    def elements(self):
        for k in range(self.n):
            for t in (0, 1):
                yield (k, t)


class CyclicGroupZn:
    """Cyclic group Z_n with generating set S_{1,k} = {+1, +k, -1, -k}.

    Elements: integers 0, 1, ..., n-1 (addition mod n).
    Generators: 'a'=+1, 'A'=-1, 'b'=+k, 'B'=-k.

    Args:
        n: modulus (n ≥ 6)
        k: step size (k ≠ 1, k ≠ 0 mod n)
    """
    GENERATORS = [('a', '+1'), ('A', '-1'), ('b', '+k'), ('B', '-k')]

    def __init__(self, n, k=2):
        if n < 6:
            raise ValueError("n must be >= 6")
        self.n = n
        self.k = k % n
        self.order = n

    def identity(self):
        return 0

    def generator_a(self):
        return 1

    def generator_b(self):
        return self.k % self.n

    def multiply(self, elem, gen):
        x = elem
        n = self.n
        if gen == 'a':
            return (x + 1) % n
        elif gen == 'A':
            return (x - 1) % n
        elif gen == 'b':
            return (x + self.k) % n
        elif gen == 'B':
            return (x - self.k) % n
        else:
            raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        return str(elem)

    def elements(self):
        for i in range(self.n):
            yield i


class SymmetricGroup:
    """Symmetric group S_n with Coxeter generators (adjacent transpositions).

    Elements: tuples (p_1, ..., p_n) where p_i is the image of i.
    Generators: σ_i = (i, i+1) for i = 1, ..., n-1.  Each σ_i is self-inverse.
    """
    def __init__(self, n):
        if n < 2:
            raise ValueError("n must be >= 2")
        self.n = n
        self.order = math.factorial(n)
        # Dynamic generators: s_i = σ_i (self-inverse, so no separate inverse symbol)
        self.GENERATORS = [(f's{i}', fr'\sigma_{{{i}}}') for i in range(1, n)]

    def identity(self):
        return tuple(range(1, self.n + 1))

    def generator_a(self):
        return self.multiply(self.identity(), 's1')

    def generator_b(self):
        return self.multiply(self.identity(), 's2') if self.n >= 3 else self.identity()

    def multiply(self, elem, gen):
        if gen[0] == 's':
            i = int(gen[1:]) - 1  # 0-indexed
            if i < 0 or i >= self.n - 1:
                raise ValueError(f"Generator {gen} out of range for n={self.n}")
            p = list(elem)
            p[i], p[i + 1] = p[i + 1], p[i]
            return tuple(p)
        raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        # Represent permutation in cycle notation (simplified)
        return ''.join(str(x) for x in elem)

    def elements(self):
        """Yield all n! permutations in lexicographic order."""
        from itertools import permutations as perms
        for p in perms(range(1, self.n + 1)):
            yield p


class AlternatingGroup:
    """Alternating group A_n (even permutations) with Carmichael generators.

    Elements: even permutations of {1,...,n}, stored as tuples.
    Generators: V_i = (1, 2, i+2) for i = 1, ..., n-2.  Each V_i has order 3.
    V_i^{-1} = V_i^2 is distinct, represented by uppercase symbol.
    """
    def __init__(self, n):
        if n < 3:
            raise ValueError("n must be >= 3")
        self.n = n
        self.order = math.factorial(n) // 2
        # Dynamic generators: v_i = V_i, V_i = V_i^{-1}
        gens = []
        for i in range(1, n - 1):
            gens.append((f'v{i}', f'V_{{{i}}}'))
            gens.append((f'V{i}', f'V_{{{i}}}^{{-1}}'))
        self.GENERATORS = gens

    def identity(self):
        return tuple(range(1, self.n + 1))

    def generator_a(self):
        return self.multiply(self.identity(), 'v1')

    def generator_b(self):
        return self.multiply(self.identity(), 'v2') if self.n >= 4 else self.identity()

    def _is_even(self, perm):
        """Check if a permutation is even."""
        seen = [False] * len(perm)
        swaps = 0
        for i in range(len(perm)):
            if not seen[i]:
                j = i
                cycle_len = 0
                while not seen[j]:
                    seen[j] = True
                    j = perm[j] - 1  # perm is 1-indexed
                    cycle_len += 1
                swaps += cycle_len - 1
        return swaps % 2 == 0

    def _apply_3cycle(self, perm, a, b, c):
        """Apply 3-cycle (a, b, c) to perm (1-indexed)."""
        p = list(perm)
        # (a b c): a→b, b→c, c→a
        pa, pb, pc = p[a - 1], p[b - 1], p[c - 1]
        p[a - 1], p[b - 1], p[c - 1] = pc, pa, pb
        return tuple(p)

    def multiply(self, elem, gen):
        # V_i = (1, 2, i+2), V_i^{-1} = (1, i+2, 2)
        if gen[0] in ('v', 'V'):
            i = int(gen[1:])  # 1-indexed generator number
            if i < 1 or i > self.n - 2:
                raise ValueError(f"Generator {gen} out of range for n={self.n}")
            if gen[0] == 'v':  # V_i = (1, 2, i+2)
                return self._apply_3cycle(elem, 1, 2, i + 2)
            else:  # V_i^{-1} = (1, i+2, 2)
                return self._apply_3cycle(elem, 1, i + 2, 2)
        raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        return ''.join(str(x) for x in elem)

    def elements(self):
        """Yield all even permutations."""
        from itertools import permutations as perms
        for p in perms(range(1, self.n + 1)):
            if self._is_even(p):
                yield p


class GeneralLinearGroupZ:
    """GL_n(Z): invertible n×n integer matrices with det = ±1.

    For n = 1:  GL_1(Z) = {[1], [-1]} ≅ Z_2 (finite, order 2).
                Only generator is D = [-1] (self-inverse).
    For n ≥ 2:  **Infinite group.**  BFS will not terminate.

    Generators (n ≥ 2):
      'a' = E_{1,2}    'A' = E_{1,2}^{-1} (same as E_{1,2} with -1)
      'b' = E_{2,1}    'B' = E_{2,1}^{-1}
      'd' = D = diag(-1, 1, ..., 1)   (self-inverse)
    """

    def __init__(self, n):
        if n < 1:
            raise ValueError("n must be >= 1")
        self.n = n
        if n == 1:
            self.order = 2                       # finite: {[1], [-1]}
            self.GENERATORS = [('d', 'D')]        # D = [-1], self-inverse
        else:
            self.order = float('inf')
            self.GENERATORS = [
                ('a', 'E_{1,2}'), ('A', 'E_{1,2}^{-1}'),
                ('b', 'E_{2,1}'), ('B', 'E_{2,1}^{-1}'),
                ('d', 'D'),
            ]

    def identity(self):
        """Return n×n identity matrix as tuple of tuples."""
        return tuple(tuple(1 if i == j else 0 for j in range(self.n)) for i in range(self.n))

    def generator_a(self):
        return self.multiply(self.identity(), 'a') if self.n >= 2 else self.identity()

    def generator_b(self):
        return self.multiply(self.identity(), 'b') if self.n >= 2 else self.identity()

    def _mat_add(self, A, B):
        return tuple(tuple(A[i][j] + B[i][j] for j in range(self.n)) for i in range(self.n))

    def _mat_mul(self, A, B):
        """A * B (both n×n tuples of tuples)."""
        n = self.n
        if n == 1:
            return ((A[0][0] * B[0][0],),)
        C = [[0] * n for _ in range(n)]
        for i in range(n):
            for k in range(n):
                if A[i][k] != 0:
                    aik = A[i][k]
                    for j in range(n):
                        C[i][j] += aik * B[k][j]
        return tuple(tuple(row) for row in C)

    def _elem_matrix(self, i, j, val):
        """Return matrix with `val` at (i,j) and identity elsewhere."""
        I = list(list(row) for row in self.identity())
        I[i][j] = val
        return tuple(tuple(row) for row in I)

    def multiply(self, elem, gen):
        n = self.n
        if gen == 'd' or gen == 'D':
            # D = diag(-1, 1, ..., 1), self-inverse
            D = list(list(row) for row in self.identity())
            D[0][0] = -1
            return self._mat_mul(elem, tuple(tuple(row) for row in D))
        elif n == 1:
            raise ValueError(f"Unknown generator for n=1: {gen}")
        elif gen == 'a':
            E = self._elem_matrix(0, 1, 1)   # E_{1,2}
            return self._mat_mul(elem, E)
        elif gen == 'A':
            E = self._elem_matrix(0, 1, -1)  # E_{1,2}^{-1}
            return self._mat_mul(elem, E)
        elif gen == 'b':
            E = self._elem_matrix(1, 0, 1)   # E_{2,1}
            return self._mat_mul(elem, E)
        elif gen == 'B':
            E = self._elem_matrix(1, 0, -1)  # E_{2,1}^{-1}
            return self._mat_mul(elem, E)
        else:
            raise ValueError(f"Unknown generator: {gen}")

    def neighborhood(self, elem, S):
        return [self.multiply(elem, s) for s in S]

    def elem_name(self, elem):
        """Compact string representation of a matrix."""
        if self.n == 1:
            return str(elem[0][0])
        rows = [''.join(f'{x:2d}' for x in row) for row in elem]
        return '|' + ';'.join(rows) + '|'

    def elements(self):
        """Enumerate group elements (finite only for n=1)."""
        if self.n == 1:
            yield ((1,),)   # [1]
            yield ((-1,),)  # [-1]
        else:
            raise NotImplementedError(
                "GL_n(Z) for n≥2 is infinite; cannot enumerate all elements for BFS.\n"
                "Consider using a finite quotient (e.g. GL_n(F_p)) or a depth-limited BFS."
            )


# ============================================================
# 2. Core Algorithms (group-agnostic)
# ============================================================

def all_distances(group, S):
    """BFS all-pairs shortest directed path distances."""
    dist = {}
    for start in group.elements():
        visited = {start: 0}
        queue = deque([start])
        while queue:
            curr = queue.popleft()
            d = visited[curr]
            for s in S:
                nxt = group.multiply(curr, s)
                if nxt not in visited:
                    visited[nxt] = d + 1
                    queue.append(nxt)
        for target, d in visited.items():
            dist[(start, target)] = d
    return dist


def min_matching_cost(sources, targets, dist):
    """Minimum-weight perfect matching via brute-force permutations."""
    n = len(sources)
    if len(targets) != n:
        return float('inf')
    if n == 0:
        return 0
    dmat = [[dist.get((s, t), float('inf')) for t in targets] for s in sources]
    best = float('inf')
    for perm in permutations(range(n)):
        cost = sum(dmat[i][perm[i]] for i in range(n))
        if cost < best:
            best = cost
    return best


def ricci_curvature(group, dist, x, y, S):
    """Compute κ(x, y) using the Ollivier-type formula."""
    r = len(S)
    Nx = group.neighborhood(x, S)
    Ny = group.neighborhood(y, S)

    if x in Ny:
        # Case 2: x is in the out-neighborhood of y
        sources = [v for v in Nx if v != y]
        targets = [v for v in Ny if v != x]
        cost = min_matching_cost(sources, targets, dist)
        kappa = 1 - (cost - 1) / r
    else:
        # Case 1
        cost = min_matching_cost(Nx, Ny, dist)
        kappa = 1 - cost / r

    return kappa


# ============================================================
# 3. Computation & LaTeX Output
# ============================================================

def generator_label(group_cls, symbol, instance=None):
    """Get LaTeX label for a generator symbol from the group's GENERATORS."""
    # Try instance-level GENERATORS first (for groups with dynamic generators)
    if instance is not None and hasattr(instance, 'GENERATORS'):
        for sym, lbl in instance.GENERATORS:
            if sym == symbol:
                return lbl
    # Fall back to class-level GENERATORS
    for sym, lbl in group_cls.GENERATORS:
        if sym == symbol:
            return lbl
    return symbol  # fallback


def compute_results(group_cls, param_values, S_list, param_name,
                    auto_extend=True, max_param=10, stable_steps=3,
                    verbose=False, init_kwargs=None):
    """Compute κ(e, g) for ALL generator symbols g, for every S in S_list.

    Args:
        init_kwargs: extra kwargs passed to group constructor (e.g. k for CyclicGroupZn)
    """
    if init_kwargs is None:
        init_kwargs = {}
    results = []
    for S, latex_label in S_list:
        # Determine ALL generator symbols (not just those in S)
        # Create a temp instance to discover dynamic generators
        all_gen_symbols = []
        if hasattr(group_cls, 'GENERATORS'):
            all_gen_symbols = [sym for sym, _ in group_cls.GENERATORS]
        else:
            temp = group_cls(param_values[0], **init_kwargs)
            all_gen_symbols = [sym for sym, _ in temp.GENERATORS]

        kappa_by_sym = {g: [] for g in all_gen_symbols}
        stable_at = None

        for p in param_values:
            g = group_cls(p, **init_kwargs)
            e = g.identity()
            dist = all_distances(g, S)
            if verbose:
                print(f"\n  [{param_name}={p}] S={S}")
                N_e = g.neighborhood(e, S)
                print(f"    N_out(e) = {[g.elem_name(x) for x in N_e]}")
            for gen_sym in all_gen_symbols:
                gen_elem = g.multiply(e, gen_sym)
                if verbose:
                    N_s = g.neighborhood(gen_elem, S)
                    print(f"    N_out({g.elem_name(gen_elem)}) = {[g.elem_name(x) for x in N_s]}")
                kappa_by_sym[gen_sym].append(ricci_curvature(g, dist, e, gen_elem, S))

        all_params = list(param_values)

        # Auto-extend
        if auto_extend and param_values:
            p = param_values[-1] + 1
            while p <= max_param:
                g = group_cls(p, **init_kwargs)
                e = g.identity()
                dist = all_distances(g, S)
                for gen_sym in all_gen_symbols:
                    gen_elem = g.multiply(e, gen_sym)
                    kappa_by_sym[gen_sym].append(ricci_curvature(g, dist, e, gen_elem, S))
                all_params.append(p)

                # Check if ALL rows have stabilized
                if len(kappa_by_sym[all_gen_symbols[0]]) >= stable_steps:
                    all_stable = True
                    for gen_sym in all_gen_symbols:
                        last_vals = kappa_by_sym[gen_sym][-stable_steps:]
                        if not all(v == last_vals[0] for v in last_vals):
                            all_stable = False
                            break
                    if all_stable:
                        stable_at = p - stable_steps + 1
                        trim_point = len(all_params) - stable_steps + 1
                        all_params = all_params[:trim_point]
                        for gen_sym in all_gen_symbols:
                            kappa_by_sym[gen_sym] = kappa_by_sym[gen_sym][:trim_point]
                        break
                p += 1

        # Build rows for ALL generator symbols
        label_instance = group_cls(all_params[0], **init_kwargs) if all_params else group_cls(param_values[0], **init_kwargs)
        rows = []
        for gen_sym in all_gen_symbols:
            lbl = generator_label(group_cls, gen_sym, instance=label_instance)
            rows.append((gen_sym, lbl, kappa_by_sym[gen_sym]))

        results.append({
            'latex_label': latex_label,
            'param_values': all_params,
            'rows': rows,
            'stable_at': stable_at,
            'S_size': len(S),          # original |S| for grouping
        })
    return results


def to_fraction(val):
    """Convert float to plain-text fraction string."""
    if math.isinf(val):
        return "-∞" if val < 0 else "∞"
    frac = Fraction(val).limit_denominator(1000)
    if frac.denominator == 1:
        return str(frac.numerator)
    else:
        return f"{frac.numerator}/{frac.denominator}"


def print_text_table(result, group_name, param_name, table_num=None):
    """Print a plain-text table to console with dynamic rows."""
    vals = result['param_values']
    stable_at = result.get('stable_at')
    rows = result['rows']
    n_cols = len(vals)
    width = 70
    label = result['latex_label']

    print()
    print("=" * width)
    num_str = f"Table {table_num}: " if table_num else ""
    if stable_at is not None:
        print(f"  {num_str}Γ({group_name}, {label})  [stabilizes at {param_name}≥{stable_at}]")
    else:
        print(f"  {num_str}Γ({group_name}, {label})  [up to {param_name}={vals[-1]}]")
    print("=" * width)

    # Header
    header = f"{'κ(x,y)':<14}"
    for i, v in enumerate(vals):
        if stable_at is not None and i == n_cols - 1:
            header += f" {param_name}≥{v:<6}"
        else:
            header += f" {param_name}={v:<6}"
    print(header)
    print("-" * width)

    # One row per generator in S
    for sym, lbl, kappa_vals in rows:
        row_label = f"κ(e,{lbl})"
        row_str = f"{row_label:<14}"
        for val in kappa_vals:
            row_str += f"{to_fraction(val):<8}"
        print(row_str)
    print("-" * width)


def to_latex(val):
    """Convert float to LaTeX math string: $\\frac{a}{b}$, $\\infty$, etc."""
    if math.isinf(val):
        return r"$-\infty$" if val < 0 else r"$\infty$"
    frac = Fraction(val).limit_denominator(1000)
    num, den = abs(frac.numerator), frac.denominator
    sign = "-" if frac.numerator < 0 else ""
    if den == 1:
        return f"${sign}{num}$"
    else:
        return f"${sign}\\frac{{{num}}}{{{den}}}$"


def print_latex_table(result, group_name, param_name, file=None):
    """Write a single result dict as a LaTeX table with dynamic rows."""
    vals = result['param_values']
    stable_at = result.get('stable_at')
    rows = result['rows']
    n_cols = len(vals)

    col_spec = "|l|" + "c|" * n_cols
    print(r"\begin{center}", file=file)
    print(r"\begin{tabular}{" + col_spec + "}", file=file)
    print(r"\hline", file=file)

    title = rf"\Gamma({group_name},{result['latex_label']})"
    print(r" & \multicolumn{" + str(n_cols) + r"}{|c|}{$" + title + r"$} \\", file=file)
    print(r"\hline", file=file)

    headers = r"$\kappa(x,y)$"
    for i, v in enumerate(vals):
        if stable_at is not None and i == n_cols - 1:
            headers += f" & ${param_name}\\geq {v}$"
        else:
            headers += f" & ${param_name}={v}$"
    print(headers + r" \\", file=file)
    print(r"\hline", file=file)

    for sym, lbl, kappa_vals in rows:
        row_label = rf"$\kappa(e,{lbl})$"
        row_str = row_label
        for val in kappa_vals:
            row_str += " & " + to_latex(val)
        print(row_str + r" \\", file=file)
        print(r"\hline", file=file)

    print(r"\end{tabular}", file=file)
    print(r"\end{center}", file=file)


# ============================================================
# 4. Group Registry & Generating Set Configurations
# ============================================================

GROUP_REGISTRY = {
    'Q_{4m}': {
        'class': Q4m,
        'param_name': 'm',
        'param_values': [2, 3, 4, 5],
        'max_param': 10,
    },
    'D_n': {
        'class': DihedralGroup,
        'param_name': 'n',
        'param_values': [3, 4, 5],
        'max_param': 10,
    },
    'Z_n (k=2)': {
        'class': CyclicGroupZn,
        'param_name': 'n',
        'param_values': [6, 7, 8, 9],
        'max_param': 15,
        'init_kwargs': {'k': 2},
    },
    'S_n': {
        'class': SymmetricGroup,
        'param_name': 'n',
        'param_values': [3, 4],
        'max_param': 5,
    },
    'A_n': {
        'class': AlternatingGroup,
        'param_name': 'n',
        'param_values': [4, 5],
        'max_param': 6,
    },
    # GL_n(Z): finite only for n=1 (GL_1(Z) ≅ Z_2).
    # For n≥2 the group is infinite — disabled by default.
    'GL_n(Z)': {
        'class': GeneralLinearGroupZ,
        'param_name': 'n',
        'param_values': [1],
        'max_param': 1,
        'auto_extend': False,
    },
}

# S configurations: (S_array, LaTeX label)
# Groups listed here use manually curated subsets; all others auto-generate
# ALL non-empty subsets of their GENERATORS via auto_generate_S_configs().
S_CONFIGS = {
    # All groups now auto-generate.  Keep this dict for future manual overrides.
}


def auto_generate_S_configs(group_cls, instance=None):
    """Auto-generate all non-empty subsets of a group's generators.

    Args:
        group_cls: group class with GENERATORS = [(symbol, latex_name), ...]
        instance: optional instance for groups with dynamic (per-parameter) GENERATORS

    Returns:
        list of (S_array, latex_label) for all non-empty subsets
    """
    from itertools import combinations as comb

    # Use instance-level GENERATORS if available, else class-level
    if instance is not None and hasattr(instance, 'GENERATORS'):
        generators = instance.GENERATORS
    else:
        generators = group_cls.GENERATORS
    n = len(generators)
    configs = []
    for r in range(1, n + 1):
        for indices in comb(range(n), r):
            symbols = [generators[i][0] for i in indices]
            labels = [generators[i][1] for i in indices]
            latex_label = r"\{" + ",".join(labels) + r"\}"
            configs.append((symbols, latex_label))
    return configs


# ============================================================
# 5. Combined Table Output (LaTeX + Markdown)
# ============================================================

def param_header(param_name, v, is_last, stable_at, idx, n_cols):
    """Format a parameter column header: 'm=2' or 'm≥6'."""
    if stable_at is not None and idx == n_cols - 1:
        return f"${param_name}\\geq {v}$"
    else:
        return f"${param_name}={v}$"


def param_header_md(param_name, v, is_last, stable_at, idx, n_cols):
    """Markdown parameter header."""
    if stable_at is not None and idx == n_cols - 1:
        return f"${param_name}\\ge {v}$"
    else:
        return f"${param_name}={v}$"


def build_combined_tables(results, param_name):
    """Group results by |S| and collect all generator symbols with labels.

    Returns: list of dicts, one per size group:
        { 'size': int,
          'configs': [ { 'label': str, 'param_vals': [...], 'stable_at': int|None,
                          'data': {gen_sym: [κ,...]} }, ... ],
          'all_gens': [gen_sym, ...],
          'gen_labels': {gen_sym: latex_label, ...} }
    """
    from collections import OrderedDict

    # Group by |S|
    by_size = OrderedDict()
    for r in results:
        sz = r.get('S_size', len(r['rows']))
        by_size.setdefault(sz, []).append(r)

    groups = []
    for sz, configs in by_size.items():
        # Collect all generator symbols and their LaTeX labels
        gen_order = []
        gen_labels = {}
        seen_gens = set()
        for cfg in configs:
            for sym, lbl, kvals in cfg['rows']:
                if sym not in seen_gens:
                    seen_gens.add(sym)
                    gen_order.append(sym)
                    gen_labels[sym] = lbl

        # Build per-config data dict
        cfgs_out = []
        for cfg in configs:
            data = {}
            for sym, lbl, kvals in cfg['rows']:
                data[sym] = list(kvals)
            cfgs_out.append({
                'label': cfg['latex_label'],
                'param_vals': list(cfg['param_values']),
                'stable_at': cfg.get('stable_at'),
                'data': data,
            })

        groups.append({
            'size': sz,
            'configs': cfgs_out,
            'all_gens': gen_order,
            'gen_labels': gen_labels,
        })
    return groups


def write_combined_latex(group_name, param_name, size_groups, file):
    """Write one combined LaTeX table per size group."""
    SIZE_NAMES = {1: "singletons", 2: "pairs", 3: "triples", 4: "all four"}

    for grp in size_groups:
        sz = grp['size']
        configs = grp['configs']
        all_gens = grp['all_gens']
        size_label = SIZE_NAMES.get(sz, f"|S|={sz}")

        # Calculate total columns: 1 (row label) + sum of param cols per config
        total_cols = 1 + sum(len(c['param_vals']) for c in configs)
        col_spec = "|c|" + "c|" * (total_cols - 1)

        print(r"\begin{center}", file=file)
        print(r"\begin{tabular}{" + col_spec + "}", file=file)
        print(r"\hline", file=file)

        # Title row
        title = rf"\Gamma({group_name}, -)\quad (|S|={sz}: {size_label})"
        print(r"\multicolumn{" + str(total_cols) + r"}{|c|}{$" + title + r"$} \\", file=file)
        print(r"\hline", file=file)

        # Header row 1: κ(x,y) (multirow=2) | S labels (multicolumn)
        header1 = r"\multirow{2}{*}{$\kappa(x,y)$}"
        for c in configs:
            n = len(c['param_vals'])
            header1 += r" & \multicolumn{" + str(n) + r"}{c|}{$" + c['label'] + r"$}"
        print(header1 + r" \\", file=file)

        # Header row 2 (cline from col 2): param values
        header2 = ""
        for c in configs:
            n = len(c['param_vals'])
            for idx, v in enumerate(c['param_vals']):
                header2 += " & " + param_header(param_name, v, (idx == n - 1),
                                                 c['stable_at'], idx, n)
        print(r"\cline{2-" + str(total_cols) + "}", file=file)
        print(header2 + r" \\", file=file)
        print(r"\hline", file=file)

        # Data rows
        for gen_sym in all_gens:
            gen_lbl = grp['gen_labels'].get(gen_sym, gen_sym)
            row = rf"$\kappa(e,{gen_lbl})$"
            for c in configs:
                if gen_sym in c['data']:
                    for val in c['data'][gen_sym]:
                        row += " & " + to_latex(val)
                else:
                    for _ in c['param_vals']:
                        row += r" & $-$"
            print(row + r" \\", file=file)
            print(r"\hline", file=file)

        print(r"\end{tabular}", file=file)
        print(r"\end{center}", file=file)
        print(file=file)


def write_combined_markdown(group_name, param_name, size_groups, file):
    """Write one combined Markdown table per size group."""
    SIZE_NAMES = {1: "singletons", 2: "pairs", 3: "triples", 4: "all four"}

    for grp in size_groups:
        sz = grp['size']
        configs = grp['configs']
        all_gens = grp['all_gens']
        size_label = SIZE_NAMES.get(sz, f"|S|={sz}")

        total_cols = 1 + sum(len(c['param_vals']) for c in configs)

        file.write(f"\n### Γ({group_name}, −)  — |S|={sz} ({size_label})\n\n")

        # Header row 1: S labels (shown only once per S group; rest blank)
        header1 = r"| $\kappa(x,y)$ |"
        for c in configs:
            n = len(c['param_vals'])
            header1 += f" {c['label']} |"
            for _ in range(n - 1):
                header1 += " |"
        file.write(header1 + "\n")

        # Header row 2: parameter values
        header2 = "| |"
        for c in configs:
            n = len(c['param_vals'])
            for idx, v in enumerate(c['param_vals']):
                p = param_header_md(param_name, v, (idx == n - 1),
                                    c['stable_at'], idx, n)
                header2 += f" {p} |"
        file.write(header2 + "\n")

        # Separator
        sep = "|:---:|" + ":---:|" * (total_cols - 1)
        file.write(sep + "\n")

        # Data rows
        for gen_sym in all_gens:
            gen_lbl = grp['gen_labels'].get(gen_sym, gen_sym)
            row = f"| **$\\kappa(e,{gen_lbl})$** |"
            for c in configs:
                if gen_sym in c['data']:
                    for val in c['data'][gen_sym]:
                        row += f" {to_latex(val)} |"
                else:
                    for _ in c['param_vals']:
                        row += " $-$ |"
            file.write(row + "\n")

        file.write("\n")


# ============================================================
# 6. Main
# ============================================================

if __name__ == "__main__":
    latex_path = r"c:\Users\User\group\ricci_tables.tex"
    md_path = r"c:\Users\User\group\ricci_tables.md"

    # Collect all results first (compute once)
    all_results = {}
    for group_name, info in GROUP_REGISTRY.items():
        group_cls = info['class']
        param_name = info['param_name']
        param_values = info['param_values']
        max_param = info.get('max_param', 10)
        init_kwargs = info.get('init_kwargs', {})

        # Use manual S_CONFIGS if provided, else auto-generate all subsets
        if group_name in S_CONFIGS:
            s_configs = S_CONFIGS[group_name]
            print(f"\n[{group_name}] Using {len(s_configs)} manually curated S-configs")
        else:
            # For groups with dynamic (instance-level) GENERATORS, create a
            # temp instance at the smallest parameter value to discover them.
            if hasattr(group_cls, 'GENERATORS'):
                temp_inst = None
            else:
                temp_inst = group_cls(param_values[0], **init_kwargs)
            s_configs = auto_generate_S_configs(group_cls, instance=temp_inst)
            print(f"\n[{group_name}] Auto-generated {len(s_configs)} S-configs from GENERATORS")

        all_results[group_name] = {
            'param_name': param_name,
            'results': compute_results(group_cls, param_values, s_configs,
                                       param_name, max_param=max_param,
                                       auto_extend=info.get('auto_extend', True),
                                       init_kwargs=init_kwargs),
        }

    # Print summary to console
    for group_name, data in all_results.items():
        print(f"\n{'#'*70}")
        print(f"#  {group_name}  —  {len(data['results'])} tables")
        print(f"{'#'*70}")
        size_groups = build_combined_tables(data['results'], data['param_name'])
        for grp in size_groups:
            print(f"  |S|={grp['size']}: {len(grp['configs'])} configs × "
                  f"{len(grp['all_gens'])} generators")

    # --- Write LaTeX ---
    with open(latex_path, 'w', encoding='utf-8') as f:
        f.write("%" * 70 + "\n")
        f.write("%  Ricci Curvature — Combined LaTeX Tables\n")
        f.write("%  Generated by ricci_curvature_all.py\n")
        f.write("%" * 70 + "\n\n")
        f.write(r"\documentclass{article}" + "\n")
        f.write(r"\usepackage{amsmath}" + "\n")
        f.write(r"\usepackage{multirow}" + "\n")
        f.write(r"\begin{document}" + "\n\n")

        for group_name, data in all_results.items():
            f.write(f"\n% {'='*60}\n")
            f.write(f"%  {group_name}\n")
            f.write(f"% {'='*60}\n\n")
            size_groups = build_combined_tables(data['results'], data['param_name'])
            write_combined_latex(group_name, data['param_name'], size_groups, file=f)

        f.write(r"\end{document}" + "\n")

    # --- Write Markdown ---
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write("# Ricci Curvature — Combined Tables\n\n")
        f.write("*Generated by `ricci_curvature_all.py`*\n\n")

        for group_name, data in all_results.items():
            f.write(f"## {group_name}\n\n")
            size_groups = build_combined_tables(data['results'], data['param_name'])
            write_combined_markdown(group_name, data['param_name'], size_groups, file=f)

    print(f"\nDone!  LaTeX → {latex_path}")
    print(f"       Markdown → {md_path}")
