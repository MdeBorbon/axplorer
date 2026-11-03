import math
from collections import Counter
from itertools import combinations
from logging import getLogger

import numpy as np

from src.envs.environment import BaseEnvironment, DataPoint
from src.envs.tokenizers import SparseTokenizerSingleInteger

logger = getLogger()

_PROJECTIVE_TOL = 1e-9
_PROJECTIVE_DECIMALS = 7


def _normalize_projective(v):
    """Return a stable affine representative of a nonzero vector in CP^2."""
    v = np.asarray(v, dtype=np.complex128).copy()
    norm = np.linalg.norm(v)
    if norm < _PROJECTIVE_TOL:
        return None

    v[np.abs(v) < _PROJECTIVE_TOL * norm] = 0
    pivot = next((i for i, z in enumerate(v) if abs(z) > _PROJECTIVE_TOL * norm), None)
    if pivot is None:
        return None

    v /= v[pivot]
    v[np.abs(v) < _PROJECTIVE_TOL] = 0
    return v


def _projective_key(v):
    """
    Hashable approximate projective key.

    Complex representatives remain unrounded for subsequent calculations;
    rounding is used only to identify equal projective points or lines.
    """
    v = _normalize_projective(v)
    if v is None:
        return None

    cutoff = 10 ** (-_PROJECTIVE_DECIMALS)
    key = []
    for z in v:
        re = 0.0 if abs(z.real) < cutoff else round(float(z.real), _PROJECTIVE_DECIMALS)
        im = 0.0 if abs(z.imag) < cutoff else round(float(z.imag), _PROJECTIVE_DECIMALS)
        key.append((re, im))
    return tuple(key)


def reflection_arrangement(m, family):
    """
    Rank-three imprimitive complex reflection arrangements.

    ``gmm`` is G(m,m,3):
        x = zeta^a y, y = zeta^a z, z = zeta^a x.

    ``gm1`` is G(m,1,3), obtained by adding x=0, y=0, z=0.
    """
    if m < 2:
        raise ValueError("reflection_m must be at least 2")
    if family not in {"gmm", "gm1"}:
        raise ValueError("seed_family must be 'gmm' or 'gm1'")

    zeta = np.exp(2j * np.pi / m)
    lines = []

    if family == "gm1":
        lines.extend(
            [
                np.array([1, 0, 0], dtype=np.complex128),
                np.array([0, 1, 0], dtype=np.complex128),
                np.array([0, 0, 1], dtype=np.complex128),
            ]
        )

    for a in range(m):
        w = zeta**a
        lines.extend(
            [
                np.array([1, -w, 0], dtype=np.complex128),
                np.array([0, 1, -w], dtype=np.complex128),
                np.array([-w, 0, 1], dtype=np.complex128),
            ]
        )

    return [_normalize_projective(line) for line in lines]


def arrangement_incidence(lines):
    """
    Compute vertices of a concrete projective line arrangement.

    Returns the number of vertices on each line and a dictionary whose values
    are ``{"point": representative, "lines": incident-line-indices}``.
    """
    vertices = {}

    for i, j in combinations(range(len(lines)), 2):
        point = _normalize_projective(np.cross(lines[i], lines[j]))
        key = _projective_key(point)
        if key is None:
            continue
        rec = vertices.setdefault(key, {"point": point, "lines": set()})
        rec["lines"].update((i, j))

    line_vertices = [set() for _ in lines]
    for key, rec in vertices.items():
        for i in rec["lines"]:
            line_vertices[i].add(key)

    counts = np.array([len(v) for v in line_vertices], dtype=np.int64)
    return counts, vertices


def hirzebruch_defect(lines, min_lines=6):
    """
    Sum of squared defects delta_L = 3*v(L)-N-3.

    For a nontrivial arrangement this vanishes exactly when N=3k and every
    line contains k+1 vertices.
    """
    n = len(lines)
    if n < min_lines:
        return 1_000_000 + (min_lines - n) ** 2

    counts, _ = arrangement_incidence(lines)
    delta = 3 * counts - n - 3
    return int(delta @ delta)


def multiplicity_profile(lines):
    """Return the intersection multiplicities t_r."""
    _, vertices = arrangement_incidence(lines)
    return Counter(len(rec["lines"]) for rec in vertices.values())


def build_candidate_pool(seed_lines):
    """
    Seed lines plus every line through a pair of vertices of the seed.

    The seed lines occupy the first coordinates of the resulting pool.
    """
    _, vertices = arrangement_incidence(seed_lines)
    points = [rec["point"] for rec in vertices.values()]

    pool = {}
    for line in seed_lines:
        pool[_projective_key(line)] = line

    for p, q in combinations(points, 2):
        line = _normalize_projective(np.cross(p, q))
        key = _projective_key(line)
        if key is not None and key not in pool:
            pool[key] = line

    return list(pool.values())


def build_intersection_table(lines):
    """
    Precompute the projective intersection ID of every pair of candidate lines.

    Local search then evaluates a proposed arrangement using only small integer
    arrays, rather than repeating complex projective cross-products.
    """
    n = len(lines)
    table = np.full((n, n), -1, dtype=np.int32)
    vertex_ids = {}

    for i, j in combinations(range(n), 2):
        key = _projective_key(np.cross(lines[i], lines[j]))
        if key not in vertex_ids:
            vertex_ids[key] = len(vertex_ids)
        vid = vertex_ids[key]
        table[i, j] = vid
        table[j, i] = vid

    return table


class HirzebruchMutationDataPoint(DataPoint):
    """
    Variable-size arrangements obtained by mutating a reflection arrangement.

    The candidate pool is the seed together with lines through pairs of seed
    vertices. Moves add, remove, or replace triples of lines.
    """

    REFLECTION_M = 4
    SEED_FAMILY = "gmm"
    SIZE_RADIUS = 3
    MIN_LINES = 6
    LOCAL_SEARCH_STEPS = 500
    INITIAL_MUTATIONS = 2
    ANNEAL_TEMPERATURE = 20.0
    CANDIDATE_TOP_K = 12
    GREEDY_MOVE_PROBABILITY = 0.35

    CANDIDATE_LINES = None
    INTERSECTION_TABLE = None
    SEED_COUNT = None
    POOL_SIZE = None
    SIZE_MIN = None
    SIZE_MAX = None
    ADD_PRIORITY = None
    DELETE_PRIORITY = None

    @classmethod
    def configure(
        cls,
        reflection_m,
        seed_family,
        size_radius,
        min_lines,
        local_search_steps,
        initial_mutations,
        anneal_temperature,
        candidate_top_k,
        greedy_move_probability,
    ):
        cls.REFLECTION_M = reflection_m
        cls.SEED_FAMILY = seed_family
        cls.SIZE_RADIUS = size_radius
        cls.MIN_LINES = min_lines
        cls.LOCAL_SEARCH_STEPS = local_search_steps
        cls.INITIAL_MUTATIONS = initial_mutations
        cls.ANNEAL_TEMPERATURE = anneal_temperature
        cls.CANDIDATE_TOP_K = candidate_top_k
        cls.GREEDY_MOVE_PROBABILITY = greedy_move_probability

        seed = reflection_arrangement(reflection_m, seed_family)
        cls.CANDIDATE_LINES = build_candidate_pool(seed)
        cls.SEED_COUNT = len(seed)
        cls.POOL_SIZE = len(cls.CANDIDATE_LINES)
        cls.INTERSECTION_TABLE = build_intersection_table(cls.CANDIDATE_LINES)

        raw_min = max(min_lines, cls.SEED_COUNT - size_radius)
        raw_max = min(cls.POOL_SIZE, cls.SEED_COUNT + size_radius)
        cls.SIZE_MIN = raw_min + ((cls.SEED_COUNT - raw_min) % 3)
        cls.SIZE_MAX = raw_max - ((raw_max - cls.SEED_COUNT) % 3)

        cls._build_seed_guidance()

    @classmethod
    def _vertex_line_sets(cls, selected):
        """Map intersection IDs to the selected candidate-line indices through them."""
        result = {}
        selected = list(selected)
        for a, i in enumerate(selected):
            for j in selected[a + 1 :]:
                vid = int(cls.INTERSECTION_TABLE[i, j])
                result.setdefault(vid, set()).update((i, j))
        return result

    @classmethod
    def _build_seed_guidance(cls):
        """
        Static geometric prior computed from the starting reflection arrangement.

        Addition priority rewards a candidate whose intersections with many seed
        lines are already seed vertices. Deletion priority rewards seed lines
        containing many double vertices. For G(m,m,3) and G(m,1,3), respectively,
        these criteria single out the three coordinate lines and make the known
        three-line extension/deletion an explicit validation target.
        """
        seed = list(range(cls.SEED_COUNT))
        seed_vertices = cls._vertex_line_sets(seed)
        seed_vertex_ids = set(seed_vertices)

        add_scores = np.zeros(cls.POOL_SIZE, dtype=np.int32)
        for candidate in range(cls.POOL_SIZE):
            reused = 0
            for i in seed:
                if i == candidate:
                    continue
                if int(cls.INTERSECTION_TABLE[candidate, i]) in seed_vertex_ids:
                    reused += 1
            add_scores[candidate] = reused

        delete_scores = np.zeros(cls.POOL_SIZE, dtype=np.int32)
        double_vertex_ids = {vid for vid, incident in seed_vertices.items() if len(incident) == 2}
        for i in seed:
            vids = {int(cls.INTERSECTION_TABLE[i, j]) for j in seed if j != i}
            delete_scores[i] = sum(vid in double_vertex_ids for vid in vids)

        cls.ADD_PRIORITY = sorted(range(cls.POOL_SIZE), key=lambda i: (-int(add_scores[i]), i))
        cls.DELETE_PRIORITY = sorted(range(cls.POOL_SIZE), key=lambda i: (-int(delete_scores[i]), i))

    def __init__(self, N, init=False):
        super().__init__()

        if self.CANDIDATE_LINES is None:
            raise RuntimeError("HirzebruchMutationDataPoint has not been configured")
        if N != self.POOL_SIZE:
            raise ValueError(f"N={N} does not match mutation pool size {self.POOL_SIZE}")

        self.N = N
        self.data = np.zeros(self.N, dtype=np.uint8)

        if init:
            self.data[: self.SEED_COUNT] = 1
            self._randomize_from_seed()
            self.calc_features()
            self.calc_score()

    def selected_indices(self):
        return np.flatnonzero(self.data).astype(int).tolist()

    def novelty(self):
        seed = np.zeros(self.N, dtype=np.uint8)
        seed[: self.SEED_COUNT] = 1
        return int(np.count_nonzero(self.data != seed))

    def _vertex_counts(self, selected):
        selected = np.asarray(selected, dtype=np.int32)
        if selected.size == 0:
            return np.array([], dtype=np.int64)

        sub = self.INTERSECTION_TABLE[np.ix_(selected, selected)]
        return np.fromiter(
            (np.unique(row[row >= 0]).size for row in sub),
            dtype=np.int64,
            count=len(selected),
        )

    def defect(self):
        selected = self.selected_indices()
        n = len(selected)
        if n < self.MIN_LINES:
            return 1_000_000 + (self.MIN_LINES - n) ** 2

        counts = self._vertex_counts(selected)
        delta = 3 * counts - n - 3
        return int(delta @ delta)

    def multiplicity_profile(self):
        selected = self.selected_indices()
        vertices = self._vertex_line_sets(selected)
        return Counter(len(incident) for incident in vertices.values())

    def _defect_bound(self):
        n = max(self.SIZE_MAX, self.MIN_LINES)
        max_abs_delta = 2 * n + 3
        return n * max_abs_delta * max_abs_delta + 1_000_000

    def calc_score(self):
        defect = self.defect()
        novelty = self.novelty()

        # Defect is lexicographically primary. Among equal defects, reward
        # distance from the starting reflection arrangement.
        defect_weight = self.POOL_SIZE + 1
        offset = defect_weight * (self._defect_bound() + 1) + self.POOL_SIZE
        self.score = int(offset - defect_weight * defect + novelty)

    def calc_features(self):
        self.features = ",".join(map(str, self.selected_indices()))

    def _allowed_sizes(self):
        return list(range(self.SIZE_MIN, self.SIZE_MAX + 1, 3))

    def _project_size(self):
        """Repair arbitrary transformer output to the nearest allowed line count."""
        rng = np.random.default_rng()
        selected = set(self.selected_indices())
        target = min(
            self._allowed_sizes(),
            key=lambda s: (abs(s - len(selected)), abs(s - self.SEED_COUNT)),
        )

        if len(selected) > target:
            nonseed = [i for i in selected if i >= self.SEED_COUNT]
            seed = [i for i in selected if i < self.SEED_COUNT]
            rng.shuffle(nonseed)
            rng.shuffle(seed)
            selected.difference_update((nonseed + seed)[: len(selected) - target])

        elif len(selected) < target:
            missing_seed = [i for i in range(self.SEED_COUNT) if i not in selected]
            other = [i for i in range(self.SEED_COUNT, self.N) if i not in selected]
            rng.shuffle(missing_seed)
            rng.shuffle(other)
            selected.update((missing_seed + other)[: target - len(selected)])

        self.data[:] = 0
        if selected:
            self.data[list(selected)] = 1

    def _choose_three(self, priority, eligible, rng, guided):
        eligible = set(eligible)
        ranked = [i for i in priority if i in eligible]
        if len(ranked) < 3:
            return None

        top_k = max(3, min(self.CANDIDATE_TOP_K, len(ranked)))
        candidates = ranked[:top_k]

        if guided and rng.random() < self.GREEDY_MOVE_PROBABILITY:
            return candidates[:3]
        return list(map(int, rng.choice(candidates, size=3, replace=False)))

    def _propose_move(self, rng, guided=True):
        selected = set(self.selected_indices())
        move_types = []

        if len(selected) + 3 <= self.SIZE_MAX:
            move_types.append("add")
        if len(selected) - 3 >= self.SIZE_MIN:
            move_types.append("remove")
        if len(selected) >= 3 and self.N - len(selected) >= 3:
            move_types.append("swap")

        if not move_types:
            return None

        move_type = str(rng.choice(move_types))
        proposal = set(selected)

        if move_type in {"remove", "swap"}:
            priority = self.DELETE_PRIORITY if guided else list(proposal)
            if not guided:
                rng.shuffle(priority)
            remove = self._choose_three(priority, proposal, rng, guided)
            if remove is None:
                return None
            proposal.difference_update(remove)

        if move_type in {"add", "swap"}:
            eligible = set(range(self.N)) - proposal
            priority = self.ADD_PRIORITY if guided else list(eligible)
            if not guided:
                rng.shuffle(priority)
            add = self._choose_three(priority, eligible, rng, guided)
            if add is None:
                return None
            proposal.update(add)

        data = np.zeros(self.N, dtype=np.uint8)
        if proposal:
            data[list(proposal)] = 1
        return data

    def _randomize_from_seed(self):
        rng = np.random.default_rng()
        n_moves = int(rng.integers(0, self.INITIAL_MUTATIONS + 1))
        for _ in range(n_moves):
            proposal = self._propose_move(rng, guided=False)
            if proposal is not None:
                self.data = proposal

    def _local_energy(self):
        # Novelty/(pool+1) is <1, so an entire unit of defect always dominates.
        return self.defect() - self.novelty() / (self.POOL_SIZE + 1)

    def _anneal(self):
        if self.LOCAL_SEARCH_STEPS <= 0:
            return

        rng = np.random.default_rng()
        current_energy = self._local_energy()
        best_data = self.data.copy()
        best_energy = current_energy

        for step in range(self.LOCAL_SEARCH_STEPS):
            old_data = self.data.copy()
            proposal = self._propose_move(rng, guided=True)
            if proposal is None:
                break

            self.data = proposal
            new_energy = self._local_energy()

            fraction = 1.0 - step / max(1, self.LOCAL_SEARCH_STEPS)
            temperature = max(0.05, self.ANNEAL_TEMPERATURE * fraction)

            accept = new_energy <= current_energy
            if not accept:
                log_accept = (current_energy - new_energy) / temperature
                accept = math.log(max(rng.random(), 1e-300)) < log_accept

            if accept:
                current_energy = new_energy
            else:
                self.data = old_data

            if current_energy < best_energy:
                best_energy = current_energy
                best_data = self.data.copy()

        self.data = best_data

    def local_search(self, improve_with_local_search):
        self._project_size()
        if improve_with_local_search:
            self._anneal()
        self.calc_features()
        self.calc_score()

    @classmethod
    def _update_class_params(cls, pars):
        cls.configure(*pars)

    @classmethod
    def _save_class_params(cls):
        return (
            cls.REFLECTION_M,
            cls.SEED_FAMILY,
            cls.SIZE_RADIUS,
            cls.MIN_LINES,
            cls.LOCAL_SEARCH_STEPS,
            cls.INITIAL_MUTATIONS,
            cls.ANNEAL_TEMPERATURE,
            cls.CANDIDATE_TOP_K,
            cls.GREEDY_MOVE_PROBABILITY,
        )


class HirzebruchMutationEnvironment(BaseEnvironment):
    k = 1
    are_coordinates_symmetric = False
    data_class = HirzebruchMutationDataPoint

    def __init__(self, params):
        super().__init__(params)

        self.data_class.configure(
            params.reflection_m,
            params.seed_family,
            params.size_radius,
            params.min_lines,
            params.local_search_steps,
            params.initial_mutations,
            params.anneal_temperature,
            params.candidate_top_k,
            params.greedy_move_probability,
        )

        # In this environment N is the candidate-pool size. The actual number
        # of selected arrangement lines is variable.
        params.N = self.data_class.POOL_SIZE

        self.tokenizer = SparseTokenizerSingleInteger(
            self.data_class,
            params.N,
            self.k,
            self.are_coordinates_symmetric,
            self.SPECIAL_SYMBOLS,
        )

        seed = reflection_arrangement(params.reflection_m, params.seed_family)
        logger.info(
            "Mutation environment: %s(m=%d), seed lines=%d, candidate pool=%d, "
            "allowed line counts=%d..%d",
            params.seed_family,
            params.reflection_m,
            len(seed),
            self.data_class.POOL_SIZE,
            self.data_class.SIZE_MIN,
            self.data_class.SIZE_MAX,
        )
        logger.info(
            "Seed Hirzebruch defect=%d, multiplicities=%s",
            hirzebruch_defect(seed, min_lines=params.min_lines),
            dict(sorted(multiplicity_profile(seed).items())),
        )

    @staticmethod
    def register_args(parser):
        parser.add_argument("--reflection_m", type=int, default=4, help="m in G(m,m,3) or G(m,1,3)")
        parser.add_argument("--seed_family", type=str, default="gmm", choices=["gmm", "gm1"])
        parser.add_argument(
            "--size_radius",
            type=int,
            default=3,
            help="Allowed line-count change from the seed; mutation moves occur in triples",
        )
        parser.add_argument("--min_lines", type=int, default=6, help="Exclude tiny/trivial arrangements")
        parser.add_argument("--local_search_steps", type=int, default=500)
        parser.add_argument(
            "--initial_mutations",
            type=int,
            default=2,
            help="Random triple moves used to diversify the initial population",
        )
        parser.add_argument(
            "--anneal_temperature",
            type=float,
            default=20.0,
            help="Initial simulated-annealing temperature in defect units",
        )
        parser.add_argument(
            "--candidate_top_k",
            type=int,
            default=12,
            help="Guided moves sample among this many top-ranked candidate lines",
        )
        parser.add_argument(
            "--greedy_move_probability",
            type=float,
            default=0.35,
            help="Chance of taking the three highest-ranked lines in a guided move",
        )
