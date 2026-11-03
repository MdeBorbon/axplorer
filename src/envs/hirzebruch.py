import math

import numpy as np

from src.envs.environment import BaseEnvironment, DataPoint
from src.envs.tokenizers import DenseTokenizer, SparseTokenizerSequenceKTokens, SparseTokenizerSingleInteger


def count_c4(data):
    """Count 4-cycles in the bipartite graph with incidence matrix ``data``."""
    a = data.astype(np.int16, copy=False)
    common = a @ a.T
    upper = common[np.triu_indices(data.shape[0], k=1)].astype(np.int64)
    return int(np.sum(upper * (upper - 1) // 2))


class HirzebruchDataPoint(DataPoint):
    """
    Candidate incidence structure for a symmetric (N_d) configuration.

    Rows represent lines, columns represent multiple points, and data[i, j] = 1
    means that line i contains point j. The target for the Hirzebruch 39-line
    problem is N=39 and d=6.
    """

    DEGREE = 6
    LOCAL_SEARCH_STEPS = 2000

    def __init__(self, N, init=False):
        super().__init__()
        self.N = N
        self.data = np.zeros((self.N, self.N), dtype=np.uint8)

        if self.DEGREE > self.N:
            raise ValueError("degree cannot exceed N")

        if init:
            self._project_to_regular()
            self._improve_c4()
            self.calc_features()
            self.calc_score()

    def _score_constants(self):
        # Absolute upper bounds valid for every N x N 0/1 matrix.
        c4_bound = math.comb(self.N, 2) ** 2
        max_degree_defect = 2 * self.N * max(self.DEGREE, self.N - self.DEGREE)

        # One unit of degree defect is more expensive than every possible C4.
        degree_weight = c4_bound + 1
        max_penalty = degree_weight * max_degree_defect + c4_bound
        return c4_bound, degree_weight, max_penalty

    def degree_defect(self):
        row_degrees = self.data.sum(axis=1).astype(np.int64)
        col_degrees = self.data.sum(axis=0).astype(np.int64)
        return int(np.abs(row_degrees - self.DEGREE).sum() + np.abs(col_degrees - self.DEGREE).sum())

    def calc_score(self):
        c4 = count_c4(self.data)
        defect = self.degree_defect()
        _, degree_weight, max_penalty = self._score_constants()
        penalty = degree_weight * defect + c4
        self.score = max_penalty - penalty

    def calc_features(self):
        self.features = ",".join(map(str, self.data.ravel().tolist()))

    def _project_to_regular(self):
        """
        Replace the current matrix by a simple d-regular bipartite graph.

        We construct d edge-disjoint perfect matchings. When the current matrix
        came from the transformer, its 1-entries are tried first, so the repair
        retains as much of the model proposal as the greedy matching permits.
        """
        preferred = self.data.copy()
        result = np.zeros_like(self.data)
        rng = np.random.default_rng()

        for _ in range(self.DEGREE):
            match_col = np.full(self.N, -1, dtype=np.int32)
            row_order = rng.permutation(self.N)
            choices = []

            for r in range(self.N):
                available = np.flatnonzero(result[r] == 0)
                preferred_cols = available[preferred[r, available] == 1].copy()
                other_cols = available[preferred[r, available] == 0].copy()
                rng.shuffle(preferred_cols)
                rng.shuffle(other_cols)
                choices.append(np.concatenate((preferred_cols, other_cols)))

            def augment(row, seen):
                for col_value in choices[row]:
                    col = int(col_value)
                    if seen[col]:
                        continue
                    seen[col] = True
                    matched_row = int(match_col[col])
                    if matched_row == -1 or augment(matched_row, seen):
                        match_col[col] = row
                        return True
                return False

            for row_value in row_order:
                row = int(row_value)
                seen = np.zeros(self.N, dtype=bool)
                if not augment(row, seen):
                    # After r<d matchings have been chosen, the unused graph is
                    # (N-r)-regular and therefore has a perfect matching.
                    raise RuntimeError("failed to construct a perfect matching")

            for col, row_value in enumerate(match_col):
                result[int(row_value), col] = 1

        self.data = result

    def _c4_delta_for_switch(self, r1, c1, r2, c2, common):
        """
        Change in the number of C4s under the bipartite 2-switch

            (r1,c1),(r2,c2) -> (r1,c2),(r2,c1).

        The overlap of rows r1 and r2 is unchanged; only their overlaps with
        the other rows contribute to the delta.
        """
        delta = 0

        for k in range(self.N):
            if k == r1 or k == r2:
                continue

            old = int(common[r1, k])
            change = int(self.data[k, c2]) - int(self.data[k, c1])
            new = old + change
            delta += new * (new - 1) // 2 - old * (old - 1) // 2

            old = int(common[r2, k])
            change = int(self.data[k, c1]) - int(self.data[k, c2])
            new = old + change
            delta += new * (new - 1) // 2 - old * (old - 1) // 2

        return delta

    def _improve_c4(self):
        """Greedy random 2-switch search, preserving all row/column degrees."""
        if self.LOCAL_SEARCH_STEPS <= 0:
            return

        rng = np.random.default_rng()
        a = self.data.astype(np.int16, copy=False)
        common = a @ a.T

        for _ in range(self.LOCAL_SEARCH_STEPS):
            r1, r2 = map(int, rng.choice(self.N, size=2, replace=False))
            c1 = int(rng.choice(np.flatnonzero(self.data[r1])))
            c2 = int(rng.choice(np.flatnonzero(self.data[r2])))

            if c1 == c2 or self.data[r1, c2] or self.data[r2, c1]:
                continue

            delta = self._c4_delta_for_switch(r1, c1, r2, c2, common)

            if delta < 0:
                self.data[r1, c1] = 0
                self.data[r2, c2] = 0
                self.data[r1, c2] = 1
                self.data[r2, c1] = 1
                a = self.data.astype(np.int16, copy=False)
                common = a @ a.T

    def local_search(self, improve_with_local_search):
        # Model samples need not satisfy the degree constraints. Project them
        # back to a d-regular incidence matrix before searching for fewer C4s.
        if self.degree_defect() != 0:
            self._project_to_regular()

        if improve_with_local_search:
            self._improve_c4()

        self.calc_features()
        self.calc_score()

    @classmethod
    def _update_class_params(cls, pars):
        cls.DEGREE, cls.LOCAL_SEARCH_STEPS = pars

    @classmethod
    def _save_class_params(cls):
        return cls.DEGREE, cls.LOCAL_SEARCH_STEPS


class HirzebruchEnvironment(BaseEnvironment):
    k = 2
    are_coordinates_symmetric = False
    data_class = HirzebruchDataPoint

    def __init__(self, params):
        super().__init__(params)
        self.data_class.DEGREE = params.degree
        self.data_class.LOCAL_SEARCH_STEPS = params.local_search_steps

        if params.encoding_tokens == "single_integer":
            self.tokenizer = SparseTokenizerSingleInteger(
                self.data_class, params.N, self.k, self.are_coordinates_symmetric, self.SPECIAL_SYMBOLS
            )
        elif params.encoding_tokens == "sequence_k_tokens":
            self.tokenizer = SparseTokenizerSequenceKTokens(
                self.data_class, params.N, self.k, self.are_coordinates_symmetric, self.SPECIAL_SYMBOLS
            )
        elif params.encoding_tokens == "adjacency":
            self.tokenizer = DenseTokenizer(
                self.data_class,
                params.N,
                self.k,
                self.are_coordinates_symmetric,
                self.SPECIAL_SYMBOLS,
                pow2base=params.pow2base,
            )
        else:
            raise ValueError(f"Invalid encoding: {params.encoding_tokens}")

    @staticmethod
    def register_args(parser):
        parser.add_argument("--N", type=int, default=39, help="Number of lines and number of multiple points")
        parser.add_argument("--degree", type=int, default=6, help="Required number of incidences at every line and point")
        parser.add_argument("--local_search_steps", type=int, default=2000, help="Number of random degree-preserving 2-switch attempts")
        parser.add_argument("--encoding_tokens", type=str, default="single_integer", help="single_integer/sequence_k_tokens/adjacency")
        parser.add_argument("--pow2base", type=int, default=1, help="Number of adjacency entries to code together")
