from logging import getLogger

import numpy as np

from src.envs.mutations_tabu import (
    HirzebruchTabuMutationDataPoint,
    HirzebruchTabuMutationEnvironment,
)

logger = getLogger()


class HirzebruchExploreDataPoint(HirzebruchTabuMutationDataPoint):
    """Exploratory mutation search with a size-neutral objective.

    Raw Hirzebruch defect D=sum(3v-N-3)^2 strongly prefers some line counts
    over others even when a larger arrangement has a greater fraction of lines
    already satisfying the Hirzebruch condition.  Here the primary objective is
    the fraction of defective lines.  Mean squared vertex error is secondary.

    Fixed-size exploration also uses one-line-for-one-line swaps.  Changes of
    line count remain in triples, as required by N=3k.
    """

    DEFECT_MAGNITUDE_WEIGHT = 0.02
    SINGLE_SWAP_PROBABILITY = 0.60

    def line_error_stats(self):
        selected = self.selected_indices()
        n = len(selected)
        if n == 0 or n % 3 != 0:
            return 1.0, float("inf"), np.array([], dtype=np.int64)

        counts = self._vertex_counts(selected)
        target = n // 3 + 1
        errors = counts - target
        bad_fraction = float(np.count_nonzero(errors)) / n
        mse = float(np.mean(errors.astype(float) ** 2))
        return bad_fraction, mse, errors

    def _exploration_energy_without_tabu(self):
        bad_fraction, mse, _ = self.line_error_stats()
        if not np.isfinite(mse):
            return 1_000_000.0

        # Novelty is only a tiny tie-breaker.  Any improvement in the
        # Hirzebruch statistics dominates it.
        novelty_bonus = self.novelty() / (10000.0 * (self.POOL_SIZE + 1))
        return bad_fraction + self.DEFECT_MAGNITUDE_WEIGHT * mse - novelty_bonus

    def _local_energy(self):
        energy = self._exploration_energy_without_tabu()
        if self.is_tabu():
            energy += self.TABU_PENALTY
        return energy

    def calc_score(self):
        # Larger Axplorer score is better.  Keep every ordinary candidate
        # positive while sending the exact reflection endpoints far down the
        # ranking through the inherited tabu energy.
        energy = self._local_energy()
        self.score = float(2_000_000_000.0 - 1_000_000.0 * energy)

    def _energy_for_data(self, data):
        old = self.data
        self.data = data
        try:
            return self._local_energy()
        finally:
            self.data = old

    def _single_swap(self, rng, guided):
        selected = self.selected_indices()
        outside = [i for i in range(self.N) if not self.data[i]]
        if not selected or not outside:
            return None

        remove_candidates = selected
        if guided:
            _, _, errors = self.line_error_stats()
            bad = [selected[i] for i, err in enumerate(errors) if err != 0]
            if bad and rng.random() < 0.75:
                remove_candidates = bad

        remove = int(rng.choice(remove_candidates))

        if not guided or rng.random() >= self.GREEDY_MOVE_PROBABILITY:
            add = int(rng.choice(outside))
            proposal = self.data.copy()
            proposal[remove] = 0
            proposal[add] = 1
            return proposal

        sample_size = min(self.CANDIDATE_TOP_K, len(outside))
        additions = list(map(int, rng.choice(outside, size=sample_size, replace=False)))
        best = None
        best_energy = None
        for add in additions:
            proposal = self.data.copy()
            proposal[remove] = 0
            proposal[add] = 1
            energy = self._energy_for_data(proposal)
            if best_energy is None or energy < best_energy:
                best_energy = energy
                best = proposal
        return best

    def _propose_move(self, rng, guided=True):
        # The fixed-N geometry is naturally connected by one-line swaps.  Use
        # them frequently.  For size-changing moves deliberately turn off the
        # old seed-specific priority: the validation experiment showed that it
        # simply hugs the known G(m,m,3) <-> G(m,1,3) transition.
        if rng.random() < self.SINGLE_SWAP_PROBABILITY:
            proposal = self._single_swap(rng, guided)
            if proposal is not None:
                return proposal

        return super()._propose_move(rng, guided=False)


class HirzebruchExploreEnvironment(HirzebruchTabuMutationEnvironment):
    data_class = HirzebruchExploreDataPoint

    def __init__(self, params):
        super().__init__(params)
        logger.info(
            "Exploratory objective: defective-line fraction + %.3f * mean-square vertex error",
            self.data_class.DEFECT_MAGNITUDE_WEIGHT,
        )
