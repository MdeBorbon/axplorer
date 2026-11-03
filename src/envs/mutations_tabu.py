from logging import getLogger

from src.envs.mutations import (
    HirzebruchMutationDataPoint,
    HirzebruchMutationEnvironment,
    _projective_key,
    reflection_arrangement,
)

logger = getLogger()


class HirzebruchTabuMutationDataPoint(HirzebruchMutationDataPoint):
    """Mutation search that makes the known reflection endpoints repulsive.

    The first mutation environment deliberately rewards the known transition
    G(m,m,3) <-> G(m,1,3) as a validation target. Once that has been recovered,
    this subclass places both exact coordinate realizations on a tabu list so
    that local search and Axplorer must move away from them.
    """

    TABU_PENALTY = 1000.0
    TABU_SETS = ()
    TABU_LABELS = {}

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
        super().configure(
            reflection_m,
            seed_family,
            size_radius,
            min_lines,
            local_search_steps,
            initial_mutations,
            anneal_temperature,
            candidate_top_k,
            greedy_move_probability,
        )
        cls._build_tabu_sets()

    @classmethod
    def _indices_for_lines(cls, lines):
        pool_index = {_projective_key(line): i for i, line in enumerate(cls.CANDIDATE_LINES)}
        indices = []
        for line in lines:
            key = _projective_key(line)
            if key not in pool_index:
                return None
            indices.append(pool_index[key])
        return frozenset(indices)

    @classmethod
    def _build_tabu_sets(cls):
        tabu = []
        labels = {}
        for family in ("gmm", "gm1"):
            indices = cls._indices_for_lines(reflection_arrangement(cls.REFLECTION_M, family))
            if indices is not None:
                tabu.append(indices)
                labels[indices] = f"{family}(m={cls.REFLECTION_M})"
        cls.TABU_SETS = tuple(tabu)
        cls.TABU_LABELS = labels

    def selected_set(self):
        return frozenset(self.selected_indices())

    def tabu_label(self):
        return self.TABU_LABELS.get(self.selected_set())

    def is_tabu(self):
        return self.selected_set() in self.TABU_SETS

    def _local_energy(self):
        energy = super()._local_energy()
        if self.is_tabu():
            energy += self.TABU_PENALTY
        return energy

    def calc_score(self):
        super().calc_score()
        if self.is_tabu():
            # Use the same scale as a defect penalty, so known exact solutions
            # rank far below genuinely exploratory near-zero configurations.
            defect_weight = self.POOL_SIZE + 1
            self.score -= int(defect_weight * self.TABU_PENALTY)


class HirzebruchTabuMutationEnvironment(HirzebruchMutationEnvironment):
    data_class = HirzebruchTabuMutationDataPoint

    def __init__(self, params):
        super().__init__(params)
        logger.info(
            "Tabu reflection configurations: %s",
            sorted(self.data_class.TABU_LABELS.values()),
        )
