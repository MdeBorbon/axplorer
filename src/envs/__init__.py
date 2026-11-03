from src.envs.cycle import SquareEnvironment
from src.envs.hirzebruch import HirzebruchEnvironment
from src.envs.isosceles import IsoscelesEnvironment
from src.envs.mutations import HirzebruchMutationEnvironment
from src.envs.mutations_explore import HirzebruchExploreEnvironment
from src.envs.mutations_tabu import HirzebruchTabuMutationEnvironment
from src.envs.sphere import SphereEnvironment

ENVS = {
    "square": SquareEnvironment,
    "hirzebruch39": HirzebruchEnvironment,
    "hirzebruch_mutations": HirzebruchMutationEnvironment,
    "hirzebruch_mutations_tabu": HirzebruchTabuMutationEnvironment,
    "hirzebruch_mutations_explore": HirzebruchExploreEnvironment,
    "isosceles": IsoscelesEnvironment,
    "sphere": SphereEnvironment,
}


def build_env(params):
    """
    Build environment.
    """
    env = ENVS[params.env_name](params)
    return env
