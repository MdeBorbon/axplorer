from src.envs.cycle import SquareEnvironment
from src.envs.hirzebruch import HirzebruchEnvironment
from src.envs.isosceles import IsoscelesEnvironment
from src.envs.mutations import HirzebruchMutationEnvironment
from src.envs.sphere import SphereEnvironment

ENVS = {
    "square": SquareEnvironment,
    "hirzebruch39": HirzebruchEnvironment,
    "hirzebruch_mutations": HirzebruchMutationEnvironment,
    "isosceles": IsoscelesEnvironment,
    "sphere": SphereEnvironment,
}


def build_env(params):
    """
    Build environment.
    """
    env = ENVS[params.env_name](params)
    return env
