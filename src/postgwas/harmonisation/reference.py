import os

def load_reference_build(build: str, base_path: str = "/refs") -> str:
    path = os.path.join(base_path, build.lower())
    if not os.path.exists(path):
        raise FileNotFoundError(f"Reference path not found for build {build}: {path}")
    print(f"Loaded reference for {build} from {path}")
    return path
