from pathlib import Path

from dotenv import load_dotenv
from loguru import logger
from importlib_resources import files
import yaml

executables_file = [
    i for i in files("a3mtools").iterdir() if i.name == "executables.yaml"
][0]
if executables_file is None:
    raise FileNotFoundError("executables.yaml file not found")
else:
    with executables_file.open() as f:
        EXECUTABLES = yaml.safe_load(f)
    print(EXECUTABLES)


colabfold_batch = EXECUTABLES["colabfold_batch"]
colabfold_data = EXECUTABLES["colabfold_data"]
# Load environment variables from .env file if it exists
load_dotenv()

# Paths
PROJ_ROOT = Path(__file__).resolve().parents[1]
logger.info(f"PROJ_ROOT path is: {PROJ_ROOT}")

DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"

MODELS_DIR = PROJ_ROOT / "models"

REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"

# If tqdm is installed, configure loguru with tqdm.write
# https://github.com/Delgan/loguru/issues/135
# try:
#     from tqdm import tqdm

#     logger.remove(0)
#     logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True)
# except ModuleNotFoundError:
#     pass
