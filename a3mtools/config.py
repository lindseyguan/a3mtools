from pathlib import Path

from dotenv import load_dotenv
from importlib_resources import files
import yaml

example_a3m_files = [
    i for i in files("a3mtools.test_data").iterdir() if i.stem != "__init__"
]


