from setuptools import setup, find_packages
from pathlib import Path


setup(
    name="LiftingForce_cfg",
    version="0.1",
    packages=["config"], # find_packages(),  # Автоматически найдет, если есть __init__.py
    install_requires=Path("requirements.txt").read_text().splitlines()
)
