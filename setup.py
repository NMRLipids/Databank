from setuptools import setup, find_packages
from os import path


def parse_requirements(filename):
    with open(filename, "r") as f:
        return f.read().splitlines()


setup(
    name="DatabankLib",
    version="1.1.0",                          # should agree with git tag!
    package_dir={"": "Scripts"},
    packages=find_packages(where="Scripts"),  # Automatically list packages
    install_requires=parse_requirements(
        path.join('.', 'Scripts', 'DatabankLib', 'requirements.txt')),
    author="NMRlipids open collaboration",  # Your name or organization
    author_email="samuli.ollila@helsinki.fi",
    description="NMRLipids Databank main package",
    url="https://github.com/NMRlipids/Databank",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9"
)
