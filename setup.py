from setuptools import setup
from os import path
import re


def get_variable_from_file(vname, filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    match = re.search(rf"^__{vname}__\s*=\s*['\"]([^'\"]+)['\"]", content, re.MULTILINE)
    if match:
        return match.group(1)
    raise RuntimeError(f"{vname} string not found")


def extract_multis_var(vname, filepath):
    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except OSError as e:
        raise RuntimeError("Could not read source file for extracting versioning"
                           f" information {filepath}: {e}") from e
    pattern = re.compile(
        rf"__{vname}__\s*=\s*\((\s+[^)]+)\)", re.MULTILINE
    )
    match = pattern.search(content)
    if not match:
        raise RuntimeError(f"{vname} string not found")
    strings = re.findall(r'"([^"]*)"', match.group(1))
    return ''.join(strings)


def parse_requirements(filename):
    try:
        with open(filename, "r") as f:
            return f.read().splitlines()
    except OSError as e:
        raise RuntimeError(f"Could not read requirements file {filename}: {e}") from e


init_fpath = path.join('.', 'Scripts', 'DatabankLib', '__init__.py')

setup(
    name="DatabankLib",
    version=get_variable_from_file('version', init_fpath),
    author=get_variable_from_file('author', init_fpath),
    author_email=get_variable_from_file('author_email', init_fpath),
    description=extract_multis_var('description', init_fpath),
    url=get_variable_from_file('url', init_fpath),
    package_dir={"": "Scripts"},
    package_data={"": ["*.yaml"]},
    packages=["DatabankLib"],
    install_requires=parse_requirements(
        path.join('.', 'Scripts', 'DatabankLib', 'requirements.txt')),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9"
)
