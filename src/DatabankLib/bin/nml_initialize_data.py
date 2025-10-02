#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
import sys
import zipfile
from io import BytesIO
from urllib.request import urlopen, urlretrieve

DBREPO = "NMRlipids/BilayerData"
CODEREPO = "NMRlipids/Databank"


def git_exists() -> bool:
    """Check if git is installed."""
    return shutil.which("git") is not None


def get_latest_release_tag() -> str:
    """Fetch the latest release tag from the GitHub API."""
    url = f"https://api.github.com/repos/{DBREPO}/releases/latest"
    with urlopen(url) as response:
        data = json.load(response)
        return data["tag_name"]


def clone_repo(tag: str) -> None:
    """Clone the repository at the specified tag."""
    subprocess.run(["git", "clone", "--branch", tag, f"https://github.com/{DBREPO}.git"], check=True)


def download_and_extract_zip(tag: str) -> None:
    """Download and extract the ZIP file of the specified tag (release)."""
    url = f"https://github.com/{DBREPO}/archive/refs/tags/{tag}.zip"
    with urlopen(url) as response:
        with zipfile.ZipFile(BytesIO(response.read())) as z:
            z.extractall()


def download_toy_folder(dest: str = "ToyData", folder_path: str = "Scripts/tests/Data") -> None:
    """Download a folder from a GitHub repository using the GitHub API."""
    branch = "main"
    api_url = f"https://api.github.com/repos/{CODEREPO}/contents/{folder_path}?ref={branch}"
    os.makedirs(dest, exist_ok=True)

    with urlopen(api_url) as response:
        if response.status != 200:
            raise Exception(f"Feil ved hente data: {response.status} {response.text}")
        data_bytes = response.read()
        data_text = data_bytes.decode("utf-8")
        items = json.loads(data_text)

    if isinstance(items, dict) and items.get("message"):
        raise Exception(f"GitHub API feil: {items['message']}")

    for item in items:
        if item["type"] == "file":
            download_url = item["download_url"]
            filename = item["name"]
            dest_path = os.path.join(dest, filename)
            urlretrieve(download_url, dest_path)
            print(".", end="", flush=True)
        else:
            download_toy_folder(os.path.join(dest, item["name"]), os.path.join(folder_path, item["name"]))


def main() -> None:
    """Code program logic."""
    print("== NMRlipids Databank Data Initializer ==")
    # can be prg toy | prg stable | prg dev
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: nml_initialize_data.py [toy|stable|dev]\n")
        sys.exit(1)

    prg = sys.argv[1].lower()
    if prg not in ["toy", "stable", "dev"]:
        sys.stderr.write("Invalid argument. Use 'toy', 'stable', or 'dev'.")
        sys.exit(1)

    if prg == "toy":
        print("Downloading toy data: [", end="", flush=True)
        download_toy_folder()
        print("] Done.")
        data_path = os.path.join(os.getcwd(), "ToyData")
        sim_path = os.path.join(data_path, "Simulations.1")
    elif prg == "stable":
        if git_exists():
            tag = get_latest_release_tag()
            print(f"Git found. Cloning {DBREPO} at {tag} ...")
            clone_repo(tag)
        else:
            tag = get_latest_release_tag()
            print(f"Git not found. Downloading ZIP for {DBREPO} at {tag} ...")
            download_and_extract_zip(tag)
            print("Done.")
        data_path = os.path.join(os.getcwd(), "BilayerData")
    elif prg == "dev":
        if not git_exists():
            sys.stderr.write("Git is required for 'dev' option. Please install git.\n")
            sys.exit(1)
        print(f"Cloning {DBREPO} at main branch ...")
        clone_repo("main")
        data_path = os.path.join(os.getcwd(), "BilayerData")
    
    # Write environment setup file
    with open("databank_env.rc", "w") as f:
        f.write(f"export NMLDB_DATA_PATH={data_path}\n")
        if sim_path is not None:
            f.write(f"export NMLDB_SIMU_PATH={sim_path}\n")

    print(f""""Data initialized into {data_path}. Please do 

   $ source databank_env.rc

to set environment variables.""")


if __name__ == "__main__":
    main()
