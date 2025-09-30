#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
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


def download_github_folder(dest: str = "ToyData", folder_path: str = "Scripts/tests/Data") -> None:
    """Download a folder from a GitHub repository using the GitHub API."""
    branch = "main"
    api_url = f"https://api.github.com/repos/{CODEREPO}/contents/{folder_path}?ref={branch}"
    print(api_url)
    os.makedirs(dest, exist_ok=True)

    with urlopen(api_url) as response:
        if response.status != 200:
            raise Exception(f"Feil ved hente data: {response.status} {response.text}")
        data_bytes = response.read()
        data_text = data_bytes.decode('utf-8')
        items = json.loads(data_text)

    if isinstance(items, dict) and items.get("message"):
        raise Exception(f"GitHub API feil: {items['message']}")

    for item in items:
        if item["type"] == "file":
            download_url = item["download_url"]
            filename = item["name"]
            dest_path = os.path.join(dest, filename)
            urlretrieve(download_url, dest_path)
        else:
            download_github_folder(os.path.join(dest, item["name"]), os.path.join(folder_path, item["name"]))


def main() -> None:
    """Code program logic."""

    download_github_folder()
    import sys
    sys.exit(1)

    if git_exists():
        tag = get_latest_release_tag()
        print(f"Git found. Cloning {DBREPO} at {tag} ...")
        clone_repo(tag)
    else:
        tag = get_latest_release_tag()
        print(f"Git not found. Downloading ZIP for {DBREPO} at {tag} ...")
        download_and_extract_zip(tag)
        print("Done.")


if __name__ == "__main__":
    main()
