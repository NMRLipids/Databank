"""
:module databankio: Inut/Output auxilary module
:description:
    Input/Output auxilary module with some small usefull functions. It includes:
    - Network communication.
    - Downloading files.
    - Checking links.
    - Resolving DOIs.
    - Calculating file hashes.
"""

import logging
import math
import hashlib
import os
import time
import socket
from typing import Mapping
import urllib.error
from tqdm import tqdm
import urllib.request

logger = logging.getLogger(__name__)
MAX_DRYRUN_SIZE = 50 * 1024 * 1024  # 50 MB, max size for dry-run download


def download_resource_from_uri(
    uri: str, dest: str, override_if_exists: bool = False,
    dry_run_mode: bool = False
) -> int:
    """
    Download file resource [from uri] to given file destination using urllib

    :param uri: (str) file URL
    :param dest: (str) file destination path
    :param override_if_exists: (bool, optional)
            Override dest. file if exists. Defaults to False.

    :raises HTTPException: An error occured during download

    :return: code (int) 
             0 - OK, 1 - skipped, 2 - redownloaded
    """
    # TODO verify file size before skipping already existing download!

    class RetrieveProgressBar(tqdm):
        # uses tqdm.update(), see docs https://github.com/tqdm/tqdm#hooks-and-callbacks
        def update_retrieve(self, b=1, bsize=1, tsize=None):
            if tsize is not None:
                self.total = tsize
            return self.update(b * bsize - self.n)

    fi_name = uri.split("/")[-1]

    # check if dest path already exists
    if not override_if_exists and os.path.isfile(dest):
        socket.setdefaulttimeout(10)  # seconds

        # compare filesize
        fi_size = urllib.request.urlopen(uri).length  # download size
        if fi_size == os.path.getsize(dest):
            logger.info(f"{dest}: file already exists, skipping")
            return 1
        else:
            logger.warning(
                f"{fi_name} filesize mismatch of local "
                f"file '{fi_name}', redownloading ..."
            )
            return 2

    # download
    socket.setdefaulttimeout(10)  # seconds
    url_size = urllib.request.urlopen(uri).length  # download size

    if dry_run_mode:
        # Download only up to MAX_DRYRUN_SIZE bytes, no tqdm
        with urllib.request.urlopen(uri) as response, open(dest, "wb") as out_file:
            total = min(url_size, MAX_DRYRUN_SIZE) if url_size else MAX_DRYRUN_SIZE
            downloaded = 0
            chunk_size = 8192
            next_report = 10 * 1024 * 1024  # print every 10 MB
            logger.info(
                "Dry-run: Downloading up to"
                f" {total // (1024*1024)} MB of {fi_name} ...")
            while downloaded < total:
                to_read = min(chunk_size, total - downloaded)
                chunk = response.read(to_read)
                if not chunk:
                    break
                out_file.write(chunk)
                downloaded += len(chunk)
                if downloaded >= next_report:
                    print(f"  Downloaded {downloaded // (1024*1024)} MB ...")
                    next_report += 10 * 1024 * 1024
            logger.info(
                "Dry-run: Finished, downloaded"
                f" {downloaded // (1024*1024)} MB of {fi_name}")
        return 0

    with RetrieveProgressBar(
        unit="B", unit_scale=True, unit_divisor=1024, miniters=1, desc=fi_name
    ) as u:
        _ = urllib.request.urlretrieve(uri, dest, reporthook=u.update_retrieve)

    # check if the file is fully downloaded
    size = os.path.getsize(dest)

    if url_size != size:
        raise Exception(f"downloaded filsize mismatch ({size}/{url_size} B)")

    return 0


def resolve_doi_url(doi: str, validate_uri: bool = True) -> str:
    """
    :meta private:
    Returns full doi link of given ressource, also checks if URL is valid.

    Args:
        doi (str): [doi] part from config
        validate_uri (bool, optional): Check if URL is valid. Defaults to True.

    Returns:
        str: full doi link
    """
    res = "https://doi.org/" + doi

    if validate_uri:
        socket.setdefaulttimeout(10)  # seconds
        _ = urllib.request.urlopen(res)
    return res


def resolve_download_file_url(
        doi: str, fi_name: str, validate_uri: bool = True,
        sleep429=5) -> str:
    """
    :meta private:
    Resolve file URI from supported DOI with given filename

    Args:
        doi (str): DOI string
        fi_name (str): name of the file to resolve from source
        validate_uri (bool, optional): Check if URI exists. Defaults to True.
        sleep429 (int, optional): Sleep in seconds if 429 HTTP code returned

    Raises:
        NotImplementedError: Unsupported DOI repository
        HTTPError: HTTP Error Status Code
        URLError: Failed to reach the server

    Returns:
        str: file URI
    """
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        uri = "https://zenodo.org/record/" + zenodo_entry_number + "/files/" + fi_name

        # check if ressource exists, may throw exception
        if validate_uri:
            try:
                socket.setdefaulttimeout(10)  # seconds
                _ = urllib.request.urlopen(uri, timeout=10)
            except TimeoutError:
                raise RuntimeError(f"Cannot open {uri}. Timeout error.")
            except urllib.error.HTTPError as hte:
                if hte.code == 429:
                    if sleep429/5 > 10:
                        raise TimeoutError(
                            "Too many iteration of increasing waiting time!")
                    logger.warning(f"HTTP error returned from URI: {uri}")
                    logger.warning(f"Site returns 429 code."
                                   f" Try to sleep {sleep429} seconds and repeat!")
                    time.sleep(sleep429)
                    return resolve_download_file_url(doi, fi_name, validate_uri,
                                                     sleep429=sleep429+5)
                else:
                    raise hte
        return uri
    else:
        raise NotImplementedError(
            "Repository not validated. Please upload the data for example to zenodo.org"
        )


def calc_file_sha1_hash(fi: str, step: int = 67108864, one_block: bool = True) -> str:
    """
    Calculates SHA1 hash of given file using hashlib.

    :param fi: (str) path to file
    :param step: (int, optional) file read bytes step. Defaults to 64MB.
    :param one_block: (bool, optional) read just a single block. Defaults to True.

    :returns str: sha1 filehash of 40 char length
    """
    sha1_hash = hashlib.sha1()
    n_tot_steps = math.ceil(os.path.getsize(fi) / step)
    with open(fi, "rb") as f:
        if one_block:
            block = f.read(step)
            sha1_hash.update(block)
        else:
            # we don't need tqdm from one-block SHA1
            with tqdm(total=n_tot_steps) as pbar:
                # Read and update hash string value in blocks of 4K
                for byte_block in iter(lambda: f.read(step), b""):
                    sha1_hash.update(byte_block)
                    pbar.update(1)
    return sha1_hash.hexdigest()


def create_databank_directories(
        sim: Mapping, 
        sim_hashes: Mapping, 
        out: str,
        dry_run_mode: bool = False
        ) -> str:
    """
    Creates nested output directory structure to save results.

    :param sim: Processed simulation entries.
    :param sim_hashes: File hashes needed for directory structure.
    :param out: Output base path (str).
    :param dry_run_mode: If True, do not create directories, just return the path.

    :returns: Output directory (str).

    :raises NotImplementedError: If the simulation software is unsupported.
    :raises OSError: If an error occurs while creating the output directory.
    :raises FileExistsError: If the output directory already exists and is not empty.
    """
    # resolve output dir naming
    if sim["SOFTWARE"] == "gromacs":
        head_dir = sim_hashes.get("TPR")[0][1][0:3]
        sub_dir1 = sim_hashes.get("TPR")[0][1][3:6]
        sub_dir2 = sim_hashes.get("TPR")[0][1]
        sub_dir3 = sim_hashes.get("TRJ")[0][1]
    elif sim["SOFTWARE"] == "openMM" or sim["SOFTWARE"] == "NAMD":
        head_dir = sim_hashes.get("TRJ")[0][1][0:3]
        sub_dir1 = sim_hashes.get("TRJ")[0][1][3:6]
        sub_dir2 = sim_hashes.get("TRJ")[0][1]
        sub_dir3 = sim_hashes.get("TRJ")[0][1]
    else:
        raise NotImplementedError(f"sim software '{sim['SOFTWARE']}' not supported")

    directory_path = os.path.join(out, head_dir, sub_dir1, sub_dir2, sub_dir3)

    logger.debug(f"output_dir = {directory_path}")

    # destination directory is not empty
    if os.path.exists(directory_path) and os.listdir(directory_path) != 0:
        raise FileExistsError(
            f"Output directory '{directory_path}' is not empty. Delete it if you wish."
        )

    # create directories
    if not dry_run_mode:
        os.makedirs(directory_path, exist_ok=True)

    return directory_path
