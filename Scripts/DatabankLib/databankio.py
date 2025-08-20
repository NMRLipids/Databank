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
import functools

logger = logging.getLogger(__name__)
MAX_DRYRUN_SIZE = 50 * 1024 * 1024  # 50 MB, max size for dry-run download

def retry_with_exponential_backoff(max_attempts=3, delay_seconds=1):
    """
    Decorator that retries a function with exponential backoff for transient network errors.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            current_delay = delay_seconds
            while attempts < max_attempts:
                try:
                    return func(*args, **kwargs)
                
                # --- New logic to handle non-retriable HTTP client errors ---
                except urllib.error.HTTPError as e:
                    # Check if the error is a client error (4xx) which is not likely to be resolved by a retry.
                    if 400 <= e.code < 500:
                        logger.error(
                            f"Function {func.__name__} failed with non-retriable client error {e.code}: {e.reason}"
                        )
                        raise e  # Re-raise the HTTPError immediately
                    
                    # For other errors (like 5xx server errors), proceed with retry logic.
                    logger.warning(f"Caught retriable HTTPError {e.code}. Proceeding with retry...")
                    # Fall through to the generic exception handling below.
                
                except (urllib.error.URLError, socket.timeout) as e:
                    # This block now primarily handles non-HTTP errors or retriable HTTP errors.
                    pass # Fall through to the retry logic below

                # --- Existing retry logic ---
                attempts += 1
                if attempts < max_attempts:
                    logger.warning(
                        f"Attempt {attempts}/{max_attempts} for {func.__name__} failed. "
                        f"Retrying in {current_delay:.1f} seconds..."
                    )
                    time.sleep(current_delay)
                    current_delay *= 2
                else:
                    logger.error(
                        f"Function {func.__name__} failed after {max_attempts} attempts."
                    )
                    # Re-raise the last exception caught to be handled by the caller
                    raise ConnectionError(
                        f"Function {func.__name__} failed after {max_attempts} attempts."
                    )
        return wrapper
    return decorator

# --- Decorated Helper Functions for Network Requests ---

@retry_with_exponential_backoff(max_attempts=5, delay_seconds=2)
def _open_url_with_retry(uri: str, timeout: int = 10):
    """
    A private helper to open a URL with a timeout and retry logic.
    Returns the response object.
    """
    return urllib.request.urlopen(uri, timeout=timeout)

@retry_with_exponential_backoff(max_attempts=5, delay_seconds=2)
def get_file_size_with_retry(uri: str) -> int:
    """Fetches the size of a file from a URI with retry logic."""
    with _open_url_with_retry(uri) as response:
        content_length = response.getheader('Content-Length')
        return int(content_length) if content_length else 0

@retry_with_exponential_backoff(max_attempts=5, delay_seconds=2)
def download_with_progress_with_retry(uri: str, dest: str, fi_name: str) -> None:
    """Downloads a file with a progress bar and retry logic."""
    class RetrieveProgressBar(tqdm):
        def update_retrieve(self, b=1, bsize=1, tsize=None):
            if tsize is not None:
                self.total = tsize
            return self.update(b * bsize - self.n)

    with RetrieveProgressBar(
        unit="B", unit_scale=True, unit_divisor=1024, miniters=1, desc=fi_name
    ) as u:
        urllib.request.urlretrieve(uri, dest, reporthook=u.update_retrieve)

# --- Main Functions ---

def download_resource_from_uri(
    uri: str, dest: str, override_if_exists: bool = False,
    dry_run_mode: bool = False
) -> int:
    """
    Download file resource from uri to given file destination using urllib.
    
    :param uri: (str) file URL
    :param dest: (str) file destination path
    :param override_if_exists: (bool, optional) Override dest. file if exists. Defaults to False.
    
    :raises ConnectionError: An error occurred after multiple download attempts.
    
    :return: code (int) 0 - OK, 1 - skipped, 2 - redownloaded
    """
    fi_name = uri.split("/")[-1]

    # Check if dest path already exists and compare file size
    if not override_if_exists and os.path.isfile(dest):
        try:
            fi_size = get_file_size_with_retry(uri)
            if fi_size == os.path.getsize(dest):
                logger.info(f"{dest}: file already exists, skipping")
                return 1
            else:
                logger.warning(
                    f"{fi_name} filesize mismatch of local file '{fi_name}', redownloading ..."
                )
                return 2
        except (ConnectionError, FileNotFoundError) as e:
            logger.error(
                f"Failed to verify file size for {fi_name}: {e}. Proceeding with redownload."
            )
            return 2
    
    # Download file in dry run mode
    if dry_run_mode:
        url_size = get_file_size_with_retry(uri)
        # Use the decorated helper for opening the URL
        with _open_url_with_retry(uri) as response, open(dest, "wb") as out_file:
            total = min(url_size, MAX_DRYRUN_SIZE) if url_size else MAX_DRYRUN_SIZE
            downloaded = 0
            chunk_size = 8192
            next_report = 10 * 1024 * 1024
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

    # Download with progress bar and check for final size match
    url_size = get_file_size_with_retry(uri)
    download_with_progress_with_retry(uri, dest, fi_name)

    size = os.path.getsize(dest)
    if url_size != 0 and url_size != size:
        raise OSError(f"Downloaded filesize mismatch ({size}/{url_size} B)")
    
    return 0

def resolve_doi_url(doi: str, validate_uri: bool = True) -> str:
    """
    Returns the full DOI URL and validates that it is a reachable address.

    :param doi: (str) The DOI identifier (e.g., "10.5281/zenodo.1234")
    :param validate_uri: (bool, optional) If True, checks if the URL is valid. Defaults to True.

    :raises urllib.error.HTTPError: If the DOI resolves to a URL, but the server returns an
                                    HTTP error code (e.g., 404 Not Found).
    :raises ConnectionError: If the server cannot be reached after multiple retries.
    
    :return: (str) The full, validated DOI link.
    """
    res = "https://doi.org/" + doi

    if validate_uri:
        try:
            # The 'with' statement ensures the connection is closed
            with _open_url_with_retry(res):
                pass
        except urllib.error.HTTPError as e:
            # This specifically catches HTTP errors like 404, 500, etc.
            logger.error(f"Validation failed for DOI {doi}. URL <{res}> returned HTTP {e.code}: {e.reason}")
            # Re-raise the specific HTTPError so the caller can handle it
            raise e
        except ConnectionError as e:
            # This catches the final error from the retry decorator for other issues
            # (e.g., DNS failure, connection refused).
            logger.error(f"Could not connect to <{res}> after multiple attempts.")
            raise e
            
    return res

def resolve_download_file_url(
        doi: str, fi_name: str, validate_uri: bool = True,
        sleep429=5) -> str:
    """
    :meta private:
    Resolve file URI from supported DOI with given filename
    """
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        uri = "https://zenodo.org/record/" + zenodo_entry_number + "/files/" + fi_name

        if validate_uri:
            try:
                # Use the decorated helper to check if the URI exists
                with _open_url_with_retry(uri):
                    pass
            except urllib.error.HTTPError as hte:
                # Special handling for 429 after retries have failed
                if hte.code == 429:
                    if sleep429 / 5 > 10:
                        raise TimeoutError(
                            "Too many iterations of increasing waiting time for 429 error!"
                        )
                    logger.warning(
                        f"Site returned 429 (Too Many Requests) for URI: {uri}. "
                        f"Sleeping {sleep429} seconds before retrying..."
                    )
                    time.sleep(sleep429)
                    # Recursive call with increased sleep time
                    return resolve_download_file_url(
                        doi, fi_name, validate_uri, sleep429=sleep429 + 5
                    )
                else:
                    # Other persistent HTTP errors are re-raised
                    raise
            except ConnectionError as ce:
                # Catch the final failure from our decorator for other network issues
                raise RuntimeError(f"Cannot open {uri}. Failed after multiple retries.") from ce
        return uri
    else:
        raise NotImplementedError(
            "Repository not validated. Please upload the data for example to zenodo.org"
        )

def calc_file_sha1_hash(fi: str, step: int = 67108864, one_block: bool = True) -> str:
    """
    Calculates SHA1 hash of given file using hashlib.
    """
    sha1_hash = hashlib.sha1()
    n_tot_steps = math.ceil(os.path.getsize(fi) / step)
    with open(fi, "rb") as f:
        if one_block:
            block = f.read(step)
            sha1_hash.update(block)
        else:
            with tqdm(total=n_tot_steps, desc="Calculating SHA1") as pbar:
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

    if os.path.exists(directory_path) and os.listdir(directory_path):
        raise FileExistsError(
            f"Output directory '{directory_path}' is not empty. Delete it if you wish."
        )

    if not dry_run_mode:
        os.makedirs(directory_path, exist_ok=True)

    return directory_path