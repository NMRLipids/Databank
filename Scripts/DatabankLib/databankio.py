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
    """Decorator that retries a function with exponential backoff.

    Args:
        max_attempts (int): The maximum number of attempts. Defaults to 3.
        delay_seconds (int): The initial delay between retries in seconds.
            The delay doubles after each failed attempt. Defaults to 1.
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
    """A private helper to open a URL with a timeout and retry logic.

    Args:
        uri (str): The URL to open.
        timeout (int): The timeout for the request in seconds. Defaults to 10.

    Returns:
        The response object from urllib.request.urlopen.
    """
    return urllib.request.urlopen(uri, timeout=timeout)

@retry_with_exponential_backoff(max_attempts=5, delay_seconds=2)
def get_file_size_with_retry(uri: str) -> int:
    """Fetches the size of a file from a URI with retry logic.

    Args:
        uri (str): The URL of the file.

    Returns:
        int: The size of the file in bytes, or 0 if the 'Content-Length'
            header is not present.
    """
    with _open_url_with_retry(uri) as response:
        content_length = response.getheader('Content-Length')
        return int(content_length) if content_length else 0

@retry_with_exponential_backoff(max_attempts=5, delay_seconds=2)
def download_with_progress_with_retry(uri: str, dest: str, fi_name: str) -> None:
    """Downloads a file with a progress bar and retry logic.

    Uses tqdm to display a progress bar during the download.

    Args:
        uri (str): The URL of the file to download.
        dest (str): The local destination path to save the file.
        fi_name (str): The name of the file, used for the progress bar
            description.
    """
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
    """Download file resource from a URI to a local destination.

    Checks if the file already exists and has the same size before downloading.
    Can also perform a partial "dry-run" download.

    Args:
        uri (str): The URL of the file resource.
        dest (str): The local destination path to save the file.
        override_if_exists (bool): If True, the file will be re-downloaded
            even if it already exists. Defaults to False.
        dry_run_mode (bool): If True, only a partial download is performed
            (up to MAX_DRYRUN_SIZE). Defaults to False.

    Returns:
        int: A status code indicating the result.
            0: Download was successful.
            1: Download was skipped because the file already exists.
            2: File was re-downloaded due to a size mismatch.

    Raises:
        ConnectionError: An error occurred after multiple download attempts.
        OSError: The downloaded file size does not match the expected size.
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
                    logger.info(f"  Downloaded {downloaded // (1024*1024)} MB ...")
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
    """Resolves a DOI to a full URL and validates that it is reachable.

    Args:
        doi (str): The DOI identifier (e.g., "10.5281/zenodo.1234").
        validate_uri (bool): If True, checks if the resolved URL is a valid
            and reachable address. Defaults to True.

    Returns:
        str: The full, validated DOI link (e.g., "https://doi.org/...").

    Raises:
        urllib.error.HTTPError: If the DOI resolves to a URL, but the server
            returns an HTTP error code (e.g., 404 Not Found).
        ConnectionError: If the server cannot be reached after multiple retries.
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
        doi: str, fi_name: str, validate_uri: bool = True) -> str:
    """:meta private:
    Resolves a download file URI from a supported DOI and filename.

    Currently supports Zenodo DOIs.

    Args:
        doi (str): The DOI identifier for the repository (e.g.,
            "10.5281/zenodo.1234").
        fi_name (str): The name of the file within the repository.
        validate_uri (bool): If True, checks if the resolved URL is a valid
            and reachable address. Defaults to True.

    Returns:
        str: The full, direct download URL for the file.

    Raises:
        RuntimeError: If the URL cannot be opened after multiple retries.
        NotImplementedError: If the DOI provider is not supported.
    """
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        uri = f"https://zenodo.org/record/{zenodo_entry_number}/files/{fi_name}"

        if validate_uri:
            try:
                # Use the decorated helper to check if the URI exists
                with _open_url_with_retry(uri):
                    pass
            except urllib.error.HTTPError:
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
    """Calculates the SHA1 hash of a file.

    Reads the file in chunks to handle large files efficiently if specified.

    Args:
        fi (str): The path to the file.
        step (int): The chunk size in bytes for reading the file.
            Defaults to 64MB. Only used if `one_block` is False.
        one_block (bool): If True, reads the first `step` bytes of the file.
            If False, reads the entire file in chunks of `step` bytes.
            Defaults to True.

    Returns:
        str: The hexadecimal SHA1 hash of the file content.
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
    """Creates a nested output directory structure to save simulation results.

    The directory structure is generated based on the hashes of the simulation
    input files.

    Args:
        sim (Mapping): A dictionary containing simulation metadata, including
            the "SOFTWARE" key.
        sim_hashes (Mapping): A dictionary mapping file types (e.g., "TPR",
            "TRJ") to their hash information. The structure is expected to be
            `{'TYPE': [('filename', 'hash')]}`.
        out (str): The root output directory where the nested structure
            will be created.
        dry_run_mode (bool): If True, the directory path is resolved but
            not created. Defaults to False.

    Returns:
        str: The full path to the created output directory.

    Raises:
        FileExistsError: If the target output directory already exists and is
            not empty.
        NotImplementedError: If the simulation software is not supported.
        RuntimeError: If the target output directory could not be created.
    """
    # resolve output dir naming
    software = sim.get("SOFTWARE")

    SOFTWARE_CONFIG = {
        "gromacs": {"primary": "TPR", "secondary": "TRJ"},
        "openMM":  {"primary": "TRJ", "secondary": "TRJ"},
        "NAMD":    {"primary": "TRJ", "secondary": "TRJ"},
    }

    config = SOFTWARE_CONFIG.get(software)
    if not config:
        raise NotImplementedError(f"sim software '{software}' not supported")

    primary_hash = sim_hashes.get(config["primary"])[0][1]
    secondary_hash = sim_hashes.get(config["secondary"])[0][1]

    head_dir = primary_hash[:3]
    sub_dir1 = primary_hash[3:6]
    sub_dir2 = primary_hash
    sub_dir3 = secondary_hash

    directory_path = os.path.join(out, head_dir, sub_dir1, sub_dir2, sub_dir3)

    logger.debug(f"output_dir = {directory_path}")

    if os.path.exists(directory_path) and os.listdir(directory_path):
        raise FileExistsError(
            f"Output directory '{directory_path}' is not empty. Delete it if you wish."
        )

    if not dry_run_mode:
        try:
            os.makedirs(directory_path, exist_ok=True)
        except Exception as e:
            raise RuntimeError(f"Could not create the output directory at {directory_path}: {e}")


    return directory_path