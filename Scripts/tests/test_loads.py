"""
`test_loads` tests ONLY functions related to downloading files and/or resolving links.

NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import os
from urllib.error import HTTPError
import pytest

# run only without mocking data
pytestmark = pytest.mark.nodata

class TestDownloadResourceFromUri:

    TESTFILENAME = 't.tpr'

    def test_justdl__download_resource_from_uri(self):
        import DatabankLib.databankio as dio
        if os.path.exists(self.TESTFILENAME):
            os.remove(self.TESTFILENAME)  # just for sure
        # first-time download
        assert dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                './' + self.TESTFILENAME) == 0
        # repeat download
        assert dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                './' + self.TESTFILENAME) == 1
        os.remove(self.TESTFILENAME)

    def test_corrupted__download_resource_from_uri(self):
        import DatabankLib.databankio as dio
        # redownload corrupted file
        with open(self.TESTFILENAME, 'w') as f:
            f.write('BABABA')
        assert dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                './' + self.TESTFILENAME) == 2
        os.remove(self.TESTFILENAME)

    def test_errs__download_resource_from_uri(self):
        import DatabankLib.databankio as dio
        # put directory instead of filename
        with pytest.raises(IsADirectoryError) as _:
            dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                './')
        # ask to write to file which you don't have an access
        with pytest.raises(PermissionError) as _:
            dio.download_resource_from_uri(
                "https://zenodo.org/records/8435138/files/pope-md313rfz.tpr",
                '/no-rights.tpr')


# resolve_doi_url
class TestResolveDoiUrl:

    def test_badDOI__resolve_doi_url(self):
        import DatabankLib.databankio as dio
        # test if bad DOI fails
        with pytest.raises(HTTPError, match='404') as _:
            dio.resolve_doi_url('10.5281/zenodo.8435a', True)
        # bad DOI doesn't fail if not to check
        assert ("https://doi.org/10.5281/zenodo.8435a" ==
                dio.resolve_doi_url('10.5281/zenodo.8435a', False))

    def test_goodDOI__resolve_doi_url(self):
        import DatabankLib.databankio as dio
        # good DOI works properly
        assert ("https://doi.org/10.5281/zenodo.8435138" ==
                dio.resolve_doi_url('10.5281/zenodo.8435138', True))

    # resolve_download_file_url

    def test__resolve_download_file_url(self):
        pass
