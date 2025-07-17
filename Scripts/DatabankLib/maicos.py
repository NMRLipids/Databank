import json

import maicos
from maicos.core import ProfilePlanarBase
from maicos.lib.weights import density_weights
import numpy as np
import MDAnalysis as mda

from .jsonEncoders import CompactJSONEncoder


class NumpyArrayEncoder(CompactJSONEncoder):
    """Encoder for 2xN numpy arrays to be used with json.dump."""

    def encode(self, o):
        if isinstance(o, np.ndarray):
            return CompactJSONEncoder.encode(self, o.tolist())
        else:
            return CompactJSONEncoder.encode(self, o)



class FormFactorPlanar(ProfilePlanarBase):
    """Form factor of a planar system based on the linear electron density profile."""

    def __init__(
        self,
        atomgroup: mda.AtomGroup,
        unwrap: bool = True,
        dim: int = 2,
        zmin: float | None = None,
        zmax: float | None = None,
        bin_width: float = 1,
        refgroup: mda.AtomGroup | None = None,
        pack: bool = True,
        output: str = "form_factor.dat",
        concfreq: int = 0,
        jitter: float = 0.0,
    ) -> None:
        self._locals = locals()
        super().__init__(
            atomgroup=atomgroup,
            unwrap=unwrap,
            pack=pack,
            jitter=jitter,
            concfreq=concfreq,
            dim=dim,
            zmin=zmin,
            zmax=zmax,
            bin_width=bin_width,
            refgroup=refgroup,
            sym=False,
            sym_odd=False,
            grouping="atoms",
            bin_method="com",
            output=output,
            weighting_function=density_weights,
            weighting_function_kwargs={"dens": "electron"},
            normalization="volume",
        )

        self.results.scattering_vectors = np.linspace(0, 1, 1000)

    def _single_frame(self) -> float:
        super()._single_frame()

        bin_pos = self._obs.bin_pos - self._obs.box_center[self.dim]

        # Define bulk region as the first 3.3 Å (two water layers) from the box edges
        # `bin_pos[-1]` is the bin corresponding to the box edge
        bulk_mask = np.abs(bin_pos) > bin_pos[-1] - 3.3
        bulk = self._obs.profile[bulk_mask].mean()

        delta = (self._obs.profile - bulk)[:, np.newaxis]
        angles = self.results.scattering_vectors * bin_pos[:, np.newaxis]

        # Here delta is e/A^3 and _bin_width in A,
        # so the result is 100 times lower than if we calculate in nm
        self._obs.ff_real = np.sum(delta * np.cos(angles) * self._bin_width, axis=0)
        self._obs.ff_imag = np.sum(delta * np.sin(angles) * self._bin_width, axis=0)

        # This value at q=0 will be used for a correlation analysis and error estimate
        return np.sqrt(self._obs.ff_real[0] ** 2 + self._obs.ff_imag[0] ** 2)

    def _conclude(self) -> None:
        super()._conclude()

        self.results.form_factor = np.sqrt(
            self.means.ff_real**2 + self.means.ff_imag**2
        )

        # error from error propagation of the form factor
        self.results.dform_factor = np.sqrt(
            (self.sems.ff_real * self.means.ff_real / self.results.form_factor) ** 2
            + (self.sems.ff_imag * self.means.ff_imag / self.results.form_factor) ** 2
        )

    def save(self):
        # perform unit conversion from Å to nm
        # (see comments in _single_frame)
        output = np.vstack(
            [
                self.results.scattering_vectors,
                self.results.form_factor * 1e2,
                self.results.dform_factor * 1e2,
            ]
        ).T
        with open(self.output, "w") as f:
            json.dump(output, f, cls=NumpyArrayEncoder)


class DensityPlanar(maicos.DensityPlanar):
    def save(self):
        # perform unit conversion from Å to nm and e/Å^3 to e/nm^3
        outdata = np.vstack(
            [
                self.results.bin_pos / 10,
                self.results.profile * 1e3,
                self.results.dprofile * 1e3,
            ]
        ).T
        with open(self.output, "w") as f:
            json.dump(outdata, f, cls=NumpyArrayEncoder)


class DielectricPlanar(maicos.DielectricPlanar):
    def save(self):
        outdata_perp = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.eps_perp,
                self.results.deps_perp,
            ]
        ).T

        with open(f"{self.output_prefix}_perp.json", "w") as f:
            json.dump(outdata_perp, f, cls=NumpyArrayEncoder)

        outdata_par = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.eps_par,
                self.results.deps_par,
            ]
        ).T

        with open(f"{self.output_prefix}_par.json", "w") as f:
            json.dump(outdata_par, f, cls=NumpyArrayEncoder)


class DiporderPlanar(maicos.DiporderPlanar):
    def save(self):
        outdata = np.vstack(
            [
                self.results.bin_pos / 10,  # Convert from Å to nm
                self.results.profile,
                self.results.dprofile,
            ]
        ).T
        with open(self.output, "w") as f:
            json.dump(outdata, f, cls=NumpyArrayEncoder)
