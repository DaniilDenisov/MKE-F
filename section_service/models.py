from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, FiniteFloat, model_validator

from section_solver import Options, Profile


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class Limits(Model):
    levels: int = Field(default=5, ge=1, le=8)
    segments: int = Field(default=64, ge=8, le=1024, multiple_of=4)
    area_fraction: float = Field(default=1/32, ge=1e-5, le=1)
    max_elements: int = Field(default=50000, ge=1, le=50000)
    max_bytes: int = Field(default=Options().max_bytes, ge=1, le=Options().max_bytes)
    seconds: float = Field(default=120, gt=0, le=120)
    max_cells: int = Field(default=1000000, ge=1, le=1000000)
    search_tolerance: float = Field(default=.001, gt=0, le=.005)


class Section(Model):
    type: Literal["rectangle","rectangular_tube","circle","circular_tube","square_tube"]
    dimensions: dict[str,FiniteFloat] = Field(max_length=4)

    @model_validator(mode="after")
    def valid_geometry(self):
        try:
            Profile.from_section(self.model_dump())
        except KeyError as err:
            raise ValueError(f"Missing profile dimension {err}") from err
        return self


class Solve(Model):
    section: Section
    nu: float = Field(gt=-1,lt=.5)
    options: Limits = Field(default_factory=Limits)
    J_model: float | None = Field(default=None,gt=0)


class Maximum(Solve):
    forces: list[FiniteFloat] = Field(min_length=6,max_length=6)


class Point(Maximum):
    y: FiniteFloat
    z: FiniteFloat
    element: int | None = Field(default=None,ge=0)


Polynomial = Annotated[list[FiniteFloat],Field(min_length=1,max_length=4)]


class Member(Solve):
    coefficients: list[Polynomial] = Field(min_length=6,max_length=6)
    include_mesh_peak: bool = False
    max_samples: int = Field(default=257,ge=2,le=257)


class Truss(Model):
    axial: Polynomial
    area: float = Field(gt=0)

