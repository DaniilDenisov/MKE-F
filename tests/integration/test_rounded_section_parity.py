"""Run after Octave's test_spatial_sections, inside the numeric Docker image."""
import json
from pathlib import Path

import numpy as np

from section_solver import Profile

rows = json.loads(Path("output/spatial-section-parity.json").read_text())
count = 0
for row in rows:
    if not row["section"]["dimensions"].get("inner_radius"):
        continue
    profile = Profile.from_section(row["section"])
    expected = [row["properties"][k] for k in ("area", "Iy", "Iz")]
    np.testing.assert_allclose(profile.properties, expected, rtol=1e-12)
    # Independent polygon boundary integrals also catch a shared formula error.
    shape = profile.polygon(4096)
    moments = np.zeros(3)
    for ring in [shape.exterior, *shape.interiors]:
        a, b = np.asarray(ring.coords)[:-1], np.asarray(ring.coords)[1:]
        cross = a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
        iz, iy = cross @ (a*a+a*b+b*b)/12
        moments += [cross.sum()/2, iy, iz]
    np.testing.assert_allclose(abs(moments), expected, rtol=1e-6)
    count += 1
assert count == 2
print("Rounded section Octave/Python parity and independent boundary integrals PASS")
