from dataclasses import dataclass
from math import isfinite, pi
from numbers import Real

from shapely.geometry import Point, box


def positive(x):
    return isinstance(x, Real) and not isinstance(x, bool) and isfinite(x) and x > 0


def validate_nu(nu):
    if isinstance(nu, bool) or not isinstance(nu, Real) or not -1 < nu < .5:
        raise ValueError("Poisson ratio must satisfy -1 < nu < 0.5")
    return float(nu)


@dataclass(frozen=True)
class Profile:
    kind: str
    h: float
    b: float
    t: float = 0
    inner_radius: float = 0

    def __post_init__(self):
        if self.kind not in ("circle", "rectangle", "square_tube", "rectangular_tube", "circular_tube"):
            raise ValueError("Unsupported profile")
        if not all(map(positive, (self.h, self.b))):
            raise ValueError("Dimensions must be finite and positive")
        if self.kind in ("square_tube", "circle", "circular_tube") and self.h != self.b:
            raise ValueError("Equal outside dimensions required")
        if self.kind.endswith("_tube"):
            if not positive(self.t) or 2*self.t >= min(self.h, self.b):
                raise ValueError("Invalid wall thickness")
        elif self.t != 0:
            raise ValueError("Solid profile has no wall thickness")
        if self.inner_radius != 0 and (self.kind not in ("square_tube","rectangular_tube")
                or not positive(self.inner_radius) or self.inner_radius+self.t >= min(self.h,self.b)/2*(1-1e-14)):
            raise ValueError("Tube corner radii must fit between the straight walls")
        try:
            if not all(map(positive,self.properties)):
                raise ValueError("Dimensions exceed the numerical range")
        except OverflowError as err:
            raise ValueError("Dimensions exceed the numerical range") from err

    @classmethod
    def from_section(cls, section):
        kind, dims = section["type"], section["dimensions"]
        if kind in ("circle", "circular_tube"):
            h = b = dims["d"]
        else:
            b = dims["b"]
            h = b if kind == "square_tube" else dims["h"]
        return cls(kind, h, b, dims["t"] if kind.endswith("_tube") else 0,dims.get("inner_radius",0))

    @property
    def round(self):
        return self.kind in ("circle", "circular_tube")

    @property
    def sharp(self):
        return self.kind in ("rectangular_tube", "square_tube") and not self.inner_radius

    @property
    def properties(self):
        """Analytic A, Iy, Iz."""
        h, b = self.h, self.b
        hi, bi = (h-2*self.t, b-2*self.t) if self.t else (0, 0)
        if self.round:
            inertia = pi*(h**4-hi**4)/64
            return pi*(h**2-hi**2)/4, inertia, inertia
        if self.inner_radius:
            # Subtract the square-minus-quadrant at each of the four corners.
            def rounded(h,b,r):
                cut = (1-pi/4)*r**2
                def inertia(span,width):
                    d = span/2-r
                    return width*span**3/12-4*(d*d*cut+d*r**3/3+(1/3-pi/16)*r**4)
                return h*b-4*cut,inertia(b,h),inertia(h,b)
            outer = rounded(h,b,self.inner_radius+self.t)
            inner = rounded(hi,bi,self.inner_radius)
            return tuple(a-b for a,b in zip(outer,inner))
        return h*b-hi*bi, (h*b**3-hi*bi**3)/12, (b*h**3-bi*hi**3)/12

    def polygon(self, segments):
        if not isinstance(segments, int) or segments < 8 or segments % 4:
            raise ValueError("Use a multiple of four circular segments")
        def boundary(h, b, r):
            if self.round:
                return Point(0, 0).buffer(h/2, quad_segs=segments//4)
            if r:
                return box(-h/2+r,-b/2+r,h/2-r,b/2-r).buffer(r,quad_segs=segments//4)
            return box(-h/2, -b/2, h/2, b/2)
        outer = boundary(self.h, self.b,self.inner_radius+self.t if self.inner_radius else 0)
        return outer.difference(boundary(self.h-2*self.t,self.b-2*self.t,self.inner_radius)) if self.t else outer

