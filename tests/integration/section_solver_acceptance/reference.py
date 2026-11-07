"""Independent Saint-Venant references in physical coordinates, used only in tests."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def tube_j(h,b,t,n=320):
    """Cell-centred Prandtl problem. All hole walls share one unknown potential."""
    step = min(h,b)/n
    ny,nz,nt = np.rint(np.array([h,b,t])/step).astype(int)
    assert np.allclose([ny*step,nz*step,nt*step],[h,b,t])
    grid = np.full((ny+2,nz+2),-1)
    material = np.ones((ny,nz),bool)
    material[nt:-nt,nt:-nt] = False
    y,z = np.where(material)
    y,z = y+1,z+1
    ids = np.arange(len(y))
    hole = len(ids)
    grid[1:-1,1:-1][~material] = hole
    grid[y,z] = ids
    rr,cc,vv = [],[],[]
    def entries(r,c,v):
        rr.extend(r)
        cc.extend(c)
        vv.extend(v)
    for dy,dz in ((1,0),(-1,0),(0,1),(0,-1)):
        other = grid[y+dy,z+dz]
        wall = (other==-1)|(other==hole)
        w = np.where(wall,2.,1.)
        entries(ids,ids,w)
        keep = other>=0
        entries(ids[keep],other[keep],-w[keep])
        inner = other==hole
        entries(other[inner],ids[inner],-w[inner])
        entries(other[inner],other[inner],w[inner])
    rhs = np.full(hole+1,2*step**2)
    rhs[-1] = 2*(h-2*t)*(b-2*t)
    matrix = coo_matrix((vv,(rr,cc)),shape=(len(rhs),len(rhs))).tocsc()
    potential = spsolve(matrix,rhs)
    assert np.linalg.norm(matrix@potential-rhs)<1e-9*np.linalg.norm(rhs)
    return float(rhs@potential)


def circle_shear(y,z,R,r,I,nu):
    a = (3+2*nu)/(8*(1+nu)*I)
    b = (2*nu-1)/(8*(1+nu)*I)
    c = -(1+2*nu)/(4*(1+nu)*I)
    dipole = a*(R*r/(y*y+z*z))**2 if r else 0.
    return a*(R*R+r*r-y*y)+b*z*z+dipole*(z*z-y*y), c*y*z-2*dipole*y*z


def hyper(x,a):
    # Stable cosh(x)/cosh(a) and sinh(x)/cosh(a) for |x| <= a.
    u = np.abs(x)
    scale = np.exp(u-a)/(1+np.exp(-2*a))
    e = np.exp(-2*u)
    return scale*(1+e),scale*np.sign(x)*(1-e)


def rectangle_shear(y,z,a,b,I,nu):
    k = nu/((1+nu)*I)
    base = (a*a-y*y)/(2*I)
    if b>a:
        ty,tz = base-k*(a*a-y*y)/2,-k*z*y
        for i in range(192):
            m = 2*i+1
            w = m*np.pi/(2*a)
            fac = 8*k*a*b*(-1)**i/(np.pi*m)**2/np.tanh(w*b)
            ch,sh = hyper(w*z,w*b)
            ty = ty+fac*ch*np.cos(w*y)
            tz = tz+fac*sh*np.sin(w*y)
    else:
        ty,tz = base+k*(3*z*z-b*b)/6,np.zeros_like(y)
        for m in range(1,193):
            w = m*np.pi/b
            fac = 2*k*b*b*(-1)**m/(m*np.pi)**2
            ch,sh = hyper(w*y,w*a)
            ty = ty-fac*ch*np.cos(w*z)
            tz = tz+fac*sh*np.sin(w*z)
    return np.where(abs(y)==a,0.,ty),np.where(abs(z)==b,0.,tz)


class Reference:
    def __init__(self,section):
        kind,d = section["type"],section["dimensions"]
        self.round = "d" in d
        self.sharp = "t" in d and not self.round
        self.h = d.get("h",d.get("b",d.get("d")))
        self.b = d.get("b",self.h)
        self.t = d.get("t",0.)
        h,b,t = self.h,self.b,self.t
        if self.round:
            R,r = h/2,h/2-t if t else 0.
            self.A = np.pi*(R**2-r**2)
            self.Iy = self.Iz = np.pi*(R**4-r**4)/4
            self.J = 2*self.Iy
        else:
            hi,bi = (h-2*t,b-2*t) if t else (0.,0.)
            self.A = h*b-hi*bi
            self.Iy,self.Iz = (h*b**3-hi*bi**3)/12,(b*h**3-bi*hi**3)/12
            if t:
                self.J = tube_j(h,b,t)
            else:
                a,c = max(h,b),min(h,b)
                m = np.arange(1,1000,2,dtype=float)
                self.J = a*c**3/3*(1-192*c/(np.pi**5*a)*sum(np.tanh(m*np.pi*a/(2*c))/m**5))

    def fields(self,y,z,nu):
        if self.round:
            R,r = self.h/2,self.h/2-self.t if self.t else 0.
            vy = circle_shear(y,z,R,r,self.Iz,nu)
            vz = circle_shear(z,y,R,r,self.Iy,nu)[::-1]
            tor = -z/self.J,y/self.J
        else:
            h,b = self.h,self.b
            vy = rectangle_shear(y,z,h/2,b/2,self.Iz,nu)
            vz = rectangle_shear(z,y,b/2,h/2,self.Iy,nu)[::-1]
            swapped = b>h
            u,v = (z,y) if swapped else (y,z)
            a,c = max(h,b)/2,min(h,b)/2
            dv,du = -2*v,np.zeros_like(v)
            for i in range(192):
                m = 2*i+1
                w = m*np.pi/(2*c)
                fac = 16*c*(-1)**i/(np.pi*m)**2
                ch,sh = hyper(w*u,w*a)
                dv = dv+fac*ch*np.sin(w*v)
                du = du+fac*sh*np.cos(w*v)
            dv,du = np.where(abs(u)==a,0.,dv),np.where(abs(v)==c,0.,du)
            tor = (-du/self.J,-dv/self.J) if swapped else (dv/self.J,du/self.J)
        return np.stack([np.stack(f,axis=-1) for f in (vy,vz,tor)],axis=-2)

    def grid(self):
        if self.round:
            r,theta = np.meshgrid(np.linspace(self.h/2-self.t if self.t else 0,self.h/2,129),
                                  np.linspace(0,2*np.pi,513))
            return (r*np.cos(theta)).ravel(),(r*np.sin(theta)).ravel()
        y,z = np.meshgrid(np.linspace(-self.h/2,self.h/2,161),np.linspace(-self.b/2,self.b/2,161))
        return y.ravel(),z.ravel()
