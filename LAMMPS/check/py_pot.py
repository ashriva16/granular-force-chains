from __future__ import print_function
import math

class LAMMPSPairPotential(object):
    def __init__(self):
        self.pmap=dict()
        self.units='lj'
    def map_coeff(self,name,ltype):
        self.pmap[ltype]=name
    def check_units(self,units):
        if (units != self.units):
           raise Exception("Conflicting units: %s vs. %s" % (self.units,units))

class LJCut(LAMMPSPairPotential):
    def __init__(self):
        super(LJCut, self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = 'lj'
        self.coeff = {'lj'  : {'lj'  : (48.0,24.0,4.0,4.0)}}

    def compute_force(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self,rsq,itype,jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv  = 1.0/rsq
        r6inv  = r2inv*r2inv*r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return (r6inv * (lj3*r6inv - lj4))


class LJCutshift(LAMMPSPairPotential):
    def __init__(self):
        super(LJCutshift, self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = 'lj'
        self.coeff = {'lj': {'lj': (48.0, 24.0, 4.0, 4.0)}}

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0/rsq
        r6inv = r2inv*r2inv*r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r_c = 2.5
        if(rsq < r_c**2):
            r2inv = 1.0/rsq
            r2inv_c = 1.0/r_c**2
            r6inv = r2inv*r2inv*r2inv
            r6inv_c = r2inv_c*r2inv_c*r2inv_c
            lj3 = coeff[2]
            lj4 = coeff[3]
            return (r6inv * (lj3*r6inv - lj4)) - (r6inv_c * (lj3*r6inv_c - lj4))
        else:
            return 0

class customhooke(LAMMPSPairPotential):
    def __init__(self):
        super(customhooke, self).__init__()
        self.units = 'lj'
        disorder = 0
        r1 = 1+disorder
        r2 = 1-disorder
        Krep = 1e5
        Katt = 1e5
        self.coeff = {
            'lj1': {
                'lj1': (Krep, Katt, r1, r1),
                'lj2': (Krep, Katt, r1, r2),
                'lj3': (Krep, Katt, r1, 1)
            },
            'lj2': {
                'lj1': (Krep, Katt, r2, r1),
                'lj2': (Krep, Katt, r2, r2),
                'lj3': (Krep, Katt,r2, 1)
            },
            'lj3': {
                'lj1': (Krep, Katt, 1, r1),
                'lj2': (Krep, Katt, 1, r2),
                'lj3': (Krep, Katt, 1, 1)
            },
        }
        self.cohesion = 0.1

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        Krep = coeff[0]
        Katt = coeff[1]
        r1 = coeff[2]*.5
        r2 = coeff[3]*.5
        radsum = r1+r2
        radcut = radsum+self.cohesion*(radsum)
        # if(rsq <= radcut*radcut):
        #     F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
        # else:
        #     F = 0

        if(rsq <= radsum*radsum):
            F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
        elif(rsq > radsum*radsum and rsq <= radcut*radcut):
            F = Katt*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
        else:
            F = 0

        return F

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        Krep = coeff[0]
        Katt = coeff[1]
        r1 = coeff[2]*.5
        r2 = coeff[3]*.5
        radsum = r1+r2
        radcut = radsum+self.cohesion*(radsum)
        # if(rsq <= radcut*radcut):
        #     U = .5*Krep*((radsum-math.sqrt(rsq))**2 -
        #                  (radsum-radcut)**2)  # end 0
        # else:
        #     U = 0
        # return U

        if(rsq <= radsum*radsum):
            U = .5*Krep*(radsum-math.sqrt(rsq))**2-.5 * Katt*(radsum-radcut)**2  # end 0
        elif(rsq > radsum*radsum and rsq <= radcut*radcut):
            U = .5*Katt*(radsum-math.sqrt(rsq))**2-.5 * Katt*(radsum-radcut)**2  # start 0
        else:
            U = 0
        return U
