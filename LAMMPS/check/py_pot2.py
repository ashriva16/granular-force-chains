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

# class customhooke(LAMMPSPairPotential):
#     def __init__(self):
#         super(customhooke, self).__init__()
#         self.units = 'lj'
        disorder = 0
#         r1 = 1+disorder
#         r2 = 1-disorder
#         Krep = 1e5
#         Katt = 1e5
#         self.coeff = {
#             'lj1': {
#                 'lj1': (Krep, Katt, r1, r1),
#                 'lj2': (Krep, Katt, r1, r2),
#                 'lj3': (Krep, Katt, r1, 1)
#             },
#             'lj2': {
#                 'lj1': (Krep, Katt, r2, r1),
#                 'lj2': (Krep, Katt, r2, r2),
#                 'lj3': (Krep, Katt, r2, 1)
#             },
#             'lj3': {
#                 'lj1': (Krep, Katt, 1, r1),
#                 'lj2': (Krep, Katt, 1, r2),
#                 'lj3': (Krep, Katt, 1, 1)
#             },
#         }
#         self.cohesion = 0.01

#     def compute_force(self, rsq, itype, jtype):
#         coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
#         Krep = coeff[0]
#         Katt = coeff[1]
#         r1 = coeff[2]*.5
#         r2 = coeff[3]*.5
#         radsum = r1+r2
#         radcut = radsum+self.cohesion*(radsum)
#         # if(rsq <= radcut*radcut):
#         #     F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#         # else:
#         #     F = 0

#         if(rsq <= radsum*radsum):
#             F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#         elif(rsq > radsum*radsum and rsq <= radcut*radcut):
#             F = Katt*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#         else:
#             F = 0

#         # print(math.sqrt(rsq), itype, jtype, r1, r2, F)
#         return F

#     def compute_energy(self, rsq, itype, jtype):
#         coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
#         Krep = coeff[0]
#         Katt = coeff[1]
#         r1 = coeff[2]*.5
#         r2 = coeff[3]*.5
#         radsum = r1+r2
#         radcut = radsum+self.cohesion*(radsum)
#         # if(rsq <= radcut*radcut):
#         #     U = .5*Krep*((radsum-math.sqrt(rsq))**2 -
#         #                  (radsum-radcut)**2)  # end 0
#         # else:
#         #     U = 0
#         # return U

#         if(rsq <= radsum*radsum):
#             U = .5*Krep*(radsum-math.sqrt(rsq))**2-.5 * \
#                 Katt*(radsum-radcut)**2  # end 0
#         elif(rsq > radsum*radsum and rsq <= radcut*radcut):
#             U = .5*Katt*(radsum-math.sqrt(rsq))**2-.5 * \
#                 Katt*(radsum-radcut)**2  # start 0
#         else:
#             U = 0
#         return U


# class customhooke2(LAMMPSPairPotential):
#     def __init__(self):
#         super(customhooke, self).__init__()
#         self.units = 'lj'
#         disorder = 0.0
#         # ri is actually a diameter
#         r1 = 1+disorder
#         r2 = 1-disorder
#         R = 25
#         Krep = 1e5
#         Katt = 1e5
#         self.coeff = {
#             'lj1': {
#                 'lj1': (Krep, Katt, r1, r1),
#                 'lj2': (Krep, Katt, r1, r2),
#                 'lj3': (Krep, 0, R, 0)
#             },
#             'lj2': {
#                 'lj1': (Krep, Katt, r2, r1),
#                 'lj2': (Krep, Katt, r2, r2),
#                 'lj3': (Krep, 0, R, 0)
#             },
#             'lj3': {
#                 'lj1': (Krep, 0, R, 0),
#                 'lj2': (Krep, 0, R, 0),
#                 'lj3': (0, 0, R, 0)
#             },
#         }
#         self.cohesion = .01

#     def compute_force(self, rsq, itype, jtype):
#         coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
#         Krep = coeff[0]
#         Katt = coeff[1]
#         r1 = coeff[2]*.5
#         r2 = coeff[3]*.5
#         radsum = r1+r2
#         radcut = radsum+self.cohesion*(radsum)
#         # if(rsq <= radcut*radcut):
#         #     F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#         # else:
#         #     F = 0
#         if(itype == 3 or jtype == 3):
#             if(rsq <= radsum*radsum):
#                 F = Krep*(radsum-math.sqrt(rsq))**2/math.sqrt(rsq)
#             else:
#                 F = 0
#         else:
#             if(rsq <= radsum*radsum):
#                 F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#             elif(rsq > radsum*radsum and rsq <= radcut*radcut):
#                 F = Katt*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
#             else:
#                 F = 0

#         return F

#     def compute_energy(self, rsq, itype, jtype):
#         coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
#         Krep = coeff[0]
#         Katt = coeff[1]
#         r1 = coeff[2]*.5
#         r2 = coeff[3]*.5
#         radsum = r1+r2
#         radcut = radsum+self.cohesion*(radsum)
#         # if(rsq <= radcut*radcut):
#         #     U = .5*Krep*((radsum-math.sqrt(rsq))**2 -
#         #                  (radsum-radcut)**2)  # end 0
#         # else:
#         #     U = 0
#         # return U

#         if(itype == 3 or jtype == 3):
#             if(rsq <= radsum*radsum):
#                 U = (1/3)*Krep*(radsum-math.sqrt(rsq))**3
#             else:
#                 U = 0
#         else:
#             if(rsq <= radsum*radsum):
#                 U = .5*Krep*(radsum-math.sqrt(rsq))**2  # end 0
#             elif(rsq > radsum*radsum and rsq <= radcut*radcut):
#                 U = .5*Katt*(radsum-math.sqrt(rsq))**2  # start 0
#             else:
#                 U = .5 * Katt*(radsum-radcut)**2

#         return U


class customhooke(LAMMPSPairPotential):
    def __init__(self):
        super(customhooke, self).__init__()
        self.units = 'lj'
        disorder = 0
        att_coeff = 1
        rep_coeff1 = 1
        rep_coeff2 = 1

        #Krep = 1e5
        #Katt = 0

        r1 = 1+disorder
        r2 = 1-disorder

        Krep1 = 1*rep_coeff1
        Krep2 = 1*rep_coeff2
        Krepmid = (Krep1+Krep2)/2
        Katt1 = att_coeff*Krep1
        Katt2 = att_coeff*Krep2
        Kattmid = att_coeff*Krepmid

        self.coeff = {
            'lj1': {
                'lj1': (Krep1, Katt1, r1, r1),
                'lj2': (Krepmid, Kattmid, r1, r2),
                'lj3': (Krepmid, 0, r1, 1)
            },
            'lj2': {
                'lj1': (Krepmid, Kattmid, r2, r1),
                'lj2': (Krep2, Katt2, r2, r2),
                'lj3': (Krepmid, 0, r2, 1)
            },
            'lj3': {
                'lj1': (Krepmid, 0, 1, r1),
                'lj2': (Krepmid, 0, 1, r2),
                'lj3': (0, 0, 1, 1)
            },
        }
        self.cohesion = .01

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        Krep = coeff[0]
        Katt = coeff[1]
        r1 = coeff[2]*.5
        r2 = coeff[3]*.5
        radsum = r1+r2
        radcutmid = radsum+.5*self.cohesion*(radsum)
        radcut = radsum+self.cohesion*(radsum)
 
        if(rsq <= radsum*radsum):
            F = Krep*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
        elif(rsq > radsum*radsum and rsq <= radcutmid*radcutmid):
            F = Katt*(radsum-math.sqrt(rsq))/math.sqrt(rsq)
        elif(rsq > radcutmid*radcutmid and rsq <= radcut*radcut):
            F = (-Katt*(radsum-math.sqrt(rsq))- Katt*self.cohesion*(radsum))/math.sqrt(rsq)
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
        radcutmid = radsum+.5*self.cohesion*(radsum)
        radcut = radsum+self.cohesion*(radsum)
        d = self.cohesion

        U_C = -Katt*radsum**2*(d**2/4+d)

        if(rsq <= radsum*radsum):
            U = .5*Krep*(radsum-math.sqrt(rsq))**2  # end 0
        elif(rsq > radsum*radsum and rsq <= radcutmid*radcutmid):
            U = .5*Katt*(radsum-math.sqrt(rsq))**2  # start 0
        elif(rsq > radcutmid*radcutmid and rsq <= radcut*radcut):
            U = -.5 * Katt*(radsum-math.sqrt(rsq))**2 + Katt * \
                d*radsum*math.sqrt(rsq) + U_C
        else:
            U = .25 * Katt*radsum**2*d**2
        return U
