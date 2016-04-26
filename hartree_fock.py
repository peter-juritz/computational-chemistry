from basis import Basis
from math import sqrt, pi, exp, factorial
from gaussian import Primative, ContrPrimative
from scipy import linalg, real, dot, transpose
from scipy.special import erf
import time
import numpy
from util import (my_dot, scal_dot, vec_add, vec_minus, vec_cross, vec_dot, vec_times, dist)


class HartreeFock:

    def __init__(self,mol,bas_file,do_print=0):
        print 'Hf module started'
        self.mol = mol
        bas = Basis()
        atom_functions = bas.read_basis(bas_file)
        self.atom_functions = atom_functions

        self.do_print=do_print

    def converge(self,F, H_core, G, P,C,S_arb,ERI,aof,mol,lim=0.00001):
        print '----------- Converging Electronic Energy -----------'
        E = None
        for i in range(3000):
            t = self.iterate(F, H_core, G, P,C,S_arb,ERI,aof)
            print 'Iteration  ', i, ' Electronic energy: ', t

            if E and (abs(E-t) < lim):
                nuc_eng = self.nuclear_energy(mol)
                print 'Converged at ', t, ', Including Nuclear energy : ', t + nuc_eng
                return t + self.nuclear_energy(mol)
            E = t

        print 'No convergance'
        return E + h.nuclear_energy(mol)

    def build_atomic_orbitals(self,ao,mol):
        functions = []
        for v in xrange(len(mol.atoms)):
            for func in ao[mol.atoms[v].name.lower()][1:]: # TODO check this
                cp = ContrPrimative()
                off = len(functions)
                functions.append(cp)
                mol.atoms[v].aof_functions.append(off)
                if (func[0] == 'S'):
                    cp.type = 0
                else:
                    cp.type = 1

                for j in func[1:]:
                    if (cp.type==0):
                        #print 'Adding ',at.name,' coef=',j[0]
                        cp.add(Primative(j[1],j[0],0,0,0,mol.atoms[v].xyz))
                    else:
                        #print 'Sp function!'
                        cp.add(Primative(j[1],j[0],0,0,0,mol.atoms[v].xyz)) #2S
                        cp.add(Primative(j[2],j[0],1,0,0,mol.atoms[v].xyz)) #2px
                        cp.add(Primative(j[2],j[0],0,1,0,mol.atoms[v].xyz)) #2py
                        cp.add(Primative(j[2],j[0],0,0,1,mol.atoms[v].xyz)) #2pz

        return functions

    def norm_const(self,cp): # s type functions only
        norm = 0
        for n in xrange(cp.n):
            for m in xrange(cp.n):
                norm += (cp.prim(n).coef*cp.prim(m).coef)/((cp.prim(n).exp + cp.prim(m).exp )**1.5)
        norm= (norm**-0.5)
        norm*= (pi)**-0.75
        return norm

    def sum_cs(self,cu,cv): # s s contracted integral
        sum = 0
        for u in cu.primatives:
            for v in cv.primatives:
                AA = my_dot(u.centre,u.centre)
                BB = my_dot(v.centre,v.centre)
                AB = my_dot(u.centre,v.centre)

                A_B = vec_minus(u.centre,v.centre)
                A_B= my_dot(A_B,A_B)
                sum += u.coef * v.coef *  exp((-1*u.exp * v.exp*(AA -2*AB + BB))/(u.exp + v.exp)) * (pi/(u.exp + v.exp))**1.5 * self.norm_prim(u)*self.norm_prim(v)
        return sum

    def choose(self,n,k):
        c = factorial(n)/(factorial(k)*factorial(n-k))
        return (c)

    def f_k(self,u,v,k,PA,PB,l1,l2):
        f = 0
        q = max(-1*k,k - 2*l2)
        end = min(k,2*l1-k) + 1
        y= u.exp + v.exp
        ##print 'f_k: looping ',q,' to ',end
        for i in xrange(q,end,2):
    #        #print 'exp is ',l1-i
            f+=self.choose(l1,(k+q)/2) * self.choose(l2,(k-q)/2) * PA**(l1-(k+q)/2) * PB**(l2 - (k-q)/2)
    #        f+=choose(l1,i) * PA**(l1-i) * PB**(l2-k+i) * choose(l2,k-i)
    #    #print
        return f

    def double_fac(self,n):
        r = 1
        for i in xrange(1,n,2):
            r*=i
        return r

    def I(self,u,v,var): #var = x,y,z
        ##print '--------------- I',var,' ----------------'
        I_sum = 0
        y = u.exp + v.exp
        P = vec_add(scal_dot(u.centre,u.exp/y),scal_dot(v.centre,v.exp/y))
        if (var == 'x'):
            l1=u.l
            l2=v.l
            PA = P[0]-u.centre[0]
            PB = P[0]-v.centre[0]
        elif(var=='y'):
            l1=u.m
            l2=v.m
            PA = P[1]-u.centre[1]
            PB = P[1]-v.centre[1]
        elif(var=='z'):
            l1=u.n
            l2=v.n
            PA = P[2]-u.centre[2]
            PB = P[2]-v.centre[2]
    #    PA = dot(z.vec_minus(P,u.centre),z.vec_minus(P,u.centre))
    #    PB=dot(z.vec_minus(P,v.centre),z.vec_minus(P,v.centre))
        ##print 'I',var,' looping to',(l1+l2)/2
        for i in xrange((l1+l2)/2 + 1) :
        #    #print 'i=',i
            I_sum+=self.f_k(u,v,2*i,PA,PB,l1,l2) * (self.double_fac(2*i-1))/((2*y)**i) * (pi/y)**0.5
        return I_sum

    def kinetic_integral(self,functions): # TODO Check for symmetry problems
        N = len(functions)
        T = numpy.zeros((N,N))
        for n in xrange(N):
            for m in xrange(N):
                cu = functions[n]
                cv = functions[m]
                for u in cu.primatives:
                    for v in cv.primatives:
                            tix = v.exp*(2*v.l + 1)*(self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,0,0)))
                            tix+=-2*v.exp**2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(+2,0,0)))
                            tix+= -1*v.l*(v.l-1)/2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(-2,0,0)))

                            tiy = v.exp*(2*v.m + 1)*(self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,0,0)))
                            tiy+=-2*v.exp**2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,+2,0)))
                            tiy+= -1*v.m*(v.m-1)/2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,-2,0)))

                            tiz = v.exp*(2*v.n + 1)*(self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,0,0)))
                            tiz+=-2*v.exp**2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,0,+2)))
                            tiz+= -1*v.n*(v.n-1)/2 * (self.arb_overlap_prim(u.prim_copy(0,0,0),v.prim_copy(0,0,-2)))

                            T[n][m] += (tix + tiy + tiz )* self.norm_prim(u)*self.norm_prim(v)*u.coef*v.coef
      # for n in xrange(N):
      #     for m in xrange(N):
      #         print 'checking T [',n,'][',m,'] = ',T[n][m],' T[',m,'][',n,'] = ',T[m][n]
              #  if (round(T[n][m],6) != round(T[m][n],6)):
              #     print 'BAD T [',n,'][',m,'] = ',T[n][m],' T[',m,'][',n,'] = ',T[m][n]
                #    raw_input('continue:')

        return T

    def arb_overlap(self,cu,cv):
        S = 0
        for u in cu.primatives:
            chunk = 0
            for v in cv.primatives:
                chunk+= self.arb_overlap_prim(u,v)*u.coef*v.coef*self.norm_prim(u)*self.norm_prim(v)
            S+=chunk

        return S

    def arb_overlap_prim(self,u,v): # so far only for primative
        S = 0
        AA = my_dot(u.centre,u.centre)
        BB = my_dot(v.centre,v.centre)
        AB = my_dot(u.centre,v.centre)
        Ix = self.I(u,v,'x')
        Iy = self.I(u,v,'y')
        Iz = self.I(u,v,'z')
        S = exp((-1*u.exp * v.exp*(AA -2*AB + BB))/(u.exp + v.exp)) *Ix*Iy*Iz

        return S


    def overlap(self,functions):
        N = len(functions)
        S = numpy.zeros((N,N))
        for u in xrange(N):
            for v in xrange(N):
                nc1 = self.norm_const(functions[u])
                nc2 = self.norm_const(functions[v])
                cu = functions[u]
                cv = functions[v]
                if (cu.n == 10 and cv.n==10): # TODO prim S S FIX
                    y = cu.prim(0).exp + cv.prim(0).exp
                    S[u][v] = nc1*nc2* cu.prim(0).coef*cv.prim(0).coef* exp(-1*cu.prim(0).exp * cv.prim(0).exp * dist(cu.centre,cv.centre)**2/y) * (pi/y)**1.5
                    continue
                else:
                    sum = self.sum_cs(cu,cv)
                    S[u][v] = sum
                    continue
    #    #print S
        return S

    def prim_nuclear(self,u,v,C,name):
        charges = {'h':1,'he':2,'c':6,'o':8}
        y = u.exp + v.exp
        #AA = my_dot(u.centre,u.centre)
        #BB = my_dot(v.centre,v.centre)
        #AB = my_dot(u.centre,v.centre)
        R1 = vec_minus(u.centre,v.centre)
        R2 = my_dot(R1,R1)


        P = vec_add(scal_dot(u.centre,u.exp/y),scal_dot(v.centre,v.exp/y))
        PC = dist(P,C)
        Z_c = charges[name]#1.2266/(0.87622069)
        if(P[0]==C[0] and P[1]==C[1] and P[2]==C[2]):
            B=Z_c
    #        #print 'P==C'
        else:
    #        #print 'P!=C'
            B =  Z_c * (pi**0.5*erf(y**0.5*PC))/(2*y**0.5 * PC)
        V =-1*(2*pi*exp((-1*u.exp * v.exp*(R2))/(y)))/(y)*B * u.coef*v.coef *((2*u.exp/pi)**0.75) * ((2*v.exp/pi)**0.75 )

        return V

    def nuclear(self,cu,cv,nucl,name):
        V = 0
        for u in cu.primatives:
            for v in cv.primatives:
                V+= self.prim_nuclear(u,v,nucl,name)
        return V#*norm_const(cu)*norm_const(cv)

    def f0(self,t):
        if (t==0):
            return 1
        else:
            return (pi**0.5)/(2*t**0.5) * erf(t**0.5)

    def norm_prim(self,u):
        return (2*u.exp/pi)**0.75

    def electron_rep_prim(self,u,v,j,k):

      yp = u.exp + v.exp
      yq = j.exp + k.exp
      P = vec_add(vec_times(u.centre,u.exp/yp),vec_times(v.centre,v.exp/yp))
      Q = vec_add(vec_times(j.centre,j.exp/yq),vec_times(k.centre,k.exp/yq))

      AA = vec_dot(u.centre,u.centre)
      BB = vec_dot(v.centre,v.centre)
      AB = vec_dot(u.centre,v.centre)

      CC = vec_dot(j.centre,j.centre)
      DD = vec_dot(k.centre,k.centre)
      CD = vec_dot(j.centre,k.centre)

      K1 =exp((-1*u.exp * v.exp*(AA -2*AB + BB))/(yp))
      K2 =exp((-1*j.exp * k.exp*(CC -2*CD + DD))/(yq))
      pq = vec_minus(P,Q)
      PQ = vec_dot(pq,pq)
      n1 = (2*u.exp/pi)**0.75
      n2 = (2*v.exp/pi)**0.75
      n3 = (2*j.exp/pi)**0.75
      n4 = (2*k.exp/pi)**0.75
      E = 2*pi**2.5*K1*K2/(yp*yq*(yp+yq)**0.5) * self.f0(PQ*yp*yq/(yp+yq)) * u.coef*v.coef*j.coef*k.coef * n1*n2*n3*n4
      return E


    def dens(self,x,y,z,aof):
        d=0
        for j in range(self.N): # FIX!
            k=0
            e=0
            for a in aof:
                for b in a.primatives:
                    t = b.eval(x,y,z)
                    e+= t*t
                # Below required to multiply in the
                #converged coefficients after SCF convergence
                # e*= self.C[k][j]
                k+=1
            d+=e
        return d

    def electron_rep(self,cu,cv,cj,ck):
        E = 0
        for u in cu.primatives:
            for v in cv.primatives:
                for j in cj.primatives:
                    for k in ck.primatives:
                        E+=self.electron_rep_prim(u,v,j,k)
        return E #* norm_const(cu)*norm_const(cv)*norm_const(cj)*norm_const(ck)

    def energy_one(self,T,V,n_basis):
        n = len(V)
        E = 0
        for i in xrange(n_basis):
            E+=T[i][i]
            for j in xrange(n):
                E+= V[j][i][i]
        return E

    def energy_two(self,ERI):
        n = len(ERI)
        E = 0
        for i in xrange(n):
            for j in xrange(n):
                E+=ERI[i][i][j][j] - ERI[i][j][j][i]
        return E

    def electronic_energy_v1(self,T,V,ERI,n_basis): # fix later , basis can be found not passed
        E = 0
        E+= self.energy_one(T,V,n_basis) + 0.5*self.energy_two(ERI)
        print 'Electronic energy: ',E
        return E

    def density(self,functions,C,N_elec): # TODO, auto calculate N_elec
        N = len(functions)
        P = numpy.zeros((N,N))
        for u in xrange(N):
            for v in xrange(N):
                s = 0
                for a in xrange(N_elec): # TODO N_elec/2
                    s+=C[u][a] * C[v][a]
                P[u][v] = 2*s
        return P

    def calc_G(self,P,ERI):
        N = len(P)
        G = numpy.zeros((N,N))
        for u in xrange(N):
            for v in xrange(N):
                sum_val = 0
                for a in xrange(N):
                    for s in xrange(N):
                        sum_val +=P[a][s]*(ERI[u][v][s][a] - 0.5*ERI[u][a][s][v])
                G[u][v] = sum_val
        return G

    def electronic_energy(self,P,H,F):
        N = len(P)
        E0 = 0
        for u in xrange(N):
            for v in xrange(N):
                E0+=P[v][u]*(H[u][v] + F[u][v])
        return 0.5*E0

    def nuclear_energy(self,mol):
        Z = {'h':1,'he':2,'o':8,'c':6}
        N = mol.n_atoms
        E=0
        for A in xrange(N):
            for B in xrange(A+1,N):
                E+=Z[mol.atoms[A].name]*Z[mol.atoms[B].name]/(dist(mol.atoms[A].xyz,mol.atoms[B].xyz))
        return E

    def matr_eq(self,a,b):
        n = len(a)
        m = len(a[0])
        for i in range(n):
            for j in range(m):
                if ((a[i][j] - b[i][j]) > 0.0000001):
                    #print 'matrix doesnt match: ',a[i][j],' and ',b[i][j]
                    return False
        return True

    def initial_calculation(self):
        eneg = {'c':2.55,'h':2.1,'o':3.44,'n':3.04}
        self.aof = self.build_atomic_orbitals(self.atom_functions,self.mol)
        self.N = len(self.aof)

        N = self.N
        self.back=[numpy.zeros((N,N)),numpy.zeros((N,N))]
        #print '----------- S(norm) -----------'
        #print self.overlap(self.aof)
        if (self.do_print ==1):
                print '----------- Overlap Matrix S(arb)  -----------'
        self.S_arb = numpy.zeros((len(self.aof),len(self.aof)))
        p=0
        q=0
        for u in xrange(self.N):
            for v in xrange(u +1):
              t = self.arb_overlap(self.aof[u],self.aof[v])
              self.S_arb[u][v] = t
              self.S_arb[v][u] = t
        if (self.do_print ==1):
                print self.S_arb
                print '-----------   Kinetic Energy T    -----------'
        self.T = self.kinetic_integral(self.aof)
        if (self.do_print ==1):
                print self.T
        self.V = numpy.zeros((self.mol.n_atoms,self.N,self.N))
        for i in xrange(len(self.mol.atoms)):
            for n in xrange(self.N):
                for m in xrange(self.N):
                    t = self.nuclear(self.aof[n],self.aof[m],self.mol.atoms[i].xyz,self.mol.atoms[i].name)
                    self.V[i][n][m] =t
            if (self.do_print ==1):
                    print '-----------  Nuclear Attraction V [',i,']  -----------'
                    print self.V[i]
        self.H_core = numpy.zeros((self.N,self.N))
        self.H_core += self.T
        for i in xrange(self.mol.n_atoms):
            self.H_core+= self.V[i]
        if (self.do_print ==1):
                print '-----------  Core Hamiltonian Hcore  -----------'
                print self.H_core
        self.ERI = numpy.zeros( (self.N,self.N,self.N,self.N))

        if (self.do_print ==1):
                print '----------- Electronic Repulsion Integrals  -----------'
        for i in xrange(self.N):
            for j in xrange(i + 1):
                for k in xrange(self.N):
                    for l in xrange(k + 1):
                        t = self.electron_rep(self.aof[i],self.aof[j],self.aof[k],self.aof[l])
                        self.ERI[i][j][k][l] = t
                        self.ERI[j][i][k][l] = t
                        self.ERI[i][j][l][k] = t
                        self.ERI[j][i][l][k] = t
                        if (self.do_print ==1):
                                print '(',i+1,' ',j+1,'|',k+1,' ',l+1,') = ', self.ERI[i][j][k][l]

        self.F = numpy.zeros((self.N,self.N))
        #for i in xrange(self.N):
        #    for j in xrange(self.N):
#                self.F[i][j] = self.H_core[i][j]
        self.F = self.H_core.copy()
        self.G = numpy.zeros((self.N,self.N))
        self.P= numpy.zeros((self.N,self.N))
        self.C = numpy.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                self.C[i][j] = 1

    def set_equal(self,A,B,N): # must find a work around
      for i in range(N):
        for j in range(N):
          A[i][j]=B[i][j]
      # make a workaround for set_equal

    def iterate(self,F, H_core, G, P,C,S_arb,ERI,aof): # F, H_core, G, P,C,ERI,aof

        N = self.N
        self.set_equal(F,H_core+G,N)
        E0 = self.electronic_energy(P,H_core,F)

        X = real(linalg.sqrtm(linalg.inv(S_arb)))
        Fd =  dot(dot(transpose(X),F),X)
        self.set_equal(C ,dot(X,linalg.eig(Fd )[1]),N) # I think this needs fixing
        # Note, the number of electrons is currently specified as a constant
        # below - this needs to be treated properly to perform
        # correct calculations on different molecules
        self.set_equal(P , self.density(aof,C,1),N) #

        self.set_equal(G , (self.calc_G(P,ERI)),N)
        return E0
