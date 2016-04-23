import math
class Primative:
    def __init__(self,coef,exp,l,m,n,centre): # XXX check if we can have <0 l,m or n
        self.coef = coef
        self.exp = exp
        self.l = l
        self.m = m
        self.n = n
        self.centre = centre
        self.norm =1

    def dist(self,a,b):
        return (math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2))

    def eval(self,x,y,z):
        d =((self.centre[0]-x)**2 + (self.centre[1]-y)**2 + (self.centre[2]-z)**2)**0.5
        t = self.coef * math.e**(-1*self.exp*self.dist(self.centre,[x,y,z])*d) * x**self.l * y**self.m *z**self.n
        return t

    def pr(self):
        return self.coef,'e^(-',self.exp,'r^2) x^',self.l,'y^',self.m,'z^',self.n

    def prim_copy(self,l,m,n):
        if (l == m == n == 0):
            return self;
        temp = Primative(self.coef,self.exp,self.l+l,self.m+m,self.n+n,self.centre)
        return temp

class ContrPrimative:
    def __init__(self):
        self.type = -1 # not set , 0 = S , 1 = SP ...
        self.n = 0
        self.primatives = []
    def add(self,prim): # fix add all the time later
        if (prim==[]):
            return
        self.n+=1
        self.primatives.append(prim)
        self.centre = prim.centre # HACK HACK HACK HACK

    def prim(self,n):
        return self.primatives[n]

