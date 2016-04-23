import math

def my_dot(a,b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def scal_dot(a,b):
        return (a[0]*b, a[1]*b , a[2]*b)

def vec_add(a,b):
        return [a[0] + b[0],a[1]+b[1],a[2]+b[2]]

def vec_minus(a,b):
        return [a[0] - b[0],a[1]-b[1],a[2]-b[2]]

def vec_cross(a,b):
        return [a[1]*b[2] -b[1]*a[2],-1*(a[0]*b[2]-a[2]*b[0]),a[0]*b[1]-a[1]*b[0]]

def vec_dot(a,b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def normalise(a):
        k = math.sqrt((a[0])**2.0 + (a[1])**2.0 + (a[2])**2.0)
        return [a[0]/k , a[1]/k , a[2]/k]

def vec_times(a,k):
        return [a[0]*k,a[1]*k,a[2]*k]

def dist(a,b):
    return (math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2))

