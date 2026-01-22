import math
import numpy as np
import random
import numbers
from mpmath import mp

import numpy as np

count = 10000
    

mp.dps = 20
mp.prec = 150

class myarray:
    def __init__(self,a):
        self.data = []
        for row in a:
            if not hasattr(row, '__len__'):
                self.data.append(row)
            else:
                self.data.append(myarray(row))
    
    def __getitem__(self, i):
        return self.data[i]
        
    def __setitem__(self, i, x):
        self.data[i] = x
    
    def __len__(self):
        return len(self.data)

    def __add__(self, a):
        x = myarray(self.data)
        for i in range(len(x)):
            x[i] += a[i]
        return x
        
    def __sub__(self, a):
        x = myarray(self.data)
        for i in range(len(x)):
            x[i] -= a[i]
        return x
        
    def __mul__(self, a):
        x = myarray(self.data)
        for i in range(len(x)):
            if hasattr(a,'__len__'):
                x[i] *= a[i]
            else:
                x[i] *= a
        return x
        
    def __truediv__(self, a):
        x = myarray(self.data)
        for i in range(len(self.data)):
            if hasattr(a,'__len__'):
                x[i] /= a[i]
            else:
                x[i] /= a
        return x
        
    def dot(self, a): #a must be a vector
        if hasattr(self.data[0],'__len__'):
            d = []
            for i in range(len(self.data)):
                d.append(self.data[i].dot(a))
            return d
        else:
            d = 0
            for i in range(len(self.data)):
                d += self.data[i] * a[i]
            return d
            
    def __str__(self):
        out = "[ "
        for b in self:
            out += str(b) + " "
        return out+"]"

array = myarray

def norm(v):
    s2 = 0
    for a in v:
        s2 += a*a
    return mth.sqrt(s2)
    
def cross(a,b):
    return array( ( a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] ) )

# author: Alexander Pruss
# license: MIT

#
# No attempts at efficiency are made in this code. The main point is to ensure readability
# for the sake of verification of formulae.
#

# inputs: H = m23-m12 = altitude of landmark 2 above landmark 1
#         d = horizontal distance between landmarks
#         rho[0], rho[1] = observed tilt angles
#         beta           = observed horizontal landmark angle from landmark 1 to 2
#         eps            = tolerance for comparison with zero
#
# outputs: array of:
#     (h1,d1,h2,d2)
#         where: hi is altitude of camera above landmark i
#                di is horizontal distance from landmark i to camera
#
# algorithm: follow the steps of the proof of Theorem 1 in
#            https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5228196

def cot(x):
    return 1/mth.tan(x)

def solve(H, d, rho, beta, eps=1e-12, assumeOnPlane=False):
    scaledEps = eps * mth.hypot(H, d)
    
    if abs(rho[0]) <= eps or abs(rho[0]-mth.pi) <= eps:
        assert not( abs(rho[1]) < eps or abs(rho[1]-mth.pi) < eps )
        d1 = 0.
        d2 = d
        h2 = -d * cot(rho[1])
        h1 = h2 + H
        return [(h1,d1,h2,d2),]
    elif abs(rho[1]) <= eps or abs(rho[1]-mth.pi) <= eps:
        d2 = 0.
        d1 = d
        h1 = -d * cot(rho[0])
        h2 = h1 - H
        return [(h1,d1,h2,d2),]

    if abs(rho[0]-mth.pi/2) <= eps:
        assert not abs(rho[1]-mth.pi/2) <= eps
        i = 1
        j = 0
    else:
        i = 0
        j = 1
    
    delta = H if i == 0 else -H
    
    a = 1 - 2 * mth.cos(beta) * cot(rho[j]) * mth.tan(rho[i]) + (cot(rho[j]) * mth.tan(rho[i])) ** 2
    b = 2 * delta * mth.tan(rho[i]) * ( mth.cos(beta) - cot(rho[j]) * mth.tan(rho[i]) )
    c = (delta * mth.tan(rho[i]))**2 - d**2
    
    assert a > 0
    
    def finishSolution(dj):
        if dj < 0:
            return None
        di = (dj * cot(rho[j]) - delta) * mth.tan(rho[1-j])
        if di < 0:
            return None
        hj = -dj * cot(rho[j])
        hi = hj + delta
        if i == 0:
            return (hi,di,hj,dj)
        else:
            return (hj,dj,hi,di)
    
    discriminant = b**2 - 4 * a * c
    
    if assumeOnPlane or abs(discriminant) <= scaledEps ** 2:
        dj = -b / (2 * a)
        solution = finishSolution(dj)
        if not solution:
            return []
        else:
            return [solution,]
            
    assert discriminant > 0, "discriminant is "+str(discriminant)
        
    solutions = []
    
    for s in [-1,1]:
        dj = ( -b + s*mth.sqrt(discriminant) ) / (2 * a)
        if dj > 0:
            solution = finishSolution(dj)
            if solution:
                solutions.append(solution)
            
    return solutions
 
# 
# find points at intersection of circles centered on A and B with radii rA and RB
#
def intersectionOfCircles(A,rA,B,rB,scaledEps):
    diff = B-A
    d = norm(diff)
    if abs(rA-rB-d) <= scaledEps:
        # one solution
        # return appropriate affine combination of points A and B
        return [ (A * rB + B * rA) / (rA+rB), ]
        
    # rotate problem so A lies at origin and B at (d,0)
    # to go from there to original problem, rotate by matrix and add A
    matrix = array( [ [diff[0],-diff[1]], [diff[1],diff[0]] ] ) / d
    # cf. https://mathworld.wolfram.com/Circle-CircleIntersection.html
    x = (rA**2 - rB**2)/(2*d) + 0.5*d
    
    assert(rA**2>=x**2)

    out = []
    for s in (-1,1):
        xy = array([ x, s*mth.sqrt(rA**2-x**2) ])
        out.append( A + array( matrix.dot( xy ) ) )
    return out
    
# inputs: 
#      solution:       (h1,d1,h2,d2)
#      landmarks[i]:   position (x,y,z) of landmark i (z=altitude)
#      beta:           observed horizontal landmark angle from landmark 1 to 2
# output: array((x,y,z))
#
# algorithm: follow ideas in proof of Lemma 2

def getCameraPosition(solution, beta, landmarks, eps=1e-15):
    h1,d1,h2,d2 = solution
    m = array([ landmarks[0], landmarks[1] ])
    scale = norm( m[1]-m[0] )
    scaledEps = scale * eps
    if d1 <= scaledEps:
        return m[0] + [0,0,h1]
    if d2 <= scaledEps:
        return m[1] + [0,0,h2]
    assert beta != 0
    A = array((m[0][0], m[0][1]))
    B = array((m[1][0], m[1][1]))
    cc = intersectionOfCircles(A,d1,B,d2,scaledEps)
    for C in cc:
        # the sign of s is the same as that of the signed angle ACB
        AC = A-C
        BC = B-C
        s = cross( (AC[0], AC[1], 0), (BC[0], BC[1], 0))[2] 
        if len(cc) == 1 or (s>0) == (beta>0):
            # camera (x,y)=C and camera has altitude m[0]+[0,0,h1]
            return array((C[0],C[1],m[0][2]+h1))
    assert False    
    
def det3x3(a):
    assert(len(a) == 3)
    assert(len(a[0]) == 3)
    return (a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1]) 
          - a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
          + a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]))          
          
#
# inputs: 
#      camera:       position (x,y,z) of camera (z=altitude)
#      landmarks[i]: position (x,y,z) of landmark i (z=altitude)
#      eps:          tolerance for comparisons with zero
# outputs: integer
# algorithm: use characterization in Theorem 1
# assuming landmarks are not at the same point
#

def howManySolutions(camera, landmarks, eps=1e-15):
    # for convenience move camera to origin
    l = array( [ landmarks[0]-camera, landmarks[1]-camera ] )
    scale = norm(l[1]-l[0])
    assert scale != 0
    scaledEps = scale * eps
    # type 1 degeneracy: camera and landmarks are colinear
    if norm(cross( l[0], l[1] )) <= scaledEps ** 2:
        return mth.inf
    # type 2 degeneracy: landmarks are on same vertical line, but not type 1
    d = mth.hypot( l[0][0]-l[1][0], l[0][1]-l[1][1] )
    if d <= scaledEps:
        return mth.inf
    # type 3 degeneracy: landmarks and camera on same horizontal plane, but not type 1
    if abs(l[0][2]) <= scaledEps and abs(l[1][2]) <= scaledEps:
        return mth.inf
    # if max(|s1|,|s2|) >= |s|, one solution
    d1 = mth.hypot( l[0][0], l[0][1])
    d2 = mth.hypot( l[1][0], l[1][1])    
    if d1 <= scaledEps or d2 <= scaledEps:
        # si would be infinite if di were zero
        return 1
    # slope from landmark 1 to landmark 2
    s = (l[1][2] - l[0][2]) / d
    # slopes from cameras to landmarks
    s1 = l[0][2] / d1
    s2 = l[1][2] / d2
    if max(abs(s1),abs(s2)) >= abs(s):
        return 1

    if abs(s) > eps:
        # if camera lies on plane defined by landmarks and horizontal line through
        # landmark 1 and perpendicular to line through landmarks, have 1 solution
        horizontalVector = cross(l[1]-l[0], array([0,0,1]))
        p3 = l[0] + horizontalVector * scale
        # check if origin lies on plane defined by l0, l1, p3
        # i.e., check if the three vectors are linearly dependent
        d = det( array( (p3,l[0],l[1]) ) )
        if abs(d) <= scaledEps**3:
            return 1
    
    return 2
    
#
# inputs:
#      camera:       position (x,y,z) of camera (z=altitude)
#      landmarks[i]: position (x,y,z) of landmark i (z=altitude)
# output:
#      (rho1,rho2),beta
#
def getObservations(camera, landmarks):
    # move camera to origin for simplicity
    m = ( landmarks[0] - camera, landmarks[1] - camera )
    v1 = m[0] / norm(m[0])
    v2 = m[1] / norm(m[1])
    rho = (mth.acos(v1[2]),mth.acos(v2[2])) 
    q1 = array((v1[0],v1[1],0.))
    q2 = array((v2[0],v2[1],0.))
    beta = mth.atan2( cross( q1,q2 )[2], q1.dot(q2) )
    return rho,beta
    
#
# inputs: 
#   camera:     ground truth position (x,y,z) with z vertical
#   landmarks:  array of two (x,y,z) positions
#
# outputs:
#   None:       if degenerate
#  or:
#   errNumber,errPosition,errObservation:
#       errNumber:      actual number of solutions minus expected number 
#       errPosition:    closest distance from a solution's camera position to ground truth
#       errObservation: furthest distance from the simulated observations of a solution to the 
#                       simulated observations of the ground truth; Euclidean distance on the
#                       triple (rho[0],rho[1],beta)
  
def computeErrors(camera, landmarks, verbose=False, assumeOnPlane=False):
    H = landmarks[1][2]-landmarks[0][2]
    d = mth.hypot( landmarks[1][0]-landmarks[0][0], landmarks[1][1]-landmarks[0][1] )
    if verbose:
        gd1 = mth.hypot( camera[0]-landmarks[0][0], camera[1]-landmarks[0][1] )
        gd2 = mth.hypot( camera[0]-landmarks[1][0], camera[1]-landmarks[1][1] )
        gh1 = camera[2]-landmarks[0][2]
        gh2 = camera[2]-landmarks[1][2]
        print("ground truth h1,d1,h2,d2: ",gh1,gd1,gh2,gd2)
    rho,beta = getObservations(camera, landmarks)
    expectedSolutions = howManySolutions(camera,landmarks) 
    if expectedSolutions == mth.inf:
        return None
    solutions = solve(H,d,rho,beta,assumeOnPlane=assumeOnPlane)
    errPosition = mth.inf
    errObservation = 0
    if not solutions:
        return 0,mth.inf,0
    else:
        for solution in solutions:
            solvedC = getCameraPosition(solution, beta, landmarks)
            if verbose:
                print(solution)
                print(solvedC)
            e = norm(solvedC-camera)
            errPosition = min(e, errPosition)
            solvedRho,solvedBeta = getObservations(solvedC, landmarks)
            e = norm(array((solvedRho[0],solvedRho[1],solvedBeta))-(rho[0],rho[1],beta))
            errObservation = max(e, errObservation)
        return len(solutions)-expectedSolutions,errPosition,errObservation
        
def generateFullyRandom():
    landmarks = ( array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ), 
                    array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ) )
    camera = array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) )
    return camera,landmarks

def generateCameraOnPlane():
    landmarks = ( array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ), 
                    array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ) )
    horiz = cross(landmarks[1]-landmarks[0], array([0,0,1]))
    u = random.uniform(-5,5)
    v = random.uniform(-5,5)
    return landmarks[0] + (landmarks[1]-landmarks[0])*u + horiz*v,landmarks
    
def generateCameraOnPlane_mp():
    def r(x):
        return mp.mpmathify(random.uniform(-x,x))
    landmarks = ( array( (r(20),r(20),r(20)) ), 
                    array( (r(20),r(20),r(20)) ) )
    horiz = cross(landmarks[1]-landmarks[0], array([0,0,1]))
    u = r(5)
    v = r(5)
    return landmarks[0] + (landmarks[1]-landmarks[0])*u + horiz*v,landmarks
    
def generateWithCohorizonatalLandmarks():
    camera, landmarks = generateFullyRandom()
    return camera, (landmarks[0],(landmarks[1][0],landmarks[1][2],landmarks[0][2]))
    
generator = generateWithCohorizonatalLandmarks
mth = math

if __name__ == '__main__':
    landmarks = array( ( (0,10,0), (4,10,3) ) )
    random.seed(1)
    errNumber = 0
    errPosition = 0
    errObservation = 0
    totalNumberError = 0
    totalPositionError = 0
    totalObservationError = 0
    
    for i in range(count):
        camera,landmarks = generator()
        try:
            n,p,o = computeErrors(camera,landmarks,assumeOnPlane=(generator==generateCameraOnPlane))
        except:
            p = 1
        if p > 1e-5 and mth != mp:
            mth = mp
            mpcamera = array([mp.mpmathify(c) for(c) in camera])
            mplandmarks = array([array([mp.mpmathify(l) for(l) in landmark]) for landmark in landmarks])
            n,p,o = computeErrors(mpcamera,mplandmarks,verbose=False,assumeOnPlane=(generator==generateCameraOnPlane))
            #print("mperrors ",n,p,o)
            mth = math
        errNumber = max(n,errNumber)
        totalNumberError += n
        errPosition = max(p,errPosition)
        totalPositionError += p
        errObservation = max(o,errObservation)
        totalObservationError += o

    print("n=",count)
    print("Worst errors: %d %g %g" % (errNumber,errPosition,errObservation))
    print("Average errors: %d %g %g" % (totalNumberError/count,totalPositionError/count,totalObservationError/count))
