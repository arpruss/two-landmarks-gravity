import math
import random
import numbers
import sys
from mpmath import mp


# author: Alexander Pruss
# license: MIT

#
# No attempts at efficiency are made in this code. The main point is to ensure readability
# for the sake of verification of formulae.
#


mth = math

count = int(1e7)
seed = 1001
eps = 1e-14
minimalLandmarkDistance = 0
cameraPositionEps = 1e-15
switchToMP = 1e-5

mp.dps = 15
mp.prec = 150

class array:
    def __init__(self,a):
        self.data = []
        for row in a:
            if not hasattr(row, '__len__'):
                self.data.append(row)
            else:
                self.data.append(array(row))
    
    def __getitem__(self, i):
        return self.data[i]
        
    def __setitem__(self, i, x):
        self.data[i] = x
    
    def __len__(self):
        return len(self.data)

    def __add__(self, a):
        x = array(self.data)
        for i in range(len(x)):
            x[i] += a[i]
        return x
        
    def __sub__(self, a):
        x = array(self.data)
        for i in range(len(x)):
            x[i] -= a[i]
        return x
        
    def __mul__(self, a):
        x = array(self.data)
        for i in range(len(x)):
            if hasattr(a,'__len__'):
                x[i] *= a[i]
            else:
                x[i] *= a
        return x
        
    def __truediv__(self, a):
        x = array(self.data)
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
            
    def __repr__(self):
        return str(self)
        
    def __str__(self):
        out = "[ "
        for b in self:
            out += str(b) + " "
        return out+"]"
        
    def mpf(self):
        out = []
        for x in self:
            if hasattr(x,'__len__'):
                out.append(x.mpf())
            else:
                out.append(mp.mpmathify(x))
        return array(out)
        
    def normSquared(self):
        s2 = 0
        for x in self:
            if hasattr(x,'__len__'):
                s2 += x.normSquared()
            else:
                s2 += x*x
        return s2
        
    def norm(self):
        return mth.sqrt(self.normSquared())
        
    def normalized(self):
        return array(self / self.norm())
        
    def cross(self,b):
        return array( ( self[1]*b[2]-self[2]*b[1], self[2]*b[0]-self[0]*b[2], self[0]*b[1]-self[1]*b[0] ) )

# inputs: H = m23-m12 = altitude of landmark 2 above landmark 1
#         d = horizontal distance between landmarks
#         rho[0], rho[1] = observed tilt angles
#         beta           = observed horizontal landmark angle from landmark 1 to 2
#         assumeOnPlane  = assume that the camera is on the "critical plane", defined by the landmarks
#                          and a horizontal line through one of them and perpendicular to the other
#
#     (h1,d1,h2,d2)
#         where: hi is altitude of camera above landmark i
#                di is horizontal distance from landmark i to camera
#
# algorithm: follow the steps of the proof of Theorem 1 in
#            https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5228196

def cot(x):
    return 1/mth.tan(x)

def solve(H, d, rho, beta, assumeOnPlane=False, verbose=False):
    if abs(rho[0]) <= eps or abs(rho[0]-mth.pi) <= eps:
        assert not( abs(rho[1]) < eps or abs(rho[1]-mth.pi) < eps )
        if verbose:
            print("solved d1=0")
        d1 = 0.
        d2 = d
        h2 = -d * cot(rho[1])
        h1 = h2 + H
        return [(h1,d1,h2,d2),]
    elif abs(rho[1]) <= eps or abs(rho[1]-mth.pi) <= eps:
        if verbose:
            print("solved d2=0")
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
    
    if assumeOnPlane or abs(discriminant) <= eps:
        dj = -b / (2 * a)
        if verbose:
            print("solved, zeroish disriminant, dj=",dj,a,c)
        solution = finishSolution(dj)
        if not solution:
            return []
        else:
            return [solution,]
            
    assert discriminant > 0, "discriminant "+str(discriminant)
        
    solutions = []
    
    for s in [-1,1]:
        dj = ( -b + s*mth.sqrt(discriminant) ) / (2 * a)
        if verbose:
            print("dj is",dj)
        if dj > 0:
            solution = finishSolution(dj)
            if solution:
                solutions.append(solution)
            
    return solutions
 
# 
# find points at intersection of circles centered on A and B with radii rA and RB
#
def intersectionOfCircles(A,rA,B,rB):
    diff = B-A
    d = diff.norm()
    
    # one solution
    singleIntersection = abs(rA-rB-d) <= eps
        
    # rotate problem so A lies at origin and B at (d,0)
    # to go from there to original problem, rotate by matrix and add A
    matrix = array( [ [diff[0],-diff[1]], [diff[1],diff[0]] ] ) / d
    if singleIntersection:
        xy = array([ rA, 0. ])
        return (  A + array( matrix.dot( xy ) ), )
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

def getCameraPosition(solution, beta, landmarks, verbose=False):
    h1,d1,h2,d2 = solution
    m = array([ landmarks[0], landmarks[1] ])
    if d1 <= cameraPositionEps:
        if verbose:
            print("small d1",d1)
        return m[0] + [0,0,h1]
    if d2 <= cameraPositionEps:
        if verbose:
            print("small d2",d2)
        return m[1] + [0,0,h2]
    assert beta != 0
    A = array((m[0][0], m[0][1]))
    B = array((m[1][0], m[1][1]))
    cc = intersectionOfCircles(A,d1,B,d2)
    if verbose:
        print("d=",(B-A).norm())
        print("circle intersections",cc)
    for C in cc:
        # the sign of s is the same as that of the signed angle ACB
        AC = A-C
        BC = B-C
        s = array((AC[0], AC[1], 0)). cross( (BC[0], BC[1], 0) )[2] 
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
#      camera:        position (x,y,z) of camera (z=altitude)
#      landmarks[i]:  position (x,y,z) of landmark i (z=altitude)
#      assumeOnPlane: assume that the camera is on the "critical plane"
#

# outputs: integer
# algorithm: use characterization in Theorem 1
# assuming landmarks are not at the same point
#

def howManySolutions(camera, landmarks, assumeOnPlane=False):
    # for convenience move camera to origin
    l = array( [ landmarks[0]-camera, landmarks[1]-camera ] )
    # type 1 degeneracy: camera and landmarks are colinear
    if l[0].cross(l[1]).norm() <= eps:
        return mth.inf
    # type 2 degeneracy: landmarks are on same vertical line, but not type 1
    d = mth.hypot( l[0][0]-l[1][0], l[0][1]-l[1][1] )
    if d <= eps:
        return mth.inf
    # type 3 degeneracy: landmarks and camera on same horizontal plane, but not type 1
    if abs(l[0][2]) <= eps and abs(l[1][2]) <= eps:
        return mth.inf
    # if max(|s1|,|s2|) >= |s|, one solution
    d1 = mth.hypot( l[0][0], l[0][1])
    d2 = mth.hypot( l[1][0], l[1][1])    
    if d1 <= eps or d2 <= eps:
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
        horizontalVector = (l[1]-l[0]).cross((0,0,1))
        p3 = l[0] + horizontalVector
        # check if origin lies on plane defined by l0, l1, p3
        # i.e., check if the three vectors are linearly dependent
        d = det3x3( (p3,l[0],l[1]) )
        if assumeOnPlane or abs(d) <= eps:
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
    v1 = m[0].normalized()
    v2 = m[1].normalized()
    rho = (mth.acos(v1[2]),mth.acos(v2[2])) 
    q1 = array((v1[0],v1[1],0.))
    q2 = array((v2[0],v2[1],0.))
    beta = mth.atan2( q1.cross( q2 )[2], q1.dot(q2) )
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
  
def compute(camera, landmarks, verbose=False, assumeOnPlaneForSolutions=False, assumeOnPlaneForPredictions=False,msg=""):
    H = landmarks[1][2]-landmarks[0][2]
    d = mth.hypot( landmarks[1][0]-landmarks[0][0], landmarks[1][1]-landmarks[0][1] )
    expectedSolutions = howManySolutions(camera,landmarks, assumeOnPlane=assumeOnPlaneForPredictions) 
    rho,beta = getObservations(camera, landmarks)
    if expectedSolutions == mth.inf:
        return None
    solutions = solve(H,d,rho,beta,assumeOnPlane=assumeOnPlaneForSolutions,verbose=verbose)
    errPosition = mth.inf
    errObservation = 0
    if verbose and len(solutions) != expectedSolutions:
        print("Actual solutions: %d\nExpected solutions: %d" % (len(solutions),expectedSolutions))
    if not solutions:
        return 0,mth.inf,0,0
    else:
        for solution in solutions:
            solvedC = getCameraPosition(solution, beta, landmarks)
            if verbose or solvedC[0]==landmarks[0][0]:
                print("message:",msg)
                print("camera",camera)
                solve(H,d,rho,beta,assumeOnPlane=assumeOnPlaneForSolutions,verbose=True)
                print("landmarks",landmarks)
                print("solutions",solutions)
                gd1 = mth.hypot( camera[0]-landmarks[0][0], camera[1]-landmarks[0][1] )
                gd2 = mth.hypot( camera[0]-landmarks[1][0], camera[1]-landmarks[1][1] )
                gh1 = camera[2]-landmarks[0][2]
                gh2 = camera[2]-landmarks[1][2]
                print("ground truth h1,d1,h2,d2: ",gh1,gd1,gh2,gd2)
                print("solution",solution)
                print("solved camera",solvedC)
                gcc = getCameraPosition([gh1,gd1,gh2,gd2],beta,landmarks,verbose=True)
                print("camera according to ground truth solution",gcc)
                print("delta",solvedC-gcc)
            e = (solvedC-camera).norm()
            errPosition = min(e, errPosition)
            # Problem here!
            solvedRho,solvedBeta = getObservations(solvedC, landmarks)
            e = (array((solvedRho[0],solvedRho[1],solvedBeta))-(rho[0],rho[1],beta)).norm()
            errObservation = max(e, errObservation)
        return len(solutions)-expectedSolutions,errPosition,errObservation,len(solutions)
        
def generateFullyRandom():
    landmarks = array( ( ( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ), 
                       ( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ) ) )
    camera = array( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) )
    return camera,landmarks

def generateCameraOnPlane():
    landmarks = array( ( ( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ), 
                    ( (random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20)) ) ) )
    joining = (landmarks[1]-landmarks[0]).normalized()
    horiz = joining.cross((0,0,1))
    u = random.uniform(-20,20)
    v = random.uniform(-20,20)
    return landmarks[0] + joining*u + horiz*v,landmarks
    
def generateCameraOnPlane_mp():
    def r(x):
        return mp.mpmathify(random.uniform(-x,x))
    landmarks = array( ( (r(20),r(20),r(20)) ), 
                    ( (r(20),r(20),r(20)) ) )
    joining = (landmarks[1]-landmarks[0]).normalized()
    horiz = joining.cross((0,0,1))
    u = r(20)
    v = r(20)
    return landmarks[0] + joining*u + horiz*v,landmarks
    
def generateWithCohorizontalLandmarks():
    camera, landmarks = generateFullyRandom()
    return camera, array((landmarks[0],((landmarks[1][0],landmarks[1][2],landmarks[0][2]))))

def setMode(mode):  
    global generator,assumeOnPlaneForSolutions,assumeOnPlaneForPredictions
    if mode.lower() == "fullyrandom":
        generator = generateFullyRandom
        assumeOnPlaneForSolutions = False
        assumeOnPlaneForPredictions = False
    elif mode.lower() == "onplane":
        generator = generateCameraOnPlane
        assumeOnPlaneForSolutions = False
        assumeOnPlaneForPredictions = True
    elif mode.lower() == "cohorizontal":
        generator = generateWithCohorizontalLandmarks
        assumeOnPlaneForSolutions = False
        assumeOnPlaneForPredictions = False
    else:
        assert False, "Unknown mode"

if __name__ == '__main__':
    if len(sys.argv)>1:
        mode = sys.argv[1]
    else:
        mode = "fullyRandom"
    setMode(mode)
    if len(sys.argv)>2:
        seed = int(sys.argv[2])
    print("mode=",mode)
    print("seed=",seed)
    print("n=",count)

    random.seed(seed)
    errNumber = 0
    errPosition = 0
    errObservation = 0
    totalNumberError = 0
    totalPositionError = 0
    totalObservationError = 0
    totalOneSolution = 0
    switchedToMP = 0
    
    for i in range(count):
        while True:
            camera,landmarks = generator()
            if (landmarks[0]-landmarks[1]).norm() < minimalLandmarkDistance:
                continue
            s = howManySolutions(camera,landmarks, assumeOnPlane=assumeOnPlaneForPredictions) 
            if s != mth.inf: # degenerate case
                break
        try:
            n = 0
            n,p,o,sols = compute(camera,landmarks,assumeOnPlaneForPredictions=assumeOnPlaneForPredictions,assumeOnPlaneForSolutions=assumeOnPlaneForSolutions,msg=str(i))
        except AssertionError as e:
            p = math.inf
        if n > 0:
            p = math.inf
        if p > switchToMP and mth != mp:
            mth = mp
            mpcamera = camera.mpf()
            mplandmarks = landmarks.mpf()
            n,p,o,sols = compute(mpcamera,mplandmarks,verbose=False,assumeOnPlaneForPredictions=assumeOnPlaneForPredictions,assumeOnPlaneForSolutions=assumeOnPlaneForSolutions)
            if p == math.inf:
                print("really bad",p)
                print("camera",mpcamera)
                print("landmarks",mplandmarks)
                compute(mpcamera,mplandmarks,verbose=True,assumeOnPlaneForPredictions=assumeOnPlaneForPredictions,assumeOnPlaneForSolutions=assumeOnPlaneForSolutions)
            switchedToMP += 1
            mth = math
        errNumber = max(errNumber,abs(n))
        totalNumberError += abs(n)
        errPosition = max(p,errPosition)
        totalPositionError += p
        errObservation = max(o,errObservation)
        totalObservationError += o
        if sols == 1:
            totalOneSolution += 1

    print("one solution ratio=",totalOneSolution/count)
    print("Worst errors: %d %g %g" % (errNumber,errPosition,errObservation))
    print("Average errors: %g %g %g" % (totalNumberError/count,totalPositionError/count,totalObservationError/count))
    print("Ratio switched to MP: %g" % (switchedToMP/count))
     