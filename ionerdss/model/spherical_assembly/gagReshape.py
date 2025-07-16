
# %%
import numpy as np  # matrix, array, math calculation
import sys          # 
import os          
import pandas as pd # data frame, input and output

# %%
###########################################################################################
###########################################################################################
# function definition
def calculateAngles(c1, c2, p1, p2, n1=None, n2=None):
    """
    Determine the angles of the reaction (theta1, theta2, phi1, phi2, omega) given the coordinates of the two Center of Mass (c1 and c2) and two reaction sites (p1 and p2), and two norm vectors (n1 and n2).

    Parameters
    ----------
    c1 : numpy.array
        Center of Mass vector for the first molecule.
    c2 : numpy.array
        Center of Mass vector for the second molecule.
    p1 : numpy.array
        Reaction site vector for the first molecule.
    p2 : numpy.array
        Reaction site vector for the second molecule.
    n1 : numpy.array
        Norm vector for the first molecule.
    n2 : numpy.array
        Norm vector for the second molecule.

    Returns
    -------
    tuple
        The tuple (theta1, theta2, phi1, phi2, omega), where theta1, theta2, phi1, phi2, omega are the angles in radians.
    """
    if n1 is None:
        n1 = np.array([0, 0, 1])
    if n2 is None:
        n2 = np.array([0, 0, 1])
    v1 = p1 - c1
    v2 = p2 - c2
    sigma1 = p1 - p2
    sigma2 = -sigma1

    theta1 = np.arccos(
        np.dot(v1, sigma1) / (np.linalg.norm(v1) * np.linalg.norm(sigma1))
    )
    theta2 = np.arccos(
        np.dot(v2, sigma2) / (np.linalg.norm(v2) * np.linalg.norm(sigma2))
    )

    t1 = np.cross(v1, sigma1)
    t2 = np.cross(v1, n1)
    norm_t1 = t1 / np.linalg.norm(t1)
    norm_t2 = t2 / np.linalg.norm(t2)
    phi1 = np.arccos(np.dot(norm_t1, norm_t2))

    # the sign of phi1 is determined by the direction of t2 relative to the right-hand rule of cross product of v1 and t1
    if np.dot(np.cross(v1, t1), t2) > 0:
        phi1 = -phi1

    t1 = np.cross(v2, sigma2)
    t2 = np.cross(v2, n2)
    norm_t1 = t1 / np.linalg.norm(t1)
    norm_t2 = t2 / np.linalg.norm(t2)
    phi2 = np.arccos(np.dot(norm_t1, norm_t2))

    # the sign of phi2 is determined by the direction of t2 relative to the right-hand rule of cross product of v2 and t1
    if np.dot(np.cross(v2, t1), t2) > 0:
        phi2 = -phi2

    if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0) and not np.isclose(
        np.linalg.norm(np.cross(v2, sigma2)), 0
    ):
        t1 = np.cross(sigma1, v1)
        t2 = np.cross(sigma1, v2)
    else:
        t1 = np.cross(sigma1, n1)
        t2 = np.cross(sigma1, n2)

    omega = np.arccos(np.dot(t1, t2) / (np.linalg.norm(t1) * np.linalg.norm(t2)))

    # the sign of omega is determined by the direction of t2 relative to the right-hand rule of cross product of sigma1 and t1
    if np.dot(np.cross(sigma1, t1), t2) > 0:
        omega = -omega

    print("n1: ", n1)
    print("n2: ", n2)
    print("theta1, theta2, phi1, phi2, omega: ", theta1, theta2, phi1, phi2, omega)
    print("sigma: ", np.linalg.norm(sigma1))

    return theta1, theta2, phi1, phi2, omega

def calculate_rmsd(centers, xyzR):  # rmsd: to describe how gags' centers are biased from a perfect sphere surface
    x0 = xyzR[0]
    y0 = xyzR[1]
    z0 = xyzR[2]
    r0 = xyzR[3]
    s = 0.0
    for i in range(0,centers.shape[0]):
        xi = centers[i,0]
        yi = centers[i,1]
        zi = centers[i,2]
        ri = np.sqrt( (xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2 )
        s = s + (ri-r0)**2
    return s

def calculate_gradient(centers,xyzR): 
    x0 = xyzR[0]
    y0 = xyzR[1]
    z0 = xyzR[2]
    r0 = xyzR[3]
    dsdx = 0
    dsdy = 0
    dsdz = 0
    dsdr = 0
    for i in range(0,centers.shape[0]):
        xi = centers[i,0]
        yi = centers[i,1]
        zi = centers[i,2]
        ri = np.sqrt( (xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2 )
        dsdx = dsdx + (-2.0/ri)*(ri-r0)*(xi-x0)
        dsdy = dsdy + (-2.0/ri)*(ri-r0)*(yi-y0)
        dsdz = dsdz + (-2.0/ri)*(ri-r0)*(zi-z0)
        dsdr = dsdr + (-2.0)*(ri-r0)
    gradient = np.array([dsdx, dsdy, dsdz, dsdr])
    return gradient

def determine_gagTemplate_structure(numGag, positionsVec):
    internalBasis = np.zeros([3,3,numGag]) # 18 gags, each gag has 3 vectors of internal basises, rowvec
    coefficients = np.zeros([5,3,numGag])  # internal coords of interfaces in the internal basis system
    # set up the internal coord system for each gag: basis vec1, vec2, vec3 
    # and then calculate the coords of each interface in this internal system: internal coords
    # then the mean value of the internal coords gives the template structure 
    for i in range (0,numGag):
        center = positionsVec[6*i,:]
        interfaces = positionsVec[1+6*i:5+1+6*i,:]
        # determine the three internal basis: vec1, vec2, vec3
        vec1 = center
        vec1 = vec1/np.linalg.norm(vec1) 
        vec2 = interfaces[0,:] - center
        vec3 = np.cross(vec1,vec2)
        vec3 = vec3/np.linalg.norm(vec3)
        vec2 = np.cross(vec3,vec1)
        vec2 = vec2/np.linalg.norm(vec2)
        internalBasis[0,:,i] = vec1
        internalBasis[1,:,i] = vec2
        internalBasis[2,:,i] = vec3
        # calculate the interal coords for the 5 interfaces
        for j in range (0,5):
            p = interfaces[j,:] - center
            A = np.array([vec1, vec2, vec3])
            coeff = np.dot(p, np.linalg.inv(A))
            if (np.linalg.norm(np.dot(coeff,A)-p) > 1e-12) : # check whether correctly calculated!
                print('Wrong calculation of coefficients\n',np.linalg.norm(np.dot(coeff,A)-p))
            coefficients[j,:,i] = coeff
    
    # regularize the gags internal coords
    coeffReg = np.zeros([5,3]) # five sites, each site has 3 coefficients of internal coords
    for i in range(0,5) :
        coeffReg[i,0] = np.mean(coefficients[i,0,:])
        coeffReg[i,1] = np.mean(coefficients[i,1,:])
        coeffReg[i,2] = np.mean(coefficients[i,2,:])

    # using the mean coefficients to calculate the structure of the first gag, and take it as the template
    chosenGagIndex = 0
    vec1 = internalBasis[0,:,chosenGagIndex]
    vec2 = internalBasis[1,:,chosenGagIndex]
    vec3 = internalBasis[2,:,chosenGagIndex]
    center1 = positionsVec[0+6*chosenGagIndex,:]
    template = np.zeros([6,3])
    template[0,:] = center1
    for i in range(0,5):
        template[i+1,:] = coeffReg[i,0] * vec1 + coeffReg[i,1] * vec2 + coeffReg[i,2] * vec3 + center1
    
    return template

def xyz_to_sphere_coordinates(position): # translate x-y-z coords to spherical coords
    x = position[0]
    y = position[1]
    z = position[2]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arccos( x/r/np.sin(theta) )
    if y < 0 :
        phi = 2.0*np.pi - phi
    spherecoordinates = [theta, phi, r]
    return spherecoordinates

def translate_gags_on_sphere(hexmer, center1, center2): # move a hexamer gags on the sphere surface
    # express the hexmer in the internal coord system, which is based on the point 'from' (center1)
    # set up the internal coordinate basis for hexmer: vec1, vec2, vec3
    vec1 = center1 / np.linalg.norm(center1) # vec1 is along the radius direction
    vec2 = center2 - center1                 # vec2 considers the translational direction
    vec3 = np.cross(vec1,vec2)               # vec3 is along the tangent line of the hexmer
    vec3 = vec3/np.linalg.norm(vec3)
    vec2 = np.cross(vec3,vec1)
    vec2 = vec2/np.linalg.norm(vec2)
    coeffs = np.zeros([hexmer.shape[0],3])
    for i in range (0,hexmer.shape[0]) :
        b = hexmer[i,:] - center1
        A = np.array([vec1, vec2, vec3])
        coeff = [0.0, 0.0, 0.0]
        if (np.linalg.norm(b) > 1e-10) :    
            coeff = np.dot(b, np.linalg.inv(A))
        coeffs[i,:] = coeff
    # move the internal coord system to the point 'to' (center2)
    # find the internal coord system based on the point 'to'
    vec1 = center2 / np.linalg.norm(center2)
    vec2 = center2 - center1 
    vec3 = np.cross(vec1,vec2)
    vec3 = vec3/np.linalg.norm(vec3)
    vec2 = np.cross(vec3,vec1)
    vec2 = vec2/np.linalg.norm(vec2)
    # rebuid the hexmer at this new point 'to'
    newhexmer = np.zeros([hexmer.shape[0],3])
    for i in range (0,hexmer.shape[0]) :
        newhexmer[i,:] = coeffs[i,0]*vec1 + coeffs[i,1]*vec2 + coeffs[i,2]*vec3 + center2
    return newhexmer
###########################################################################################
###########################################################################################

# %%

R0 = 25.0           # the target radius of the gag capsid, nm
distanceCC = 10.0   # the distance between two hexamers, center-to-center distance, nm
# read gag positions
# positions = pd.read_excel('gagspositions.txt', header=None)
positions = pd.read_csv('gagpositions.txt', header=None, sep='\s+')
positions.columns = ['x','y','z']
# change angstrom to nm
positions['x'] = positions['x'] / 10.0
positions['y'] = positions['y'] / 10.0
positions['z'] = positions['z'] / 10.0
positionsVec = positions.to_numpy()
##############################################
# find the sphere radius and the sphere center
# 18 gags, center + 5 nodes' positions, so each gag has 6 positions. I will
# add the membrane-bind and RNA-bind sites later in this code
##############################################

# first, using the centers of gags to calculate the sphere radius and sphere center 
numGag = 18
centersVec = np.zeros([18,3])
for i in range(0,numGag):
    centersVec[i] = positionsVec[6*i]

sphereXYZR = [0,0,0,70] # initial trial values for sphere x, y, z, and R respectively 
rmsdOld = calculate_rmsd(centersVec,sphereXYZR)
isForceSmallEnough = False
iteration_count = 0
while isForceSmallEnough == False :
    force = calculate_gradient(centersVec,sphereXYZR)
    stepSize = 1.0
    tempXYZR = sphereXYZR - force * stepSize
    rmsdNew = calculate_rmsd(centersVec,tempXYZR)
    while rmsdNew > rmsdOld :
        stepSize = stepSize * 0.8
        tempXYZR = sphereXYZR - force * stepSize
        rmsdNew = calculate_rmsd(centersVec,tempXYZR)
    sphereXYZR = tempXYZR
    rmsdOld = calculate_rmsd(centersVec,sphereXYZR)
    if ( np.linalg.norm(force) < 0.01 ):
        isForceSmallEnough = True
    iteration_count = iteration_count + 1

print('Sphere center position [x,y,z] and radius [R] are, respectively: \n',sphereXYZR) 
x0 = sphereXYZR[0] 
y0 = sphereXYZR[1] 
z0 = sphereXYZR[2] 
r0 = sphereXYZR[3] 

##############################################
# Second, move the center of the sphere to the axis origin, equivalently move the gags 
positionsVec[:,0] = positionsVec[:,0] - x0
positionsVec[:,1] = positionsVec[:,1] - y0
positionsVec[:,2] = positionsVec[:,2] - z0
centersVec[:,0] = centersVec[:,0] - x0
centersVec[:,1] = centersVec[:,1] - y0
centersVec[:,2] = centersVec[:,2] - z0

##############################################
# Third, move the centers of gag to the sphere surface
for i in range (0,numGag):
    center = centersVec[i,:]
    move = center/np.linalg.norm(center) * r0 - center
    centersVec[i,:] = centersVec[i,:] + move
    for j in range (0,6):
        positionsVec[j+6*i,:] = positionsVec[j+6*i,:] + move

##############################################
# Fourth, determine the template of the gag (automatically as the first gag). Other gags will be got by translation and rotation of this template gag
# gagTemplate is the positions of the gag center and five interfaces
# gagTemplateInterCoeffs is the coefficients of the gag 5 interfaces in the internal basis system
gagTemplate = determine_gagTemplate_structure(numGag, positionsVec)

##############################################
# Fiveth, adjust the hexmerGags 0, 3, 6, 9, 12, 15 
center0 = centersVec[0,:]   
center3 = centersVec[3,:]   
center6 = centersVec[6,:]   
center9 = centersVec[9,:]
center12 = centersVec[12,:]
center15 = centersVec[15,:]
hexmerCenter = (center0 + center9) / 2.0 
# the hexamerCenter is almost along the Z axis, just set as along Z axis, then it would be easier for the following rotation and translation of gags
hexmerCenter[0] = 0.0
hexmerCenter[1] = 0.0
hexmerCenter[2] = center0[2]
# set up the internal coordinate system of the first gag: 3 basis vecs: interBaseVec0, interBaseVec1, interBaseVec2
interBaseVec0 = center0 / np.linalg.norm(center0)                   # along the radius diraction
interBaseVec1 = center0 - hexmerCenter
interBaseVec2 = np.cross(interBaseVec0,interBaseVec1)               # along the tangent line of the hexamer circumference
interBaseVec2 = interBaseVec2 / np.linalg.norm(interBaseVec2) 
interBaseVec1 = np.cross(interBaseVec2,interBaseVec0)
interBaseVec1 = interBaseVec1 / np.linalg.norm(interBaseVec1)
# calculate gag1's (also gagTemplate) coordinates in its internal coordinate-system
interCoords = np.zeros([5,3]) # 5 sites, each needs 3 coordinates
for i in range (0,5) :
    p = gagTemplate[i+1,:] - gagTemplate[0,:] 
    A = np.array([interBaseVec0, interBaseVec1, interBaseVec2])
    interCoords[i,:] = np.dot(p, np.linalg.inv(A))
sphereCoords0 = xyz_to_sphere_coordinates(center0)
theta = sphereCoords0[0]
phi   = sphereCoords0[1]
r     = sphereCoords0[2]
r     = R0     # change the radius as the one inputted
theta = np.arcsin(r0*np.sin(theta)/r) # recaculate the theta angle according to the target radius of the sphere
# the postions of gags 3, 6, 9, 12, 15 are generated by rotating gag0 along the z-axis
deltaAngle = 2.0 * np.pi / 6.0
sixCenters = np.zeros([6,3]) # to store the centers of the hexamer gags
for i in range (0,6) :
    thetai = theta
    phii = phi + i*deltaAngle
    ri = r
    sixCenters[i,:] = [ri*np.sin(thetai)*np.cos(phii), ri*np.sin(thetai)*np.sin(phii), ri*np.cos(thetai)]

hexmerCenter[0] = 0  # since the radius is rebuilt by the target radius R0, then the center of the hexamer needs also updated
hexmerCenter[1] = 0
hexmerCenter[2] = sixCenters[0,2]
hexmerPositionsVec = np.zeros([6*6,3]) # to store the postions of the hexamer gags, both their center position and interface position
for i in range (0,6):
    center = sixCenters[i,:]
    vec1 = center/np.linalg.norm(center) # set up the internal coord system, the three basis vectors. The way is similar as for the template 
    vec2 = center - hexmerCenter
    vec3 = np.cross(vec1,vec2)
    vec3 = vec3/np.linalg.norm(vec3)
    vec2 = np.cross(vec3,vec1)
    vec2 = vec2/np.linalg.norm(vec2)
    interfaces = np.zeros([5,3])
    for j in range (0,5):   # using the template internal coords to determine the five interfaces of each gag
        interfaces[j,:] = interCoords[j,0]*vec1 + interCoords[j,1]*vec2 + interCoords[j,2]*vec3 + center
    hexmerPositionsVec[6*i+0,:] = center
    hexmerPositionsVec[6*i+1:6*i+6,:] = interfaces
    
# store the hexmerPositionsVec to the new positions
newPositionsVec = np.zeros([numGag*6,3])

newPositionsVec[6*0:6*0+6,:]   = hexmerPositionsVec[0:6,:]    # gag0
newPositionsVec[6*3:6*3+6,:]   = hexmerPositionsVec[6:12,:]   # gag3      
newPositionsVec[6*6:6*6+6,:]   = hexmerPositionsVec[12:18,:]  # gag6
newPositionsVec[6*9:6*9+6,:]   = hexmerPositionsVec[18:24,:]  # gag9   
newPositionsVec[6*12:6*12+6,:] = hexmerPositionsVec[24:30,:]  # gag12
newPositionsVec[6*15:6*15+6,:] = hexmerPositionsVec[30:36,:]  # gag15

##############################################
# Sixth, find the positions of gag 1,2,4,5,7,8,10,11,13,14,16,17, by moving the hexamer on the sphere surface
# determine the target position that the original hexamer will move to. I use spherical coordinates
moveVec = positionsVec[6*2,:] - positionsVec[6*12,:]
toPosition = moveVec + np.array([0,0,r])
sphereCrds = xyz_to_sphere_coordinates(toPosition) 
phi = sphereCrds[1] 
theta = distanceCC/r # the hexamer center moves 10 nm on the sphere

fromPosition = np.array([0,0,r]) # center of the original hexamer
deltaAngle = 2.0 * np.pi / 6.0
hexmer = hexmerPositionsVec

# gag 1, 2
phi = phi
fromPosition = fromPosition 
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)]) # center of the newhexamer
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*1:6+6*1,:] = newhexmer[6*3:6+6*3,:]     # gag 1 comes from the new position of hexmer gag3 
newPositionsVec[6*2:6+6*2,:] = newhexmer[6*4:6+6*4,:]     # gag 2 comes from the new position of hexmer gag4

# gag 4, 5
phi = phi + deltaAngle
fromPosition = fromPosition
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*4:6+6*4,:] = newhexmer[6*4:6+6*4,:]     # gag 4 comes from the new position of hexmer gag4 
newPositionsVec[6*5:6+6*5,:] = newhexmer[6*5:6+6*5,:]     # gag 5 comes from the new position of hexmer gag5

# gag 7, 8
phi = phi + deltaAngle
fromPosition = fromPosition
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*7:6+6*7,:] = newhexmer[6*5:6+6*5,:]     # gag 7 comes from the new position of hexmer gag5 
newPositionsVec[6*8:6+6*8,:] = newhexmer[6*0:6+6*0,:]     # gag 8 comes from the new position of hexmer gag0

# gag 10, 11
phi = phi + deltaAngle
fromPosition = fromPosition
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*10:6+6*10,:] = newhexmer[6*0:6+6*0,:]     # gag 10 comes from the new position of hexmer gag0 
newPositionsVec[6*11:6+6*11,:] = newhexmer[6*1:6+6*1,:]     # gag 11 comes from the new position of hexmer gag1

# gag 13, 14
phi = phi + deltaAngle
fromPosition = fromPosition
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*13:6+6*13,:] = newhexmer[6*1:6+6*1,:]     # gag 13 comes from the new position of hexmer gag1 
newPositionsVec[6*14:6+6*14,:] = newhexmer[6*2:6+6*2,:]     # gag 14 comes from the new position of hexmer gag2

# gag 16, 17
phi = phi + deltaAngle
fromPosition = fromPosition
toPosition = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])
newhexmer = translate_gags_on_sphere(hexmerPositionsVec, fromPosition, toPosition)
newPositionsVec[6*16:6+6*16,:] = newhexmer[6*2:6+6*2,:]     # gag 16 comes from the new position of hexmer gag2 
newPositionsVec[6*17:6+6*17,:] = newhexmer[6*3:6+6*3,:]     # gag 17 comes from the new position of hexmer gag3

##############################################
# seventh, add the membrane-bind and RNA-bind sites, then each gag has 8 points
finalPositionsVec = np.zeros([8*18,3])
for i in range(0,numGag) :
    center = newPositionsVec[6*i,:]
    surfaceSite = center + 1.0*center/np.linalg.norm(center) # the distance between surfacesite and the center is set as 1nm
    rnaSite = center - 1.0*center/np.linalg.norm(center) # the distance between rnasite and the center is set as 1nm
    finalPositionsVec[0+8*i,:] = center
    finalPositionsVec[1+8*i,:] = surfaceSite
    finalPositionsVec[2+8*i,:] = rnaSite
    finalPositionsVec[3+8*i:8+8*i,:] = newPositionsVec[1+6*i:6+6*i,:]


# %%
###########################################################################################
###########################################################################################
# output the coordinates of each gag
gagNames = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"]
for i in range (0,numGag):
    print(gagNames[i],'\n')
    positions = finalPositionsVec[8*i:8+8*i,:]
    print(positions,'\n')

# ouput the coord.txt file
coordFileName = f'coordR{int(R0)}d{int(distanceCC)}.txt'
with open(coordFileName, 'w') as f:
    for i in range(0,numGag):
        f.write(f'{gagNames[i]}\n')
        positions = finalPositionsVec[8*i:8+8*i,:]
        for j in range(0,8):
            f.write(f'{positions[j,0]:.8f} {positions[j,1]:.8f} {positions[j,2]:.8f}\n')


# %%
# we will output the relative positions of gag centers and interfaces with the center at [0,0,0]
relativePositionsGag = finalPositionsVec[0:8, :]
relativePositionsGag = relativePositionsGag - relativePositionsGag[0, :]
print("COM", relativePositionsGag[0, :])
print("mem", relativePositionsGag[1, :])
print("rna", relativePositionsGag[2, :])
print("dim", relativePositionsGag[3, :])#IF1
print("tr1", relativePositionsGag[4, :])#IF2
print("tr2", relativePositionsGag[7, :])#IF5
print("hx1", relativePositionsGag[6, :])#IF4
print("hx2", relativePositionsGag[5, :])#IF3


# we will calculate the binding parameters for reactions

# dimer
# A(IF1)-B(IF1) is dimerization
# A is GAG 0, B is GAG 1
c1 = np.array(finalPositionsVec[0, :])
c2 = np.array(finalPositionsVec[8, :])
p1 = np.array(finalPositionsVec[3, :])
p2 = np.array(finalPositionsVec[11, :])
print("dimer:")
calculateAngles(c1, c2, p1, p2)

# trimer
# A(IF2)-C(IF5) is trimmerization
# A is GAG 0, C is GAG 2
c1 = np.array(finalPositionsVec[0, :])
c2 = np.array(finalPositionsVec[16, :])
p1 = np.array(finalPositionsVec[4, :])
p2 = np.array(finalPositionsVec[23, :])
print("trimer:")
calculateAngles(c1, c2, p1, p2)

# hexamer
# A(IF3)-D(IF4) is hexamerization
# A is GAG 0, D is GAG 3
c1 = np.array(finalPositionsVec[0, :])
c2 = np.array(finalPositionsVec[24, :])
p1 = np.array(finalPositionsVec[5, :])
p2 = np.array(finalPositionsVec[30, :])
print("hexamer:")
calculateAngles(c1, c2, p1, p2)



