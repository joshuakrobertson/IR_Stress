## All Python Code Necessary for Analyses in Tabh et al (2021), Physiological Reports.

# First loading in necessary packages.

import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd 

# Defining required functions. 

def isRotationMatrix(R) :
	Rt = np.transpose(R)
	shouldBeIdentity = np.dot(Rt, R)
	I = np.identity(3, dtype = R.dtype)
	n = np.linalg.norm(I - shouldBeIdentity)
	return n < 1e-6
	
def rotationMatrixToEulerAngles(R) :
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    return np.array([x, y, z])

# Now testing out ePnP algorithm on a single bird before proceeding to bulk analysis. To do so, beginning by loading in an image of a given pigeon.

P1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ6_Handling00010001.jpg")
P1_Size = P1.shape

# Next, setting camera calibration (assuming focal point at center of image).

Lfocus = P1_Size[1]
Centre = (P1_Size[1],P1_Size[0])
Camera = np.array([
    [Lfocus, 0, Centre[0]], 
    [0, Lfocus, Centre[1]], 
    [0, 0, 1]], 
    dtype = "double")

print("Camera Matrix :\n {0}".format(Camera))

Distortion_Coef = np.zeros((4,1))

# Using custom 2D landmarks and general 3D model of pigeon to test our PnP model. Below are 2D landmarks from our pigeon image loaded in above.

Land = np.array([
    (465.833,118.167), # Beak tip
    (501.500,130.833), # Beak base
    (505.833,103.5), # Cyr base
    (535.17,120.167), # Left eyeball
    (548.5,134.833), # Left eye back
    (516.167,115.5)# Left eye front
], dtype = "double")

# Now adding model pigeon landmarks. Note that landmarks that cannot be seen in our imaged individual are commented out below.

Ref_Pig = np.array([
(0,0,0), # Beak tip
(19.991,5.712,0), # Beak base
(13.439,18.816,0), # Cyr base
(34.439,26.04,12.2), # Left eyeball
(42.671,25.032,9.7), # Left eye back
(26.543,20.496,9.7) # Left eye front
# (34.439,26.04,-12.2), # Right eyeball
# (42.671,25.032,-9.7), # Right eye back
# (26.543,20.496,-9.7) # Right eye front
], dtype = "double")

# Setting up camera matrix with 0 distortion (typical of a thermographic camera).

Lfocus = P1_Size[1]
Centre = (P1_Size[1],P1_Size[0])
Camera = np.array([
    [Lfocus, 0, Centre[0]], 
    [0, Lfocus, Centre[1]], 
    [0, 0, 1]], 
    dtype = "double")

print("Camera Matrix :\n {0}".format(Camera))

# Here, distortion is set.

Distortion_Coef = np.zeros((4,1))

# Now calling ePnP function in cv2 package. Note that random sample consensus ("ransac") is also used below.

(success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Pig, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

# Good. Now projecting a line ~ 1 cm directly in front of the bill and drawing it in our image. This projection is completed using our estimated rotation and translation matrices and will help us visually evaluate the accuracy of those matrices.

(bill_end_point2D, jacobian) = cv2.projectPoints(np.array([(-11.307,-8.922,0)]), rotation_vector, translation_vector, Camera, Distortion_Coef)

# Printing the line on our image.

for p in Land:
	    cv2.circle(P1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)

p1 = ( int(Land[0][0]), int(Land[0][1]))
p2 = ( int(bill_end_point2D[0][0][0]), int(bill_end_point2D[0][0][1]))

cv2.line(P1, p1, p2, (0,0,0), 2)

# And now showing image.

plt.imshow(cv2.cvtColor(P1, cv2.COLOR_BGR2RGB))
plt.show()

# Visual assessment is helpful, however, somewhat subjective. So, below, we calculate Euclidean residuals describing the summed and averaged distance between predicted 2D landmarks, and true 2D landmarks. We can use these Euclidean residuals to, more objectively, assess the accuracy of our epnp predictions.  

(yhat, jacobian) = cv2.projectPoints(Ref_Pig, rotation_vector, translation_vector, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0)

# Good. Now proceeding to analysing all data in bulk. Note that we randomly sample images per individual to visually assess the fit of our epnp projections. Nevertheless, Euclidean residuals are calculated for each image of each bird and are later used to filter the quality of rotation estimates for analysis.  

######################################################################################
################################## Bird 1 ############################################
######################################################################################

# Loading in landmark data and assessing random sample

TZ6 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ6/TZ6_Handling.csv") 
TZ6.head()

S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ6/Images/TZ6_Handling00010248.jpg")
S1Land = np.array([
    (TZ6.iloc[247,3],TZ6.iloc[247,4]), # Beak tip
    (TZ6.iloc[247,5],TZ6.iloc[247,6]), # Beak base
    (TZ6.iloc[247,7],TZ6.iloc[247,8]), # Cyr base
    (TZ6.iloc[247,9],TZ6.iloc[247,10]), # Left eyeball
    (TZ6.iloc[247,11],TZ6.iloc[247,12]), # Left eye back
    (TZ6.iloc[247,13],TZ6.iloc[247,14]), # Left eye front
    (TZ6.iloc[247,15],TZ6.iloc[247,16]), # Right eyeball
    (TZ6.iloc[247,17],TZ6.iloc[247,18]), # Right eye back
    (TZ6.iloc[247,19],TZ6.iloc[247,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ6.iloc[247,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good.

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ6/Images/TZ6_Handling00010481.jpg")
S2Land = np.array([
    (TZ6.iloc[480,3],TZ6.iloc[480,4]), # Beak tip
    (TZ6.iloc[480,5],TZ6.iloc[480,6]), # Beak base
    (TZ6.iloc[480,7],TZ6.iloc[480,8]), # Cyr base
    (TZ6.iloc[480,9],TZ6.iloc[480,10]), # Left eyeball
    (TZ6.iloc[480,11],TZ6.iloc[480,12]), # Left eye back
    (TZ6.iloc[480,13],TZ6.iloc[480,14]), # Left eye front
    (TZ6.iloc[480,15],TZ6.iloc[480,16]), # Right eyeball
    (TZ6.iloc[480,17],TZ6.iloc[480,18]), # Right eye back
    (TZ6.iloc[480,19],TZ6.iloc[480,20]) # Right eye front
    ], dtype = "double")    
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ6.iloc[480,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Again, good

S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ6/Images/TZ6_Handling00020045.jpg")
S3Land = np.array([
    (TZ6.iloc[534,3],TZ6.iloc[534,4]), # Beak tip
    (TZ6.iloc[534,5],TZ6.iloc[534,6]), # Beak base
    (TZ6.iloc[534,7],TZ6.iloc[534,8]), # Cyr base
    (TZ6.iloc[534,9],TZ6.iloc[534,10]), # Left eyeball
    (TZ6.iloc[534,11],TZ6.iloc[534,12]), # Left eye back
    (TZ6.iloc[534,13],TZ6.iloc[534,14]), # Left eye front
    (TZ6.iloc[534,15],TZ6.iloc[534,16]), # Right eyeball
    (TZ6.iloc[534,17],TZ6.iloc[534,18]), # Right eye back
    (TZ6.iloc[534,19],TZ6.iloc[534,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ6.iloc[534,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Checking rotation matrix

rmat, jac = cv2.Rodrigues(rotation_vector_S3)
angles, mtxR, mtxQ, Qx, Qy, Qz = cv2.RQDecomp3x3(rmat) # angles = pitch, yaw, roll; we are largely interested in the yaw
rotationMatrixToEulerAngles(rmat)*(180/math.pi)

# Great. Proceeding to looping through all images.

for i in range(len(TZ6.index)):
    Land = np.array([
    (TZ6.iloc[i,3],TZ6.iloc[i,4]), # Beak tip
    (TZ6.iloc[i,5],TZ6.iloc[i,6]), # Beak base
    (TZ6.iloc[i,7],TZ6.iloc[i,8]), # Cyr base
    (TZ6.iloc[i,9],TZ6.iloc[i,10]), # Left eyeball
    (TZ6.iloc[i,11],TZ6.iloc[i,12]), # Left eye back
    (TZ6.iloc[i,13],TZ6.iloc[i,14]), # Left eye front
    (TZ6.iloc[i,15],TZ6.iloc[i,16]), # Right eyeball
    (TZ6.iloc[i,17],TZ6.iloc[i,18]), # Right eye back
    (TZ6.iloc[i,19],TZ6.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ6.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ6.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ6.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ6/Rotation_Out.csv')

######################################################################################
################################## Bird 2 ############################################
######################################################################################

TZ18 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ18/TZ18_Control.csv") 

S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ18/Images/TZ18_Control00010001.jpg")
S1Land = np.array([
    (TZ18.iloc[0,3],TZ18.iloc[0,4]), # Beak tip
    (TZ18.iloc[0,5],TZ18.iloc[0,6]), # Beak base
    (TZ18.iloc[0,7],TZ18.iloc[0,8]), # Cyr base
    (TZ18.iloc[0,9],TZ18.iloc[0,10]), # Left eyeball
    (TZ18.iloc[0,11],TZ18.iloc[0,12]), # Left eye back
    (TZ18.iloc[0,13],TZ18.iloc[0,14]), # Left eye front
    (TZ18.iloc[0,15],TZ18.iloc[0,16]), # Right eyeball
    (TZ18.iloc[0,17],TZ18.iloc[0,18]), # Right eye back
    (TZ18.iloc[0,19],TZ18.iloc[0,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ18.iloc[0,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Mean euclidean distance is rather large, but projection appears okay. Look into this.

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ18/Images/TZ18_Control00010093.jpg")
S2Land = np.array([
    (TZ18.iloc[92,3],TZ18.iloc[92,4]), # Beak tip
    (TZ18.iloc[92,5],TZ18.iloc[92,6]), # Beak base
    (TZ18.iloc[92,7],TZ18.iloc[92,8]), # Cyr base
    (TZ18.iloc[92,9],TZ18.iloc[92,10]), # Left eyeball
    (TZ18.iloc[92,11],TZ18.iloc[92,12]), # Left eye back
    (TZ18.iloc[92,13],TZ18.iloc[92,14]), # Left eye front
    (TZ18.iloc[92,15],TZ18.iloc[92,16]), # Right eyeball
    (TZ18.iloc[92,17],TZ18.iloc[92,18]), # Right eye back
    (TZ18.iloc[92,19],TZ18.iloc[92,20]) # Right eye front
    ], dtype = "double")    
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ18.iloc[92,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Again, mean euclidean distance is high, but projection looks great.

S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ18/Images/TZ18_Control00010215.jpg")
S3Land = np.array([
    (TZ18.iloc[214,3],TZ18.iloc[214,4]), # Beak tip
    (TZ18.iloc[214,5],TZ18.iloc[214,6]), # Beak base
    (TZ18.iloc[214,7],TZ18.iloc[214,8]), # Cyr base
    (TZ18.iloc[214,9],TZ18.iloc[214,10]), # Left eyeball
    (TZ18.iloc[214,11],TZ18.iloc[214,12]), # Left eye back
    (TZ18.iloc[214,13],TZ18.iloc[214,14]), # Left eye front
    (TZ18.iloc[214,15],TZ18.iloc[214,16]), # Right eyeball
    (TZ18.iloc[214,17],TZ18.iloc[214,18]), # Right eye back
    (TZ18.iloc[214,19],TZ18.iloc[214,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ18.iloc[214,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# One more;

S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ18/Images/TZ18_Control00010238.jpg")
S4Land = np.array([
    (TZ18.iloc[237,3],TZ18.iloc[237,4]), # Beak tip
    (TZ18.iloc[237,5],TZ18.iloc[237,6]), # Beak base
    (TZ18.iloc[237,7],TZ18.iloc[237,8]), # Cyr base
    (TZ18.iloc[237,9],TZ18.iloc[237,10]), # Left eyeball
    (TZ18.iloc[237,11],TZ18.iloc[237,12]), # Left eye back
    (TZ18.iloc[237,13],TZ18.iloc[237,14]), # Left eye front
    (TZ18.iloc[237,15],TZ18.iloc[237,16]), # Right eyeball
    (TZ18.iloc[237,17],TZ18.iloc[237,18]), # Right eye back
    (TZ18.iloc[237,19],TZ18.iloc[237,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ18.iloc[237,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great; looping through this individual.

TZ18 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ18/TZ18_Control_Adj.csv") 

for i in range(len(TZ18.index)):
    Land = np.array([
    (TZ18.iloc[i,3],TZ18.iloc[i,4]), # Beak tip
    (TZ18.iloc[i,5],TZ18.iloc[i,6]), # Beak base
    (TZ18.iloc[i,7],TZ18.iloc[i,8]), # Cyr base
    (TZ18.iloc[i,9],TZ18.iloc[i,10]), # Left eyeball
    (TZ18.iloc[i,11],TZ18.iloc[i,12]), # Left eye back
    (TZ18.iloc[i,13],TZ18.iloc[i,14]), # Left eye front
    (TZ18.iloc[i,15],TZ18.iloc[i,16]), # Right eyeball
    (TZ18.iloc[i,17],TZ18.iloc[i,18]), # Right eye back
    (TZ18.iloc[i,19],TZ18.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ18.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ18.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ18.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ18/Rotation_Out.csv')

######################################################################################
################################## Bird 3 ############################################
######################################################################################

TZ11 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ11/TZ11_Handling.csv") 

S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010001.jpg")
S1Land = np.array([
    (TZ11.iloc[0,3],TZ11.iloc[0,4]), # Beak tip
    (TZ11.iloc[0,5],TZ11.iloc[0,6]), # Beak base
    (TZ11.iloc[0,7],TZ11.iloc[0,8]), # Cyr base
    (TZ11.iloc[0,9],TZ11.iloc[0,10]), # Left eyeball
    (TZ11.iloc[0,11],TZ11.iloc[0,12]), # Left eye back
    (TZ11.iloc[0,13],TZ11.iloc[0,14]), # Left eye front
    (TZ11.iloc[0,15],TZ11.iloc[0,16]), # Right eyeball
    (TZ11.iloc[0,17],TZ11.iloc[0,18]), # Right eye back
    (TZ11.iloc[0,19],TZ11.iloc[0,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ11.iloc[0,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Second image 

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010002.jpg")
S2Land = np.array([
    (TZ11.iloc[1,3],TZ11.iloc[1,4]), # Beak tip
    (TZ11.iloc[1,5],TZ11.iloc[1,6]), # Beak base
    (TZ11.iloc[1,7],TZ11.iloc[1,8]), # Cyr base
    (TZ11.iloc[1,9],TZ11.iloc[1,10]), # Left eyeball
    #(247.833,228.5), # Left eye back
    (TZ11.iloc[1,11],TZ11.iloc[1,12]), # Left eye back
    (TZ11.iloc[1,13],TZ11.iloc[1,14]), # Left eye front
    (TZ11.iloc[1,15],TZ11.iloc[1,16]), # Right eyeball
    (TZ11.iloc[1,17],TZ11.iloc[1,18]), # Right eye back
    (TZ11.iloc[1,19],TZ11.iloc[1,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ11.iloc[1,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Third image

S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010296.jpg")
S3Land = np.array([
    (TZ11.iloc[295,3],TZ11.iloc[295,4]), # Beak tip
    (TZ11.iloc[295,5],TZ11.iloc[295,6]), # Beak base
    #(370.75,285), # Cyr base
    (TZ11.iloc[295,7],TZ11.iloc[295,8]), # Cyr base
    (TZ11.iloc[295,9],TZ11.iloc[295,10]), # Left eyeball
    (TZ11.iloc[295,11],TZ11.iloc[295,12]), # Left eye back
    (TZ11.iloc[295,13],TZ11.iloc[295,14]), # Left eye front
    (TZ11.iloc[295,15],TZ11.iloc[295,16]), # Right eyeball
    (TZ11.iloc[295,17],TZ11.iloc[295,18]), # Right eye back
    #(307.75,240.25), # Right eye back
    (TZ11.iloc[295,19],TZ11.iloc[295,20]) # Right eye front
    #(337.5,270.25) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ11.iloc[295,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Fourth image

S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010246.jpg")
S4Land = np.array([
    (TZ11.iloc[245,3],TZ11.iloc[245,4]), # Beak tip
    (TZ11.iloc[245,5],TZ11.iloc[245,6]), # Beak base
    #(370.75,285), # Cyr base
    (TZ11.iloc[245,7],TZ11.iloc[245,8]), # Cyr base
    (TZ11.iloc[245,9],TZ11.iloc[245,10]), # Left eyeball
    (TZ11.iloc[245,11],TZ11.iloc[245,12]), # Left eye back
    (TZ11.iloc[245,13],TZ11.iloc[245,14]), # Left eye front
    (TZ11.iloc[245,15],TZ11.iloc[245,16]), # Right eyeball
    (TZ11.iloc[245,17],TZ11.iloc[245,18]), # Right eye back
    #(307.75,240.25), # Right eye back
    (TZ11.iloc[245,19],TZ11.iloc[245,20]) # Right eye front
    #(337.5,270.25) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ11.iloc[245,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S4, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Trying another

S5 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010423.jpg")
S5Land = np.array([
    (TZ11.iloc[422,3],TZ11.iloc[422,4]), # Beak tip
    (TZ11.iloc[422,5],TZ11.iloc[422,6]), # Beak base
    (TZ11.iloc[422,7],TZ11.iloc[422,8]), # Cyr base
    (TZ11.iloc[422,9],TZ11.iloc[422,10]), # Left eyeball
    (TZ11.iloc[422,11],TZ11.iloc[422,12]), # Left eye back
    (TZ11.iloc[422,13],TZ11.iloc[422,14]), # Left eye front
    (TZ11.iloc[422,15],TZ11.iloc[422,16]), # Right eyeball
    (TZ11.iloc[422,17],TZ11.iloc[422,18]), # Right eye back
    (TZ11.iloc[422,19],TZ11.iloc[422,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I5 = (TZ11.iloc[422,]==-1).values
row_vals = [t for t, x in enumerate(I5) if x]
if len(I5) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S5Land = np.delete(S5Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S5, rotation_vector_S5, translation_vector_S5, inliers_S5) = cv2.solvePnPRansac(Ref_Bird, S5Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S5, jacobian_S5) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
for p in S5Land:
	    cv2.circle(S5, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S5Land[0][0]), int(S5Land[0][1]))
p2 = ( int(bill_end_point2D_S5[0][0][0]), int(bill_end_point2D_S5[0][0][1]))
cv2.line(S5, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S5, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# One more 
S6 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11_New/TZ11_handling00010629.jpg")
#S6 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ11/Images/TZ11_handling00010629.jpg")
S6Land = np.array([
    (TZ11.iloc[628,3],TZ11.iloc[628,4]), # Beak tip
    (TZ11.iloc[628,5],TZ11.iloc[628,6]), # Beak base
    (TZ11.iloc[628,7],TZ11.iloc[628,8]), # Cyr base
    (TZ11.iloc[628,9],TZ11.iloc[628,10]), # Left eyeball
    (TZ11.iloc[628,11],TZ11.iloc[628,12]), # Left eye back
    (TZ11.iloc[628,13],TZ11.iloc[628,14]), # Left eye front
    (TZ11.iloc[628,15],TZ11.iloc[628,16]), # Right eyeball
    (TZ11.iloc[628,17],TZ11.iloc[628,18]), # Right eye back
    (TZ11.iloc[628,19],TZ11.iloc[628,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I6 = (TZ11.iloc[628,]==-1).values
row_vals = [t for t, x in enumerate(I6) if x]
if len(I6) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S6Land = np.delete(S6Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S6, rotation_vector_S6, translation_vector_S6, inliers_S6) = cv2.solvePnPRansac(Ref_Bird, S6Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S6, jacobian_S6) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
for p in S6Land:
	    cv2.circle(S6, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S6Land[0][0]), int(S6Land[0][1]))
p2 = ( int(bill_end_point2D_S6[0][0][0]), int(bill_end_point2D_S6[0][0][1]))
cv2.line(S6, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S6, cv2.COLOR_BGR2RGB))
plt.show()
plt.savefig('/home/joshk/Desktop/Figure_2.jpeg', dpi=1000)

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. Looping through.

TZ11 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ11/TZ11_Handling_Adj.csv") 

for i in range(len(TZ11.index)):
    Land = np.array([
    (TZ11.iloc[i,3],TZ11.iloc[i,4]), # Beak tip
    (TZ11.iloc[i,5],TZ11.iloc[i,6]), # Beak base
    (TZ11.iloc[i,7],TZ11.iloc[i,8]), # Cyr base
    (TZ11.iloc[i,9],TZ11.iloc[i,10]), # Left eyeball
    (TZ11.iloc[i,11],TZ11.iloc[i,12]), # Left eye back
    (TZ11.iloc[i,13],TZ11.iloc[i,14]), # Left eye front
    (TZ11.iloc[i,15],TZ11.iloc[i,16]), # Right eyeball
    (TZ11.iloc[i,17],TZ11.iloc[i,18]), # Right eye back
    (TZ11.iloc[i,19],TZ11.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ11.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ11.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ11.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ11/Rotation_Out.csv')

######################################################################################
################################## Bird 4 ############################################
######################################################################################

# Continuing to TZ09
# First image

TZ9 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ9/TZ9_Control_Adj.csv") 
im = 1
S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ9/Images/TZ9_control00010002.jpg")
S1Land = np.array([
    (TZ9.iloc[im,3],TZ9.iloc[im,4]), # Beak tip
    (TZ9.iloc[im,5],TZ9.iloc[im,6]), # Beak base
    (TZ9.iloc[im,7],TZ9.iloc[im,8]), # Cyr base
    (TZ9.iloc[im,9],TZ9.iloc[im,10]), # Left eyeball
    (TZ9.iloc[im,11],TZ9.iloc[im,12]), # Left eye back
    (TZ9.iloc[im,13],TZ9.iloc[im,14]), # Left eye front
    (TZ9.iloc[im,15],TZ9.iloc[im,16]), # Right eyeball
    (TZ9.iloc[im,17],TZ9.iloc[im,18]), # Right eye back
    (TZ9.iloc[im,19],TZ9.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ9.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Second image 
#S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ9_New/TZ9_control00010264.jpg")
S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ9/Images/TZ9_control00010264.jpg")
im = 263
S2Land = np.array([
    (TZ9.iloc[im,3],TZ9.iloc[im,4]), # Beak tip
    (TZ9.iloc[im,5],TZ9.iloc[im,6]), # Beak base
    (TZ9.iloc[im,7],TZ9.iloc[im,8]), # Cyr base
    (TZ9.iloc[im,9],TZ9.iloc[im,10]), # Left eyeball
    (TZ9.iloc[im,11],TZ9.iloc[im,12]), # Left eye back
    (TZ9.iloc[im,13],TZ9.iloc[im,14]), # Left eye front
    (TZ9.iloc[im,15],TZ9.iloc[im,16]), # Right eyeball
    (TZ9.iloc[im,17],TZ9.iloc[im,18]), # Right eye back
    (TZ9.iloc[im,19],TZ9.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ9.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()
plt.savefig('/home/joshk/Desktop/Figure_1.jpeg', dpi=1000)

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Third image

im = 322
S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ9/Images/TZ9_control00010323.jpg")
S3Land = np.array([
    (TZ9.iloc[im,3],TZ9.iloc[im,4]), # Beak tip
    (TZ9.iloc[im,5],TZ9.iloc[im,6]), # Beak base
    (TZ9.iloc[im,7],TZ9.iloc[im,8]), # Cyr base
    (TZ9.iloc[im,9],TZ9.iloc[im,10]), # Left eyeball
    (TZ9.iloc[im,11],TZ9.iloc[im,12]), # Left eye back
    (TZ9.iloc[im,13],TZ9.iloc[im,14]), # Left eye front
    (TZ9.iloc[im,15],TZ9.iloc[im,16]), # Right eyeball
    (TZ9.iloc[im,17],TZ9.iloc[im,18]), # Right eye back
    (TZ9.iloc[im,19],TZ9.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ9.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. Looping again.

TZ9 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ9/TZ9_Control_Adj.csv") 

for i in range(len(TZ9.index)):
    Land = np.array([
    (TZ9.iloc[i,3],TZ9.iloc[i,4]), # Beak tip
    (TZ9.iloc[i,5],TZ9.iloc[i,6]), # Beak base
    (TZ9.iloc[i,7],TZ9.iloc[i,8]), # Cyr base
    (TZ9.iloc[i,9],TZ9.iloc[i,10]), # Left eyeball
    (TZ9.iloc[i,11],TZ9.iloc[i,12]), # Left eye back
    (TZ9.iloc[i,13],TZ9.iloc[i,14]), # Left eye front
    (TZ9.iloc[i,15],TZ9.iloc[i,16]), # Right eyeball
    (TZ9.iloc[i,17],TZ9.iloc[i,18]), # Right eye back
    (TZ9.iloc[i,19],TZ9.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ9.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ9.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ9.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ9/Rotation_Out.csv')

######################################################################################
################################## Bird 5 ############################################
######################################################################################

# Loading in and sampling TZ22

TZ22 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ22/TZ22_Control.csv") 
im = 0
S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ22/Images/TZ22_Control0002-10001.jpg")
S1Land = np.array([
    (TZ22.iloc[im,3],TZ22.iloc[im,4]), # Beak tip
    (TZ22.iloc[im,5],TZ22.iloc[im,6]), # Beak base
    (TZ22.iloc[im,7],TZ22.iloc[im,8]), # Cyr base
    (TZ22.iloc[im,9],TZ22.iloc[im,10]), # Left eyeball
    (TZ22.iloc[im,11],TZ22.iloc[im,12]), # Left eye back
    (TZ22.iloc[im,13],TZ22.iloc[im,14]), # Left eye front
    (TZ22.iloc[im,15],TZ22.iloc[im,16]), # Right eyeball
    (TZ22.iloc[im,17],TZ22.iloc[im,18]), # Right eye back
    (TZ22.iloc[im,19],TZ22.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ22.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Second image 

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ22/Images/TZ22_Control0002-10002.jpg")
im = 1
S2Land = np.array([
    (TZ22.iloc[im,3],TZ22.iloc[im,4]), # Beak tip
    (TZ22.iloc[im,5],TZ22.iloc[im,6]), # Beak base
    (TZ22.iloc[im,7],TZ22.iloc[im,8]), # Cyr base
    (TZ22.iloc[im,9],TZ22.iloc[im,10]), # Left eyeball
    (TZ22.iloc[im,11],TZ22.iloc[im,12]), # Left eye back
    (TZ22.iloc[im,13],TZ22.iloc[im,14]), # Left eye front
    (TZ22.iloc[im,15],TZ22.iloc[im,16]), # Right eyeball
    (TZ22.iloc[im,17],TZ22.iloc[im,18]), # Right eye back
    (TZ22.iloc[im,19],TZ22.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ22.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Third image

im = 535
S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ22/Images/TZ22_Control0002-10536.jpg")
S3Land = np.array([
    (TZ22.iloc[im,3],TZ22.iloc[im,4]), # Beak tip
    (TZ22.iloc[im,5],TZ22.iloc[im,6]), # Beak base
    (TZ22.iloc[im,7],TZ22.iloc[im,8]), # Cyr base
    (TZ22.iloc[im,9],TZ22.iloc[im,10]), # Left eyeball
    (TZ22.iloc[im,11],TZ22.iloc[im,12]), # Left eye back
    (TZ22.iloc[im,13],TZ22.iloc[im,14]), # Left eye front
    (TZ22.iloc[im,15],TZ22.iloc[im,16]), # Right eyeball
    (TZ22.iloc[im,17],TZ22.iloc[im,18]), # Right eye back
    (TZ22.iloc[im,19],TZ22.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ22.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Running another 

im = 93
S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ22/Images/TZ22_Control0002-10094.jpg")
S4Land = np.array([
    (TZ22.iloc[im,3],TZ22.iloc[im,4]), # Beak tip
    (TZ22.iloc[im,5],TZ22.iloc[im,6]), # Beak base
    (TZ22.iloc[im,7],TZ22.iloc[im,8]), # Cyr base
    (TZ22.iloc[im,9],TZ22.iloc[im,10]), # Left eyeball
    (TZ22.iloc[im,11],TZ22.iloc[im,12]), # Left eye back
    (TZ22.iloc[im,13],TZ22.iloc[im,14]), # Left eye front
    (TZ22.iloc[im,15],TZ22.iloc[im,16]), # Right eyeball
    (TZ22.iloc[im,17],TZ22.iloc[im,18]), # Right eye back
    (TZ22.iloc[im,19],TZ22.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ22.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S4, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

## Good. One more

im = 200
S5 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ22/Images/TZ22_Control0002-10201.jpg")
S5Land = np.array([
    (TZ22.iloc[im,3],TZ22.iloc[im,4]), # Beak tip
    (TZ22.iloc[im,5],TZ22.iloc[im,6]), # Beak base
    (TZ22.iloc[im,7],TZ22.iloc[im,8]), # Cyr base
    (TZ22.iloc[im,9],TZ22.iloc[im,10]), # Left eyeball
    (TZ22.iloc[im,11],TZ22.iloc[im,12]), # Left eye back
    (TZ22.iloc[im,13],TZ22.iloc[im,14]), # Left eye front
    (TZ22.iloc[im,15],TZ22.iloc[im,16]), # Right eyeball
    (TZ22.iloc[im,17],TZ22.iloc[im,18]), # Right eye back
    (TZ22.iloc[im,19],TZ22.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I5 = (TZ22.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I5) if x]
if len(I5) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S5Land = np.delete(S5Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S5, rotation_vector_S5, translation_vector_S5, inliers_S5) = cv2.solvePnPRansac(Ref_Bird, S5Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S5, jacobian_S5) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
for p in S5Land:
	    cv2.circle(S5, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S5Land[0][0]), int(S5Land[0][1]))
p2 = ( int(bill_end_point2D_S5[0][0][0]), int(bill_end_point2D_S5[0][0][1]))
cv2.line(S5, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S5, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Quite reasonable, despite awkward positioning of bird. Good.
# Looping images.

TZ22 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ22/TZ22_Control.csv") 

for i in range(len(TZ22.index)):
    Land = np.array([
    (TZ22.iloc[i,3],TZ22.iloc[i,4]), # Beak tip
    (TZ22.iloc[i,5],TZ22.iloc[i,6]), # Beak base
    (TZ22.iloc[i,7],TZ22.iloc[i,8]), # Cyr base
    (TZ22.iloc[i,9],TZ22.iloc[i,10]), # Left eyeball
    (TZ22.iloc[i,11],TZ22.iloc[i,12]), # Left eye back
    (TZ22.iloc[i,13],TZ22.iloc[i,14]), # Left eye front
    (TZ22.iloc[i,15],TZ22.iloc[i,16]), # Right eyeball
    (TZ22.iloc[i,17],TZ22.iloc[i,18]), # Right eye back
    (TZ22.iloc[i,19],TZ22.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ22.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ22.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ22.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ22/Rotation_Out.csv')

######################################################################################
################################## Bird 6 ############################################
######################################################################################

# Testing TZ24

TZ24 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ24/TZ24_Handling.csv") 
im = 0
S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010001.jpg")
S1Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great.
# Second image 

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010035.jpg")
im = 34
S2Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    #(-1,-1),
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Third image

im = 49
S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010050.jpg")
S3Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. Next.

im = 109
S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010110.jpg")
S4Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S4, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# One more.

im = 127
S5 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010128.jpg")
S5Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I5 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I5) if x]
if len(I5) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S5Land = np.delete(S5Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S5, rotation_vector_S5, translation_vector_S5, inliers_S5) = cv2.solvePnPRansac(Ref_Bird, S5Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S5, jacobian_S5) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
for p in S5Land:
	    cv2.circle(S5, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S5Land[0][0]), int(S5Land[0][1]))
p2 = ( int(bill_end_point2D_S5[0][0][0]), int(bill_end_point2D_S5[0][0][1]))
cv2.line(S5, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S5, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good. Continuing with a few later images to be cautious.

im = 207
S6 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010208.jpg")
S6Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I6 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I6) if x]
if len(I6) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S6Land = np.delete(S6Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S6, rotation_vector_S6, translation_vector_S6, inliers_S6) = cv2.solvePnPRansac(Ref_Bird, S6Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S6, jacobian_S6) = cv2.projectPoints(np.array([(0, 24.269, 0.0)]), rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
for p in S6Land:
	    cv2.circle(S6, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S6Land[0][0]), int(S6Land[0][1]))
p2 = ( int(bill_end_point2D_S6[0][0][0]), int(bill_end_point2D_S6[0][0][1]))
cv2.line(S6, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S6, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Note that the line here projects to directly anterior to the left eye. Good.

im = 499
S7 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24/Images/TZ24_expII00010500.jpg")
S7Land = np.array([
    (TZ24.iloc[im,3],TZ24.iloc[im,4]), # Beak tip
    (TZ24.iloc[im,5],TZ24.iloc[im,6]), # Beak base
    (TZ24.iloc[im,7],TZ24.iloc[im,8]), # Cyr base
    (TZ24.iloc[im,9],TZ24.iloc[im,10]), # Left eyeball
    (TZ24.iloc[im,11],TZ24.iloc[im,12]), # Left eye back
    (TZ24.iloc[im,13],TZ24.iloc[im,14]), # Left eye front
    (TZ24.iloc[im,15],TZ24.iloc[im,16]), # Right eyeball
    (TZ24.iloc[im,17],TZ24.iloc[im,18]), # Right eye back
    (TZ24.iloc[im,19],TZ24.iloc[im,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I7 = (TZ24.iloc[im,]==-1).values
row_vals = [t for t, x in enumerate(I7) if x]
if len(I7) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S7Land = np.delete(S7Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S7, rotation_vector_S7, translation_vector_S7, inliers_S7) = cv2.solvePnPRansac(Ref_Bird, S7Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S7, jacobian_S7) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
for p in S7Land:
	    cv2.circle(S7, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S7Land[0][0]), int(S7Land[0][1]))
p2 = ( int(bill_end_point2D_S7[0][0][0]), int(bill_end_point2D_S7[0][0][1]))
cv2.line(S7, p1, p2, (255,255,255), 2)
plt.imshow(cv2.cvtColor(S7, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. Looping.

TZ24 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ24/TZ24_Handling.csv") 

for i in range(len(TZ24.index)):
    Land = np.array([
    (TZ24.iloc[i,3],TZ24.iloc[i,4]), # Beak tip
    (TZ24.iloc[i,5],TZ24.iloc[i,6]), # Beak base
    (TZ24.iloc[i,7],TZ24.iloc[i,8]), # Cyr base
    (TZ24.iloc[i,9],TZ24.iloc[i,10]), # Left eyeball
    (TZ24.iloc[i,11],TZ24.iloc[i,12]), # Left eye back
    (TZ24.iloc[i,13],TZ24.iloc[i,14]), # Left eye front
    (TZ24.iloc[i,15],TZ24.iloc[i,16]), # Right eyeball
    (TZ24.iloc[i,17],TZ24.iloc[i,18]), # Right eye back
    (TZ24.iloc[i,19],TZ24.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ24.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ24.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ24.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ24/Rotation_Out.csv')

##########################################################################
# Running through TZ8
##########################################################################

TZ8 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ8/TZ8_Handling.csv") 
TZ8.head()

S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00011734.jpg")
S1Land = np.array([
    (TZ8.iloc[17,3],TZ8.iloc[17,4]), # Beak tip
    (TZ8.iloc[17,5],TZ8.iloc[17,6]), # Beak base
    (TZ8.iloc[17,7],TZ8.iloc[17,8]), # Cyr base
    (TZ8.iloc[17,9],TZ8.iloc[17,10]), # Left eyeball
    (TZ8.iloc[17,11],TZ8.iloc[17,12]), # Left eye back
    (TZ8.iloc[17,13],TZ8.iloc[17,14]), # Left eye front
    (TZ8.iloc[17,15],TZ8.iloc[17,16]), # Right eyeball
    (TZ8.iloc[17,17],TZ8.iloc[17,18]), # Right eye back
    (TZ8.iloc[17,19],TZ8.iloc[17,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ8.iloc[17,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good.

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00011735.jpg")
S2Land = np.array([
    (TZ8.iloc[18,3],TZ8.iloc[18,4]), # Beak tip
    #(-1,-1),
    (TZ8.iloc[18,5],TZ8.iloc[18,6]), # Beak base
    (TZ8.iloc[18,7],TZ8.iloc[18,8]), # Cyr base
    (TZ8.iloc[18,9],TZ8.iloc[18,10]), # Left eyeball
    (TZ8.iloc[18,11],TZ8.iloc[18,12]), # Left eye back
    (TZ8.iloc[18,13],TZ8.iloc[18,14]), # Left eye front
    (TZ8.iloc[18,15],TZ8.iloc[18,16]), # Right eyeball
    (TZ8.iloc[18,17],TZ8.iloc[18,18]), # Right eye back
    (TZ8.iloc[18,19],TZ8.iloc[18,20]) # Right eye front
    ], dtype = "double")    
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ8.iloc[18,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Again, good

S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00011736.jpg")
S3Land = np.array([
    (TZ8.iloc[19,3],TZ8.iloc[19,4]), # Beak tip
    (TZ8.iloc[19,5],TZ8.iloc[19,6]), # Beak base
    (TZ8.iloc[19,7],TZ8.iloc[19,8]), # Cyr base
    (TZ8.iloc[19,9],TZ8.iloc[19,10]), # Left eyeball
    (TZ8.iloc[19,11],TZ8.iloc[19,12]), # Left eye back
    (TZ8.iloc[19,13],TZ8.iloc[19,14]), # Left eye front
    (TZ8.iloc[19,15],TZ8.iloc[19,16]), # Right eyeball
    (TZ8.iloc[19,17],TZ8.iloc[19,18]), # Right eye back
    (TZ8.iloc[19,19],TZ8.iloc[19,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ8.iloc[19,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good. Continuing with a few more images.

S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00011736.jpg")
S4Land = np.array([
    (TZ8.iloc[22,3],TZ8.iloc[22,4]), # Beak tip
    (TZ8.iloc[22,5],TZ8.iloc[22,6]), # Beak base
    (TZ8.iloc[22,7],TZ8.iloc[22,8]), # Cyr base
    (TZ8.iloc[22,9],TZ8.iloc[22,10]), # Left eyeball
    (TZ8.iloc[22,11],TZ8.iloc[22,12]), # Left eye back
    (TZ8.iloc[22,13],TZ8.iloc[22,14]), # Left eye front
    (TZ8.iloc[22,15],TZ8.iloc[22,16]), # Right eyeball
    (TZ8.iloc[22,17],TZ8.iloc[22,18]), # Right eye back
    (TZ8.iloc[22,19],TZ8.iloc[22,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ8.iloc[22,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S4, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good. 

S5 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00020009.jpg")
S5Land = np.array([
    (TZ8.iloc[109,3],TZ8.iloc[109,4]), # Beak tip
    (TZ8.iloc[109,5],TZ8.iloc[109,6]), # Beak base
    (TZ8.iloc[109,7],TZ8.iloc[109,8]), # Cyr base
    (TZ8.iloc[109,9],TZ8.iloc[109,10]), # Left eyeball
    (TZ8.iloc[109,11],TZ8.iloc[109,12]), # Left eye back
    (TZ8.iloc[109,13],TZ8.iloc[109,14]), # Left eye front
    (TZ8.iloc[109,15],TZ8.iloc[109,16]), # Right eyeball
    (TZ8.iloc[109,17],TZ8.iloc[109,18]), # Right eye back
    (TZ8.iloc[109,19],TZ8.iloc[109,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I5 = (TZ8.iloc[109,]==-1).values
row_vals = [t for t, x in enumerate(I5) if x]
if len(I5) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S5Land = np.delete(S5Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S5, rotation_vector_S5, translation_vector_S5, inliers_S5) = cv2.solvePnPRansac(Ref_Bird, S5Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S5, jacobian_S5) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
for p in S5Land:
	    cv2.circle(S5, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S5Land[0][0]), int(S5Land[0][1]))
p2 = ( int(bill_end_point2D_S5[0][0][0]), int(bill_end_point2D_S5[0][0][1]))
cv2.line(S5, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S5, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Fair. One more.

S6 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00020208.jpg")
S6Land = np.array([
    (TZ8.iloc[309,3],TZ8.iloc[309,4]), # Beak tip
    (TZ8.iloc[309,5],TZ8.iloc[309,6]), # Beak base
    (TZ8.iloc[309,7],TZ8.iloc[309,8]), # Cyr base
    (TZ8.iloc[309,9],TZ8.iloc[309,10]), # Left eyeball
    (TZ8.iloc[309,11],TZ8.iloc[309,12]), # Left eye back
    (TZ8.iloc[309,13],TZ8.iloc[309,14]), # Left eye front
    (TZ8.iloc[309,15],TZ8.iloc[309,16]), # Right eyeball
    (TZ8.iloc[309,17],TZ8.iloc[309,18]), # Right eye back
    (TZ8.iloc[309,19],TZ8.iloc[309,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I6 = (TZ8.iloc[309,]==-1).values
row_vals = [t for t, x in enumerate(I6) if x]
if len(I6) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S6Land = np.delete(S6Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S6, rotation_vector_S6, translation_vector_S6, inliers_S6) = cv2.solvePnPRansac(Ref_Bird, S6Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S6, jacobian_S6) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
for p in S6Land:
	    cv2.circle(S6, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S6Land[0][0]), int(S6Land[0][1]))
p2 = ( int(bill_end_point2D_S6[0][0][0]), int(bill_end_point2D_S6[0][0][1]))
cv2.line(S6, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S6, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# During handling?

S7 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ8/Images/TZ8_Handling00020548.jpg")
S7Land = np.array([
    (TZ8.iloc[648,3],TZ8.iloc[648,4]), # Beak tip
    (TZ8.iloc[648,5],TZ8.iloc[648,6]), # Beak base
    (TZ8.iloc[648,7],TZ8.iloc[648,8]), # Cyr base
    (TZ8.iloc[648,9],TZ8.iloc[648,10]), # Left eyeball
    (TZ8.iloc[648,11],TZ8.iloc[648,12]), # Left eye back
    (TZ8.iloc[648,13],TZ8.iloc[648,14]), # Left eye front
    (TZ8.iloc[648,15],TZ8.iloc[648,16]), # Right eyeball
    (TZ8.iloc[648,17],TZ8.iloc[648,18]), # Right eye back
    (TZ8.iloc[648,19],TZ8.iloc[648,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I7 = (TZ8.iloc[648,]==-1).values
row_vals = [t for t, x in enumerate(I7) if x]
if len(I7) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S7Land = np.delete(S7Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S7, rotation_vector_S7, translation_vector_S7, inliers_S7) = cv2.solvePnPRansac(Ref_Bird, S7Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S7, jacobian_S7) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
for p in S7Land:
	    cv2.circle(S7, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S7Land[0][0]), int(S7Land[0][1]))
p2 = ( int(bill_end_point2D_S7[0][0][0]), int(bill_end_point2D_S7[0][0][1]))
cv2.line(S7, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S7, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. Looping through images.

TZ8 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ8/TZ8_Handling.csv") 

for i in range(len(TZ8.index)):
    Land = np.array([
    (TZ8.iloc[i,3],TZ8.iloc[i,4]), # Beak tip
    (TZ8.iloc[i,5],TZ8.iloc[i,6]), # Beak base
    (TZ8.iloc[i,7],TZ8.iloc[i,8]), # Cyr base
    (TZ8.iloc[i,9],TZ8.iloc[i,10]), # Left eyeball
    (TZ8.iloc[i,11],TZ8.iloc[i,12]), # Left eye back
    (TZ8.iloc[i,13],TZ8.iloc[i,14]), # Left eye front
    (TZ8.iloc[i,15],TZ8.iloc[i,16]), # Right eyeball
    (TZ8.iloc[i,17],TZ8.iloc[i,18]), # Right eye back
    (TZ8.iloc[i,19],TZ8.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ8.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ8.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ8.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ8/Rotation_Out.csv')

############################################################################
## Validating TZ24_2
############################################################################

TZ24_2 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/TZ24_Control.csv") 
TZ24_2.head()

S1 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040068.jpg")
S1Land = np.array([
    (TZ24_2.iloc[67,3],TZ24_2.iloc[67,4]), # Beak tip
    (TZ24_2.iloc[67,5],TZ24_2.iloc[67,6]), # Beak base
    (TZ24_2.iloc[67,7],TZ24_2.iloc[67,8]), # Cyr base
    (TZ24_2.iloc[67,9],TZ24_2.iloc[67,10]), # Left eyeball
    (TZ24_2.iloc[67,11],TZ24_2.iloc[67,12]), # Left eye back
    (TZ24_2.iloc[67,13],TZ24_2.iloc[67,14]), # Left eye front
    (TZ24_2.iloc[67,15],TZ24_2.iloc[67,16]), # Right eyeball
    (TZ24_2.iloc[67,17],TZ24_2.iloc[67,18]), # Right eye back
    (TZ24_2.iloc[67,19],TZ24_2.iloc[67,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I1 = (TZ24_2.iloc[67,]==-1).values
row_vals = [t for t, x in enumerate(I1) if x]
if len(I1) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S1Land = np.delete(S1Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S1, rotation_vector_S1, translation_vector_S1, inliers_S1) = cv2.solvePnPRansac(Ref_Bird, S1Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S1, jacobian_S1) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
for p in S1Land:
	    cv2.circle(S1, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S1Land[0][0]), int(S1Land[0][1]))
p2 = ( int(bill_end_point2D_S1[0][0][0]), int(bill_end_point2D_S1[0][0][1]))
cv2.line(S1, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S1, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S1, translation_vector_S1, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S1Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good.

S2 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040077.jpg")
S2Land = np.array([
    (TZ24_2.iloc[76,3],TZ24_2.iloc[76,4]), # Beak tip
    (TZ24_2.iloc[76,5],TZ24_2.iloc[76,6]), # Beak base
    (TZ24_2.iloc[76,7],TZ24_2.iloc[76,8]), # Cyr base
    (TZ24_2.iloc[76,9],TZ24_2.iloc[76,10]), # Left eyeball
    (TZ24_2.iloc[76,11],TZ24_2.iloc[76,12]), # Left eye back
    (TZ24_2.iloc[76,13],TZ24_2.iloc[76,14]), # Left eye front
    (TZ24_2.iloc[76,15],TZ24_2.iloc[76,16]), # Right eyeball
    (TZ24_2.iloc[76,17],TZ24_2.iloc[76,18]), # Right eye back
    (TZ24_2.iloc[76,19],TZ24_2.iloc[76,20]) # Right eye front
    ], dtype = "double")    
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I2 = (TZ24_2.iloc[76,]==-1).values
row_vals = [t for t, x in enumerate(I2) if x]
if len(I2) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):
        S2Land = np.delete(S2Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S2, rotation_vector_S2, translation_vector_S2, inliers_S2) = cv2.solvePnPRansac(Ref_Bird, S2Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S2, jacobian_S2) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
for p in S2Land:
	    cv2.circle(S2, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S2Land[0][0]), int(S2Land[0][1]))
p2 = ( int(bill_end_point2D_S2[0][0][0]), int(bill_end_point2D_S2[0][0][1]))
cv2.line(S2, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S2, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S2, translation_vector_S2, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S2Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Again, good

S3 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040100.jpg")
S3Land = np.array([
    (TZ24_2.iloc[99,3],TZ24_2.iloc[99,4]), # Beak tip
    (TZ24_2.iloc[99,5],TZ24_2.iloc[99,6]), # Beak base
    (TZ24_2.iloc[99,7],TZ24_2.iloc[99,8]), # Cyr base
    (TZ24_2.iloc[99,9],TZ24_2.iloc[99,10]), # Left eyeball
    (TZ24_2.iloc[99,11],TZ24_2.iloc[99,12]), # Left eye back
    (TZ24_2.iloc[99,13],TZ24_2.iloc[99,14]), # Left eye front
    (TZ24_2.iloc[99,15],TZ24_2.iloc[99,16]), # Right eyeball
    (TZ24_2.iloc[99,17],TZ24_2.iloc[99,18]), # Right eye back
    (TZ24_2.iloc[99,19],TZ24_2.iloc[99,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I3 = (TZ24_2.iloc[99,]==-1).values
row_vals = [t for t, x in enumerate(I3) if x]
if len(I3) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S3Land = np.delete(S3Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S3, rotation_vector_S3, translation_vector_S3, inliers_S3) = cv2.solvePnPRansac(Ref_Bird, S3Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S3, jacobian_S3) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
for p in S3Land:
	    cv2.circle(S3, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S3Land[0][0]), int(S3Land[0][1]))
p2 = ( int(bill_end_point2D_S3[0][0][0]), int(bill_end_point2D_S3[0][0][1]))
cv2.line(S3, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S3, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S3, translation_vector_S3, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S3Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. 

S4 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040147.jpg")
S4Land = np.array([
    (TZ24_2.iloc[146,3],TZ24_2.iloc[146,4]), # Beak tip
    (TZ24_2.iloc[146,5],TZ24_2.iloc[146,6]), # Beak base
    (-1,-1),
    #(TZ24_2.iloc[146,7],TZ24_2.iloc[146,8]), # Cyr base
    (TZ24_2.iloc[146,9],TZ24_2.iloc[146,10]), # Left eyeball
    (TZ24_2.iloc[146,11],TZ24_2.iloc[146,12]), # Left eye back
    (TZ24_2.iloc[146,13],TZ24_2.iloc[146,14]), # Left eye front
    (TZ24_2.iloc[146,15],TZ24_2.iloc[146,16]), # Right eyeball
    (TZ24_2.iloc[146,17],TZ24_2.iloc[146,18]), # Right eye back
    (TZ24_2.iloc[146,19],TZ24_2.iloc[146,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I4 = (TZ24_2.iloc[146,]==-1).values
row_vals = [t for t, x in enumerate(I4) if x]
if len(I4) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S4Land = np.delete(S4Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S4, rotation_vector_S4, translation_vector_S4, inliers_S4) = cv2.solvePnPRansac(Ref_Bird, S4Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S4, jacobian_S4) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
for p in S4Land:
	    cv2.circle(S4, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S4Land[0][0]), int(S4Land[0][1]))
p2 = ( int(bill_end_point2D_S4[0][0][0]), int(bill_end_point2D_S4[0][0][1]))
cv2.line(S4, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S4, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S4, translation_vector_S4, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S4Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good. 

S5 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040148.jpg")
S5Land = np.array([
    (TZ24_2.iloc[147,3],TZ24_2.iloc[147,4]), # Beak tip
    (TZ24_2.iloc[147,5],TZ24_2.iloc[147,6]), # Beak base
    (TZ24_2.iloc[147,7],TZ24_2.iloc[147,8]), # Cyr base
    (TZ24_2.iloc[147,9],TZ24_2.iloc[147,10]), # Left eyeball
    (TZ24_2.iloc[147,11],TZ24_2.iloc[147,12]), # Left eye back
    (TZ24_2.iloc[147,13],TZ24_2.iloc[147,14]), # Left eye front
    (TZ24_2.iloc[147,15],TZ24_2.iloc[147,16]), # Right eyeball
    (TZ24_2.iloc[147,17],TZ24_2.iloc[147,18]), # Right eye back
    (TZ24_2.iloc[147,19],TZ24_2.iloc[147,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I5 = (TZ24_2.iloc[147,]==-1).values
row_vals = [t for t, x in enumerate(I5) if x]
if len(I5) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S5Land = np.delete(S5Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S5, rotation_vector_S5, translation_vector_S5, inliers_S5) = cv2.solvePnPRansac(Ref_Bird, S5Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S5, jacobian_S5) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
for p in S5Land:
	    cv2.circle(S5, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S5Land[0][0]), int(S5Land[0][1]))
p2 = ( int(bill_end_point2D_S5[0][0][0]), int(bill_end_point2D_S5[0][0][1]))
cv2.line(S5, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S5, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S5, translation_vector_S5, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S5Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Fair. One more.

S6 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040149.jpg")
S6Land = np.array([
    (TZ24_2.iloc[148,3],TZ24_2.iloc[148,4]), # Beak tip
    (TZ24_2.iloc[148,5],TZ24_2.iloc[148,6]), # Beak base
    (TZ24_2.iloc[148,7],TZ24_2.iloc[148,8]), # Cyr base
    (TZ24_2.iloc[148,9],TZ24_2.iloc[148,10]), # Left eyeball
    (TZ24_2.iloc[148,11],TZ24_2.iloc[148,12]), # Left eye back
    (TZ24_2.iloc[148,13],TZ24_2.iloc[148,14]), # Left eye front
    (TZ24_2.iloc[148,15],TZ24_2.iloc[148,16]), # Right eyeball
    (TZ24_2.iloc[148,17],TZ24_2.iloc[148,18]), # Right eye back
    (TZ24_2.iloc[148,19],TZ24_2.iloc[148,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I6 = (TZ24_2.iloc[148,]==-1).values
row_vals = [t for t, x in enumerate(I6) if x]
if len(I6) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S6Land = np.delete(S6Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S6, rotation_vector_S6, translation_vector_S6, inliers_S6) = cv2.solvePnPRansac(Ref_Bird, S6Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S6, jacobian_S6) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
for p in S6Land:
	    cv2.circle(S6, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S6Land[0][0]), int(S6Land[0][1]))
p2 = ( int(bill_end_point2D_S6[0][0][0]), int(bill_end_point2D_S6[0][0][1]))
cv2.line(S6, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S6, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S6, translation_vector_S6, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S6Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# A few more.

S7 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040191.jpg")
S7Land = np.array([
    (TZ24_2.iloc[190,3],TZ24_2.iloc[190,4]), # Beak tip
    (TZ24_2.iloc[190,5],TZ24_2.iloc[190,6]), # Beak base
    (TZ24_2.iloc[190,7],TZ24_2.iloc[190,8]), # Cyr base
    (TZ24_2.iloc[190,9],TZ24_2.iloc[190,10]), # Left eyeball
    (TZ24_2.iloc[190,11],TZ24_2.iloc[190,12]), # Left eye back
    (TZ24_2.iloc[190,13],TZ24_2.iloc[190,14]), # Left eye front
    (TZ24_2.iloc[190,15],TZ24_2.iloc[190,16]), # Right eyeball
    (TZ24_2.iloc[190,17],TZ24_2.iloc[190,18]), # Right eye back
    (TZ24_2.iloc[190,19],TZ24_2.iloc[190,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I7 = (TZ24_2.iloc[190,]==-1).values
row_vals = [t for t, x in enumerate(I7) if x]
if len(I7) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S7Land = np.delete(S7Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S7, rotation_vector_S7, translation_vector_S7, inliers_S7) = cv2.solvePnPRansac(Ref_Bird, S7Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S7, jacobian_S7) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
for p in S7Land:
	    cv2.circle(S7, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S7Land[0][0]), int(S7Land[0][1]))
p2 = ( int(bill_end_point2D_S7[0][0][0]), int(bill_end_point2D_S7[0][0][1]))
cv2.line(S7, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S7, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S7, translation_vector_S7, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S7Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Great. 

S8 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040210.jpg")
S8Land = np.array([
    (TZ24_2.iloc[209,3],TZ24_2.iloc[209,4]), # Beak 
    (TZ24_2.iloc[209,5],TZ24_2.iloc[209,6]), # Beak base
    (TZ24_2.iloc[209,7],TZ24_2.iloc[209,8]), # Cyr base
    (TZ24_2.iloc[209,9],TZ24_2.iloc[209,10]), # Left eyeball
    (TZ24_2.iloc[209,11],TZ24_2.iloc[209,12]), # Left eye back
    (TZ24_2.iloc[209,13],TZ24_2.iloc[209,14]), # Left eye front
    (TZ24_2.iloc[209,15],TZ24_2.iloc[209,16]), # Right eyeball
    (TZ24_2.iloc[209,17],TZ24_2.iloc[209,18]), # Right eye back
    (TZ24_2.iloc[209,19],TZ24_2.iloc[209,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I8 = (TZ24_2.iloc[209,]==-1).values
row_vals = [t for t, x in enumerate(I8) if x]
if len(I8) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S8Land = np.delete(S8Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S8, rotation_vector_S8, translation_vector_S8, inliers_S8) = cv2.solvePnPRansac(Ref_Bird, S8Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S8, jacobian_S8) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S8, translation_vector_S8, Camera, Distortion_Coef)
for p in S8Land:
	    cv2.circle(S8, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S8Land[0][0]), int(S8Land[0][1]))
p2 = ( int(bill_end_point2D_S8[0][0][0]), int(bill_end_point2D_S8[0][0][1]))
cv2.line(S8, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S8, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S8, translation_vector_S8, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S8Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S8Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Acceptable.

S9 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040228.jpg")
S9Land = np.array([
    (TZ24_2.iloc[227,3],TZ24_2.iloc[227,4]), # Beak 
    (TZ24_2.iloc[227,5],TZ24_2.iloc[227,6]), # Beak base
    (TZ24_2.iloc[227,7],TZ24_2.iloc[227,8]), # Cyr base
    (TZ24_2.iloc[227,9],TZ24_2.iloc[227,10]), # Left eyeball
    (TZ24_2.iloc[227,11],TZ24_2.iloc[227,12]), # Left eye back
    (TZ24_2.iloc[227,13],TZ24_2.iloc[227,14]), # Left eye front
    (TZ24_2.iloc[227,15],TZ24_2.iloc[227,16]), # Right eyeball
    (TZ24_2.iloc[227,17],TZ24_2.iloc[227,18]), # Right eye back
    (TZ24_2.iloc[227,19],TZ24_2.iloc[227,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I9 = (TZ24_2.iloc[227,]==-1).values
row_vals = [t for t, x in enumerate(I9) if x]
if len(I9) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S9Land = np.delete(S9Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S9, rotation_vector_S9, translation_vector_S9, inliers_S9) = cv2.solvePnPRansac(Ref_Bird, S9Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S9, jacobian_S9) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S9, translation_vector_S9, Camera, Distortion_Coef)
for p in S9Land:
	    cv2.circle(S9, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S9Land[0][0]), int(S9Land[0][1]))
p2 = ( int(bill_end_point2D_S9[0][0][0]), int(bill_end_point2D_S9[0][0][1]))
cv2.line(S9, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S9, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S9, translation_vector_S9, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S9Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S9Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Good. 

S10 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040420.jpg")
S10Land = np.array([
    (TZ24_2.iloc[419,3],TZ24_2.iloc[419,4]), # Beak 
    (TZ24_2.iloc[419,5],TZ24_2.iloc[419,6]), # Beak base
    #(-1,-1),
    (TZ24_2.iloc[419,7],TZ24_2.iloc[419,8]), # Cyr base
    (TZ24_2.iloc[419,9],TZ24_2.iloc[419,10]), # Left eyeball
    (TZ24_2.iloc[419,11],TZ24_2.iloc[419,12]), # Left eye back
    (TZ24_2.iloc[419,13],TZ24_2.iloc[419,14]), # Left eye front
    (TZ24_2.iloc[419,15],TZ24_2.iloc[419,16]), # Right eyeball
    (TZ24_2.iloc[419,17],TZ24_2.iloc[419,18]), # Right eye back
    (TZ24_2.iloc[419,19],TZ24_2.iloc[419,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I10 = (TZ24_2.iloc[419,]==-1).values
row_vals = [t for t, x in enumerate(I10) if x]
if len(I10) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S10Land = np.delete(S10Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S10, rotation_vector_S10, translation_vector_S10, inliers_S10) = cv2.solvePnPRansac(Ref_Bird, S10Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S10, jacobian_S10) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S10, translation_vector_S10, Camera, Distortion_Coef)
for p in S10Land:
	    cv2.circle(S10, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S10Land[0][0]), int(S10Land[0][1]))
p2 = ( int(bill_end_point2D_S10[0][0][0]), int(bill_end_point2D_S10[0][0][1]))
cv2.line(S10, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S10, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S10, translation_vector_S10, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S10Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S10Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Again.

S11 = cv2.imread("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Images/TZ24_Control00040302.jpg")
S11Land = np.array([
    (TZ24_2.iloc[301,3],TZ24_2.iloc[301,4]), # Beak 
    (TZ24_2.iloc[301,5],TZ24_2.iloc[301,6]), # Beak base
    (TZ24_2.iloc[301,7],TZ24_2.iloc[301,8]), # Cyr base
    (TZ24_2.iloc[301,9],TZ24_2.iloc[301,10]), # Left eyeball
    (TZ24_2.iloc[301,11],TZ24_2.iloc[301,12]), # Left eye back
    (TZ24_2.iloc[301,13],TZ24_2.iloc[301,14]), # Left eye front
    (TZ24_2.iloc[301,15],TZ24_2.iloc[301,16]), # Right eyeball
    (TZ24_2.iloc[301,17],TZ24_2.iloc[301,18]), # Right eye back
    (TZ24_2.iloc[301,19],TZ24_2.iloc[301,20]) # Right eye front
    ], dtype = "double")
Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
], dtype = "double")
I11 = (TZ24_2.iloc[301,]==-1).values
row_vals = [t for t, x in enumerate(I11) if x]
if len(I11) > 0:
    to_drop = np.array([row_vals], dtype = "single")
    to_drop = np.unique(np.floor((to_drop-3)/2))
    for j in reversed(range(len(to_drop))):

        S11Land = np.delete(S11Land,int(to_drop[j]),0)
        Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
else: 
    print("All landmarks visible")
(success_S11, rotation_vector_S11, translation_vector_S11, inliers_S11) = cv2.solvePnPRansac(Ref_Bird, S11Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)

(bill_end_point2D_S11, jacobian_S11) = cv2.projectPoints(np.array([(-5.544, -4.368, 0.0)]), rotation_vector_S11, translation_vector_S11, Camera, Distortion_Coef)
for p in S11Land:
	    cv2.circle(S11, (int(p[0]), int(p[1])), 3, (0,0,0), -1)
p1 = ( int(S11Land[0][0]), int(S11Land[0][1]))
p2 = ( int(bill_end_point2D_S11[0][0][0]), int(bill_end_point2D_S11[0][0][1]))
cv2.line(S11, p1, p2, (0,0,0), 2)
plt.imshow(cv2.cvtColor(S11, cv2.COLOR_BGR2RGB))
plt.show()

(yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector_S11, translation_vector_S11, Camera, Distortion_Coef)
yhat = np.concatenate(yhat)
SEuc_D = np.sqrt(np.square(yhat-S11Land).sum(axis = 1)).sum(axis = 0)
AEuc_D = np.sqrt(np.square(yhat-S11Land).sum(axis = 1)).mean(axis = 0)
print(SEuc_D, AEuc_D)

# Looping through images.

TZ24_2 = pd.read_csv("/home/joshk/Desktop/Testing_Pigeons/TZ24_2/TZ24_Control.csv") 

for i in range(len(TZ24_2.index)):
    Land = np.array([
    (TZ24_2.iloc[i,3],TZ24_2.iloc[i,4]), # Beak tip
    (TZ24_2.iloc[i,5],TZ24_2.iloc[i,6]), # Beak base
    (TZ24_2.iloc[i,7],TZ24_2.iloc[i,8]), # Cyr base
    (TZ24_2.iloc[i,9],TZ24_2.iloc[i,10]), # Left eyeball
    (TZ24_2.iloc[i,11],TZ24_2.iloc[i,12]), # Left eye back
    (TZ24_2.iloc[i,13],TZ24_2.iloc[i,14]), # Left eye front
    (TZ24_2.iloc[i,15],TZ24_2.iloc[i,16]), # Right eyeball
    (TZ24_2.iloc[i,17],TZ24_2.iloc[i,18]), # Right eye back
    (TZ24_2.iloc[i,19],TZ24_2.iloc[i,20]) # Right eye front
    ], dtype = "double")
    Ref_Bird = np.array([
    (0,0,0), # Beak tip
    (19.991,5.712,0), # Beak base
    (13.439,18.816,0), # Cyr base
    (34.439,26.04,12.2), # Left eyeball
    (42.671,25.032,9.7), # Left eye back
    (26.543,20.496,9.7), # Left eye front
    (34.439,26.04,-12.2), # Right eyeball
    (42.671,25.032,-9.7), # Right eye back
    (26.543,20.496,-9.7) # Right eye front
    ], dtype = "double")
    Invisible = (TZ24_2.iloc[i,]==-1).values
    row_vals = [t for t, x in enumerate(Invisible) if x]
    if len(Invisible) > 0:
        to_drop = np.array([row_vals], dtype = "single")
        to_drop = np.unique(np.floor((to_drop-3)/2))
        for j in reversed(range(len(to_drop))):
            Land = np.delete(Land,int(to_drop[j]),0)
            Ref_Bird = np.delete(Ref_Bird,int(to_drop[j]),0)
    else: 
        print("All landmarks visible")
    (success, rotation_vector, translation_vector, inliers) = cv2.solvePnPRansac(Ref_Bird, Land, Camera, Distortion_Coef, flags = cv2.SOLVEPNP_EPNP)
    rmat, jac = cv2.Rodrigues(rotation_vector) # New additions to identify euler angles
    pyr=rotationMatrixToEulerAngles(rmat)*(180/math.pi)
    (yhat, jacobian) = cv2.projectPoints(Ref_Bird, rotation_vector, translation_vector, Camera, Distortion_Coef)
    yhat = np.concatenate(yhat) # Pulling residuals
    SEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).sum(axis = 0) # Calculating summed euclidean distance among residuals
    AEuc_D = np.sqrt(np.square(yhat-Land).sum(axis = 1)).mean(axis = 0) # Calculating average euclidean distance among residuals
    if i == 0:
        as_named_vector = np.concatenate((np.array([TZ24_2.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        hold = pd.DataFrame({"Slice" : [as_named_vector[0]],
            "R1" : [as_named_vector[1]],
            "R2" : [as_named_vector[2]],
            "R3" : [as_named_vector[3]],
            "P" : [as_named_vector[4]],
            "Y" : [as_named_vector[5]],
            "R" : [as_named_vector[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
        )
    else:
        as_named_vector_add = np.concatenate((np.array([TZ24_2.iloc[i,2]]),np.concatenate(rotation_vector),pyr))
        in_pd_add = pd.DataFrame({"Slice" : [as_named_vector_add[0]],
            "R1" : [as_named_vector_add[1]],
            "R2" : [as_named_vector_add[2]],
            "R3" : [as_named_vector_add[3]],
            "P" : [as_named_vector_add[4]],
            "Y" : [as_named_vector_add[5]],
            "R" : [as_named_vector_add[6]],
            "SEuc" : [SEuc_D],
            "MEuc" : [AEuc_D]}
            )
        binding = [hold,in_pd_add]
        hold = pd.concat(binding)

hold.to_csv('/home/joshk/Desktop/Testing_Pigeons/TZ24_2/Rotation_Out.csv')
