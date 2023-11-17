import matplotlib.pyplot as plt
import numpy as np
import camera_functions

# Load data from .npy file
data = np.load('hw2.npy',allow_pickle=True).item()
verts3d = data['verts3d']
vcolors = data['vcolors']
faces = data['faces']
c_org = data['c_org']
c_lookat = data['c_lookat']
c_up = data['c_up']
t_1 = data['t_1']
t_2 = data['t_2']
u = data['u']
phi = data['phi']

# Specify some variables
img_h = 512
img_w = 512
cam_h = 15
cam_w = 15
f = 70

# STEP ZERO
print('Loading first image...')
# transpose vertices matrix for functionality
verts3d = verts3d.T
image = camera_functions.render_object(verts3d,faces,vcolors,img_h,img_w,cam_h,cam_w,f,c_org,c_lookat,c_up)
print('First image saved!\n\n')

fig = plt.figure(0)
plt.xticks([])
plt.yticks([])

plt.imshow(image)
plt.savefig(('0.jpg'))
# plt.show()

# STEP ONE
print('Loading second image...')
# Apply shift
verts3d = camera_functions.affine_transform(verts3d,None,None,t_1)
image = camera_functions.render_object(verts3d,faces,vcolors,img_h,img_w,cam_h,cam_w,f,c_org,c_lookat,c_up)
print('Second image saved!\n\n')

fig = plt.figure(1)
plt.xticks([])
plt.yticks([])

plt.imshow(image)
plt.savefig(('1.jpg'))
# plt.show()

# STEP TWO
print('Loading third image...')
# Apply rotation
verts3d = camera_functions.affine_transform(verts3d,phi,u,None)
image = camera_functions.render_object(verts3d,faces,vcolors,img_h,img_w,cam_h,cam_w,f,c_org,c_lookat,c_up)
print('Third image saved!\n\n')

fig = plt.figure(2)
plt.xticks([])
plt.yticks([])

plt.imshow(image)
plt.savefig(('2.jpg'))
# plt.show()

# STEP THREE
print('Loading fourth image...')
# Apply shift for second time
verts3d = camera_functions.affine_transform(verts3d,None,None,t_2)
image = camera_functions.render_object(verts3d,faces,vcolors,img_h,img_w,cam_h,cam_w,f,c_org,c_lookat,c_up)
print('Fourth image saved!\n\n')

fig = plt.figure(3)
plt.xticks([])
plt.yticks([])

plt.imshow(image)
plt.savefig(('3.jpg'))
# plt.show()
