import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
d=30.0 
Nd=30 
N=Nd*Nd  
dx=d/(Nd-1)  
dy=dx 

if Nd%2==0:
    p_par=Nd//2-0.5
else:
    p_par=Nd//2
    
R0=dx
R1=((Nd-3)//2)*dx
A=np.zeros([N,N])
B=np.zeros(N)

def nastav_rovnici_i(x,y):
    r=np.sqrt(x**2+y**2)
    i=int((y+p_par)*Nd+x+p_par)
    if r<=R1 and r>=R0:
        A[i, i] = -4.0
        A[i, i - 1] = 1.0
        A[i, i + 1] = 1.0
        A[i, i - Nd] = 1.0
        A[i, i + Nd] = 1.0
    else:
        A[i, i] = 1.0
        
    if (r<R0):
        B[i]=0.0
        return
    if (r>R1):
        B[i]=1.0
        return

# nastav koeficienty v matici A
for iy in range(int(-(Nd//2)), int((Nd//2)+1)):
    # if statement to handle even matrices!!!!
    if Nd%2==0:
        iy-=0.5
    for ix in range(int(-(Nd//2)), int((Nd//2)+1)):
        if Nd%2==0:
            ix-=0.5
        nastav_rovnici_i(ix,iy)
df = pd.DataFrame(A)
df.to_csv("C:/Users/jirka/Desktop/100 days challenge/Python VUT/Laplace/Laplace/matrix.csv", index=False)

sparse_matrix=csr_matrix(A)

# Propably most important row!!!
Phi_A = spsolve(sparse_matrix, B)
Phi=Phi_A.reshape(Nd, Nd)

x = np.arange(0, d + dx / 2, dx)
y = np.arange(0, d + dy / 2, dy)
X, Y = np.meshgrid(x, y)

# Create a single figure with subplots
fig = plt.figure(figsize=(12, 6))

# First subplot: 2D heatmap
ax1 = fig.add_subplot(1, 2, 1)  # 1 row, 2 columns, position 1
heatmap = ax1.imshow(Phi, aspect='auto', cmap='rainbow')
fig.colorbar(heatmap, ax=ax1)
ax1.set_aspect('equal')
ax1.set_title('Potential Heatmap')
ax1.set_xlabel('x')
ax1.set_ylabel('y')

# Second subplot: 3D surface plot
ax2 = fig.add_subplot(1, 2, 2, projection='3d')  # 1 row, 2 columns, position 2
surface = ax2.plot_surface(X, Y, Phi, cmap='viridis')
fig.colorbar(surface, ax=ax2, shrink=0.5, aspect=5)
ax2.set_title('3D Potential Surface')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('Phi [V]')

# Show the combined figure
plt.savefig('combined_plot.png')
plt.show()

print("konec")