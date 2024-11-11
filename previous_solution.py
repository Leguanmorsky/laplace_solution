import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
d=10.0  # velikost oblasti dxd
Nd=10  # pocet deleni ve smeru x a y
N=Nd*Nd  # pocet neznamych
dx=d/(Nd-1)  # krok site ve smeru x
dy=dx        # krok site ve smeru y

A=np.zeros([N,N])  # matice soustavy
B=np.zeros(N)    # vektor prave strany

def nastav_rovnici_i(i):
    if i>Nd-1 and i<N-Nd and i%Nd !=0 and (i-(Nd-1))%Nd !=0:
        A[i,i]=-4.0
    else: 
        A[i,i]=1.0
    if A[i,i]==-4.0:
        A[i,i-1]=1.0
        A[i,i-Nd]=1.0
        A[i,i+1]=1.0
        A[i,i+Nd]=1.0

def nastav_dirichlet_i(i,Phi_i):
    
    B[i]=Phi_i  # nastav pozadovany potencial do prave strany


def nastav_okrajovou_podminku(ix,iy):
    i=(iy)*Nd +ix
    if (ix==0):
        # leva hranice
        nastav_dirichlet_i(i,0.0)
        return
    if (ix==Nd-1):
        # prava hranice
        nastav_dirichlet_i(i,1.0)
        return
    if (iy==0):
        # dolni hranice
        nastav_dirichlet_i(i,0.0)
        return
    if (iy==Nd-1):
        # horni hranice
        nastav_dirichlet_i(i,1.0)
        return
    


# nastav koeficienty v matici A
for i in range(N):
    nastav_rovnici_i(i)

# nastav Dirichletovu podminku na leve a prave strane
for iy in range(Nd):
    for ix in range(0,Nd,Nd-1):
        nastav_okrajovou_podminku(ix,iy)

# nastav Dirichletovu podminku na dolni a horni strane
for iy in range(0,Nd,Nd-1):
    for ix in range(Nd):
        nastav_okrajovou_podminku(ix,iy)

# vyres soustavu rovnic
Phi=np.linalg.solve(A,B).reshape(Nd, Nd)  

# Create a single figure with subplots
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(1, 2, 1)  # 1 row, 2 columns, position 1

heatmap = ax1.imshow(Phi, aspect='auto', cmap='rainbow')
fig.colorbar(heatmap, ax=ax1)
ax1.set_aspect('equal')
ax1.set_title('Potential Heatmap')
ax1.set_xlabel('x')
ax1.set_ylabel('y')

x = np.arange(0, d+dx/2, dx)  # hodnoty x souradnic sitovych bodu
y = np.arange(0, d+dy/2, dy)
X, Y = np.meshgrid(x, y)

# Second subplot: 3D surface plot
ax2 = fig.add_subplot(1, 2, 2, projection='3d')  # 1 row, 2 columns, position 2
surface = ax2.plot_surface(X, Y, Phi, cmap='viridis')
fig.colorbar(surface, ax=ax2, shrink=0.5, aspect=5)
ax2.set_title('3D Potential Surface')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('Phi [V]')

plt.savefig('plot.png')
plt.show()

print("konec")