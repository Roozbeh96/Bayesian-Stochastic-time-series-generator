
import numpy as np
import matplotlib.pyplot as plt

class Velocity_Field:

    def __init__(self, case, u, w, z,
                 delta, u_tau, Delx, nu, Structure_Len):
        
        self.type = case
        self.u = u
        self.w = w
        self.z = z
        self.delta = delta
        self.u_tau = u_tau
        self.Delx = Delx
        self.Delz = z[1][0]-z[0][0]
        self.Retau = u_tau*delta/nu
        self.lambdaT = self.Delx
        self.Structurelen = Structure_Len



    def __str__(self):
        return (f"{self.type} Velocity field size is {self.u.shape[1]}.\n " 
               f"length is from {self.z[0][0]} [m] to {self.z[-1][0]} [m].\n"
               f"value of \delta is {self.delta} [m].\n"
               f"value of friction velocity is {self.u_tau} [m/s].\n"
               f"value of \Delta x is {self.Delx} [m].\n"
               f"value of \Delta z is {self.Delz} [m].\n"
               f"value of Re_tau is {self.Retau}.\n"
               f"value of \lambda_T is {self.lambdaT}.\n")

    def Lambdaci(self):
        """consider velcoity fields lasts Structurelen*\delta"""



        
        L = int(np.floor(self.Structurelen *self.delta/self.Delx))
        numerator = int(np.floor(self.u.shape[1]/L))
        self.lambdaci = np.zeros((self.u.shape[0],L,numerator))
        self.dudx = np.zeros((self.u.shape[0],L,numerator))
        self.dudz = np.zeros((self.u.shape[0],L,numerator))
        self.dwdx = np.zeros((self.u.shape[0],L,numerator))
        self.dwdz = np.zeros((self.u.shape[0],L,numerator))


        for S in range(numerator):

            u = self.u[:, S*L:(S+1)*L]
            w = self.w[:, S*L:(S+1)*L]

            for r in range(1,u.shape[0]-1):

                for c in range(1,u.shape[1]-1):

                    dudx = (u[r, c+1]-u[r, c-1])/(2*self.Delx)
                    dudz = (u[r+1, c]-u[r-1, c])/(2*self.Delz)
                    dwdx = (w[r, c+1]-w[r, c-1])/(2*self.Delx)
                    dwdz = (w[r+1, c]-w[r-1, c])/(2*self.Delz)
                    Mat = np.array([[dudx, dudz],[dwdx, dwdz]])
                    Lambda = np.linalg.eigvals(Mat)
                    omega = dwdx - dudz
                    self.lambdaci[r, c, S] = np.unique(np.abs(np.imag(Lambda)))*np.sign(omega)
                    self.dudx[r, c, S] = dudx
                    self.dudz[r, c, S] = dudz
                    self.dwdx[r, c, S] = dwdx
                    self.dwdz[r, c, S] = dwdz


    def contour(self, Field, framenumb):
        try:

            if framenumb<0 or framenumb>Field.shape[2]:
                raise ValueError(f"Frame number {framenumb} does not exist.")

            X, Y = np.meshgrid(np.array(list(range(0,Field.shape[1])))*self.Delx, self.z)
            plt.figure(figsize=(12,6))
            plt.contourf(X, Y, Field[:,:,framenumb], levels = np.linspace(-200, 200, 75), cmap = 'seismic')
            plt.colorbar()
            plt.xlabel('x [m]')
            plt.ylabel('z [m]')
            plt.show()

        except ValueError as e:
            print(f"Error: {e}")

    def remap(self):


        zsize = 0.4*(self.lambdaT)
        z_grid_numb = int(np.floor(self.u.shape[0]/(1+np.floor(zsize/self.Delz))))+1
  

        L = int(np.floor(self.Structurelen *self.delta/self.Delx))
        numerator = int(np.floor(self.u.shape[1]/L))

        self.lambdaciremap = np.zeros((z_grid_numb,L,numerator))
        self.dudxremap = np.zeros((z_grid_numb,L,numerator))
        self.dudzremap = np.zeros((z_grid_numb,L,numerator))
        self.dwdxremap = np.zeros((z_grid_numb,L,numerator))
        self.dwdzremap = np.zeros((z_grid_numb,L,numerator))
        self.zremap = np.zeros(z_grid_numb)
        for i in range(0,z_grid_numb):
            self.zremap[i] = (i)*zsize + self.z[0] 



        for S in range(numerator):
 
            u = self.u[0::int(np.floor(zsize/self.Delz))+1, S*L:(S+1)*L]
            w = self.w[0::int(np.floor(zsize/self.Delz))+1, S*L:(S+1)*L]
            print(u.shape)
            for r in range(1,z_grid_numb-1):

                for c in range(1, u.shape[1]-1):

                    dudx = (u[r, c+1]-u[r, c-1])/(2*self.Delx)
                    dudz = (u[r+1, c]-u[r-1, c])/(2*zsize)
                    dwdx = (w[r, c+1]-w[r, c-1])/(2*self.Delx)
                    dwdz = (w[r+1, c]-w[r-1, c])/(2*zsize)
                    Mat = np.array([[dudx, dudz],[dwdx, dwdz]])
                    Lambda = np.linalg.eigvals(Mat)
                    omega = dwdx - dudz
                    self.lambdaciremap[r, c, S] = np.unique(np.abs(np.imag(Lambda)))*np.sign(omega)
                    self.dudxremap[r, c, S] = dudx
                    self.dudzremap[r, c, S] = dudz
                    self.dwdxremap[r, c, S] = dwdx
                    self.dwdzremap[r, c, S] = dwdz
 




