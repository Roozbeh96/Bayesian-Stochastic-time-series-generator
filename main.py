

'''Data reading'''


import pre_proc as data

'''Packages'''


import numpy as np
import matplotlib.pyplot as plt
import h5py
from VorXfeaturedVF import Velocity_Field



'''Object Creation'''



LSM_len = 3

Velocity_FieldGenWT7 = Velocity_Field('GenWT7',data.uGenWT7, data.wGenWT7, data.zGenWT7,
                                       data.deltaWT7, data.u_tauWT7, data.DelxGenWT7,
                                        data.nuWT7, LSM_len)
Velocity_FieldGenWT10 = Velocity_Field('GenWT10',data.uGenWT10, data.wGenWT10, data.zGenWT10,
                                        data.deltaWT10, data.u_tauWT10, data.DelxGenWT10,
                                        data.nuWT10, LSM_len)
Velocity_FieldGenASL = Velocity_Field('GenASL', data.uGenASL, data.wGenASL, data.zGenASL, 
                                      data.deltaASL, data.u_tauASL, data.DelxGenASL,
                                        data.nuASL, LSM_len)


print(Velocity_FieldGenWT7)
print(Velocity_FieldGenWT10)
print(Velocity_FieldGenASL)

'''Vortex Structures(Lambda_{Ci})'''





Velocity_FieldGenWT7.Lambdaci()
Velocity_FieldGenWT10.Lambdaci()
Velocity_FieldGenASL.Lambdaci()

'''Vortex Structures in remap velocity field'''

Velocity_FieldGenWT7.remap()
Velocity_FieldGenWT10.remap()
Velocity_FieldGenASL.remap()

'''countor plot'''


Velocity_FieldGenWT7.contour(Velocity_FieldGenWT7.lambdaci, 49)
Velocity_FieldGenWT10.contour(Velocity_FieldGenWT10.lambdaci, 50)
Velocity_FieldGenASL.contour(Velocity_FieldGenASL.lambdaci, 1)


'''Saving'''




for property, value in vars(Velocity_FieldGenWT10).items():
    print(property, ":", type(value))


with h5py.File('Velocity_FieldGenWT7.h5', 'w') as f:
    for attr, value in Velocity_FieldGenWT7.__dict__.items():
        f.create_dataset(attr, data=value)

with h5py.File('Velocity_FieldGenWT10.h5', 'w') as f:
    for attr, value in Velocity_FieldGenWT10.__dict__.items():
        f.create_dataset(attr, data=value)

with h5py.File('Velocity_FieldGenASL.h5', 'w') as f:
    for attr, value in Velocity_FieldGenASL.__dict__.items():
        f.create_dataset(attr, data=value)