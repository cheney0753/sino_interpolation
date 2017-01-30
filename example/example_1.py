# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:36:54 2017

@author: zhong

This is the test file for the SINT algorithm python implementation
"""

import astra
import numpy as np
import odl

sizeX = 300
thetas = np.array(range(0, 180, 10))*np.pi/180
thetas_new = np.array(range(0, 180, 5))*np.pi/180

reco_space = odl.uniform_discr(
    min_pt=[-20, -20], max_pt=[20, 20], shape=[sizeX, sizeX],
    dtype='float32')

# Create a discrete Shepp-Logan phantom (modified version)
phantom = odl.phantom.shepp_logan(reco_space, modified=True)

phantom_arr = phantom.asarray()

vol_geom = astra.create_vol_geom(sizeX, sizeX)

proj_geom = astra.create_proj_geom('parallel', 1, sizeX, thetas)

proj_id = astra.create_projector('line', proj_geom, vol_geom )

[sino_id, sino ]= astra.create_sino(phantom_arr, proj_id)

astra.projector.delete(proj_id)
astra.data2d.delete(sino_id)

import os
os.chdir('../')
import sint
cd('./example')

sino_estimated = sint(sino, thetas_new)