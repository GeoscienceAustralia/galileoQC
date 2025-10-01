#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
In development.
"""
"""
Created on Sat Jul 18 16:43:31 2020

author: Mark Dransfield

"""

# def gravity_from_curv(whizzfile, gne_chan, guv_chan, x=None, y=None, altitude=None):
"""

- call whizz_to_xarray (expanded to read more than 3 channels) to load data into ds.
- create area mask
- grid data with minc
- pad data with zeros
- 2D FFT of Craig function
- multiply by kernel
- inverse 2D FFT
- mask results
- analyse imag part and report
- return real part as gD

"""
# my_dataset = whizz_to_xarray(whizz_file, z_chan, *, n_chan='', e_chan='', 
# 	lines=[], remove_mean=False, diff_one=False, skipcontrols=False, controls=[]):
