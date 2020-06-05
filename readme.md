# conti2d

![Download](https://img.shields.io/badge/Download-orange.svg)   [![download](https://img.shields.io/badge/binary-win-green.svg)](https://github.com/CosmicScholar/conti2d/archive/2.0.zip)   [![download](https://img.shields.io/badge/binary-OSX-green.svg)](https://github.com/CosmicScholar/conti2d)[![download](https://img.shields.io/badge/source-v2.0-green.svg)](https://github.com/CosmicScholar/conti2d)

![Download](https://img.shields.io/badge/Documents-orange.svg) [![download](https://img.shields.io/badge/Doxygen-doc-green.svg)](https://zguoch.github.io/conti2d/)

[![download](https://img.shields.io/badge/doi-10.1016/j.cageo.2020.104405-green.svg)](https://www.doi.org/10.1016/j.cageo.2020.104405)

## Citation
The bibtex entry of [10.1016/j.cageo.2020.104405](https://www.doi.org/10.1016/j.cageo.2020.104405) is listed below.

```bash
@article{Guo2020,
author = {Guo, Zhikui and Tao, Chunhui},
doi = {10.1016/j.cageo.2020.104405},
issn = {00983004},
journal = {Computers and Geosciences},
keywords = {Downward continuation,Fredholm integration,New kernel function,Potential field,Spatial domain,Surface to surface},
month = {mar},
publisher = {Elsevier Ltd},
title = {{Potential field continuation in spatial domain: A new kernel function and its numerical scheme}},
volume = {136},
year = {2020}
}

```

# Manual

[Tutorial video of GUI version](https://youtu.be/jQcnD3rnlkg)

## Hellow World

![](figures/snapshot/hello.png)

    The data and run scripts are in the example folder

## 1. Upward Continuation from Plane to Plane

```bash
field_path=../data
field_in=mag_plane_0
field_out=${field_in}_uwc_p2p.vtk
height_uwc=8
conti2d ${field_path}/${field_in}.grd -G${field_out} -H$height_uwc
```
![](figures/snapshot/uwc_p2p.png)

## 2. Upward Continuation from Plane to Surface

```bash
field_path=../data
field_in=mag_plane_0
field_out=${field_in}_uwc_p2s.vtk
topo_in=topo
conti2d ${field_path}/${field_in}.grd -G${field_out} -T0 -H${field_path}/${topo_in}.grd 
```
![](figures/snapshot/uwc_p2s.png)

## 3. Downward Continuation from Plane to Plane

```bash
field_path=../data
field_in=mag_plane_8
field_out=${field_in}_dwc_p2p.grd
height_dwc=8
conti2d ${field_path}/${field_in}.grd -G${field_out} -E5 -H0 -T$height_dwc -D+L1500  
```
![](figures/snapshot/dwc_p2p.png)

## 4. Downward Continuation from Surface to Plane

```bash
field_path=../data
field_in=mag_topo   
field_out=${field_in}_dwc_s2p.vtk
topofile=${field_path}/topo.grd
conti2d ${field_path}/${field_in}.grd -G${field_out} -E0 -H0 -T$topofile -D+L1500 
```
![](figures/snapshot/dwc_s2p.png)