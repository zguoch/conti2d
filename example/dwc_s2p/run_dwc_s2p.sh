field_path=../data
field_in=mag_topo   
field_out=${field_in}_dwc_s2p.vtk
topofile=${field_path}/topo.grd
conti2d ${field_path}/${field_in}.grd -G${field_out} -E0 -H0 -T$topofile -D+L100 
