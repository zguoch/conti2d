field_path=../data
field_in=mag_plane_0
field_out=${field_in}_uwc_p2s.vtk
topo_in=topo
conti2d ${field_path}/${field_in}.grd -G${field_out} -T0 -H${field_path}/${topo_in}.grd 
