field_path=../data
field_in=mag_plane_0
field_out=${field_in}_uwc_p2p.vtk
height_uwc=8
conti2d ${field_path}/${field_in}.grd -G${field_out} -H$height_uwc