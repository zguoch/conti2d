field_path=../data
field_in=mag_plane_8
field_out=${field_in}_dwc_p2p.grd
height_dwc=8
conti2d ${field_path}/${field_in}.grd -G${field_out} -E5 -H0 -T$height_dwc -D+L150  