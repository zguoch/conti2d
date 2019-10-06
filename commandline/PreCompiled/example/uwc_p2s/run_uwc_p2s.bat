
set bin_path=..\..\bin
set field_path=..\data
set field_in=mag_plane_0
set field_out=%field_in%_uwc_p2s.grd
set topo_in=topo
set height_uwc=8
%bin_path%\conti2d_win.exe %field_path%/%field_in%.grd -G%field_out% -H%field_path%\%topo_in%.grd 

pause