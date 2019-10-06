
set bin_path=..\..\bin
set field_path=..\data
set field_in=mag_plane_0
set field_out=%field_in%_uwc_p2p.grd
set height_uwc=8
%bin_path%\conti2d_win.exe %field_path%/%field_in%.grd -G%field_out% -H%height_uwc%

pause