;replace path in expand_path command below
; by the path where you have stored your lib_idl folder

!PATH = '.:'+Expand_Path('+$HOME/lib_idl/proLEV') + ':' + !PATH
!PATH=!PATH+':'+Expand_Path('+$HOME/lib_idl/proJO')
!PATH=!PATH+':'+Expand_Path('+$HOME/lib_idl/proMAST')
;device,retain=2,decompose=0,true_color=24
;
;now, define those programs that shall be compiled at startup
;
;load constants in cgs units
const_cgs, /HELP
;
;load globally stored line transitions
line_transitions
