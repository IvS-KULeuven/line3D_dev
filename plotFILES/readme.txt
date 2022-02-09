levpy/
   contains the python library for postprocessing data

idl_lib/
   contains the idl library for postprocessing data

gdl_lib/
   contains the gdl library for postprocessing data


#to automatically load the libraries to gdl or idl
(please double check the startup.pro files and corresponding paths)
export GDL_STARTUP=path_to_my_gdl_lib/startup.pro
export IDL_STARTUP=path_to_my_idl_lib/startup.pro

#to automatically include the python libraries
LEVPYPATH="path_to_levpy"
export PYTHONPATH="${PYTHONPATH}:$LEVPYPATH"
LEVPYPATH="path_to_levpy/version_sc3d"
export PYTHONPATH="${PYTHONPATH}:$LEVPYPATH"
