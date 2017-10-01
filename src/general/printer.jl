# ---------------------------------------------------------------------------- #
#
#   printer.jl
#
#   Type for printer
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Printer

Printer type:
Holds name of 3D printer, manufacturing constraints, etc.
"""
type Printer

  name::String     # Name of printer



end

function Printer( )

    Printer("Printer 1")

end
