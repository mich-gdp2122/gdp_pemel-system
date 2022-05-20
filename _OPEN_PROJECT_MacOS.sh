#!/bin/zsh

HERE="${0:h}"

# open matlab project
/Applications/MATLAB_R2022a.app/bin/matlab -r "open ('$HERE/GDPPEMELSystem.prj')" ; disown

exit 0