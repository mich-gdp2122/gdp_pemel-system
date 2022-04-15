#!/bin/zsh

HERE="${0:h}"

# copy repo ssh keys to working dir
if [[ ! -e ~/".ssh/" ]]; then
    mkdir ~/".ssh/"
fi
cp -f "$HERE/.ssh/id_rsa" ~/".ssh/id_rsa"
cp -f "$HERE/.ssh/id_rsa.pub" ~/".ssh/id_rsa.pub"

# open matlab project
/Applications/MATLAB_R2022a.app/bin/matlab -r "open ('$HERE/GDPPEMELSystem.prj')" ; disown

exit 0