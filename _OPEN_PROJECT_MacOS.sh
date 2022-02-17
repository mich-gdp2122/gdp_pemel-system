#!/bin/sh

# copy repo ssh keys to working dir
if [[ ! -e "~/.ssh" ]]; then
	mkdir "~/.ssh"
fi
cp -f ".ssh/id_rsa" "~/.ssh/id_rsa"
cp -f ".ssh/id_rsa.pub" "~/.ssh/id_rsa.pub"

# open matlab project
matlab -r "setenv('HOME', [getenv('USERPROFILE')]); open('GDPPEMELSystem.prj')" && disown || exit 1

exit 0