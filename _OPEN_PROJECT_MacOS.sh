#!/bin/sh

# copy repo ssh keys to working dir
if [[ ! -e "~/.ssh" ]]; then
	mkdir "~/.ssh"
fi
if [[ ! -e "~/.ssh/id_rsa" ]]; then
	cp ".ssh/id_rsa" "~/.ssh/id_rsa"
fi
if [[ ! -e "~/.ssh/id_rsa.pub" ]]; then
	cp ".ssh/id_rsa.pub" "~/.ssh/id_rsa.pub"
fi

# open matlab project
matlab -r "setenv('HOME', [getenv('USERPROFILE')]); open('GDPPEMELSystem.prj')" && disown || exit 1

exit 0