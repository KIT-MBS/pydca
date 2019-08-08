#!/usr/bin/bash

#----------------------------------------
# Installs the current version of pydca from PyPI 
# into the current user's home directory
# ---------------------------------------

MINIMUM_PYTHON_VERSION='3.5.0'
MINIMUM_PYTHON_VERSION_COMPARE=350
echo "python3 minimum version required: ${MINIMUM_PYTHON_VERSION}"
if command -v python3 &> /dev/null; then
	PYTHON_VERSION=$(python3 -c 'import sys; version=sys.version_info[:3]; print("{0}.{1}.{2}".format(*version))')
	PYTHON_VERSION_COMPARE=$(python3 -c 'import sys; version=sys.version_info[:3]; print("{0}{1}{2}".format(*version))')
	if [ "${PYTHON_VERSION_COMPARE}" -lt "${MINIMUM_PYTHON_VERSION_COMPARE}" ]; then 
			echo "You have older Python3: version=${PYTHON_VERSION}"
	else
		echo "Python3 version found: ${PYTHON_VERSION}"
		if python3 -c 'import pkgutil; exit(not pkgutil.find_loader("tkinter"))'; then
			echo 'tkinter is installed'
			echo 'installing pydca in current user'
			pip3 install pydca --user
		else
			echo 'You need to install python3-tk'
		fi
	fi
else 
	echo 'You need Python3 version 3.5 or later to install pydca'
fi
