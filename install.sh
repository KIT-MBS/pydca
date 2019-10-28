#!/usr/bin/env bash

#----------------------------------------
# Installs the current version of pydca from PyPI 
# into the current user's home directory
# ---------------------------------------

MINIMUM_PYTHON_VERSION='3.5'
MINIMUM_PYTHON_VERSION_COMPARE=35
echo "Python 3 minimum version required: ${MINIMUM_PYTHON_VERSION}"
if command -v python3 &> /dev/null; then
	PYTHON_VERSION=$(python -c 'import sys; version=sys.version_info[:3]; print("{0}.{1}".format(*version))')
	PYTHON_VERSION_COMPARE=$(python -c 'import sys; version=sys.version_info[:3]; print("{0}{1}".format(*version))')
	if [ "${PYTHON_VERSION_COMPARE}" -lt "${MINIMUM_PYTHON_VERSION_COMPARE}" ]; then 
			echo "ERROR: You have older Python version: version=${PYTHON_VERSION}"
	else
		echo "Python 3 version found: ${PYTHON_VERSION}"
		echo 'installing pydca'
		pip install pydca
	fi
else 
	echo 'You need Python 3 version 3.5 or later to install pydca'
fi
