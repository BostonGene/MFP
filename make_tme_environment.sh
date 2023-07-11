#!/usr/bin/env bash
set -e

VENV_NAME="TME_env"

if [ -d "$VENV_NAME" ]
then
    echo "Directory $VENV_NAME already exist."
    echo "If you want to reinstall environment, execute the commands below and then restart script:"
    echo "rm -rf $VENV_NAME" | tr '[:upper:]' '[:lower:]'
    exit 1
fi

echo "Create new virtual environment with name: '$VENV_NAME'"
python3.10 -m venv $VENV_NAME

echo "Enter virtual environment"
source $VENV_NAME/bin/activate

echo "Install packages"
pip install --upgrade pip wheel --no-cache-dir
pip install -r requirements.txt --no-cache-dir

echo "Create jupyter kernel with name $VENV_NAME" | tr '[:upper:]' '[:lower:]'
python -m ipykernel install --user --name="$VENV_NAME" 
