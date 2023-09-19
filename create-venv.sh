#!/bin/bash
set -eux

# sudo apt install python3-venv -y  #For Ubuntu,you have to run this command line first. 
python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r requirements.txt

