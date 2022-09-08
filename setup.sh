#!/bin/bash
conda create --name FunUniFrac python==3.10.4
activate FunUniFrac
pip install phylodm
pip install pyemd
conda install -y pandas
conda install -y networkx
