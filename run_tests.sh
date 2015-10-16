#!/bin/bash
python -m unittest discover -s tests -p "test_matrix.py" -v
python -m unittest discover -s tests -p "test_clusters.py" -v
python -m unittest discover -s tests -p "test_ebc.py" -v
