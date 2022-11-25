#!/usr/bin/env python
from distutils.core import setup

setup(name="flopro",
  version="0.1",
  description="Flow program for protein-protein interaction networks",
  author="Aaron Baker",
  author_email="abaker@cs.wisc.edu",
  packages=['flopro', 'flopro.parsers'],
  scripts=[
    'scripts/flow_sim_pipeline.py',
    'scripts/flow_sim_screens.py',
    'scripts/flow.py',
    'scripts/flow_sim_frequency.py',
    'scripts/flow_sim_signif.py',
    'scripts/flow_alt_pipeline.py',
    'scripts/get_revision_no.sh',
    'scripts/conda_env_runner.sh',
    'scripts/map_ensp_to_hgnc.py',
    'scripts/write_pvals.py'
  ]
)
