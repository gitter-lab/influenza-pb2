#!/usr/bin/env python
import sys, argparse
import os, os.path
import numpy as np
import networkx as nx
import flopro.parsers.abc
from flopro import script_utils
import flow

# TODO dont assume targets file is already pre-processed to be the same size as targets-file?
# TODO check to assume all targets are in the network

def main():
  parser = argparse.ArgumentParser(description="""
Sub-sample source nodes for flow analysis from an alternative virus and compare node 
frequency to null (flu) virus results. --sources-file and --targets-file refer to the
null virus sources and targets.
""")
  flow.add_flow_args(parser)
  parser.add_argument('--alt-sources-file', required=True)
  parser.add_argument('--alt-targets-file', required=True)
  parser.add_argument('--n-simulation', required=True, help="Number of sub-samplings to perform")
  parser.add_argument('--local', action='store_true')
  parser.add_argument('--dry-run', action='store_true')
  args = parser.parse_args()
  script_utils.log_script(sys.argv)

  job_graph = nx.DiGraph()
  job_id = 0

  # sub-sample source lists from significant ICC-MS Hep C hits 
  attrs = {
    'exe': 'flow_sim_screens.py',
    'args': ['--n-simulation', args.n_simulation, '--sources-file', args.sources_file, '--edges-file', args.edges_file, '--outdir', args.outdir, '--alt-sources-file', args.alt_sources_file],
    'out': os.path.join(args.outdir, 'flow_sim_screens.out'),
    'err': os.path.join(args.outdir, 'flow_sim_screens.err'),
    'env': 'flu'
  }
  job_graph.add_node(job_id, **attrs)
  sim_screens_id = job_id
  job_id += 1
  
  # NOTE the following is the same procedure as flow_sim_pipeline.py except for use of --alt-targets-file

  # simulated runs of flow.py
  flow_result_fps = []
  for i in range(int(args.n_simulation)):
    flow_outdir = os.path.join(args.outdir, 'flow{}'.format(i))
    if not os.path.exists(flow_outdir):
      os.mkdir(flow_outdir) # TODO mkdir_p
    sim_fp = os.path.join(args.outdir, 'sim{}.txt'.format(i))
    attrs = {
      'exe': 'flow.py',
      'args': ['--sources-file', sim_fp, '--edges-file', args.edges_file, '--targets-file', args.alt_targets_file, '--outdir', flow_outdir, '--mapping-file', args.mapping_file, '--min-sources', str(args.min_sources), '--min-targets', str(args.min_targets), '--flow-only'],
      'out': os.path.join(flow_outdir, 'flow.out'),
      'err': os.path.join(flow_outdir, 'flow.err'),
      'env': 'flu'
    }
    job_graph.add_node(job_id, **attrs)
    job_graph.add_edge(sim_screens_id, job_id)
    flow_result_fps.append(os.path.join(flow_outdir, 'flow_result.graphml'))
    job_id += 1

  # compute node frequency
  attrs = {
    'exe': 'flow_sim_frequency.py',
    'args': ['--flow-results'] + flow_result_fps + ['--edges-file', args.edges_file, '--outdir', args.outdir],
    'out': os.path.join(args.outdir, 'flow_sim_frequency.out'),
    'err': os.path.join(args.outdir, 'flow_sim_frequency.err'),
    'env': 'flu'
  }
  freq_job_id = job_id
  job_graph.add_node(freq_job_id, **attrs)
  for i in range(int(args.n_simulation)):
    sim_flow_job_id = i+1
    job_graph.add_edge(sim_flow_job_id, freq_job_id)
  job_id += 1

  # do the real flow run
  flow_outdir = os.path.join(args.outdir, 'flow_real')
  if not os.path.exists(flow_outdir):
    os.mkdir(flow_outdir)
  attrs = {
    'exe': 'flow.py',
    'args': ['--sources-file', args.sources_file, '--edges-file', args.edges_file, '--targets-file', args.targets_file, '--outdir', flow_outdir, '--mapping-file', args.mapping_file, '--min-sources', str(args.min_sources), '--min-targets', str(args.min_targets), '--node-weights', os.path.join(args.outdir, 'node_frequency.csv')],
    'out': os.path.join(flow_outdir, 'flow.out'),
    'err': os.path.join(flow_outdir, 'flow.err'),
    'env': 'flu'
  }
  job_graph.add_node(job_id, **attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # do significance test
  node_frequency_fp = os.path.join(args.outdir, 'node_frequency.csv')
  flow_result_fp = os.path.join(flow_outdir, 'flow_result.graphml')
  attrs = {
    'exe': 'flow_sim_signif.py',
    'args': ['--flow-result', flow_result_fp, '--node-frequency', node_frequency_fp, '--outdir', args.outdir],
    'out': os.path.join(args.outdir, 'flow_sim_signif.out'),
    'err': os.path.join(args.outdir, 'flow_sim_signif.err'),
    'env': 'flu'
  }
  job_graph.add_node(job_id, **attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1

  # add p-value table
  attrs = {
    'exe':'write_pvals.py',
    'args': ['-q', node_frequency_fp, '-n', args.n_simulation, '-f', flow_result_fp, '-o', os.path.join(args.outdir, 'pvals.csv')],
    'out': os.path.join(args.outdir, 'write_pvals.out'),
    'err': os.path.join(args.outdir, 'write_pvals.err'),
    'env': 'flu'
  }
  job_graph.add_node(job_id, **attrs)
  job_graph.add_edge(job_id-1, job_id)
  job_id += 1
  
  condor = (not args.local)
  script_utils.run_digraph(args.outdir, job_graph, condor=condor, dry_run=args.dry_run)

if __name__ == "__main__":
  main()
