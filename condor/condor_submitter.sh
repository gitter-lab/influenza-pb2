#!/bin/bash

if [ -z "${CONDOR_OUT_DIR}" ]; then
  echo "usage: set CONDOR_OUT_DIR environment variable" >&2
  exit 1
fi

# try cmd as a filepath and as a name of something on the PATH
cmd="$(which "$1")"
if [ ! -x "$cmd" ]; then
  echo "usage: provide command as first argument" >&2
  exit 2
fi

args=""
for arg in "$@"; do
  if [ "$arg" != "$1" ]; then
    args="${args} ${arg}"
  fi
done
echo "args: $args"

condor_submit executable="$cmd" arguments="$args" ~/.condor/submitter.sub
