#!/bin/bash

FILE_OR_PATH="$1"
if [ -z "$FILE_OR_PATH" ]; then
  echo "usage: get_revision_no.sh <FILE>" >&2
  exit 1
fi
PATHNAME="$(readlink -e "$FILE_OR_PATH")"
DIRNAME="$(dirname $PATHNAME)"
cd "$DIRNAME"

# TODO could check if dirname is actually tracked by git

echo "Revision: $(git rev-parse HEAD)"

# I considered erroring if there are uncommitted changes but that seems too strict
# instead, report uncommitted changes so that there is at least a log
echo "Uncommitted changes:"
git status -s
