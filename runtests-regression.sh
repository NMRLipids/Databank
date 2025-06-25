#!/bin/bash
# Run all the tests from here including regression tests. Must work on both Linux and MacOS

# Determine the platform (Linux or MacOS)
OS_TYPE=$(uname)

if [[ "$OS_TYPE" == "Linux" || "$OS_TYPE" == "Darwin" ]]; then
  echo "Running tests on $OS_TYPE"
else
  echo "Unsupported operating system: $OS_TYPE"
  exit 1
fi

status=0
# Run the tests using pytest
if command -v pytest &> /dev/null; then
  pytest Scripts/tests --cmdopt sim2 || status=1
  pytest Scripts/tests --cmdopt sim1 || status=1
  pytest Scripts/tests --cmdopt adddata || status=1
  pytest Scripts/tests --cmdopt nodata || status=1
else
  echo "pytest is not installed. Please install required libraries using"
  echo ">> pip install -e . -r Scripts/DatabankLib/requirements-dev.txt <<"
  status=1
fi
exit $status