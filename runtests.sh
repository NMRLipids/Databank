#!/bin/bash
# Run all the tests from here. Must work on both Linux and MacOS

# Determine the platform (Linux or MacOS)
OS_TYPE=$(uname)

if [[ "$OS_TYPE" == "Linux" || "$OS_TYPE" == "Darwin" ]]; then
  echo "Running tests on $OS_TYPE"
else
  echo "Unsupported operating system: $OS_TYPE"
  exit 1
fi

# Run the tests using pytest
if command -v pytest &> /dev/null; then
  pytest Scripts/tests --cmdopt sim2
  pytest Scripts/tests --cmdopt sim1
  pytest Scripts/tests --cmdopt nodata
else
  echo "pytest is not installed. Please install required libraries using"
  echo ">> pip install -e . -r Scripts/DatabankLib/requirements-dev.txt <<"
  exit 1
fi