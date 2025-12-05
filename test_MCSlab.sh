#!/bin/bash
# Create a temporary directory
mkdir -p test_sandbox
cd test_sandbox

pwd

# run tests
../build/run_tests

# delete sandbox directory
cd ../
rm -rf test_sandbox

