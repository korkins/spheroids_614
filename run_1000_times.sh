#!/bin/bash

echo running...
for i in {1..1000}
do
  ./run_sdls > /dev/null
done