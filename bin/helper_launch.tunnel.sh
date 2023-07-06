#!/bin/bash

PORT=$1

ssh -v -N -L "$PORT":127.0.0.1:"$PORT" \
rreggiar@courtyard.gi.ucsc.edu \
> tunnel_tmp.txt 2>&1 &

