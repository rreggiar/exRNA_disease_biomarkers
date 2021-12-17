#!/bin/bash
ssh -v -N -i /Users/romanreggiardo/.ssh/id_rsa_pl -L 3838:127.0.0.1:3838 rreggiar@plaza.gi.ucsc.edu > docker.out.tmp 2>&1 &

