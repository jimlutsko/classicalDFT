#!/bin/bash
qsub -l nodes=1:ppn=20 -l walltime=1:00:00  xRun.hydra
