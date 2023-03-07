## This script can be called from anywhere on ceci cluster machines and will point to the
## right modules. You can also directly load ~/ceci_modules from the submission script.

#! /bin/bash

host=$(cat /etc/hostname)
echo "Hostname: $host"
host=${host%-*}
echo "Truncated: $host"

echo "Loading modules in ~/ceci_modules"
source ~/ceci_modules
