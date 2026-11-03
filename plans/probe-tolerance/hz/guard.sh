#!/bin/sh
# Stops the extra runs, never the 384 one, if free memory falls below 15%.
while pgrep -f "runmore2.sh" > /dev/null; do
  free=$(memory_pressure | tail -1 | grep -o '[0-9]*')
  if [ "$free" -lt 15 ]; then echo "$(date +%H:%M) free $free%: stopping extra runs"; pkill -f "genfix3a.mjs [0-9] d256"; pkill -f "genfix3a.mjs [0-9] d128"; fi
  sleep 30
done
