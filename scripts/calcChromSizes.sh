#!/bin/sh

UCSCgenome="$1"
chromSizes="$2"
script_path=$(dirname "${BASH_SOURCE[0]}")

$script_path/UCSC_kent_commands/fetchChromSizes $UCSCgenome > $chromSizes
