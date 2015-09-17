#!/bin/sh

UCSCgenome="$1"
chromSizes="$2"

module load ucscUtils/262
fetchChromSizes $UCSCgenome > $chromSizes