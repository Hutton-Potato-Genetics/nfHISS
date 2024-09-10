#!/bin/bash
java -jar -Xmx$1 "$(dirname "$0")/AgRenSeq_CreatePresenceMatrix.jar" "${@:2}"
