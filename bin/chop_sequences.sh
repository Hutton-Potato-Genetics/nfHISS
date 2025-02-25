#!/bin/bash
java $1 -jar "$(dirname "$0")/ChopSequence.jar" "${@:2}"
