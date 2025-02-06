#!/bin/bash
java $1 -jar "$(dirname "$0")/NLR-Annotator.jar" "${@:2}"
