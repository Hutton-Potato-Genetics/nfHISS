#!/bin/bash
MEME_XML="$(dirname "$0")/meme.xml"
MAST="$(which mast)"
java -jar "$(dirname "$0")/NLR-Parser3.jar" -y $MAST -x $MEME_XML "$@"
