#!/bin/bash

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
    exit 1
}

if [ "$pam" != "NGG" ]; then
    error_json
fi

if [ "$sl" != "20" ]; then
    error_json
fi

if [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json
fi

if [ "$URL_CRISPR" = "" ]; then
    error_json
fi


sleep 10
printenv > end.log
echo "curl --noproxy '*' -X GET $URL_CRISPR/handshake" > handshake.cmd
curl --noproxy '*' -X GET $URL_CRISPR/handshake &> handshake.log

echo '{"test" : "done"}'