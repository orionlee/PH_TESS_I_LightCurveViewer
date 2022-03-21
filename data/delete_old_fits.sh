#/bin/sh

find . -type f -name "*.fits" -mtime +365 -delete
