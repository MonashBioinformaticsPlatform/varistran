#!/bin/sh
set -eu

backyard/make_pdf.sh
sudo cp varistran.pdf /data/www/doc
