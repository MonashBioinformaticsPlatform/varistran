#!/bin/sh
set -eu

development/make_pdf.sh
sudo cp varistran.pdf /data/www/doc
