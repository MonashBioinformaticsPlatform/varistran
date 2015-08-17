#!/bin/sh
set -eu

R CMD Rd2pdf . -o varistran.pdf --force --no-preview
