#!/bin/sh
set -eu

development/pre-commit

R CMD Rd2pdf . -o varistran.pdf --force --no-preview
