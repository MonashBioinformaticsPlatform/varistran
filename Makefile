

check :
	Rscript -e "devtools::check()"


document :
	Rscript -e "devtools::document()"


paper.pdf : paper.md paper.bib
	pandoc --filter pandoc-citeproc -o paper.pdf paper.md


#
# Release procedure
#
# 1. Update DESCRIPTION to new version number.
#
# 2. make check
#    To update documentation and check for problems.
#
# 3. git commit
#
# 4. git tag vX.X.X
#
# 5. git push --tags
#
