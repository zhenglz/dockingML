#!/bin/bash

# reference
# https://help.github.com/articles/adding-a-file-to-a-repository-using-the-command-line/

git add --all
git commit -m "rewrite PdbIndex module"
git push -u origin v1.2.2
