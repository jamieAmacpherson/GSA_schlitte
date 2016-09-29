#!/bin/bash
if [ $# -ne 1 ]
then
    echo "Incorrect number of arguments..."
    echo "Usage: git_startup.sh <remote host>"
    exit 1
fi


git init .
git add .
git commit -a

git remote add origin $1

git push -u origin master
