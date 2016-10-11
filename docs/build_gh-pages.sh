#!/bin/sh

# simple script to build and push to gh-pages
# designed to be run from a release tag or master on the verge of a release.

# this script will clone a new copy of the repo (f it doesn't exist)
# then checkpout the gh-pages branch, pull the latest, copy the generated
# docs to the copy, then add teh changes and push it all back to gitHub

# This script will create another copy of the repo, right next this
# one, but named "pyugrid.gh-pages"
# NOTE: the user runnign this will have to have permissions to push to teh gitHub repo
#       and have their sytem set up to not need an interactive password.

GHPAGESDIR=../../pyugrid.gh-pages/
REPO_URL=https://github.com/pyugrid/pyugrid.git

# make sure the gh-pages repo exists. If not, create it

if  [ -d $GHPAGESDIR ]
then
    echo "There is already a clone for gh-pages"
else
    git clone $REPO_URL $GHPAGESDIR
fi

# checkout the gh-pages branch
pushd $GHPAGESDIR
git checkout gh-pages
popd

# make the docs locally
make html
# copy to other repo (on the gh-pages branch)
cp -R _build/html/ $GHPAGESDIR

pushd $GHPAGESDIR
git add * # in case there are new files added
git commit -a -m "updating published docs"
git pull -s ours --no-edit
# git push
popd

