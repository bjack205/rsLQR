#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
cd $SCRIPT_DIR/..
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
echo "Fetch"
git fetch
echo "Checkout"
git checkout gh-pages || exit 1
echo "Copy"
cp -r build/docs/html/* html || exit 1
sleep 0.5
echo "Add"
git add html/* || exit 1
echo "Commit"
git pull
git commit -m "Update documentation"
echo "Push"
git push -f origin gh-pages
echo "Checkout"
git checkout $CUR_BRANCH

cd $CUR_DIR

