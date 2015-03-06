# Commands for building the api docs for the mindboggle website.

# MASTER BRANCH:

git checkout master
cd doc
make html
cd ..
git add doc/build/html/api 
git add doc/build/html/genindex.html
git add doc/build/html/documentation.html
git commit -m "Update software documentation (api) for website."
git push origin master

# GH-PAGES BRANCH:

git checkout gh-pages
git checkout master -- doc/build/html/api
git checkout master -- doc/build/html/genindex.html
git checkout master -- doc/build/html/documentation.html

rm -r api
mv doc/build/html/* .

sed -i .bk 's/_static/shared/g' *.html
sed -i .bk 's/_static/shared/g' api/generated/*.html
rm *.bk

git add api genindex.html documentation.html
git commit -m "Update api and faq for website."
git push origin gh-pages
