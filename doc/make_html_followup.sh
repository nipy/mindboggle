

# MASTER:

git checkout master
cd doc
bash make_html.sh
git add build/html
git commit -m “update html"
git push origin master


# PAGES:

git checkout gh-pages
git checkout master -- doc/build/html

git mv doc/build/html/* .
git rm -rf doc
#git mv _static shared
git mv shared/new_index.html index.html

sed -i .bk 's/_static/shared/g' *.html
rm *.bk
cd users
sed -i .bk 's/_static/shared/g' *.html
rm *.bk
cd ../faq
sed -i .bk 's/_static/shared/g' *.html
rm *.bk
cd ..

ga *
gc -m "update html"
gpop

