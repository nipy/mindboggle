git checkout master -— doc/build/html/

mv _static shared

mv shared/new_index.html index.html

sed -i .bk 's/_static/shared/g' *.html
rm *.bk
