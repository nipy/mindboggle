
make html

rsync -avz --sparse --exclude-from=/homedir/.rsync-exclude -e /usr/bin/ssh build/html/* binarybottle@binaryb\
ottle.com:/home/binarybottle/mindboggle.info
