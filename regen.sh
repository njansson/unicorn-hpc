#!/bin/sh


echo "Updating configuration..."
autoheader
echo "Running libtoolize"
libtoolize -i

echo "Running aclocal"
if which aclocal > /dev/null 2>&1; then
  aclocal -I m4 --install
else
  echo "No aclocal found on your system"
  exit 1
fi

echo "Running autoconf"
autoconf
echo "Running automake"
automake -a

echo "Deleting autom4te.cache directory"
rm -r autom4te.cache

