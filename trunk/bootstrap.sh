#!/bin/sh
aclocal
automake --add-missing
autoconf
autoheader
aclocal
automake --add-missing
autoconf
autoheader
