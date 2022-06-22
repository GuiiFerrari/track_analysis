#!/bin/bash
latex $1".tex"
bibtex  $1
latex $1".tex"
latex $1".tex"
dvips   $1".dvi" -o $1".ps"
ps2pdf -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dAutoFilterGrayImages=false -dGrayImageFilter=/FlateEncode  $1".ps"
rm $1".ps" $1".log" $1".bbl" $1".aux" $1".blg" $1".dvi" $1".out" $1".len" $1".snm" $1".toc" $1".tex.backup"  $1".nav"

