#!/bin/sh

gnuplot  $1".gnu"
latex $1".tex"
dvips  $1".dvi" -o $1".eps"
rm $1".log" $1".dvi" $1".tex" $1"-inc.eps" $1."aux"
