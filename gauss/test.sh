#!/bin/sh
skip=49

tab='	'
nl='
'
IFS=" $tab$nl"

umask=`umask`
umask 77

gztmpdir=
trap 'res=$?
  test -n "$gztmpdir" && rm -fr "$gztmpdir"
  (exit $res); exit $res
' 0 1 2 3 5 10 13 15

case $TMPDIR in
  / | /*/) ;;
  /*) TMPDIR=$TMPDIR/;;
  *) TMPDIR=/tmp/;;
esac
if type mktemp >/dev/null 2>&1; then
  gztmpdir=`mktemp -d "${TMPDIR}gztmpXXXXXXXXX"`
else
  gztmpdir=${TMPDIR}gztmp$$; mkdir $gztmpdir
fi || { (exit 127); exit 127; }

gztmp=$gztmpdir/$0
case $0 in
-* | */*'
') mkdir -p "$gztmp" && rm -r "$gztmp";;
*/*) gztmp=$gztmpdir/`basename "$0"`;;
esac || { (exit 127); exit 127; }

case `printf 'X\n' | tail -n +1 2>/dev/null` in
X) tail_n=-n;;
*) tail_n=;;
esac
if tail $tail_n +$skip <"$0" | gzip -cd > "$gztmp"; then
  umask $umask
  chmod 700 "$gztmp"
  (sleep 5; rm -fr "$gztmpdir") 2>/dev/null &
  "$gztmp" ${1+"$@"}; res=$?
else
  printf >&2 '%s\n' "Cannot decompress $0"
  (exit 127); res=127
fi; exit $res
����gtest.sh �T[KA~�_q:و�n6�b�`�Ԁ���C�6�I�e�kw'��lԦ��E�R��H���Q��e}�_��^4��󰛜��w��Ιܒ2�.ed��0i�y2 �Pb�0{O%�)�z8N�Rs0B D�B?��-[$[0 ������Q������*�(ʩ%�1a29<�F��ŰP�V0�>^=�-�������Ns����ej�7:_�O����Y�Yd}��~��x�~[o�7NV�s�W?�����pr��b��X 8}c������jt��w�o����ѱ��GC#1,M˦�iDKf��\�,I(3	������,E���ÞYXp$`�̇�<�뛝���'���Z���Y1��]�+N�|i�0O���ĄާV)�������b����r�$K��J)Q�eJ�����0>��Q��h$AD��-��"H�)1��Z�B3U#�췮ћ��>K�3Q�G�MR�{��F�4~���������b���µJ.tda9[ �m�2�;
y�Sf��/�M��akiű�����Ngm���Uk�Wj��_�
�I�u��ve�L$�6]�i�<+�+2B�Pr�0(�,�?5s~A�w��9���F���æi�Wf��L���f2�a�(��?�B٧���;3x��%=�{�>��`>
������C4B([2M��4�4WX6���c1X����`�y�4#������$T�+�v?a����58����م鮄�0�y�ݛ�WwFdľ`���e���  