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
�P�htest.sh �T]O�V�?���T��$�LA�Z�X�J���|B<9v;0)�4A�tCM���aK�Ti�ieğ��䊿���M(��ܜ��y?��y�	݊$U=���4	��)��a`����5������fa������O>53 �@�,������:\t.m�`����C�qb���Ν?w��N�a7J���/��Z#)���Ѹ8=5�� :.��XQ ����߬��-����D2@B���w��m�%w�?l��"s��q���v=أsPj��B�pw��S=�ϗY���/_�� .D�rN�4�%ҋ�y9o���]#�*��`��P�6>Z�qA���Nu)x"� �4��{�ܪ��Jws��gvm�����|o��1�����~zN_�.�`�<�t�E���(���S�:�kX+&��U�<^q˫Ne��z�<="�H)�Q�R��'�zS��S*7Ojͳ��1��-�9#GM�zN�1����5�I$��P���e,g�z\,pzh7C�����T>�Q-�*�i`Q��0>:"Oʛh�8�#�ž�v�HG����ƫ�P4��o�,5C���N��I�x���Es������dͤ�,�U���n�en�f������R��ປ�Y�j�~b__}D?Q�v�ٵ�h�776�5FV.O�G���[Ǭ��?��2좱������:�^�;?B/���X��\���]ߝ�{�����]�����Q��տă�>�2[�'��*;'[�D��7���l.��&�/H�q�U��y]���r9#wm4�D{�LՎ�b!HTzx�_�J�U?q~9��x����Z�C�lө���+xo��/��_��T~�⌰�3���n%�K(�"[�ᯥpF
+R���R
O�F�1�`<ؽ�"L�s;���4�(�g96�I-Y�@�!������B"z^��c��|h�:����Ƭ"�.%X<�	�^t3��ևu�>ru��##�8PDx��}��  