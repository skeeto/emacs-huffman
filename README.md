# Huffman Coding in Emacs Lisp

This library builds a Huffman coding tree from any sequence (list,
vector, string) of atoms and uses this tree to encode/decode a
sequence atoms to/from a run-length list of bits. The Huffman coding
tree is represented as a binary tree of cons cells.

Emacs is very, very slow at bit-twiddling and will never be effective
at compression/decompression in pure Emacs Lisp, so this library is
not written for speed.

## Examples

~~~el
(setf tree (huffman-tree "hello, world!"))
;; => (((?d . ?!) . ?l) (?o ?h . ?e) (?, . ?\s) ?w . ?r)

(huffman-encode-symbol tree ?l)
;; => (0 1)

(huffman-encode-symbol tree ?o)
;; => (1 0 0)

(huffman-encode-symbol tree ?w)
;; => (1 1 1 0)

(huffman-encode tree "how")
;; => (1 0 1 0 1 0 0 1 1 1 0)

(huffman-decode tree '(1 0 1 0 1 0 0 1 1 1 0) 'string)
;; => "how"
~~~
