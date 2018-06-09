.POSIX:
EMACS = emacs

compile: huffman.elc

clean:
	rm -f huffman.elc

.SUFFIXES: .el .elc
.el.elc:
	$(EMACS) -Q -batch -f batch-byte-compile $<
