;;; huffman.el --- huffman coding library -*- lexical-binding: t; -*-

;; This is free and unencumbered software released into the public domain.

;;; Commentary:

;; This library builds a Huffman coding tree from any sequence (list,
;; vector, string) of atoms and uses this tree to encode/decode a
;; sequence atoms to/from a run-length list of bits. The Huffman
;; coding tree is represented as a binary tree of cons cells.

;; Emacs is very, very slow at bit-twiddling and will never be
;; effective at compression/decompression in pure Emacs Lisp, so this
;; library is not written for speed.

;; Example usage:

;; (setf tree (huffman-tree "hello, world!"))
;; => (((?d . ?!) . ?l) (?o ?h . ?e) (?, . ?\s) ?w . ?r)

;; (huffman-encode-symbol tree ?l)
;; => (0 1)

;; (huffman-encode-symbol tree ?o)
;; => (1 0 0)

;; (huffman-encode-symbol tree ?w)
;; => (1 1 1 0)

;; (huffman-encode tree "how")
;; => (1 0 1 0 1 0 0 1 1 1 0)

;; (huffman-decode tree '(1 0 1 0 1 0 0 1 1 1 0) 'string)
;; => "how"

;;; Code:

(require 'cl-lib)

(defun huffman--queue-create (&optional sequence)
  "Create a new queue object from SEQUENCE."
  (let ((copy (append sequence nil)))
    (cons copy (last copy))))

(defun huffman--queue-fake-length (queue)
  "Returns 0, 1, or 2 if there are two or more elements in QUEUE."
  (cond ((null (car queue)) 0)
        ((eq (car queue) (cdr queue)) 1)
        (2)))

(defun huffman--queue-push (queue e)
  "Enqueue E at the end of QUEUE, returning E."
  (prog1 e
    (if (car queue)
        (setf (cdr (cdr queue)) (list e)
              (cdr queue) (cdr (cdr queue)))
      (setf (car queue) (cons e (car queue))
            (cdr queue) (car queue)))))

(defun huffman--queue-peek (queue)
  "Return but don't remove the front element of QUEUE, or nil if empty."
  (caar queue))

(defun huffman--queue-pop (queue)
  "Remove and return the front element of QUEUE, or nil if empty."
  (prog1 (huffman--queue-peek queue)
    (if (eq (car queue) (cdr queue))
        (setf (car queue) nil
              (cdr queue) nil)
      (setf (car queue) (cdr (car queue))))))

(defun huffman--count (sequence)
  "Return a sorted probability list of the elements of SEQUENCE."
  (let ((table (make-hash-table :test 'eql)))
    (cl-loop for symbol elements of sequence
             do (cl-incf (gethash symbol table 0)))
    (cl-loop for symbol hash-keys in table using (hash-values count)
             collect (cons symbol count) into results
             finally (cl-return (cl-sort results #'< :key #'cdr)))))

(defun huffman--next (queue-a queue-b)
  "Pop the lowest element from QUEUE-A or QUEUE-B, using QUEUE-A on ties."
  (let ((a (huffman--queue-peek queue-a))
        (b (huffman--queue-peek queue-b)))
    (cond ((null a) (huffman--queue-pop queue-b))
          ((null b) (huffman--queue-pop queue-a))
          ((< (cdr b) (cdr a)) (huffman--queue-pop queue-b))
          ((huffman--queue-pop queue-a)))))

(cl-defun huffman-tree (sequence &key canonical)
  "Compute and return the Huffman tree for SEQUENCE.
Given a sorting predicate to :canonical, return the canonical
tree. A canonical Huffman tree can be more compactly serialized."
  (let ((queue-a (huffman--queue-create (huffman--count sequence)))
        (queue-b (huffman--queue-create)))
    (while (> (+ (huffman--queue-fake-length queue-a)
                 (huffman--queue-fake-length queue-b))
              1)
      (cl-destructuring-bind (s0 . p0) (huffman--next queue-a queue-b)
        (cl-destructuring-bind (s1 . p1) (huffman--next queue-a queue-b)
          (huffman--queue-push queue-b `((,s0 . ,s1) . ,(+ p0 p1))))))
    (let ((tree (car (huffman--queue-pop queue-b))))
      (if canonical
          (huffman--canonicalize tree canonical)
        tree))))

(defun huffman-encode-symbol (tree symbol &optional prefix)
  "Uring Huffman TREE, return the bit encoding for SYMBOL."
  (cl-destructuring-bind (b0 . b1) tree
    (cond ((eql symbol b0) (nconc prefix (list 0)))
          ((eql symbol b1) (nconc prefix (list 1)))
          ((or (and (consp b0) (huffman-encode-symbol
                                b0 symbol (append prefix (list 0))))
               (and (consp b1) (huffman-encode-symbol
                                b1 symbol (append prefix (list 1)))))))))

(defun huffman-decode-symbol (tree code)
  "Using Huffman TREE, return (symbol . unused-bits) from the front of CODE."
  (cl-labels ((select (tree b) (if (= b 0) (car tree) (cdr tree))))
    (cl-destructuring-bind (b . rest) code
      (let ((next (select tree b)))
        (if (consp next)
            (huffman-decode-symbol (select tree b) rest)
          (cons next rest))))))

(defun huffman-encode (tree sequence)
  "Encode SEQUENCE into a bit sequence using TREE."
  (let ((table (make-hash-table :test 'eql)))
    (cl-loop for (node . bits) in (huffman--flatten tree)
             do (setf (gethash node table) bits))
    (cl-loop for symbol elements of sequence
             nconc (gethash symbol table))))

(cl-defun huffman-decode (tree bits &optional (type 'list))
  "Decode BITS into a TYPE sequence of symbols using TREE.
The TYPE argument is the same as the last argument to `cl-coerce'."
  (cl-loop while bits
           for (symbol . rest) = (huffman-decode-symbol tree bits)
           do (setf bits rest)
           collect symbol into symbols
           finally (cl-return (cl-coerce symbols type))))

;; Tree Transforms:

(defun huffman--listing-sort (listing predicate)
  "Sort LISTING according to code length, then PREDICATE."
  (let ((sorted (cl-sort listing predicate :key #'car)))
    (cl-stable-sort sorted #'< :key #'length)))

(defun huffman--flatten (tree &optional predicate)
  "Return listing form of TREE."
  (let ((listing ()))
    (cl-labels ((recur (node prefix)
                   (if (atom node)
                       (push (cons node (reverse prefix)) listing)
                     (recur (car node) (cons 0 prefix))
                     (recur (cdr node) (cons 1 prefix)))))
      (recur tree ())
      (if predicate
          (huffman--listing-sort listing predicate)
        (nreverse listing)))))

(defun huffman--unflatten (listing)
  "Transform LISTING back into a Huffman tree."
  (if (= 1 (length listing))
      (caar listing)
    (let ((b0 ())
          (b1 ()))
      (dolist (item listing)
        (cl-destructuring-bind (node . bits) item
          (cl-ecase (car bits)
            (0 (push (cons node (cdr bits)) b0))
            (1 (push (cons node (cdr bits)) b1)))))
      (cons (huffman--unflatten b0) (huffman--unflatten b1)))))

(defun huffman--increment (bits)
  "Increment the value of BITS, least-significant bit in first."
  (cl-case (car bits)
    (0 (cons 1 (cdr bits)))
    (1 (cons 0 (huffman--increment (cdr bits))))
    (t (list 1))))

(defun huffman--extend (bits length)
  "Add 0's to front of BITS so that it is LENGTH bits long."
  (let ((actual (length bits)))
    (if (< actual length)
        (nconc (make-list (- length actual) 0) bits)
      bits)))

(defun huffman--canonicalize (tree predicate)
  "Return canonical listing of TREE according to sort by PREDICATE."
  (let* ((flat (huffman--flatten tree  predicate)))
    (cl-loop for new-bits = () then (huffman--increment new-bits)
             for (node . bits) in flat
             do (setf new-bits (huffman--extend new-bits (length bits)))
             collect (cons node (reverse new-bits)) into listing
             finally (cl-return (huffman--unflatten listing)))))

;; Tree Analysis:

(cl-defun huffman-to-graphviz (tree file &key print-char)
  "Write description of TREE in Graphviz dot format to FILE."
  (let ((print-escape-newlines t)
        (font-node "Droid Sans")
        (font-edge "Droid Mono"))
    (cl-macrolet ((f (s &rest args) `(insert (format ,s ,@args))))
      (with-temp-file file
        (f "graph {\n")
        (f "  node[fontname=\"%s\",fontsize=\"16\"];\n" font-node)
        (f "  edge[fontname=\"%s\",fontsize=\"12\",labeldistance=2.5];\n"
           font-edge)
        (cl-labels
            ((recur (node name)
               (if (consp node)
                   (cl-destructuring-bind (b0 . b1) node
                     (let ((name0 (format "%s0" name))
                           (name1 (format "%s1" name)))
                       (f "  %s[shape=point];\n" name)
                       (f "  %s -- %s" name name0)
                       (f "[headlabel=0,labelangle=20];\n")
                       (f "  %s -- %s " name name1)
                       (f "[headlabel=1,labelangle=-20];\n")
                       (recur b0 name0)
                       (recur b1 name1)))
                 (f "  %s[label=%S,shape=square];\n" name
                    (cond ((symbolp node) (symbol-name node))
                          ((and (integerp node) print-char) (format "%c" node))
                          ((format "%s" node)))))))
          (recur tree "b"))
        (insert "}\n")))))

;;; huffman.el ends here
