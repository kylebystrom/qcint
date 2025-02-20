;;; Qcint is a general GTO integral library for computational chemistry
;;; Copyright (C) 2014- Qiming Sun <osirpt.sun@gmail.com>
;;;
;;; This file is part of Qcint.
;;;
;;; Qcint is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

(load "utility.cl")
(load "parser.cl")
(load "derivator.cl")

; TODO
;
;(defun gen-raw-output (foutput vs1 &optional vs2)
;  (format foutput "~%e1:")
;  (format foutput "~a" vs1)
;  (format foutput "~%e2:")
;  (format foutput "~a" vs2))
;
;(defun gen-tex-subscript (foutput vs1 &optional vs2)
;  (defun tex-subscript-filter (cell)
;    (format foutput "~a,~{~a~},~{~a~} "
;            (phase-of cell)
;            (scripts-of (consts-of cell))
;            (scripts-of (ops-of cell))))
;  (gen-subscript (set-cells-streamer #'tex-subscript-filter foutput)
;                 (tensor-of-cellss vs1 vs2)))
;
;(defun gen-ternary-subscript (foutput vs1 &optional vs2)
;  (defun ternary-filter (cell)
;    (format foutput "~a,~a,~a "
;            (phase-of cell)
;            (ternary-subscript (consts-of cell))
;            (ternary-subscript (ops-of cell))))
;  (gen-subscript (set-cells-streamer #'ternary-filter foutput)
;                 (tensor-of-cellss vs1 vs2)))

(defun gen-subscript (cells-streamer raw-script)
  (labels ((gen-tex-iter (raw-script)
             (cond ((null raw-script) raw-script)
                   ((vector? raw-script)
                    (list (gen-tex-iter (comp-x raw-script))
                          (gen-tex-iter (comp-y raw-script))
                          (gen-tex-iter (comp-z raw-script))))
                   ((cells? raw-script)
                    (funcall cells-streamer raw-script))
                   (t (mapcar cells-streamer raw-script)))))
    (gen-tex-iter raw-script)))
(defun flatten-raw-script (raw-script)
  (let ((terms ()))
    (gen-subscript (lambda (x) (push x terms)) raw-script)
    (reverse terms)))

(defun convert-from-n-sys (ls n)
  (reduce (lambda (x y) (+ (* x n) y)) ls
          :initial-value 0))

(defun xyz-to-ternary (xyzs)
  (cond ((eql xyzs 'x) 0)
        ((eql xyzs 'y) 1)
        ((eql xyzs 'z) 2)
        (t (error " unknown subscript ~a" xyzs))))

(defun ternary-subscript (ops)
  "convert the polynomial xyz to the ternary"
  (cond ((null ops) ops)
        (t (convert-from-n-sys (mapcar #'xyz-to-ternary 
                                       (remove-if (lambda (x) (eql x 's))
                                                  (scripts-of ops)))
                               3))))
(defun cell-converter (cell fout)
  (let ((fac (realpart (phase-of cell)))
        (const@3 (ternary-subscript (consts-of cell)))
        (op@3    (ternary-subscript (ops-of cell))))
    (cond ((equal fac 1)
           (cond ((null const@3)
                  (if (null op@3)
                    (format fout " + s\[0\]" )
                    (format fout " + s\[~a\]" op@3)))
                 ((null op@3)
                  (format fout " + c\[~a\]*s\[0\]" const@3))
                 (t (format fout " + c\[~a\]*s\[~a\]" const@3 op@3))))
          ((equal fac -1)
           (cond ((null const@3)
                  (if (null op@3)
                    (format fout " - s\[0\]" )
                    (format fout " - s\[~a\]" op@3)))
                 ((null op@3)
                  (format fout   " - c\[~a\]*s\[0\]" const@3))
                 (t (format fout " - c\[~a\]*s\[~a\]" const@3 op@3))))
          ((< fac 0)
           (cond ((null const@3)
                  (if (null op@3)
                    (format fout " ~a*s\[0\]" fac)
                    (format fout " ~a*s\[~a\]" fac op@3)))
                 ((null op@3)
                  (format fout   " ~a*c\[~a\]*s\[0\]" fac const@3))
                 (t (format fout " ~a*c\[~a\]*s\[~a\]" fac const@3 op@3))))
          (t
            (cond ((null const@3)
                   (if (null op@3)
                     (format fout " + ~a*s\[0\]" fac)
                     (format fout " + ~a*s\[~a\]" fac op@3)))
                  ((null op@3)
                   (format fout   " + ~a*c\[~a\]*s\[0\]" fac const@3))
                  (t (format fout " + ~a*c\[~a\]*s\[~a\]" fac const@3 op@3)))))))

(defun cell-converter-avx (cell fout)
  (let ((fac (realpart (phase-of cell)))
        (const@3 (ternary-subscript (consts-of cell)))
        (op@3    (ternary-subscript (ops-of cell))))
    (cond ((equal fac 1)
           (cond ((null const@3)
                  (if (null op@3)
                    (format fout " + rs\[0\]" )
                    (format fout " + rs\[~a\]" op@3)))
                 ((null op@3)
                  (format fout   " + MM_SET1(c\[~a\])*rs\[0\]" const@3))
                 (t (format fout " + MM_SET1(c\[~a\])*rs\[~a\]" const@3 op@3))))
          ((equal fac -1)
           (cond ((null const@3)
                  (if (null op@3)
                    (format fout " - rs\[0\]" )
                    (format fout " - rs\[~a\]" op@3)))
                 ((null op@3)
                  (format fout   " - MM_SET1(c\[~a\])*rs\[0\]" const@3))
                 (t (format fout " - MM_SET1(c\[~a\])*rs\[~a\]" const@3 op@3))))
          (t (cond ((null const@3)
                    (if (null op@3)
                      (format fout " + MM_SET1(~a)*rs\[0\]" fac)
                      (format fout " + MM_SET1(~a)*rs\[~a\]" fac op@3)))
                   ((null op@3)
                    (format fout   " + MM_SET1(~a)*MM_SET1(c\[~a\])*rs\[0\]" fac const@3))
                   (t (format fout " + MM_SET1(~a)*MM_SET1(c\[~a\])*rs\[~a\]" fac const@3 op@3)))))))

(defun to-c-code-string (fout c-converter flat-script)
  (flet ((c-streamer (cs)
           (with-output-to-string (tmpout)
             (cond ((null cs) (format tmpout " 0"))
                   ((cell? cs) (funcall c-converter cs tmpout))
                   (t (mapcar (lambda (c) (funcall c-converter c tmpout)) cs))))))
    (mapcar #'c-streamer flat-script)))

(defun gen-c-block (fout flat-script)
  (let ((assemb (to-c-code-string fout #'cell-converter flat-script))
        (comp (length flat-script)))
    (loop for s in assemb
          for gid from 0 do
          (format fout "gout[n*~a+~a] =~a;~%" comp gid s))))

(defun gen-c-block-avx (fout flat-script)
  (let ((assemb (to-c-code-string fout #'cell-converter-avx flat-script))
        (comp (length flat-script)))
    (loop for s in assemb
          for gid from 0 do
          (if (equal s " 0")
            (format fout "GOUT_SCATTER(gout, n*~a+~a, MM_SET1(0.));~%" comp gid)
            (format fout "r1 =~a; GOUT_SCATTER(gout, n*~a+~a, r1);~%" s comp gid)))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; effective keys are p,r,ri,...
(defun effect-keys (ops)
  (remove-if-not (lambda (x)
                   (or (member x *act-left-right*)
                       (member x *intvar*)))
                 ops))
(defun g?e-of (key)
  (case key
    ((p ip nabla px py pz p* ip* nabla* px* py* pz*) "D_")
    ((r x y z) "R_") ; the vector origin is on the center of the basis it acts on
    ((ri rj rk rl) "RC") ; the vector origin is R[ijkl]
    ((xi xj xk xl) "RC")
    ((yi yj yk yl) "RC")
    ((zi zj zk zl) "RC")
    ((r0 x0 y0 z0 g) "R0") ; R0 ~ the vector origin is (0,0,0)
    ((rc xc yc zc) "RC") ; the vector origin is set in env[PTR_COMMON_ORIG]
    ((nabla-rinv nabla-r12 breit-r1 breit-r2) "D_")
    (otherwise (error "unknown key ~a~%" key))))

(defun dump-header (fout)
  (format fout "/*
 * Qcint is a general GTO integral library for computational chemistry
 * Copyright (C) 2014- Qiming Sun <osirpt.sun@gmail.com>
 *
 * This file is part of Qcint.
 *
 * Qcint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Description: code generated by  gen-code.cl
 */
#include <stdlib.h>
#include \"cint_bas.h\"
#include \"simd.h\"
#include \"cart2sph.h\"
#include \"g1e.h\"
#include \"g2e.h\"
#include \"optimizer.h\"
#include \"cint1e.h\"
#include \"cint2e.h\"
#include \"misc.h\"
#include \"c2f.h\"
"))

; FIXME: how to declare dri if rc and r[ijkl] coexist?
(defun dump-declare-dri-for-rc (fout i-ops symb)
  (when (intersection '(rc xc yc zc) i-ops)
    (format fout "double dr~a[3];~%" symb)
    (format fout "dr~a[0] = envs->r~a[0] - envs->env[PTR_COMMON_ORIG+0];~%" symb symb)
    (format fout "dr~a[1] = envs->r~a[1] - envs->env[PTR_COMMON_ORIG+1];~%" symb symb)
    (format fout "dr~a[2] = envs->r~a[2] - envs->env[PTR_COMMON_ORIG+2];~%" symb symb))
  (when (intersection '(ri xi yi zi) i-ops)
    (if (intersection '(rc xc yc zc) i-ops)
      (error "Cannot declare dri because rc and ri coexist"))
    (format fout "double dr~a[3];~%" symb)
    (format fout "dr~a[0] = envs->r~a[0] - envs->ri[0];~%" symb symb)
    (format fout "dr~a[1] = envs->r~a[1] - envs->ri[1];~%" symb symb)
    (format fout "dr~a[2] = envs->r~a[2] - envs->ri[2];~%" symb symb))
  (when (intersection '(rj xj yj zj) i-ops)
    (if (intersection '(rc xc yc zc) i-ops)
      (error "Cannot declare drj because rc and rj coexist"))
    (format fout "double dr~a[3];~%" symb)
    (format fout "dr~a[0] = envs->r~a[0] - envs->rj[0];~%" symb symb)
    (format fout "dr~a[1] = envs->r~a[1] - envs->rj[1];~%" symb symb)
    (format fout "dr~a[2] = envs->r~a[2] - envs->rj[2];~%" symb symb))
  (when (intersection '(rk xk yk zk) i-ops)
    (if (intersection '(rc xc yc zc) i-ops)
      (error "Cannot declare drk because rc and rk coexist"))
    (format fout "double dr~a[3];~%" symb)
    (format fout "dr~a[0] = envs->r~a[0] - envs->rk[0];~%" symb symb)
    (format fout "dr~a[1] = envs->r~a[1] - envs->rk[1];~%" symb symb)
    (format fout "dr~a[2] = envs->r~a[2] - envs->rk[2];~%" symb symb))
  (when (intersection '(rl xl yl zl) i-ops)
    (if (intersection '(rc xc yc zc) i-ops)
      (error "Cannot declare drl because rc and rl coexist"))
    (format fout "double dr~a[3];~%" symb)
    (format fout "dr~a[0] = envs->r~a[0] - envs->rl[0];~%" symb symb)
    (format fout "dr~a[1] = envs->r~a[1] - envs->rl[1];~%" symb symb)
    (format fout "dr~a[2] = envs->r~a[2] - envs->rl[2];~%" symb symb)))

(defun dump-declare-giao-ij (fout bra ket)
  (let ((n-giao (count 'g (append bra ket))))
    (when (> n-giao 0)
      (format fout "double rirj[3], c[~a];~%" (expt 3 n-giao))
      (format fout "rirj[0] = envs->ri[0] - envs->rj[0];~%" )
      (format fout "rirj[1] = envs->ri[1] - envs->rj[1];~%" )
      (format fout "rirj[2] = envs->ri[2] - envs->rj[2];~%" )
      (loop
        for i upto (1- (expt 3 n-giao)) do
        (format fout "c[~a] = 1" i)
        (loop
          for j from (1- n-giao) downto 0
          and res = i then (multiple-value-bind (int res) (floor res (expt 3 j))
                             (format fout " * rirj[~a]" int)
                             res))
        (format fout ";~%")))))

; l-combo searches op_bit from left to right
;  o100 o010 o001|g...>  =>  |g...> = o100 |g0..>
; a operator can only be applied to the left of the existed ones
;  |g100,l> = o100 |g000,l+1>
;  |g101,l> = o100 |g001,l+1>
;  |g110,l> = o100 |g010,l+1>
;  |g111,l> = o100 |g011,l+1>
; as a result, g* intermediates are generated from the previous one whose
; id (in binary) can be obtained by removing the first bit 1 from current
; id (in binary), see combo-bra function, eg
;  000  g0
;  001  g1 from g0 (000)
;  010  g2 from g0 (000)
;  011  g3 from g1 (001)
;  100  g4 from g0 (000)
;  101  g5 from g1 (001)
;  110  g6 from g2 (010)
;  111  g7 from g3 (011)
; r-combo searches op_bit from right to left
;  o100 o010 o001|g...>  =>  |g...> = o001 |g..0>
; a operator can only be applied to the right of the exsited ones
;  |g100,l+2> = o100 |g000,l+3>
;  |g101,l  > = o001 |g100,l+1>
;  |g110,l+1> = o010 |g100,l+2>
;  |g111,l  > = o001 |g110,l+1>
; as a result, g* intermediates are generated from the previous one whose
; id (in binary) can be obtained by removing the last bit 1 from current
; id (in binary), see combo-ket function, eg
;  000  g0
;  001  g1 from g0 (000)
;  010  g2 from g0 (000)
;  011  g3 from g2 (010)
;  100  g4 from g0 (000)
;  101  g5 from g4 (100)
;  110  g6 from g4 (100)
;  111  g7 from g6 (110)
; [lr]-combo have no connection with <bra| or |ket>
;    def l_combinator(self, ops, ig, mask, template):
(defun first-bit1 (n)
  (loop
    for i upto 31
    thereis (if (zerop (ash n (- i))) (1- i))))
(defun last-bit1 (n)
  ; how many 0s follow the last bit 1
  (loop
    for i upto 31
    thereis (if (oddp (ash n (- i))) i)))
(defun combo-bra (fout fmt ops-rev n-ops ig mask)
  (let* ((right (first-bit1 (ash ig (- mask))))
         (left (- n-ops right 1))
         (ig0 (- ig (ash 1 (+ mask right))))
         (op (nth right ops-rev)))
    (format fout fmt (g?e-of op) ig ig0 left)))
(defun combo-ket (fout fmt ops-rev i-len ig mask)
  (let* ((right (last-bit1 (ash ig (- mask))))
         (ig0 (- ig (ash 1 (+ mask right))))
         (op (nth right ops-rev)))
    (format fout fmt (g?e-of op) ig ig0 i-len right)))
(defun combo-opj (fout fmt-op fmt-j opj-rev i-len j-len ig mask)
  (let ((right (last-bit1 (ash ig (- mask)))))
    (if (< right j-len) ; does not reach op yet
      (combo-ket fout fmt-j opj-rev i-len ig mask)
      (let ((ig0 (- ig (ash 1 (+ mask right))))
            (op (nth right opj-rev)))
        (if (member op *act-left-right*)
          (format fout fmt-op
                  (g?e-of op) ig ig0 i-len right
                  (g?e-of op) (1+ ig) ig0 i-len right
                  ig (1+ ig))
          (if (and (intersection *act-left-right* opj-rev)
                   (< (1+ right) (length opj-rev))) ; ops have *act-left-right* but the rightmost op is not
              (format fout fmt-j (g?e-of op) ig ig0 (1+ i-len) right)
              (format fout fmt-j (g?e-of op) ig ig0 i-len right)))))))

(defun power2-range (n &optional (shift 0))
  (range (+ shift (ash 1 n)) (+ shift (ash 1 (1+ n)))))
(defun dump-combo-braket (fout fmt-i fmt-op fmt-j i-rev op-rev j-rev mask)
  (let* ((i-len (length i-rev))
         (j-len (length j-rev))
         (op-len (length op-rev))
         (opj-rev (append j-rev op-rev)))
    (loop
      for right from mask to (+ mask j-len op-len -1) do
      (loop
        for ig in (power2-range right) do
        (combo-opj fout fmt-op fmt-j opj-rev i-len j-len ig mask)))
    (let ((shft (+ op-len j-len mask)))
      (loop
        for right from shft to (+ shft i-len -1) do
        (loop
          for ig in (power2-range right) do
          (combo-bra fout fmt-i i-rev i-len ig shft))))))

(defun dec-to-ybin (n)
  (parse-integer (substitute #\0 #\2 (write-to-string n :base 3))
                 :radix 2))
(defun dec-to-zbin (n)
  (parse-integer (substitute #\1 #\2
                             (substitute #\0 #\1
                                         (write-to-string n :base 3)))
                 :radix 2))

(defun name-c2sor (sfx sp sf ts)
  (cond ((eql sp 'spinor)
         (if (eql sf 'sf)
           (if (eql ts 'ts)
             (format nil "&c2s_sf_~a" sfx)
             (format nil "&c2s_sf_~ai" sfx))
           (if (eql ts 'ts)
             (format nil "&c2s_si_~a" sfx)
             (format nil "&c2s_si_~ai" sfx))))
         ((eql sp 'spheric)
          (format nil "&c2s_sph_~a" sfx))
         (t (format nil "&c2s_cart_~a" sfx))))

;;; g-intermediates are g0, g1, g2, ...
(defun num-g-intermediates (tot-bits op i-len j-len)
  (if (and (intersection *act-left-right* op)
           ; nabla-rinv is the last one in the operation list
           (<= (+ i-len (length op) j-len) 1))
      (ash 1 tot-bits)
      (1- (ash 1 tot-bits))))

;!!! Be very cautious of the reverse on i-operators and k operators!
;!!! When multiple tensor components (>= rank 2) provided by the operators
;!!! on bra functions, the ordering of the multiple tensor components are
;!!! also reversed in the generated integral code
(defun gen-code-gout1e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i)) ;<i| already in reverse order
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout1e_~a" intname)
      (format fout "(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
CINTg1e_ovlp(g, envs, count);
int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int ix, iy, iz, n;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a  + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
;;; generate g_(bin)
;;; for the operators act on the |ket>, the reversed scan order and r_combinator
;;; is required; for the operators acto on the <bra|, the normal scan order and
      (let ((fmt-i "G1E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, 0);~%")
            (fmt-op (mkstr "G1E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0);
G1E_~aI(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G1E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev 0))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];~%")
      (dump-s-for-nroots-avx fout tot-bits 1)
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-gout1e-nuc (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i)) ;<i| already in reverse order
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout1e_~a" intname)
      (format fout "(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int nfc = nf * ~a;
double *gout1 = gout + nfc*SIMDD;~%" goutinc)
      (format fout "int ix, iy, iz, ia, n, i;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a  + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
      (format fout "for (n = 0; n < nfc*SIMDD; n++) { gout1[n] = 0; }
for (ia = 0; ia < envs->natm; ia++) {
CINTg1e_nuc(g, envs, count, ia);~%")
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, 0, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);
G2E_~aI(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev 0))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "}~%")
      (let ((assemb (to-c-code-string fout #'cell-converter-avx flat-script)))
        (loop for s in assemb
              for gid from 0 do
              (unless (equal s " 0")
                (format fout "r1 = MM_LOAD(gout1+(n*~a+~a)*SIMDD) +~a; " goutinc gid s)
                (format fout "MM_STORE(gout1+(n*~a+~a)*SIMDD, r1);~%" goutinc gid))))
      (format fout "}
}
CINTsort_gout(gout, gout1, nfc, SIMDD);
}~%" goutinc)
      goutinc)))

(defun gen-code-gout1e-rinv (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i)) ;<i| already in reverse order
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout1e_~a" intname)
      (format fout "(double *gout, double *g, int *idx, CINTEnvVars *envs, int count) {
CINTg1e_nuc(g, envs, count, -1);
int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int nrys_roots = envs->nrys_roots;
int ix, iy, iz, n, i;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, 0, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);
G2E_~aI(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, 0, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev 0))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "}~%")
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-int1e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i)) ;<i| already in reverse order
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (raw-script (eval-int raw-infix))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (ts (car raw-script))
           (sf (cadr raw-script))
           (goutinc (length flat-script))
           (e1comps (if (eql sf 'sf) 1 4))
           (tensors (if (eql sf 'sf) goutinc (/ goutinc 4)))
           (ngdef (with-output-to-string (tmpout)
                    (if (or (member 'nuc raw-infix)
                            (member 'rinv raw-infix)
                            (member 'nabla-rinv raw-infix))
                      (format tmpout "int ng[] = {~d, ~d, 0, 0, ~d, ~d, 0, ~d};~%"
                              i-len (+ op-len j-len) tot-bits e1comps tensors)
                      (format tmpout "int ng[] = {~d, ~d, 0, 0, ~d, ~d, 1, ~d};~%"
                              i-len (+ op-len j-len) tot-bits e1comps tensors))))
           (envs-common (with-output-to-string (tmpout)
                          (format tmpout ngdef)
                          (format tmpout "CINTEnvVars envs;~%")
                          (format tmpout "CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);~%")
                          (format tmpout "envs.f_gout = &CINTgout1e_~a;~%" intname)
                          (unless (eql (factor-of raw-infix) 1)
                            (format tmpout "envs.common_factor *= ~a;~%" (factor-of raw-infix))))))
      (format fout "/* <~{~a ~}i|~{~a ~}|~{~a ~}j> */~%" bra-i op ket-j)
      (cond ((member 'nuc raw-infix)
             (gen-code-gout1e-nuc fout intname raw-infix flat-script))
            ((or (member 'rinv raw-infix)
                 (member 'nabla-rinv raw-infix))
             (gen-code-gout1e-rinv fout intname raw-infix flat-script))
            (t (gen-code-gout1e fout intname raw-infix flat-script)))
      (format fout "void ~a_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {~%" intname)
      (format fout ngdef)
      (format fout "CINTall_1e_optimizer(opt, ng, atm, natm, bas, nbas, env);~%}~%")
;;; _cart
      (format fout "int ~a_cart(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = envs.nfi * envs.x_ctr[0];
counts[1] = envs.nfj * envs.x_ctr[1];
counts[2] = 1;
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
} // ~a_cart~%" intname)
;;; _sph
      (format fout "int ~a_sph(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = (envs.i_l*2+1) * envs.x_ctr[0];
counts[1] = (envs.j_l*2+1) * envs.x_ctr[1];
counts[2] = 1;
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT1e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
} // ~a_sph~%" intname)
;;; _spinor
      (format fout "int ~a_spinor(double complex *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
counts[2] = 1;
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1];
for (i = 0; i < envs.ncomp_tensor; i++) {
c2s_zset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT1e_spinor_drv(out, dims, &envs, opt, cache, ~a);
} // ~a_spinor~%" (name-c2sor "1e" 'spinor sf ts) intname)))
;;; int1e -> cint1e
  (format fout "ALL_CINT1E(~a)~%" intname)
  (format fout "//ALL_CINT1E_FORTRAN_(c~a)~%" intname))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dump-declare-giao-ijkl (fout opi opj opk opl)
  (let ((n-ij (count 'g (append opi opj)))
        (n-kl (count 'g (append opk opl))))
    (when (> n-ij 0)
      (format fout "double rirj[3];~%")
      (format fout "rirj[0] = envs->ri[0] - envs->rj[0];~%" )
      (format fout "rirj[1] = envs->ri[1] - envs->rj[1];~%" )
      (format fout "rirj[2] = envs->ri[2] - envs->rj[2];~%" ))
    (when (> n-kl 0)
      (format fout "double rkrl[3];~%")
      (format fout "rkrl[0] = envs->rk[0] - envs->rl[0];~%" )
      (format fout "rkrl[1] = envs->rk[1] - envs->rl[1];~%" )
      (format fout "rkrl[2] = envs->rk[2] - envs->rl[2];~%" ))
    (when (> (+ n-ij n-kl) 0)
      (format fout "double c[~a];~%" (expt 3 (+ n-ij n-kl)))
      (loop
        for i upto (1- (expt 3 (+ n-ij n-kl))) do
        (format fout "c[~a] = 1" i)
        (loop
          for j from (+ n-ij n-kl -1) downto n-kl
          and res = i then (multiple-value-bind (int res) (floor res (expt 3 j))
                             (format fout " * rirj[~a]" int)
                             res))
        (loop
          for j from (1- n-kl) downto 0
          and res = (nth-value 1 (floor i (expt 3 n-kl)))
                    then (multiple-value-bind (int res) (floor res (expt 3 j))
                           (format fout " * rkrl[~a]" int)
                           res))
        (format fout ";~%")))))

(defun dump-s-for-nroots-avx (fout tot-bits nroots)
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "rs[~a] = MM_LOAD(g~a+ix*SIMDD) * MM_LOAD(g~a+iy*SIMDD) * MM_LOAD(g~a+iz*SIMDD);~%"
              i xbin ybin zbin)
      (loop
        for k from 1 to (1- nroots) do
        (format fout "rs[~a]+= MM_LOAD(g~a+(ix+~a)*SIMDD) * MM_LOAD(g~a+(iy+~a)*SIMDD) * MM_LOAD(g~a+(iz+~a)*SIMDD);~%"
                i xbin k ybin k zbin k i)))))
(defun dump-s-loop-avx (fout tot-bits)
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "rs[~a] += MM_LOAD(g~a+(ix+i)*SIMDD) * MM_LOAD(g~a+(iy+i)*SIMDD) * MM_LOAD(g~a+(iz+i)*SIMDD);~%"
              i xbin ybin zbin i)))) ; end do i = 1, envs->nrys_roots
(defun dump-s-2e-avx (fout tot-bits &optional (deriv-max 3))
  (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
  (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];~%")
  (if (< tot-bits deriv-max)
    (progn
      (format fout "switch (nrys_roots) {~%")
      (loop
        for i from 1 to 4 do
        (format fout "case ~a:~%" i)
        (dump-s-for-nroots-avx fout tot-bits i)
        (format fout "break;~%" ))
      (format fout "default:
for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "} break;}~%")) ; else
    (progn
      (format fout "for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "}~%"))))

(defun dump-s-for-nroots (fout tot-bits nroots)
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "s[~a] = " i)
      (loop
        for k upto (1- nroots) do
        (format fout "+ g~a[ix+~a]*g~a[iy+~a]*g~a[iz+~a]"
                xbin k ybin k zbin k))
      (format fout ";~%"))))
(defun dump-s-loop (fout tot-bits )
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "s[~a] += g~a[ix+i] * g~a[iy+i] * g~a[iz+i];~%"
              i xbin ybin zbin)))) ; end do i = 1, envs->nrys_roots
(defun dump-s-2e (fout tot-bits &optional (deriv-max 2))
  (format fout "double s[~a];~%" (expt 3 tot-bits))
  (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];~%")
  (if (< tot-bits deriv-max)
    (progn
      (format fout "switch (nrys_roots) {~%")
      (loop
        for i from 1 to 4 do
        (format fout "case ~a:~%" i)
        (dump-s-for-nroots fout tot-bits i)
        (format fout "break;~%" ))
      (format fout "default:
for (i = 0; i < ~a; i++) { s[i] = 0; }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop fout tot-bits)
      (format fout "} break;}~%")) ; else
    (progn
      (format fout "for (i = 0; i < ~a; i++) { s[i] = 0; }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop fout tot-bits)
      (format fout "}~%"))))

;;; generate function gout2e
(defun gen-code-gout4c2e-avx (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (l-rev (reverse (effect-keys ket-l)))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (l-len (length l-rev))
           (tot-bits (+ i-len j-len op-len k-len l-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout2e_~a(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout ket-j "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-dri-for-rc fout ket-l "l")
      (dump-declare-giao-ijkl fout bra-i ket-j bra-k ket-l)
;;; generate g_(bin)
;;; TODO: 2 pass to optimize redundant operations, eg
;;; G2E_D_I(g1, g0, envs->i_l+1, envs->j_l, envs->k_l, envs->l_l);
;;; G2E_D_I(g2, g0, envs->i_l+0, envs->j_l, envs->k_l, envs->l_l);
;;; can be converted to g2 = g1;
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, envs->l_l);~%"))
              (fmt-op "")
              (fmt-l (mkstr "G2E_~aL(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, envs->l_l+~a);~%")))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() l-rev 0))
;!! FIXME, if there are operators of electron 2 in op
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, envs->l_l);~%"))
              (fmt-op "")
              (fmt-l (mkstr "G2E_~aL(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, envs->l_l+~a);~%")))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() l-rev 0)))
;;; The derivative part of Breit term (nabla-r12) acts on both bra and ket.
;;; In current treatments nabla-r12 are combined to j-operators.
;;; j-operators are handled before i-operators, the nabla (of nabla-r12) on
;;; i-operators is applied the earliest, than all other operators on <i|
;;; Because nabla (of nabla-r12) on <i| is always on the left-most of i-operators,
;;; applying it first on <i| is exactly what we want.
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, envs->l_l);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, envs->l_l);
G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, envs->l_l);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, envs->l_l);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev (+ k-len l-len)))
;;; generate gout
      (dump-s-2e-avx fout tot-bits)
;;; dump result of eval-int
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

;;; generate function gout2e_simd1
(defun gen-code-gout4c2e-simd1 (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (l-rev (reverse (effect-keys ket-l)))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (l-len (length l-rev))
           (tot-bits (+ i-len j-len op-len k-len l-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout2e_~a_simd1(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout ket-j "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-dri-for-rc fout ket-l "l")
      (dump-declare-giao-ijkl fout bra-i ket-j bra-k ket-l)
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, envs->l_l);~%"))
              (fmt-op "")
              (fmt-l (mkstr "G2E_~aL_SIMD1(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, envs->l_l+~a);~%")))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() l-rev 0))
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, envs->l_l);~%"))
              (fmt-op "")
              (fmt-l (mkstr "G2E_~aL_SIMD1(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, envs->l_l+~a);~%")))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() l-rev 0)))
      (let ((fmt-i "G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, envs->l_l);~%")
            (fmt-op (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, envs->l_l);
G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, envs->l_l);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, envs->l_l);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev (+ k-len l-len)))
      (dump-s-2e fout tot-bits)
      (gen-c-block fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-int4c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (l-rev (reverse (effect-keys ket-l)))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (l-len (length l-rev))
           (tot-bits (+ i-len j-len op-len k-len l-len))
           (raw-script (eval-int raw-infix))
           (ts1 (car raw-script))
           (sf1 (cadr raw-script))
           (ts2 (caddr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (with-output-to-string (tmpout)
                    (if (intersection *act-left-right* op)
                      (format tmpout "int ng[] = {~d, ~d, ~d, ~d, ~d, ~d, ~d, ~d};~%"
                              (1+ i-len) (1+ j-len) k-len l-len tot-bits e1comps e2comps tensors)
                      (format tmpout "int ng[] = {~d, ~d, ~d, ~d, ~d, ~d, ~d, ~d};~%"
                              i-len (+ op-len j-len) k-len l-len tot-bits e1comps e2comps tensors))))
           (envs-common (with-output-to-string (tmpout)
                          (format tmpout ngdef)
                          (format tmpout "CINTEnvVars envs;~%")
                          (format tmpout "CINTinit_int2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);~%")
                          (format tmpout "envs.f_gout = &CINTgout2e_~a;~%" intname)
                          (format tmpout "envs.f_gout_simd1 = &CINTgout2e_~a_simd1;~%" intname)
                          (unless (eql (factor-of raw-infix) 1)
                            (format tmpout "envs.common_factor *= ~a;~%" (factor-of raw-infix))))))
      (format fout "/* <~{~a ~}k ~{~a ~}i|~{~a ~}|~{~a ~}j ~{~a ~}l> : i,j \\in electron 1; k,l \\in electron 2~%"
              bra-k bra-i op ket-j ket-l)
      (format fout " * = (~{~a ~}i ~{~a ~}j|~{~a ~}|~{~a ~}k ~{~a ~}l) */~%"
              bra-i ket-j op bra-k ket-l)
      (gen-code-gout4c2e-avx fout intname raw-infix flat-script)
      (gen-code-gout4c2e-simd1 fout intname raw-infix flat-script)
      (format fout "void ~a_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {~%" intname)
      (format fout ngdef)
      (format fout "CINTall_2e_optimizer(opt, ng, atm, natm, bas, nbas, env);~%}~%")
;;; _cart
      (format fout "int ~a_cart(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (and (eql sf1 'si) (eql ts1 'tas)
                 (eql sf2 'si) (eql ts2 'tas))
        (format fout "envs.common_factor *= -1;~%"))
      (when (member 'g raw-infix)
        (format fout "int i, nout;
int counts[4];~%")
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
counts[0] = envs.nfi * envs.x_ctr[0];
counts[1] = envs.nfj * envs.x_ctr[1];
counts[2] = envs.nfk * envs.x_ctr[2];
counts[3] = envs.nfl * envs.x_ctr[3];
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%"))
        (when (or (member 'g bra-k) (member 'g ket-l))
          (format fout "if (out != NULL && envs.shls[2] == envs.shls[3]) {
counts[0] = envs.nfi * envs.x_ctr[0];
counts[1] = envs.nfj * envs.x_ctr[1];
counts[2] = envs.nfk * envs.x_ctr[2];
counts[3] = envs.nfl * envs.x_ctr[3];
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT2e_cart_drv(out, dims, &envs, opt, cache);~%} // ~a_cart~%" intname)
;;; _sph
      (format fout "int ~a_sph(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (and (eql sf1 'si) (eql ts1 'tas)
                 (eql sf2 'si) (eql ts2 'tas))
        (format fout "envs.common_factor *= -1;~%"))
      (when (member 'g raw-infix)
        (format fout "int i, nout;
int counts[4];~%")
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
counts[0] = (envs.i_l*2+1) * envs.x_ctr[0];
counts[1] = (envs.j_l*2+1) * envs.x_ctr[1];
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = (envs.l_l*2+1) * envs.x_ctr[3];
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%"))
        (when (or (member 'g bra-k) (member 'g ket-l))
          (format fout "if (out != NULL && envs.shls[2] == envs.shls[3]) {
counts[0] = (envs.i_l*2+1) * envs.x_ctr[0];
counts[1] = (envs.j_l*2+1) * envs.x_ctr[1];
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = (envs.l_l*2+1) * envs.x_ctr[3];
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT2e_spheric_drv(out, dims, &envs, opt, cache);~%} // ~a_sph~%" intname)
;;; _spinor
      (format fout "int ~a_spinor(double complex *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (format fout "int i, nout;
int counts[4];~%")
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
counts[2] = CINTcgto_spinor(envs.shls[2], envs.bas);
counts[3] = CINTcgto_spinor(envs.shls[3], envs.bas);
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_tensor; i++) {
c2s_zset0(out+nout*i, dims, counts); }
return 0; }~%"))
        (when (or (member 'g bra-k) (member 'g ket-l))
          (format fout "if (out != NULL && envs.shls[2] == envs.shls[3]) {
counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
counts[2] = CINTcgto_spinor(envs.shls[2], envs.bas);
counts[3] = CINTcgto_spinor(envs.shls[3], envs.bas);
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2] * dims[3];
for (i = 0; i < envs.ncomp_tensor; i++) {
c2s_zset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT2e_spinor_drv(out, dims, &envs, opt, cache, ~a, ~a);
} // ~a_spinor~%" (name-c2sor "2e1" 'spinor sf1 ts1) (name-c2sor "2e2" 'spinor sf2 ts2) intname)))
;;; int2e -> cint2e
  (format fout "ALL_CINT(~a)~%" intname)
  (format fout "//ALL_CINT_FORTRAN_(c~a)~%" intname))


(defun gen-code-gout3c2e-avx (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "static void CINTgout2e_~a(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout ket-j "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-giao-ij fout bra-i ket-j)
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0))
;!! FIXME, if there are operators of electron 2 in op
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0)))
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);
G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev k-len))
;;; generate gout
      (dump-s-2e-avx fout tot-bits)
;;; dump result of eval-int
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

;;; generate function gout2e_simd1
(defun gen-code-gout3c2e-simd1 (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout2e_~a_simd1(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout ket-j "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-giao-ij fout bra-i ket-j)
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" (1+ i-len) ", envs->j_l+" (1+ j-len)
                            ", envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0))
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                            ", envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0)))
      (let ((fmt-i "G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);
G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev k-len))
      (dump-s-2e fout tot-bits)
      (gen-c-block fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-int3c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (ts1 (car raw-script))
           (sf1 (cadr raw-script))
           (ts2 (caddr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (with-output-to-string (tmpout)
                    (if (intersection *act-left-right* op)
                      (format tmpout "int ng[] = {~d, ~d, ~d, 0, ~d, ~d, ~d, ~d};~%"
                              (1+ i-len) (1+ j-len) k-len tot-bits e1comps e2comps tensors)
                      (format tmpout "int ng[] = {~d, ~d, ~d, 0, ~d, ~d, ~d, ~d};~%"
                              i-len (+ op-len j-len) k-len tot-bits e1comps e2comps tensors))))
           (envs-common (with-output-to-string (tmpout)
                          (format tmpout ngdef)
                          (format tmpout "CINTEnvVars envs;~%")
                          (format tmpout "CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);~%")
                          (format tmpout "envs.f_gout = &CINTgout2e_~a;~%" intname)
                          (format tmpout "envs.f_gout_simd1 = &CINTgout2e_~a_simd1;~%" intname)
                          (unless (eql (factor-of raw-infix) 1)
                            (format tmpout "envs.common_factor *= ~a;~%" (factor-of raw-infix))))))
      (format fout "/* (~{~a ~}i ~{~a ~}j|~{~a ~}|~{~a ~}k) */~%"
              bra-i ket-j op bra-k)
      (gen-code-gout3c2e-avx fout intname raw-infix flat-script)
      (gen-code-gout3c2e-simd1 fout intname raw-infix flat-script)
      (format fout "void ~a_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {~%" intname)
      (format fout ngdef)
      (format fout "CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);~%}~%")
;;; _cart
      (format fout "int ~a_cart(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = envs.nfi * envs.x_ctr[0];
counts[1] = envs.nfj * envs.x_ctr[1];
counts[2] = envs.nfk * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c2e_cart_drv(out, dims, &envs, opt, cache);~%} // ~a_cart~%" intname)
;;; _sph
      (format fout "int ~a_sph(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = (envs.i_l*2+1) * envs.x_ctr[0];
counts[1] = (envs.j_l*2+1) * envs.x_ctr[1];
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c2e_spheric_drv(out, dims, &envs, opt, cache);~%} // ~a_sph~%" intname)
;;; _spinor
      (format fout "int ~a_spinor(double complex *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_tensor; i++) {
c2s_zset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c2e_spinor_drv(out, dims, &envs, opt, cache, ~a);
} // ~a_spinor~%" (name-c2sor "3c2e1" 'spinor sf1 ts1) intname)))
;;; int2e -> cint2e
  (format fout "ALL_CINT(~a)~%" intname)
  (format fout "//ALL_CINT_FORTRAN_(c~a)~%" intname))


(defun gen-code-gout2c2e-avx (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "void CINTgout2e_~a(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len 0)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout bra-k "k")
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" (1+ i-len) ", 0, envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0))
;!! FIXME, if there are operators of electron 2 in op
        (let ((fmt-k (mkstr "G2E_~aK(g~a, g~a, envs->i_l+" i-len ", 0, envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0)))
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, 0, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, 0+~a, envs->k_l, 0);
G2E_~aI(g~a, g~a, envs->i_l+~a, 0+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, 0, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev '() k-len))
;;; generate gout
      (dump-s-2e-avx fout tot-bits)
;;; dump result of eval-int
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

;;; generate function gout2e_simd1
(defun gen-code-gout2c2e-simd1 (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "static void CINTgout2e_~a_simd1(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs) {~%" intname)
      (format fout "int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int ix, iy, iz, i, n;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len 0)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout bra-k "k")
      (if (intersection *act-left-right* op)
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" (1+ i-len) ", 0, envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0))
        (let ((fmt-k (mkstr "G2E_~aK_SIMD1(g~a, g~a, envs->i_l+" i-len ", 0, envs->k_l+~a, 0);~%"))
              (fmt-op "")
              (fmt-l ""))
          (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0)))
      (let ((fmt-i "G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, 0, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, 0+~a, envs->k_l, 0);
G2E_~aI_SIMD1(g~a, g~a, envs->i_l+~a, 0+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ_SIMD1(g~a, g~a, envs->i_l+~d, 0, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev '() k-len))
      (dump-s-2e fout tot-bits)
      (gen-c-block fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-int2c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (k-len (length k-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (ts1 (car raw-script))
           (sf1 (cadr raw-script))
           (ts2 (caddr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (with-output-to-string (tmpout)
                    (if (intersection *act-left-right* op)
                      (format tmpout "int ng[] = {~d, 1, ~d, 0, ~d, ~d, ~d, ~d};~%"
                              (1+ i-len) k-len tot-bits e1comps e2comps tensors)
                      (format tmpout "int ng[] = {~d, ~d, ~d, 0, ~d, ~d, ~d, ~d};~%"
                              i-len op-len k-len tot-bits e1comps e2comps tensors))))
           (envs-common (with-output-to-string (tmpout)
                          (format tmpout ngdef)
                          (format tmpout "CINTEnvVars envs;~%")
                          (format tmpout "CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);~%")
                          (format tmpout "envs.f_gout = &CINTgout2e_~a;~%" intname)
                          (format tmpout "envs.f_gout_simd1 = &CINTgout2e_~a_simd1;~%" intname)
                          (unless (eql (factor-of raw-infix) 1)
                            (format tmpout "envs.common_factor *= ~a;~%" (factor-of raw-infix))))))
      (format fout "/* (~{~a ~}i |~{~a ~}|~{~a ~}j) */~%"
              bra-i op bra-k)
      (gen-code-gout2c2e-avx fout intname raw-infix flat-script)
      (gen-code-gout2c2e-simd1 fout intname raw-infix flat-script)
      (format fout "void ~a_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {~%" intname)
      (format fout ngdef)
      (format fout "CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);~%}~%")
;;; _cart
      (format fout "int ~a_cart(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (format fout "return CINT2c2e_cart_drv(out, dims, &envs, opt, cache);~%} // ~a_cart~%" intname)
;;; _sph
      (format fout "int ~a_sph(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (format fout "return CINT2c2e_spheric_drv(out, dims, &envs, opt, cache);~%} // ~a_sph~%" intname)
;;; _spinor
      (format fout "int ~a_spinor(double complex *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      ; c2s_ function is incorrect if e1/e2 are spin-included operators
      (format fout "return CINT2c2e_spinor_drv(out, dims, &envs, opt, cache, ~a);
} // ~a_spinor~%" (name-c2sor "1e" 'spinor sf1 ts1) intname)))
;;; int2e -> cint2e
  (format fout "ALL_CINT(~a)~%" intname)
  (format fout "//ALL_CINT_FORTRAN_(c~a)~%" intname))


(defun gen-code-gout3c1e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "static void CINTgout1e_~a(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {~%" intname)
      (format fout "CINTg3c1e_ovlp(g, envs, count);
int nf = envs->nf;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int ix, iy, iz, n;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
      (let ((fmt-k (mkstr "G1E_~aK(g~a, g~a, envs->i_l+" i-len ", envs->j_l+" (+ op-len j-len)
                          ", envs->k_l+~a);~%"))
            (fmt-op "")
            (fmt-l ""))
        (dump-combo-braket fout fmt-k fmt-op fmt-l k-rev '() '() 0))
      (let ((fmt-i "G1E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l);~%")
            (fmt-op (mkstr "G1E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l);
G1E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G1E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev k-len))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];~%")
      (dump-s-for-nroots-avx fout tot-bits 1)
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-gout3c1e-nuc (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "static void CINTgout1e_~a" intname)
      (format fout "(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {
int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int nfc = nf * ~a;
double *gout1 = gout + nfc*SIMDD;~%" goutinc)
      (format fout "int ix, iy, iz, ia, n, i;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
      (format fout "for (n = 0; n < nfc*SIMDD; n++) { gout1[n] = 0; }
for (ia = 0; ia < envs->natm; ia++) {
CINTg3c1e_nuc(g, envs, count, ia);~%")
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);
G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev k-len))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "}~%")
      (let ((assemb (to-c-code-string fout #'cell-converter-avx flat-script)))
        (loop for s in assemb
              for gid from 0 do
              (unless (equal s " 0")
                (format fout "r1 = MM_LOAD(gout1+(n*~a+~a)*SIMDD) +~a; " goutinc gid s)
                (format fout "MM_STORE(gout1+(n*~a+~a)*SIMDD, r1);~%" goutinc gid))))
      (format fout "}
}
CINTsort_gout(gout, gout1, nfc, SIMDD);
}~%" goutinc)
      goutinc)))

(defun gen-code-gout3c1e-rinv (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (goutinc (length flat-script)))
      (format fout "static void CINTgout1e_~a" intname)
      (format fout "(double *RESTRICT gout,
double *RESTRICT g, int *RESTRICT idx, CINTEnvVars *envs, int count) {
CINTg3c1e_nuc(g, envs, count, -1);
int nf = envs->nf;
int nrys_roots = envs->nrys_roots;
int nfc = nf * ~a;~%" goutinc)
      (format fout "int ix, iy, iz, n, i;
DECLARE_GOUT;
double *RESTRICT g0 = g;~%")
      (loop
        for i in (range (num-g-intermediates tot-bits op i-len j-len)) do
        (format fout "double *RESTRICT g~a = g~a + envs->g_size * 3 * SIMDD;~%" (1+ i) i))
      (dump-declare-dri-for-rc fout bra-i "i")
      (dump-declare-dri-for-rc fout (append op ket-j) "j")
      (dump-declare-dri-for-rc fout bra-k "k")
      (dump-declare-giao-ij fout bra-i (append op ket-j))
      (format fout "__MD r1;
__MD rs[~a];~%" (expt 3 tot-bits))
      (let ((fmt-i "G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l, envs->k_l, 0);~%")
            (fmt-op (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);
G2E_~aI(g~a, g~a, envs->i_l+~a, envs->j_l+~a, envs->k_l, 0);
for (ix = 0; ix < envs->g_size * 3 * SIMDD; ix++) {g~a[ix] += g~a[ix];}~%"))
            (fmt-j (mkstr "G2E_~aJ(g~a, g~a, envs->i_l+~d, envs->j_l+~a, envs->k_l, 0);~%")))
        (dump-combo-braket fout fmt-i fmt-op fmt-j i-rev op-rev j-rev k-len))
      (format fout "for (n = 0; n < nf; n++) {
ix = idx[0+n*3];
iy = idx[1+n*3];
iz = idx[2+n*3];
for (i = 0; i < ~a; i++) { rs[i] = MM_SET1(0.); }
for (i = 0; i < nrys_roots; i++) {~%" (expt 3 tot-bits))
      (dump-s-loop-avx fout tot-bits)
      (format fout "}~%")
      (gen-c-block-avx fout flat-script)
      (format fout "}}~%")
      goutinc)))

(defun gen-code-int3c1e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (ts (car raw-script))
           (sf (cadr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf 'sf) 1 4))
           (tensors (if (eql sf 'sf) goutinc (/ goutinc 4)))
           (ngdef (with-output-to-string (tmpout)
                    (if (or (member 'nuc raw-infix)
                            (member 'rinv raw-infix)
                            (member 'nabla-rinv raw-infix))
                      (format tmpout "int ng[] = {~d, ~d, ~d, 0, ~d, ~d, 0, ~d};~%"
                              i-len (+ op-len j-len) k-len tot-bits e1comps tensors)
                      (format tmpout "int ng[] = {~d, ~d, ~d, 0, ~d, ~d, 1, ~d};~%"
                              i-len (+ op-len j-len) k-len tot-bits e1comps tensors))))
           (envs-common (with-output-to-string (tmpout)
                          (format tmpout ngdef)
                          (format tmpout "CINTEnvVars envs;~%")
                          (format tmpout "CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);~%")
                          (format tmpout "envs.f_gout = &CINTgout1e_~a;~%" intname)
                          (unless (eql (factor-of raw-infix) 1)
                            (format tmpout "envs.common_factor *= ~a;~%" (factor-of raw-infix))))))
      (cond ((member 'nuc raw-infix)
             (gen-code-gout3c1e-nuc fout intname raw-infix flat-script))
            ((or (member 'rinv raw-infix)
                 (member 'nabla-rinv raw-infix))
             (gen-code-gout3c1e-rinv fout intname raw-infix flat-script))
            (t (gen-code-gout3c1e fout intname raw-infix flat-script)))
      (format fout "void ~a_optimizer(CINTOpt **opt, int *atm, int natm, int *bas, int nbas, double *env) {~%" intname)
      (format fout ngdef)
      (format fout "CINTall_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env);~%}~%")
;;; _cart
      (format fout "int ~a_cart(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = envs.nfi * envs.x_ctr[0];
counts[1] = envs.nfj * envs.x_ctr[1];
counts[2] = envs.nfk * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c1e_cart_drv(out, dims, &envs, opt, cache);~%} // ~a_cart~%" intname)
;;; _sph
      (format fout "int ~a_sph(double *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = (envs.i_l*2+1) * envs.x_ctr[0];
counts[1] = (envs.j_l*2+1) * envs.x_ctr[1];
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor; i++) {
c2s_dset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c1e_spheric_drv(out, dims, &envs, opt, cache);~%} // ~a_sph~%" intname)
;;; _spinor
      (format fout "int ~a_spinor(double complex *out, int *dims, int *shls,
int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache) {~%" intname)
      (format fout envs-common)
      (when (member 'g raw-infix)
        (when (or (member 'g bra-i) (member 'g ket-j))
          (format fout "if (out != NULL && envs.shls[0] == envs.shls[1]) {
int i, nout;
int counts[4];
counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
counts[2] = (envs.k_l*2+1) * envs.x_ctr[2];
counts[3] = 1;
if (dims == NULL) { dims = counts; }
nout = dims[0] * dims[1] * dims[2];
for (i = 0; i < envs.ncomp_tensor; i++) {
c2s_zset0(out+nout*i, dims, counts); }
return 0; }~%")))
      (format fout "return CINT3c1e_spinor_drv(out, dims, &envs, opt, cache, ~a);
} // ~a_spinor~%" (name-c2sor "3c2e1" 'spinor sf ts) intname)))
  (format fout "ALL_CINT(~a)~%" intname)
  (format fout "//ALL_CINT_FORTRAN_(c~a)~%" intname))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun gen-cint (filename &rest items)
  "sp can be one of 'spinor 'spheric 'cart"
  (with-open-file (fout (mkstr filename)
                        :direction :output :if-exists :supersede)
    (dump-header fout)
    (flet ((gen-code (item)
             (let ((intname (mkstr (car item)))
                   (raw-infix (cadr item)))
               (cond ((int3c1e? raw-infix)
                      (gen-code-int3c1e fout intname raw-infix))
                     ((one-electron-int? raw-infix)
                      (gen-code-int1e fout intname raw-infix))
                     ((int4c2e? raw-infix)
                      (gen-code-int4c2e fout intname raw-infix))
                     ((int3c2e? raw-infix)
                      (gen-code-int3c2e fout intname raw-infix))
                     ((int2c2e? raw-infix)
                      (gen-code-int2c2e fout intname raw-infix))))))
      (mapcar #'gen-code items))))

; gcl -load sigma.o -batch -eval "( .. )"

;; vim: ft=lisp
