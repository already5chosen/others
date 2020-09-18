  .text
  .def   @feat.00;
  .scl  3;
  .type 0;
  .endef
  .globl  @feat.00
.set @feat.00, 0
  .file "chol_Factorize_a.s"
  .def   chol_Factorize_a;
  .scl  2;
  .type 32;
  .endef

  .section  .rdata,"dr"
  .p2align  3
.LDIAG_MIN:
  .quad 2.9387358770557187699218413430556e-39 # double 2**-128
.LDIAG_SUBST:
  .quad 1.0889035741470030830827987437817e+40 # double 2**133
.LONE:
  .quad 1.0

  .text
  .globl  chol_Factorize_a
  .p2align  4, 0x90
# extern "C" bool chol_Factorize_a(std::complex<double>* triang, unsigned N)
chol_Factorize_a:
# RCX - triang
# RDX - N
.seh_proc chol_Factorize_a
  pushq %rbx
  .seh_pushreg %rbx
  pushq %rsi
  .seh_pushreg %rsi
  pushq %rdi
  .seh_pushreg %rdi
  subq  $112, %rsp
  .seh_stackalloc 112
  vmovapd %xmm12,96(%rsp)
  .seh_savexmm %xmm12, 96
  vmovapd %xmm11,80(%rsp)
  .seh_savexmm %xmm11, 80
  vmovapd %xmm10,64(%rsp)
  .seh_savexmm %xmm10, 64
  vmovapd %xmm9, 48(%rsp)
  .seh_savexmm %xmm9, 48
  vmovapd %xmm8, 32(%rsp)
  .seh_savexmm %xmm8, 32
  vmovapd %xmm7, 16(%rsp)
  .seh_savexmm %xmm7, 16
  vmovapd %xmm6, (%rsp)
  .seh_savexmm %xmm6, 0
  .seh_endprologue

  xor %r11, %r11
  test $1, %edx
  jz .start_even
   # special handling for the first row of matrix with odd number of elements
   # process top row
   add        $16, %rcx                 # ++triang - skip over padding
   vmovsd     (%rcx), %xmm0             # ymm0 = aa = triang[0].re, 0,0,0 // diagonal element
   # check that we are positive definite
   vcmpltsd .LDIAG_MIN(%rip), %xmm0, %xmm1             # xmm1 = ltmsk =[aa<DIAG_MIN ...]
   vblendvpd  %xmm1, .LDIAG_SUBST(%rip), %xmm0,%xmm0   # xmm0 = aa < DIAG_MIN ? DIAG_SUBST : aa
   vmovmskpd  %xmm1, %r11d
   and        $1,    %r11d
   vmovsd     .LONE(%rip),  %xmm1       # xmm1 = 1.0,zero
   vsqrtsd    %xmm0, %xmm0, %xmm0       # xmm0 = aaSqrt=sqrt(aa) ...
   vdivsd     %xmm0, %xmm1, %xmm1       # xmm1 = aaInvSqrt=1/sqrt(aa) ...
   vunpcklpd  %xmm1, %xmm0, %xmm0       # ymm0 = aaSqrt aaInvSqrt ...
   vmovupd    %xmm0, (%rcx)             # triang[0].re = aaSqrt triang[0].im = aaInvSqrt
   vunpcklpd  %xmm1, %xmm1, %xmm0       # ymm0 = aaInvSqrt aaInvSqrt 0 0

   shr $1,    %edx                      # edx = hlen
   jz .done                             # x0 was last (and only) row

   # Multiply top row by invSqrt
   add        $16,  %rcx                # ++triang
   mov        %rcx, %rdi                # rdi  = x0 = triang

   vperm2f128 $0,%ymm0,%ymm0,%ymm0      # ymm0 = aa0InvSqrt x4

   mov        %edx, %ebx                # ebx  = cnt = hlen
   .odd_top_rows_loop:
    vmulpd  (%rdi),      %ymm0, %ymm3   # ymm3 = vx = x0[c]*aa0InvSqrt (re im re im)
    vmovupd   %ymm3,(%rdi)              # x0[c]= vx
    add     $32,  %rdi                  # c += 2
   dec      %rbx                        # --cnt
   jnz      .odd_top_rows_loop

   # Subtract outer product of top row from lower part of the matrix
   # Process two output rows together
   mov      %edx, %ebx                  # ebx = chlen = hlen
   mov      %edx, %eax
   shl      $5,   %eax                  # eax = (y0-x0)*sizeof(DComplex)
  .caxpy1x2_outer_loop:
    mov     %ebx, %esi
    shl     $5,   %esi                  # esi = (y1-y0)*sizeof(DComplex) = chlen*2*sizeof(DComplex)
    add     %eax, %esi                  # esi = (y1-x0)*sizeof(DComplex)
    mov     %rcx, %rdi                  # rdi = x0 = triang
    vbroadcastsd   (%rcx),      %ymm0   # ymm0 = f0.re = x0[0].re  x4
    vbroadcastsd  8(%rcx),      %ymm1   # ymm1 = f0.im = x0[0].im  x4
    vbroadcastsd 16(%rcx),      %ymm2   # ymm2 = f1.re = x0[1].re  x4
    vbroadcastsd 24(%rcx),      %ymm3   # ymm3 = f1.im = x0[1].im  x4
    mov     %ebx, %r9d                  # r9d = cnt = chlen
    .caxpy1x2_inner_loop:
      vmovupd (%rdi),       %ymm4       # ymm4= vx  = x0[c] (re im re im)
      vmovupd (%rdi,%rax),  %ymm5       # ymm5= vy0 = y0[c] (re im re im)
      vmovupd (%rdi,%rsi),  %ymm6       # ymm6= vy1 = y1[c] (re im re im)

      vmulpd    %ymm4, %ymm0, %ymm7
      vsubpd    %ymm7, %ymm5, %ymm5     # ymm5= vy0.re-vx0.re*f00.re vy0.im-vx0.im*f00.re
      vmulpd    %ymm4, %ymm2, %ymm7
      vsubpd    %ymm7, %ymm6, %ymm6     # ymm6= vy1.re-vx0.re*f01.re vy1.im-vx0.im*f01.re
      vpermilpd $5,  %ymm4, %ymm4       # ymm4= vx0 { im re im re }

      vmulpd    %ymm4, %ymm1, %ymm7
      vaddsubpd %ymm7, %ymm5, %ymm5     # ymm5= vy0.re-vx0.im*f00.im vy0.im+vx0.re*f00.im
      vmulpd    %ymm4, %ymm3, %ymm7
      vaddsubpd %ymm7, %ymm6, %ymm6     # ymm6= vy1.re-vx0.im*f01.im vy1.im+vx0.re*f01.im

      vmovupd %ymm5, (%rdi,%rax)        # ymm5= vy0 = y0[c] (re im re im)
      vmovupd %ymm6, (%rdi,%rsi)        # ymm6= vy1 = y1[c] (re im re im)

      add    $32, %rdi                  # x0 += 2
    dec      %r9d                       # --cnt
    jnz .caxpy1x2_inner_loop

    sub %eax, %esi                      # esi = chlen*2*sizeof(DComplex)
    lea -32(%eax,%esi,2), %rax          # eax = y0-triang += chlen*4-2
    add    $32, %rcx                    # triang += 2
  dec      %rbx                         # --chlen
  jnz .caxpy1x2_outer_loop
  jmp .even_loop_entry

.start_even:
  shr $1, %edx                          # edx = rhlen = N/2
  jnz .even_loop_entry
  jmp .done

.even_loop:
  # Handle pair of top rows
  vperm2f128 $0,%ymm0,%ymm0,%ymm0       # ymm0 = aa0InvSqrt x4
  vperm2f128 $0,%ymm1,%ymm1,%ymm1       # ymm1 = aa1InvSqrt x4
  vperm2f128 $0,%ymm2,%ymm2,%ymm2       # ymm2 = f.re       x4
  vperm2f128 $0,%ymm3,%ymm3,%ymm3       # ymm3 = f.im       x4

  lea 32(%rcx), %rdi                    # rdi  = x0 = &triang[2]
  lea -1(%edx), %r10d                   # r10d = cnt = rhlen-1
  .top_rows_loop:
    vmulpd  (%rdi),      %ymm0, %ymm8   # ymm8 = vx0 = x0[c]*aa0InvSqrt (re im re im)
    vmovupd (%rdi,%rbx), %ymm9          # ymm9 = vx1 = x1[c] (re im re im)
    vmovupd %ymm8, (%rdi)               # x0[c] = vx0
    vmulpd    %ymm8, %ymm2,%ymm11
    vsubpd    %ymm11,%ymm9,%ymm9        # ymm9 = vx1.re-vx0.re*f.re vx1.im-vx0.im*f.re
    vpermilpd $5,  %ymm8,  %ymm8        # ymm8 = vx0 { im re im re }
    vmulpd    %ymm8, %ymm3,%ymm11       # ymm11=        vx0.im*f.im        vx0.re*f.im
    vaddsubpd %ymm11,%ymm9,%ymm9        # ymm9 = vx1.re-vx0.im*f.im vx1.im+vx0.re*f.im
    vmulpd    %ymm9, %ymm1,%ymm9        # ymm9 = vx1 = (vx1 - vx0*conj(f))*aa1InvSqrt (re im re im)
    vmovupd   %ymm9, (%rdi,%rbx)        # x1[c] = vx1
    add $32, %rdi                       # c += 2
  dec %r10d
  jnz .top_rows_loop

  add $32,  %rcx                        # triang += 2

  # subtract outer product of two top rows from lower part of the matrix
  # process two output rows together
  dec %edx                              # --rhlen
  lea (,%ebx,2), %eax                   # EAX  = rhlen*4*sizeof(DComplex) = y0-triang
  mov %edx, %r10d                       # r10d = chlen = rhlen
  .caxpy2x2_outer_loop:
    sub $32,  %eax                      # RAX = (y0-triang) -= 32 -- account for increment of triang
    mov %r10d,%esi
    shl $5,   %esi                      # esi  = chlen*2*sizeof(DComplex) = (y1-y0)*sizeof(DComplex)
    add %eax, %esi                      # esi  = (y1-x0)*sizeof(DComplex)

    vbroadcastsd   (%rcx),      %ymm0   # ymm0 = f00.re = x0[0].re  x4
    vbroadcastsd 8 (%rcx),      %ymm1   # ymm1 = f00.im = x0[0].im  x4
    vbroadcastsd 16(%rcx),      %ymm2   # ymm2 = f01.re = x0[1].re  x4
    vbroadcastsd 24(%rcx),      %ymm3   # ymm3 = f01.im = x0[1].im  x4

    vbroadcastsd   (%rcx,%rbx), %ymm4   # ymm4 = f10.re = x1[0].re  x4
    vbroadcastsd 8 (%rcx,%rbx), %ymm5   # ymm5 = f10.im = x1[0].im  x4
    vbroadcastsd 16(%rcx,%rbx), %ymm6   # ymm6 = f11.re = x1[1].re  x4
    vbroadcastsd 24(%rcx,%rbx), %ymm7   # ymm7 = f11.im = x1[1].im  x4

    mov %rcx, %rdi                      # rdi = &x0[c] = x0
    mov %r10d,%r9d                      # r9d = cnt = chlen
#  .p2align  3, 0x90
    .caxpy2x2_inner_loop:
      vmovupd (%rdi),      %ymm8        # ymm8 = vx0 = x0[c] (re im re im)
      vmovupd (%rdi,%rbx), %ymm9        # ymm9 = vx1 = x1[c] (re im re im)
      vmovupd (%rdi,%rax), %ymm10       # ymm10= vy0 = y0[c] (re im re im)
      vmovupd (%rdi,%rsi), %ymm11       # ymm11= vy1 = y1[c] (re im re im)

      vmulpd  %ymm8, %ymm0, %ymm12
      vsubpd  %ymm12,%ymm10,%ymm10      # ymm10= vy0.re-vx0.re*f00.re vy0.im-vx0.im*f00.re
      vmulpd  %ymm8, %ymm2, %ymm12
      vsubpd  %ymm12,%ymm11,%ymm11      # ymm11= vy1.re-vx0.re*f01.re vy1.im-vx0.im*f01.re
      vpermilpd $5,  %ymm8,%ymm8        # ymm8 = vx0 { im re im re }
      vmulpd  %ymm9, %ymm4, %ymm12
      vsubpd  %ymm12,%ymm10,%ymm10      # ymm10= vy0.re-vx1.re*f10.re vy0.im-vx1.im*f10.re
      vmulpd  %ymm9, %ymm6, %ymm12
      vsubpd  %ymm12,%ymm11,%ymm11      # ymm11= vy1.re-vx1.re*f11.re vy1.im-vx1.im*f11.re
      vpermilpd $5,  %ymm9,%ymm9        # ymm9 = vx1 { im re im re }

      vmulpd    %ymm8, %ymm1, %ymm12
      vaddsubpd %ymm12,%ymm10,%ymm10    # ymm10= vy0.re-vx0.im*f00.im vy0.im+vx0.re*f00.im
      vmulpd    %ymm8, %ymm3, %ymm12
      vaddsubpd %ymm12,%ymm11,%ymm11    # ymm11= vy1.re-vx0.im*f01.im vy1.im+vx0.re*f01.im
      vmulpd    %ymm9, %ymm5, %ymm12
      vaddsubpd %ymm12,%ymm10,%ymm10    # ymm10= vy0.re-vx1.im*f10.im vy0.im+vx1.re*f10.im
      vmulpd    %ymm9, %ymm7, %ymm12
      vaddsubpd %ymm12,%ymm11,%ymm11    # ymm11= vy1.re-vx1.im*f11.im vy1.im+vx1.re*f11.im

      vmovupd %ymm10, (%rdi,%rax)       # ymm10= vy0 = y0[c] (re im re im)
      vmovupd %ymm11, (%rdi,%rsi)       # ymm11= vy1 = y1[c] (re im re im)
      add $32, %rdi                     # c += 2
    dec %r9d                            # --cnt
    jnz .caxpy2x2_inner_loop

    sub %eax, %esi                      # esi = chlen*2*sizeof(DComplex)
    lea (%eax,%esi,2), %rax             # RAX = y0-triang += chlen*4
    add $32,  %rcx                      # x0 += 2
  dec %r10d                             # --chlen
  jnz .caxpy2x2_outer_loop

  add %rbx, %rcx                        # triang = x0 + (x1-x0)
.even_loop_entry:
  # RDX = rhlen
  mov %edx, %ebx
  shl $5, %ebx                    # rbx  = rhlen*2*sizeof(DComplex) = (x1-x0)*sizeof(DComplex)
  vmovsd   (%rcx), %xmm0          # ymm0 = aa0 = triang[0].re, 0,0,0
  vmovsd 16(%rcx,%rbx), %xmm1     # ymm1 = aa1 = x1[1].re, 0,0,0
  vmulsd    %xmm0, %xmm1, %xmm1
  vmovupd 16(%rcx),%xmm2          # ymm2 = {f.re f.im 0 0} = triang[1],0,0
  # aa1 = aa1*aa0 - norm(f)
  vmulpd    %xmm2, %xmm2, %xmm3   # xmm3 = {f.re^2 f.im^2}
  vhaddpd   %xmm3, %xmm3, %xmm3   # xmm3 = {f.re^2+f.im^2 ...}
  vsubsd    %xmm3, %xmm1, %xmm1   # xmm1 = {aa1*aa0 - norm(f) ...}
  vunpcklpd %xmm1, %xmm0, %xmm0   # ymm0 = aa0, aa1, 0, 0

  vmovddup .LDIAG_MIN(%rip), %xmm1
  vcmpltpd  %xmm1, %xmm0, %xmm1         # xmm1 = ltmsk =[aa0<DIAG_MIN aa1<DIAG_MIN]
  vmovddup  .LDIAG_SUBST(%rip), %xmm3
  vblendvpd %xmm1, %xmm3, %xmm0,%xmm0   # xmm0 = aa[0:1] < DIAG_MIN ? DIAG_SUBST : aa[0:1]
  vmovmskpd %xmm1, %eax
  or        %eax,  %r11d
  vsqrtpd   %xmm0, %xmm0                # xmm0 = [sqrt(aa0), sqrt(aa1)]
  vmovddup  .LONE(%rip), %xmm1
  vdivpd    %xmm0, %xmm1,%xmm1          # xmm1 = [1/sqrt(aa0), 1/sqrt(aa1)]

  vunpcklpd %xmm1, %xmm0,%xmm3          # ymm3 = aa0Sqrt=sqrt(aa0), aa0InvSqrt=1/sqrt(aa0), 0, 0
  vunpckhpd %xmm0, %xmm1,%xmm4          # ymm4 = 1/sqrt(aa1), sqrt(aa1), 0, 0
  vmovupd   %xmm3, (%rcx)               # triang[0].re = sqrt(aa0) triang[0].im = 1/sqrt(aa0)
  vmovddup  %xmm1, %xmm0                # xmm0 = {aa0InvSqrt aa0InvSqrt 0 0}

  vmulpd    %xmm4, %xmm3,%xmm4          # ymm4  = aa1InvSqrt = 1/sqrt(aa1)*aa0Sqrt, aa1Sqrt=sqrt(aa1)*aa0InvSqrt, 0, 0
  vmovddup  %xmm4, %xmm1                # xmm1 = {aa1InvSqrt aa1InvSqrt 0 0}
  vpermilpd $1, %xmm4,   %xmm4          # xmm4  = aa1Sqrt, aa1InvSqrt, 0, 0
  vmovupd   %xmm4, 16(%rcx,%rbx)        # x1[1] = (aa1Sqrt, aa1InvSqrt)

  vmulpd    %xmm0, %xmm2,%xmm2          # xmm2  = {f.re*=aa0InvSqrt f.im*=aa0InvSqrt 0 0}
  vmovupd   %xmm2, 16(%rcx)             # triang[1] = {f.re f.im }

  vunpckhpd %xmm2, %xmm2, %xmm3         # ymm3 = f.im f.im 0 0
  vmovddup  %xmm2, %xmm2                # xmm2 = f.re f.re 0 0

  cmp       $1, %edx
  jg       .even_loop                   # x1 was not a last row
.done:
  xor %eax, %eax
  test $3, %r11
  sete %al

  vmovaps   (%rsp), %xmm6
  vmovaps 16(%rsp), %xmm7
  vmovaps 32(%rsp), %xmm8
  vmovaps 48(%rsp), %xmm9
  vmovaps 64(%rsp), %xmm10
  vmovaps 80(%rsp), %xmm11
  vmovaps 96(%rsp), %xmm12
  addq  $112, %rsp
  popq  %rdi
  popq  %rsi
  popq  %rbx
  vzeroupper
  retq
  .seh_handlerdata
  .text
  .seh_endproc
                                        # -- End function
  .addrsig
