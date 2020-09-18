  .text
  .def   @feat.00;
  .scl  3;
  .type 0;
  .endef
  .globl  @feat.00
.set @feat.00, 0
  .file "chol_SolveBwd_a.s"
  .def   chol_SolveBwd_a;
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
  .globl  chol_SolveBwd_a
  .p2align  4, 0x90
# extern "C" void chol_SolveBwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
chol_SolveBwd_a:
# RCX - x
# RDX - N
# R8  - triang
.seh_proc chol_SolveBwd_a
	subq	$24, %rsp
	.seh_stackalloc 24
	vmovapd	%xmm6, 0(%rsp)          # 16-byte Spill
	.seh_savexmm %xmm6, 0
	.seh_endprologue

  lea  1(%edx),%eax   # eax  = N+1
  imul %eax, %eax     # eax  = (N+1)*(N+1)
  shr  $1,   %eax     # eax  = trlen = (N+1)*(N+1)/2 = number of elements in triang array
  shl  $4,   %eax     # eax  = trlen*sizeof(DComplex)
  add  %rax, %r8      # r8   = triang = &triang[(N+1)*(N+1)/2]; // points past end of triang
  shl  $4,   %edx     # edx  = Nx = N*sizeof(DComplex)
  add  %rdx, %rcx     # rcx  = x += N;                          // points past end of x
  mov  %edx, %r11d    # r11d = Nx = N*sizeof(DComplex)

  mov  $32,  %edx     # edx = rlenx = rlen*sizeof(DComplex)=2*sizeof(DComplex)
  cmp  %edx, %r11d
  jb   .is_odd
    sub       $32,   %rcx   # rcx  = x -= 2;
    vmovupd   (%rcx),%xmm0  # ymm0 = acc0    = (x[-2].re, x[-2].im, 0, 0)
    vmovupd 16(%rcx),%xmm2  # ymm2 = acc1    = (x[-1].re, x[-1].im, 0, 0)
    sub       $64,   %r8    # r8   = triang -= rlenx*2
    jmp       .rloop_entry

  # process two rows per iteration
  .rloop:
    # RDX=EDX=rlenx=rlen*sizeof(DComplex)=y1-y0
    # rxlen > 2*sizeof(DComplex)
    sub %rdx, %r8
    sub %rdx, %r8           # r8   = triang -= rlenx*2
    vmovupd -32(%rcx),%xmm0 # ymm0 = acc0    = (x[-2].re, x[-2].im, 0, 0)
    vmovupd -16(%rcx),%xmm2 # ymm2 = acc1    = (x[-1].re, x[-1].im, 0, 0)
    lea 32(%r8),%rax        # rax  = y0 = &triang[2]
    mov        %rcx,  %r10          # r10  = x
    sub        %rax,  %rcx          # rcx  = x-y0
    vxorpd     %ymm4, %ymm4, %ymm4  # ymm4 = zeros
    vunpckhpd  %xmm4, %xmm0, %xmm1  # ymm1 = acc0_im = (x[-2].im, 0, 0, 0)
    vunpckhpd  %xmm4, %xmm2, %xmm3  # ymm3 = acc1_im = (x[-1].im, 0, 0, 0)
    vmovq      %xmm0, %xmm0         # ymm0 = acc0_re = (x[-2].re, 0, 0, 0)
    vmovq      %xmm2, %xmm2         # ymm2 = acc1_re = (x[-1].re, 0, 0, 0)
    lea      -32(%edx),%r9d         # r9d  = cntx = (rlen-2)*sizeof(DComplex)
    .cloop:
      # acc0 -= x[c] * conj(y0[c]);
      # acc1 -= x[c] * conj(y1[c]);
      vmovupd  (%rax, %rcx), %ymm4     # ymm4 = vx  = x[c] { re im re im }
      vmovupd  (%rax),       %ymm5     # ymm5 = y0[c]      { re im re im }
      vpermilpd    $5,    %ymm4, %ymm6 # ymm5 = vxx = x[c] { im re im re }
      vfnmadd231pd %ymm4, %ymm5, %ymm0 # ymm0 = acc0_re -= x[c]*y0[c] { re*re im*im re*re im*im }
      vfnmadd231pd %ymm6, %ymm5, %ymm1 # ymm1 = acc0_im -= x[c]*y0[c] { im*re re*im im*re re*im }

      vmovupd  (%rax, %rdx), %ymm5     # ymm5 = y1[c]      { re im re im }
      add      $32,   %rax             # eax  = y0 += 2
      vfnmadd231pd %ymm5, %ymm4, %ymm2 # ymm2 = acc1_re -= x[c]*y1[c] { re*re im*im re*re im*im }
      vfnmadd231pd %ymm5, %ymm6, %ymm3 # ymm3 = acc1_im -= x[c]*y1[c] { im*re re*im im*re re*im }
    sub        $32,   %r9d             # r9d  = cntx -= 2*sizeof(DComplex)
    jnz .cloop
    # reduce accumulators and pack re/im
    vhaddpd    %ymm2, %ymm0, %ymm2         # ymm2 = { acc0a_re acc1a_re acc0b_re acc1b_re }
    vhsubpd    %ymm3, %ymm1, %ymm1         # ymm1 = { acc0a_im acc1a_im acc0b_im acc1b_im }
    vunpcklpd  %ymm1, %ymm2, %ymm0         # ymm0 = { acc0a_re acc0a_im acc0b_re acc0b_im }
    vunpckhpd  %ymm1, %ymm2, %ymm1         # ymm1 = { acc1a_re acc1a_im acc1b_re acc1b_im }
    vperm2f128 $1,    %ymm0, %ymm0, %ymm2  # ymm2 = { acc0b_re acc0b_im x x }
    vaddpd     %xmm2, %xmm0, %xmm0         # ymm0 = { acc0_re acc0_im 0 0 }
    vperm2f128 $1,    %ymm1, %ymm1, %ymm2  # ymm2 = { acc1b_re acc1b_im x x }
    vaddpd     %xmm2, %xmm1, %xmm2         # ymm2 = { acc1_re acc1_im 0 0 }
    lea        -32(%r10), %rcx             # rcx = x -= 2

    .rloop_entry:
    # R8   = y0 = triang
    # RDX=EDX=rlenx=rlen*sizeof(DComplex)=y1-y0
    # RCX  = &x[-2]
    # YMM0 = acc0 = (x[-2].re, x[-2].im, 0, 0)
    # YMM2 = acc1 = (x[-1].re, x[-1].im, 0, 0)
    vmovddup 24(%r8,%rdx),  %xmm1     # ymm1 = { y1[1].im y1[1].im 0 0 } // imag() of diag element contains inverse of it's real()
    vmulpd   %xmm2, %xmm1,  %xmm1     # acc1 *= y1[1].im
    # acc0 -= acc1*conj(y0[1]);
    vmovddup 16(%r8),       %xmm2     # ymm2 = { y0[1].re y0[1].re 0 0 }
    vmovddup 24(%r8),       %xmm3     # ymm3 = { y0[1].im y0[1].im 0 0 }
    vmovddup  8(%r8),       %xmm4     # ymm4 = { y0[0].im y0[0].im 0 0 }
    vmulpd    %xmm1,%xmm2,  %xmm2     # ymm2 = { acc1.re*y0[1].re  acc1.im*y0[1].re 0 0 ]
    vsubpd    %xmm2,%xmm0,  %xmm0     # ymm0 = acc0.re -= acc1.re*y0[1].re acc0.im -= acc1.im*y0[1].re
    vpermilpd $1,   %xmm1,  %xmm2     # ymm2 = { acc1.im acc1.re 0 0 }
    vmulpd    %xmm3,%xmm2,  %xmm2     # ymm2 = { acc1.im*y0[1].im  acc1.re*y0[1].im 0 0 ]
    vaddsubpd %xmm2,%xmm0,  %xmm0     # ymm0 = acc0.re -= acc1.im*y0[1].im acc0.im += acc1.re*y0[1].im
    vmulpd    %xmm4,%xmm0,  %xmm0     # ymm0 = acc0 *= y0[0].imag();

    vmovupd   %xmm0,  (%rcx)          # x[0] = acc0
    vmovupd   %xmm1,16(%rcx)          # x[1] = acc1
    add       $32,    %edx            # rlenx += 2*sizeof(DComplex);
  cmp %r11d,  %edx
  jbe .rloop          # rlenx <= Nx


  .is_odd:
  test $16, %r11d
  jz .done
    # handle the first row of matrix with odd number of elements
    vmovupd -16(%rcx),%xmm0 # ymm0 = acc1    = (x[-1].re, x[-1].im, 0, 0)
    sub     %rdx, %r8       # R8   = triang -= N+1
    sub     $32,  %edx      # RDX  = rlenx -= 2*sizeof(DComplex)
    jz      .last_r_loop_done
      lea   32(%r8), %rax   # RAX = y = &triang[2]
      mov   %rcx,  %r10     # r10  = x
      sub   %rax,  %rcx     # rcx  = x - y
      vxorpd     %ymm1, %ymm1, %ymm1  # ymm1 = zeros
      vunpckhpd  %xmm1, %xmm0, %xmm1  # ymm1 = acc_im = (x[-1].im, 0, 0, 0)
      vmovq      %xmm0, %xmm0         # ymm0 = acc_re = (x[-1].re, 0, 0, 0)
      .last_r_loop:
        # acc -= x[c] * conj(y[c]);
        vmovupd  (%rax, %rcx), %ymm2  # ymm2 = vx  = x[c] { re im re im }
        vmovupd  (%rax),       %ymm3  # ymm3 = y[c]       { re im re im }
        add       $32,  %rax          # eax  = y += 2
        vpermilpd $5,   %ymm2, %ymm4  # ymm4 = vxx = x[c] { im re im re }
        vfnmadd231pd %ymm2, %ymm3, %ymm0 # ymm0 = acc0_re -= x[c]*y[c] { re*re im*im re*re im*im }
        vfnmadd231pd %ymm4, %ymm3, %ymm1 # ymm1 = acc0_im -= x[c]*y[c] { im*re re*im im*re re*im }
      sub   $32, %edx                 # edx = rlenx -= 2*sizeof(DComplex)
      jnz  .last_r_loop
      # reduce accumulators and pack re/im
      vhaddpd    %ymm0, %ymm0, %ymm0         # ymm0 = { acca_re acca_re accb_re accb_re }
      vhsubpd    %ymm1, %ymm1, %ymm1         # ymm1 = { acca_im acca_im accb_im accb_im }
      vblendpd   $0xA,  %ymm1, %ymm0, %ymm0  # ymm0 = { acca_re acca_im accb_re accb_im }
      vperm2f128 $1,    %ymm0, %ymm0, %ymm1  # ymm1 = { accb_re accb_im x x }
      vaddpd     %xmm1, %xmm0, %xmm0         # ymm0 = { acc_re  acc_im 0 0 }
      mov %r10, %rcx                         # rcx  = x
    .last_r_loop_done:
    vmovddup 24(%r8),       %xmm1     # ymm1 = { y1[1].im y1[1].im 0 0 } // imag() of diag element contains inverse of it's real()
    vmulpd   %xmm1, %xmm0,  %xmm0     # acc *= y1[1].im
    vmovupd  %xmm0, -16(%rcx)         # x[-1] = acc
  .done:

  vmovaps	0(%rsp), %xmm6
  add     $24, %rsp
  vzeroupper
  retq
.seh_handlerdata
  .text
  .seh_endproc
                                        # -- End function
  .addrsig
