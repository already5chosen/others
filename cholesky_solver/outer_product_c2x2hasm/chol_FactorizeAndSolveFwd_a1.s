  .text
  .def   @feat.00;
  .scl  3;
  .type 0;
  .endef
  .globl  @feat.00
.set @feat.00, 0
  .file "chol_FactorizeAndSolveFwd_a.s"
  .def   chol_FactorizeAndSolveFwd_a;
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
  .globl  chol_FactorizeAndSolveFwd_a
  .p2align  4, 0x90
# extern "C" bool chol_FactorizeAndSolveFwd_a(double* triang, unsigned N, double* result)
chol_FactorizeAndSolveFwd_a:
# RCX - triang
# RDX - N
# R8  - result
.seh_proc chol_FactorizeAndSolveFwd_a
  pushq %rbx
  .seh_pushreg %rbx
  pushq %rsi
  .seh_pushreg %rsi
  pushq %rdi
  .seh_pushreg %rdi
  pushq %rbp
  .seh_pushreg %rbp
  pushq %r12
  .seh_pushreg %r12
  subq  $16*9, %rsp
  .seh_stackalloc 16*9
  vmovapd %xmm14,      16*8(%rsp)
  .seh_savexmm %xmm14, 16*8
  vmovapd %xmm13,      16*7(%rsp)
  .seh_savexmm %xmm13, 16*7
  vmovapd %xmm12,      16*6(%rsp)
  .seh_savexmm %xmm12, 16*6
  vmovapd %xmm11,      16*5(%rsp)
  .seh_savexmm %xmm11, 16*5
  vmovapd %xmm10,      16*4(%rsp)
  .seh_savexmm %xmm10, 16*4
  vmovapd %xmm9,       16*3(%rsp)
  .seh_savexmm %xmm9,  16*3
  vmovapd %xmm8,       16*2(%rsp)
  .seh_savexmm %xmm8,  16*2
  vmovapd %xmm7,       16*1(%rsp)
  .seh_savexmm %xmm7,  16*1
  vmovapd %xmm6,       16*0(%rsp)
  .seh_savexmm %xmm6,  16*0
  .seh_endprologue

  xor %r11d, %r11d
  shl $3,  %edx                         # edx = rlenx = rlen*sizeof(double) = N*sizeof(double)
  jz  .done
  test $8, %edx
  jz .even_loop_entry
   # special handling for the first row of matrix with odd number of elements
   # process top row
   add        $8,  %edx                 # edx = rlenx += 1*sizeof(double) = (N+1)*sizeof(double)
   mov        %edx,%ebp
   and        $16, %ebp                 # RBP  = xix = (rlen & 2) * sizeof(double)
   vmovsd    8(%rcx,%rbp), %xmm0        # ymm0 = aa = triang[xi+1], 0,0,0 // diagonal element
   # check that we are positive definite
   vcmpltsd .LDIAG_MIN(%rip), %xmm0, %xmm1             # xmm1 = ltmsk =[aa<DIAG_MIN ...]
   vblendvpd  %xmm1, .LDIAG_SUBST(%rip), %xmm0,%xmm0   # xmm0 = aa < DIAG_MIN ? DIAG_SUBST : aa
   vmovmskpd  %xmm1, %r11d
   and        $1,    %r11d
   vmovsd     .LONE(%rip),  %xmm1       # xmm1 = 1.0,zero
   mov        %r8,   %rax               # rax  = result
   vsqrtsd    %xmm0, %xmm0, %xmm14      # xmm14= aaSqrt=sqrt(aa) ...
   vdivsd     %xmm14,%xmm1, %xmm0       # xmm0 = aaInvSqrt=1/sqrt(aa) ...
   vmulsd    8(%r8,%rbp),%xmm0, %xmm1   # ymm1 = r0_re = result[xix+1+0] * aaInvSqrt, 0, 0, 0
   vmulsd   40(%r8,%rbp),%xmm0, %xmm2   # ymm2 = r0_im = result[xix+1+4] * aaInvSqrt, 0, 0, 0
   # store first values
   vmovsd     %xmm14, 8(%rcx,%rbp)      # triang[xi+1+0] = aaSqrt
   vmovsd     %xmm0, 40(%rcx,%rbp)      # triang[xi+1+4] = aaInvSqrt - store inverse of diagonal at imag()
   vmovsd     %xmm1,  8(%r8, %rbp)      # result[xix+1+0] = r0_re
   vmovsd     %xmm2, 40(%r8, %rbp)      # result[xix+1+4] = r0_im

   sub $16,   %edx                      # edx = rlenx -= 2*sizeof(double) = (N-1)*sizeof(double)
   jz .done                             # N==1, nothing more to do

   # Combine multiplication by invSqrt with Forward propagation
   mov        %r8,  %rsi
   sub        %rcx, %rsi                # rsi  = (result-triang)*sizeof(double)
   lea (%rcx, %rbp,4),%rbx              # rbx  = x = &triang[xi*4]

   vmovddup   %xmm0, %xmm0              # ymm0 = aa0InvSqrt x4
   vmovddup   %xmm1, %xmm1              # ymm1 = r0_re      x2
   vmovddup   %xmm2, %xmm2              # ymm2 = r0_im      x2
   vperm2f128 $0,%ymm0,%ymm0,%ymm0      # ymm0 = aa0InvSqrt x4
   vperm2f128 $0,%ymm1,%ymm1,%ymm1      # ymm1 = r0.re      x4
   vperm2f128 $0,%ymm2,%ymm2,%ymm2      # ymm2 = r0.im      x4

   mov        %edx, %eax                # eax  = cnt = rlenx
   .odd_top_rows_loop:
    vmulpd    (%rbx), %ymm0, %ymm3      # ymm3 = x_re = x[c+0]*aa0InvSqrt
    vmulpd  32(%rbx), %ymm0, %ymm4      # ymm4 = x_im = x[c+4]*aa0InvSqrt
    vmovupd   (%rbx,%rsi),   %ymm5      # ymm5 = r_re = pr[c+0]
    vmovupd 32(%rbx,%rsi),   %ymm6      # ymm6 = r_im = pr[c+4]
    vmovupd   %ymm3,  (%rbx)            # x[c+0] = x_re
    vmovupd   %ymm4,32(%rbx)            # x[c+4] = x_im

    vmulpd    %ymm3,%ymm1,   %ymm7      # ymm7 = x_re*r0_re
    vsubpd    %ymm7,%ymm5,   %ymm5      # ymm5 = r_re -= x_re*r0_re
    vmulpd    %ymm3,%ymm2,   %ymm7      # ymm7 = x_re*r0_im
    vsubpd    %ymm7,%ymm6,   %ymm6      # ymm6 = r_im -= x_re*r0_im

    vmulpd    %ymm4,%ymm2,   %ymm7      # ymm7 = x_im*r0_im
    vaddpd    %ymm7,%ymm5,   %ymm5      # ymm5 = r_re += x_im*r0_im
    vmulpd    %ymm4,%ymm1,   %ymm7      # ymm7 = x_im*r0_re
    vsubpd    %ymm7,%ymm6,   %ymm6      # ymm6 = r_im -= x_im*r0_re

    vmovupd   %ymm5,  (%rbx,%rsi)       # pr[c+0] = r_re
    vmovupd   %ymm6,32(%rbx,%rsi)       # pr[c+4] = r_im
    add     $64,  %rbx                  # c += 8
   sub      $32,  %rax                  # cnt -= 4
   jg      .odd_top_rows_loop

   # repeat store of first values since when xix==0 originals are destroyed
   vmovsd     %xmm14, 8(%rcx,%rbp)      # triang[xi+1+0] = aaSqrt
   vmovsd     %xmm0, 40(%rcx,%rbp)      # triang[xi+1+4] = aaInvSqrt - store inverse of diagonal at imag()
   vmovsd     %xmm1,  8(%r8, %rbp)      # result[xix+1+0] = r0_re
   vmovsd     %xmm2, 40(%r8, %rbp)      # result[xix+1+4] = r0_im
   lea (%r8,  %rbp,4),%r8               # r8   = result += xi*4

   # Subtract outer product of top row from lower part of the matrix
   # Process two output rows together
   lea (%rcx, %rbp,4),%r9               # r9  = x = &triang[xi*4]
   mov      %rbx, %rcx                  # rcx = triang = &triang[((N+3) & -4)*2]
   sub      %r9,  %rbx                  # rbx = yix = (y0-x)*sizeof(double)
   mov      %edx, %r10d                 # r10 = clenx = clen*sizeof(double) = rlenx = (N-1)*sizeof(double)
  .caxpy1x2_outer_loop:
    xor      $16,  %ebp                 # rbp = xix = ^= 1*sizeof(double)
    vbroadcastsd   (%r9 ,%rbp), %ymm0   # ymm0 = f0_re = x[xi+0+0]  x4
    vbroadcastsd  8(%r9 ,%rbp), %ymm2   # ymm2 = f1_re = x[xi+1+0]  x4
    vbroadcastsd 32(%r9 ,%rbp), %ymm1   # ymm1 = f0_im = x[xi+0+4]  x4
    vbroadcastsd 40(%r9 ,%rbp), %ymm3   # ymm3 = f1_im = x[xi+1+4]  x4
    lea     (%r10d,%ebp),%eax           # rax = ylenx = (clen+xi)*sizeof(double) = (y1-y0)/2*sizeof(double)
    mov     %r9,  %rsi                  # rsi = &x[c]
    lea    (%ebx,%eax,2),%edi           # edi = (y1-x)*sizeof(double) = ((y0-x)+(y1-y0))*sizeof(double)
    .caxpy1x2_inner_loop:
      # y0[c] = y0[c] - x[c]*conj(f0)
      # y1[c] = y1[c] - x[c]*conj(f1)
      vmovupd   (%rsi),      %ymm4      # ymm4 = x_re  = x[c*8+0]
      vmovupd   (%rsi,%rbx), %ymm6      # ymm6 = y0_re = y0[c*8+0]
      vmovupd 32(%rsi,%rbx), %ymm7      # ymm7 = y0_im = y0[c*8+4]
      vmovupd   (%rsi,%rdi), %ymm8      # ymm8 = y1_re = y1[c*8+0]
      vmovupd 32(%rsi,%rdi), %ymm9      # ymm9 = y1_im = y1[c*8+4]

      vmulpd  %ymm4 ,%ymm0,  %ymm5
      vsubpd  %ymm5 ,%ymm6,  %ymm6      # ymm6 = y0_re -= x_re*f0_re
      vmulpd  %ymm4 ,%ymm1,  %ymm5
      vaddpd  %ymm5 ,%ymm7,  %ymm7      # ymm7 = y0_im += x_re*f0_im
      vmulpd  %ymm4 ,%ymm2,  %ymm5
      vsubpd  %ymm5 ,%ymm8 , %ymm8      # ymm8 = y1_re -= x_re*f1_re
      vmulpd  %ymm4 ,%ymm3,  %ymm5
      vaddpd  %ymm5 ,%ymm9 , %ymm9      # ymm9 = y1_im += x_re*f1_im

      vmovupd 32(%rsi),      %ymm4      # ymm4 = x_im = x[c*8+4]
      vmulpd  %ymm4 ,%ymm1,  %ymm5
      vsubpd  %ymm5 ,%ymm6,  %ymm6      # ymm6 = y0_re -= x_im*f0_im
      vmulpd  %ymm4 ,%ymm0,  %ymm5
      vsubpd  %ymm5 ,%ymm7,  %ymm7      # ymm7 = y0_im -= x_im*f0_re
      vmulpd  %ymm4 ,%ymm3,  %ymm5
      vsubpd  %ymm5 ,%ymm8 , %ymm8      # ymm8 = y1_re -= x_im*f1_im
      vmulpd  %ymm4 ,%ymm2,  %ymm5
      vsubpd  %ymm5 ,%ymm9 , %ymm9      # ymm9 = y1_im -= x_im*f1_re

      vmovupd %ymm6,   (%rsi,%rbx)      # y0[c*8+0] = y0_re
      vmovupd %ymm7, 32(%rsi,%rbx)      # y0[c*8+4] = y0_im
      vmovupd %ymm8 ,  (%rsi,%rdi)      # y1[c*8+0] = y1_re
      vmovupd %ymm9 ,32(%rsi,%rdi)      # y1[c*8+4] = y1_im
      add    $64, %rsi                  # &x[c] += 8
    sub      $32, %eax                  # ylenx -= 4*sizeof(double)
    jnz .caxpy1x2_inner_loop

    lea (%r9, %rbp, 4),%r9              # R9   = x   += xi*4
    lea (%ebx,%r10d,4),%ebx             # RBX  = yix += clenx*4

  sub      $16, %r10d                   # clenx -= 2*sizeof(double)
  jnz .caxpy1x2_outer_loop
  jmp .even_loop_entry

.even_loop:
# RAX  - xlenx = (x1-x0)*sizeof(double)
# RCX  - triang
# RDX  - rlenx = rlen*sizeof(double) >= 4*sizeof(double)
# R8   - result
# R9   = xix = index of diagonal element of x0 in triang[] * sizeof(double)
# YMM0 - aa0InvSqrt x2, 0,0
# YMM1 - aa1InvSqrt x2, 0,0
# YMM2 - f_re       x2, 0,0
# YMM3 - f_im       x2, 0,0
# YMM4 - r0_re      x2, 0,0
# YMM5 - r0_im      x2, 0,0
# YMM6 - r1_re      x2, 0,0
# YMM7 - r1_im      x2, 0,0
# YMM13- x aa1Sqrt,     0,0
# YMM14- aa0Sqrt x,     0,0
  # Combine handling of pair of top rows with forward propagation of pair of results
  vperm2f128 $0,%ymm0,%ymm0,%ymm0       # ymm0 = aa0InvSqrt x4
  vperm2f128 $0,%ymm1,%ymm1,%ymm1       # ymm1 = aa1InvSqrt x4
  vperm2f128 $0,%ymm2,%ymm2,%ymm2       # ymm2 = f_re       x4
  vperm2f128 $0,%ymm3,%ymm3,%ymm3       # ymm3 = f_im       x4
  vperm2f128 $0,%ymm4,%ymm4,%ymm4       # ymm4 = r0_re      x4
  vperm2f128 $0,%ymm5,%ymm5,%ymm5       # ymm5 = r0_im      x4
  vperm2f128 $0,%ymm6,%ymm6,%ymm6       # ymm6 = r1_re      x4
  vperm2f128 $0,%ymm7,%ymm7,%ymm7       # ymm7 = r1_im      x4
  vblendpd   $1,%xmm14,%xmm13,%xmm14    # ymm14= aa0Sqrt aa1Sqrt 0 0

  lea (%rcx,%r9,4),%rbx                 # rbx  = x0 = &triang[(rlen & 2)*4]
  mov  %r8,        %rsi
  sub  %rcx,       %rsi                 # rsi  = drx = result-triang
  mov  %edx,       %edi
  shr  $5,         %edi                 # edi  = cqlen = rlen/4==rlenx/32
  .top_rows_loop:
    ## ax0 = x0[k] * aa0InvSqrt;
    ## ax1 = x1[k];
    ## x0[k] = ax0;
    ## x1[k] = (ax1 - ax0*conj(f))*aa1InvSqrt;
    ## pr[k] -= x0[c]*r0 + x1[c]*r1;
    vmulpd    (%rbx),    %ymm0, %ymm8   # ymm8 = x0_re = x0[c*8+0]*aa0InvSqrt
    vmulpd  32(%rbx),    %ymm0, %ymm9   # ymm9 = x0_im = x0[c*8+4]*aa0InvSqrt
    vmovupd   (%rbx,%rax),      %ymm10  # ymm10= x1_re = x1[c*8+0]
    vmovupd 32(%rbx,%rax),      %ymm11  # ymm11= x1_im = x1[c*8+4]
    vmovupd   (%rbx,%rsi),      %ymm12  # ymm12= r_re  = pr[c*8+0]

    vmulpd    %ymm8, %ymm2, %ymm13
    vsubpd    %ymm13,%ymm10,%ymm10      # ymm10= x1_re -= x0_re*f_re
    vmulpd    %ymm8, %ymm3, %ymm13
    vaddpd    %ymm13,%ymm11,%ymm11      # ymm11= x1_im += x0_re*f_im
    vmulpd    %ymm8, %ymm4, %ymm13
    vsubpd    %ymm13,%ymm12,%ymm12      # ymm12= r_re  -= x0_re*r0_re

    vmulpd    %ymm9, %ymm3, %ymm13
    vsubpd    %ymm13,%ymm10,%ymm10      # ymm10= x1_re -= x0_im*f_im
    vmulpd    %ymm9, %ymm2, %ymm13
    vsubpd    %ymm13,%ymm11,%ymm11      # ymm11= x1_im -= x0_im*f_re
    vmulpd    %ymm9, %ymm5, %ymm13
    vaddpd    %ymm13,%ymm12,%ymm12      # ymm12= r_re  += x0_im*r0_im

    vmovupd 32(%rbx,%rsi),      %ymm13  # ymm13= r_im  = pr[c*8+4]

    vmovupd   %ymm8,  (%rbx)            # x0[c*8+0] = x0_re
    vmovupd   %ymm9,32(%rbx)            # x0[c*8+4] = x0_im

    vmulpd    %ymm1, %ymm10,%ymm10      # ymm10= x1_re *= aa1InvSqrt
    vmulpd    %ymm1, %ymm11,%ymm11      # ymm11= x1_im *= aa1InvSqrt

    vmovupd   %ymm10,  (%rbx,%rax)      # x1[c*8+0] = x1_re
    vmovupd   %ymm11,32(%rbx,%rax)      # x1[c*8+4] = x1_im

    vmulpd    %ymm8, %ymm5, %ymm8
    vsubpd    %ymm8, %ymm13,%ymm13      # ymm13= r_im  -= x0_re*r0_im
    vmulpd    %ymm9, %ymm4, %ymm8
    vsubpd    %ymm8, %ymm13,%ymm13      # ymm13= r_im  -= x0_im*r0_re

    vmulpd    %ymm10,%ymm6, %ymm8
    vsubpd    %ymm8, %ymm12,%ymm12      # ymm12= r_re  -= x1_re*r1_re
    vmulpd    %ymm10,%ymm7, %ymm8
    vsubpd    %ymm8, %ymm13,%ymm13      # ymm13= r_im  -= x1_re*r1_im

    vmulpd    %ymm11,%ymm7, %ymm8
    vaddpd    %ymm8, %ymm12,%ymm12      # ymm12= r_re  += x1_im*r1_im
    vmulpd    %ymm11,%ymm6, %ymm8
    vsubpd    %ymm8, %ymm13,%ymm13      # ymm13= r_im  -= x1_im*r1_re

    vmovupd   %ymm12,  (%rbx,%rsi)      # pr[c*8+0] = r_re
    vmovupd   %ymm13,32(%rbx,%rsi)      # pr[c*8+4] = r_im

    add $64, %rbx                       # c += 1
  dec %edi
  jnz .top_rows_loop

  # store last diag and near diag
  lea (%rcx,%r9),   %rsi                # rsi  = x0 = &triang[xi]
  vblendpd $1, %xmm14,%xmm2, %xmm2      # ymm2 = aa0Sqrt,   f_re,0,0
  vblendpd $1, %xmm0, %xmm3, %xmm3      # ymm3 = aa0InvSqrt,f_im,0,0
  vmovupd  %xmm2,   (%rsi)              # x0[0:1] = [aa0Sqrt,    f_re]
  vmovupd  %xmm3, 32(%rsi)              # x0[4:5] = [aa0InvSqrt, f_im]
  vmovupd  %xmm14,  (%rsi,%rax)         # x1[0:1] = [x        aa1Sqrt]
  vmovsd   %xmm1, 40(%rsi,%rax)         # x1[5]   = aa1InvSqrt
  # store last results
  lea (%r8, %r9),   %rsi                # rsi  = pr = &result[xi] points to result[0].real()
  vblendpd $1, %xmm4, %xmm6, %xmm6      # ymm6 = r0_re,r1_re,0,0
  vblendpd $1, %xmm5, %xmm7, %xmm7      # ymm7 = r0_im,r1_im,0,0
  vmovupd  %xmm6,   (%rsi)              # pr[0:1] = [r0_re r1_re]
  vmovupd  %xmm7, 32(%rsi)              # pr[4:5] = [r0_im r1_im]

  lea (%rcx,%r9, 4), %r10               # r10  = x0 = &triang[(rlen & 2)*4]
  lea (%r8, %r9, 4), %r8                # r8   = result += (rlen & 2)*4
  lea (%rcx,%rax,2), %rcx               # rcx  = triang += xlenx*2

  # subtract outer product of two top rows from lower part of the matrix
  # process two output rows together
# RAX  - xlenx = (x1-x0)*sizeof(double)
# RCX  - triang
# RDX  - rlenx = rlen*sizeof(double) >= 4*sizeof(double)
# R8   - result
# R10  - x0
# R11  - ltMask
# RBX, RSI, RDI, RBP, R9, R12 - free
  lea (,%edx,4),%edi                    # RDI  = yix == (y0-x0)*sizeof(double) = rlenx*4
  sub $16, %edx                         # RDX  = rlenx -= 2*sizeof(double)
  jz .caxpy2x2_outer_loop_done
  mov %edx,%esi                         # RSI  = clenx = clen*sizeof(double) = rlenx
  mov %edx,%r9d
  and $16, %r9d                         # R9   = xix = (clen & 2)*sizeof(double)
  .caxpy2x2_outer_loop:
    lea (%esi,%r9d),   %r12d            # R12  = ylenx = ((clen+2) & -4)*sizeof(double) = clenx+xix = ((y1-y0)/2)*sizeof(double)
    lea (%edi,%r12d,2),%ebp             # RBP  = y1-x0 = (y0-x0)+(y1-y0) = (y0-x0)+ylenx*2
    lea (%r10,%r9),    %rbx             # RBX  = &x0[xi]
    vbroadcastsd   (%rbx),      %ymm0   # ymm0 = f00_re = x0[xi+0+0] x4
    vbroadcastsd 32(%rbx),      %ymm1   # ymm1 = f00_im = x0[xi+0+4] x4
    vbroadcastsd  8(%rbx),      %ymm2   # ymm2 = f01_re = x0[xi+1+0] x4
    vbroadcastsd 40(%rbx),      %ymm3   # ymm3 = f01_im = x0[xi+1+4] x4

    vbroadcastsd   (%rbx,%rax), %ymm4   # ymm4 = f10_re = x1[xi+0+0] x4
    vbroadcastsd 32(%rbx,%rax), %ymm5   # ymm5 = f10_im = x1[xi+0+4] x4
    vbroadcastsd  8(%rbx,%rax), %ymm6   # ymm6 = f11_re = x1[xi+1+0] x4
    vbroadcastsd 40(%rbx,%rax), %ymm7   # ymm7 = f11_im = x1[xi+1+4] x4

    mov %r10, %rbx                      # RBX  = &x0[c] = x1

#  .p2align  3, 0x90
    .caxpy2x2_inner_loop:
      # y0[c] = y0[c] - x0[c]*conj(f00) - x1[c]*conj(f10);
      # y1[c] = y1[c] - x0[c]*conj(f01) - x1[c]*conj(f11);
      vmovupd (%rbx),        %ymm12     # ymm12= x0_re = x0[c*8+0]
      vmovupd (%rbx,%rdi),   %ymm8      # ymm8 = y0_re = y0[c*8+0]
      vmovupd 32(%rbx,%rdi), %ymm9      # ymm9 = y0_im = y0[c*8+4]
      vmovupd (%rbx,%rbp),   %ymm10     # ymm10= y1_re = y1[c*8+0]
      vmovupd 32(%rbx,%rbp), %ymm11     # ymm11= y1_im = y1[c*8+4]

      vmulpd  %ymm12,%ymm0,  %ymm13
      vsubpd  %ymm13,%ymm8,  %ymm8      # ymm8 = y0_re -= x0_re*f00_re
      vmulpd  %ymm12,%ymm1,  %ymm13
      vaddpd  %ymm13,%ymm9,  %ymm9      # ymm9 = y0_im += x0_re*f00_im
      vmulpd  %ymm12,%ymm2,  %ymm13
      vsubpd  %ymm13,%ymm10, %ymm10     # ymm10= y1_re -= x0_re*f01_re
      vmulpd  %ymm12,%ymm3,  %ymm13
      vaddpd  %ymm13,%ymm11, %ymm11     # ymm11= y1_im += x0_re*f01_im

      vmovupd (%rbx,%rax),   %ymm12     # ymm12= x1_re = x1[c*8+0]
      vmulpd  %ymm12,%ymm4,  %ymm13
      vsubpd  %ymm13,%ymm8,  %ymm8      # ymm8 = y0_re -= x1_re*f10_re
      vmulpd  %ymm12,%ymm5,  %ymm13
      vaddpd  %ymm13,%ymm9,  %ymm9      # ymm9 = y0_im += x1_re*f10_im
      vmulpd  %ymm12,%ymm6,  %ymm13
      vsubpd  %ymm13,%ymm10, %ymm10     # ymm10= y1_re -= x1_re*f11_re
      vmulpd  %ymm12,%ymm7,  %ymm13
      vaddpd  %ymm13,%ymm11, %ymm11     # ymm11= y1_im += x1_re*f11_im

      vmovupd 32(%rbx),      %ymm12     # ymm12= x0_im = x0[c*8+4]
      vmulpd  %ymm12,%ymm1,  %ymm13
      vsubpd  %ymm13,%ymm8,  %ymm8      # ymm8 = y0_re -= x0_im*f00_im
      vmulpd  %ymm12,%ymm0,  %ymm13
      vsubpd  %ymm13,%ymm9,  %ymm9      # ymm9 = y0_im -= x0_im*f00_re
      vmulpd  %ymm12,%ymm3,  %ymm13
      vsubpd  %ymm13,%ymm10, %ymm10     # ymm10= y1_re -= x0_im*f01_im
      vmulpd  %ymm12,%ymm2,  %ymm13
      vsubpd  %ymm13,%ymm11, %ymm11     # ymm11= y1_im -= x0_im*f01_re

      vmovupd 32(%rbx,%rax), %ymm12     # ymm12= x1_im = x1[c*8+4]
      vmulpd  %ymm12,%ymm5,  %ymm13
      vsubpd  %ymm13,%ymm8,  %ymm8      # ymm8 = y0_re -= x1_im*f10_im
      vmulpd  %ymm12,%ymm4,  %ymm13
      vsubpd  %ymm13,%ymm9,  %ymm9      # ymm9 = y0_im -= x1_im*f10_re
      vmulpd  %ymm12,%ymm7,  %ymm13
      vsubpd  %ymm13,%ymm10, %ymm10     # ymm10= y1_re -= x1_im*f11_im
      vmulpd  %ymm12,%ymm6,  %ymm13
      vsubpd  %ymm13,%ymm11, %ymm11     # ymm11= y1_im -= x1_im*f11_re

      vmovupd %ymm8,   (%rbx,%rdi)      # y0[c*8+0] = y0_re
      vmovupd %ymm9, 32(%rbx,%rdi)      # y0[c*8+4] = y0_im
      vmovupd %ymm10,  (%rbx,%rbp)      # y1[c*8+0] = y1_re
      vmovupd %ymm11,32(%rbx,%rbp)      # y1[c*8+4] = y1_im
      add $64, %rbx                     # c += 1
    sub $32,   %r12d                    # ylenx -= 4*sizeof(double)
    jnz .caxpy2x2_inner_loop

    lea (%r10,%r9, 4),%r10              # R10  = x0    += xi*4
    lea (%edi,%esi,4),%edi              # RDI  = yix   += clenx*4

    xor $16,        %r9d                # R9   = xix = ((xi+2)%4)*sizeof(double)
  sub $16, %esi                         # RSI  = clenx -= 2*sizeof(double)
  jnz .caxpy2x2_outer_loop

  .caxpy2x2_outer_loop_done:
.even_loop_entry:
# RCX - triang
# RDX - rlenx = rlen*sizeof(double)
# R8  - result
# R11 - ltMask
# RAX, RBX, RDI, RSI, RBP, R9, R10, R12 - free
  # process diagonal and near-diagonal of 2 top rows
  mov  %edx,        %r9d
  and  $16,         %r9d                # r9   = xix = index of diagonal element of x0 in triang[] * sizeof(double)
  lea (%rcx,%r9),   %rbx                # rbx  = x0  = &triang[xi]
  lea (%edx,%r9d),  %eax                # eax  = (rlen+xi)*sizeof(double)
  add  %eax,        %eax                # rax  = xlenx = (rlen+xi)*2*sizeof(double) = (x1-x0)*sizeof(double)
  vmovsd   (%rbx),           %xmm0      # ymm0 = x0[0]        = aa0, 0,0,0
  vmulsd  8(%rbx,%rax),%xmm0,%xmm1      # ymm1 = x1[1].re*aa0 = aa1*aa0, 0,0,0
  vmovsd  8(%rbx),           %xmm2      # ymm2 = x0[1]        = f_re, 0,0,0
  vmovsd 40(%rbx),           %xmm3      # ymm3 = x0[5]        = f_im, 0,0,0
  # aa1 = aa1*aa0 - norm(f)
  vmulsd    %xmm2, %xmm2, %xmm4
  vsubsd    %xmm4, %xmm1, %xmm1         # ymm1 = aa1 -= f_re*f_re
  vmulsd    %xmm3, %xmm3, %xmm4
  vsubsd    %xmm4, %xmm1, %xmm1         # ymm1 = aa1 -= f_im*f_im
  vunpcklpd %xmm1, %xmm0, %xmm0         # ymm0 = aa0, aa1, 0, 0

  # check that we are positive definite
  vmovddup .LDIAG_MIN(%rip), %xmm1
  vcmpltpd  %xmm1, %xmm0, %xmm1         # xmm1 = ltmsk =[aa0<DIAG_MIN aa1<DIAG_MIN]
  vmovddup  .LDIAG_SUBST(%rip), %xmm4
  vblendvpd %xmm1, %xmm4, %xmm0,%xmm0   # xmm0 = aa[0:1] < DIAG_MIN ? DIAG_SUBST : aa[0:1]
  vmovmskpd %xmm1, %esi
  or        %esi,  %r11d
  vsqrtpd   %xmm0, %xmm14               # ymm14= sqrt(aa0), sqrt(aa1), 0, 0
  vmovddup  .LONE(%rip), %xmm1
  vdivpd    %xmm14,%xmm1,%xmm1          # xmm1 = [1/sqrt(aa0), 1/sqrt(aa1)]

  vmulsd    %xmm1, %xmm2,%xmm2          # xmm2  = f_re*=aa0InvSqrt 0 0 0
  vmulsd    %xmm1, %xmm3,%xmm3          # xmm3  = f_im*=aa0InvSqrt 0 0 0

  vmovddup  %xmm1, %xmm0                # ymm0 = 1/sqrt(aa0), 1/sqrt(aa0) = aa0InvSqrt x2, 0,0
  vunpckhpd %xmm1, %xmm1, %xmm1         # ymm1 = 1/sqrt(aa1), 1/sqrt(aa1) 0,0
  vmulsd    %xmm14,%xmm1, %xmm1         # ymm1 = aa1InvSqrt = 1/sqrt(aa1)*aa0Sqrt, 0,0,0
  vmulpd    %xmm0, %xmm14,%xmm13        # ymm13= x aa1Sqrt=sqrt(aa1)*aa0InvSqrt, 0,0
  vmovddup  %xmm1, %xmm1                # ymm1 = aa1InvSqrt x2, 0,0
  vmovddup  %xmm2, %xmm2                # ymm2 = f_re       x2, 0,0
  vmovddup  %xmm3, %xmm3                # ymm3 = f_im       x2, 0,0

  # process 2 top results
  lea (%r8, %r9),   %rsi                # rsi  = pr = points to result[0].real()
  vmovupd   (%rsi),       %xmm5         # ymm5 = r0_re,r1_re,0,0 = pr[0:1]
  vmovupd 32(%rsi),       %xmm6         # ymm6 = r0_im,r1_im,0,0 = pr[4:5]
  vunpcklpd %xmm6, %xmm5, %xmm4         # ymm4 = r0_re,r0_im,0,0
  vunpckhpd %xmm6, %xmm5, %xmm5         # ymm5 = r1_re,r1_im,0,0
  vmulpd    %xmm0, %xmm4, %xmm4         # ymm4 = r0_re,r0_im,0,0 *= aa0InvSqrt
  # r1 -= r0*f
  vmulpd    %xmm2, %xmm4, %xmm6         # ymm6 = r0_re*f_re,r0_im*f_re,0,0
  vpermilpd $1,    %xmm4, %xmm7         # ymm7 = r0_im,r0_re,0,0
  vmulpd    %xmm3, %xmm7, %xmm7         # ymm7 = r0_im*f_im,r0_re*f_im,0,0
  vaddsubpd %xmm7, %xmm6, %xmm6         # ymm6 = r0*f = r0_re*f_re-r0_im*f_im,r0_im*f_re+r0_re*f_im,0,0
  vsubpd    %xmm6, %xmm5, %xmm5         # ymm5 = r1_re,r1_im,0,0 -= r0*f
  # r1 *= aa1InvSqrt
  vmulpd    %xmm1, %xmm5, %xmm5         # ymm5 = r1_re,r1_im,0,0 *= aa1InvSqrt

  vmovddup  %xmm5,        %xmm6         # ymm6 = r1_re       x2, 0,0
  vunpckhpd %xmm5, %xmm5, %xmm7         # ymm7 = r1_im       x2, 0,0
  vunpckhpd %xmm4, %xmm4, %xmm5         # ymm5 = r0_im       x2, 0,0
  vmovddup  %xmm4,        %xmm4         # ymm4 = r0_re       x2, 0,0

  cmp       $16, %edx
  jg       .even_loop                   # x1 was not a last row

  # store last diag and near diag
  vblendpd $1, %xmm14,%xmm2, %xmm2      # ymm2 = aa0Sqrt,   f_re,0,0
  vblendpd $1, %xmm0, %xmm3, %xmm3      # ymm3 = aa0InvSqrt,f_im,0,0
  vmovupd  %xmm2,   (%rbx)              # x0[0:1] = [aa0Sqrt,    f_re]
  vmovupd  %xmm3, 32(%rbx)              # x0[4:5] = [aa0InvSqrt, f_im]
  vmovupd  %xmm13,  (%rbx,%rax)         # x1[0:1] = [x        aa1Sqrt]
  vmovsd   %xmm1, 40(%rbx,%rax)         # x1[5]   = aa1InvSqrt
  # store last results
  vblendpd $1, %xmm4, %xmm6, %xmm6      # ymm6 = r0_re,r1_re,0,0
  vblendpd $1, %xmm5, %xmm7, %xmm7      # ymm7 = r0_im,r1_im,0,0
  vmovupd  %xmm6,   (%rsi)              # pr[0:1] = [r0_re r1_re]
  vmovupd  %xmm7, 32(%rsi)              # pr[4:5] = [r0_im r1_im]

.done:
  xor %eax, %eax
  and  $3,  %r11b
  sete %al

  vmovaps 0*16(%rsp), %xmm6
  vmovaps 1*16(%rsp), %xmm7
  vmovaps 2*16(%rsp), %xmm8
  vmovaps 3*16(%rsp), %xmm9
  vmovaps 4*16(%rsp), %xmm10
  vmovaps 5*16(%rsp), %xmm11
  vmovaps 6*16(%rsp), %xmm12
  vmovaps 7*16(%rsp), %xmm13
  vmovaps 8*16(%rsp), %xmm14
  addq  $16*9, %rsp
  popq  %r12
  popq  %rbp
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
