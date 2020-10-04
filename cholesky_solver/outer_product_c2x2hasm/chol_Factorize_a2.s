  .text
  .def   @feat.00;
  .scl  3;
  .type 0;
  .endef
  .globl  @feat.00
.set @feat.00, 0
  .file "chol_Factorize_a2.s"
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
# extern "C" bool chol_Factorize_a(double* triang, unsigned N)
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
  pushq %rbp
  .seh_pushreg %rbp
  subq  $16*8+8, %rsp
  .seh_stackalloc 16*8+8
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
   vsqrtsd    %xmm0, %xmm0, %xmm7       # xmm7 = aaSqrt=sqrt(aa) ...
   vdivsd     %xmm7, %xmm1, %xmm0       # xmm0 = aaInvSqrt=1/sqrt(aa) ...

   sub $16,   %edx                      # edx = rlenx -= 2*sizeof(double) = (N-1)*sizeof(double)
   jz .odd_done                         # N==1, nothing more to do

   # Multiply top row by invSqrt
   lea (%rcx, %rbp,4),%rbx              # rbx  = x = &triang[xi*4]
   vmovddup   %xmm0, %xmm0              # ymm0 = aa0InvSqrt x2
   vperm2f128 $0,%ymm0,%ymm0,%ymm0      # ymm0 = aa0InvSqrt x4
   mov        %edx, %eax                # eax  = cnt = rlenx
   .odd_top_rows_loop:
    vmulpd    (%rbx), %ymm0, %ymm3      # ymm3 = x_re = x[c+0]*aa0InvSqrt
    vmulpd  32(%rbx), %ymm0, %ymm4      # ymm4 = x_im = x[c+4]*aa0InvSqrt
    vmovupd   %ymm3,  (%rbx)            # x[c+0] = x_re
    vmovupd   %ymm4,32(%rbx)            # x[c+4] = x_im
    add     $64,  %rbx                  # c += 8
   sub      $32,  %rax                  # cnt -= 4
   jg      .odd_top_rows_loop

   # store first values
   vmovsd     %xmm7,  8(%rcx,%rbp)      # triang[xi+1+0] = aaSqrt
   vmovsd     %xmm0, 40(%rcx,%rbp)      # triang[xi+1+4] = aaInvSqrt - store inverse of diagonal at imag()

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
      vmovupd 32(%rsi),      %ymm5      # ymm5 = x_im = x[c*8+4]

      vfnmadd231pd %ymm4 ,%ymm0, %ymm6  # ymm6 = y0_re -= x_re*f0_re
      vfmadd231pd  %ymm4 ,%ymm1, %ymm7  # ymm7 = y0_im += x_re*f0_im
      vfnmadd231pd %ymm4 ,%ymm2, %ymm8  # ymm8 = y1_re -= x_re*f1_re
      vfmadd231pd  %ymm4 ,%ymm3, %ymm9  # ymm9 = y1_im += x_re*f1_im

      vfnmadd231pd %ymm5 ,%ymm1, %ymm6  # ymm6 = y0_re -= x_im*f0_im
      vfnmadd231pd %ymm5 ,%ymm0, %ymm7  # ymm7 = y0_im -= x_im*f0_re
      vfnmadd231pd %ymm5 ,%ymm3, %ymm8  # ymm8 = y1_re -= x_im*f1_im
      vfnmadd231pd %ymm5 ,%ymm2, %ymm9  # ymm9 = y1_im -= x_im*f1_re

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

  .odd_done:
   # store first values
   vmovsd     %xmm7,  8(%rcx,%rbp)      # triang[xi+1+0] = aaSqrt
   vmovsd     %xmm0, 40(%rcx,%rbp)      # triang[xi+1+4] = aaInvSqrt - store inverse of diagonal at imag()
  jmp .done

.even_loop:
# RAX  - xlenx = (x1-x0)*sizeof(double)
# RCX  - triang
# RDX  - rlenx = rlen*sizeof(double) >= 4*sizeof(double)
# R9   = xix = index of diagonal element of x0 in triang[] * sizeof(double)
# YMM0 - aa0InvSqrt x4
# YMM1 - aa1InvSqrt x4
# YMM2 - f_re       x4
# YMM3 - f_im       x4
# YMM10- x aa1Sqrt,     0,0
# YMM11- aa0Sqrt x,     0,0
  # Handle a pair of top rows

  lea (%rcx,%r9,4),%rbx                 # rbx  = x0 = &triang[(rlen & 2)*4]
  mov  %edx,       %edi
  shr  $5,         %edi                 # edi  = cqlen = rlen/4==rlenx/32
  .top_rows_loop:
    ## ax0 = x0[k] * aa0InvSqrt;
    ## ax1 = x1[k];
    ## x0[k] = ax0;
    ## x1[k] = (ax1 - ax0*conj(f))*aa1InvSqrt;
    vmulpd    (%rbx),    %ymm0, %ymm4   # ymm4 = x0_re = x0[c*8+0]*aa0InvSqrt
    vmulpd  32(%rbx),    %ymm0, %ymm5   # ymm5 = x0_im = x0[c*8+4]*aa0InvSqrt
    vmovupd   (%rbx,%rax),      %ymm6   # ymm6= x1_re = x1[c*8+0]
    vmovupd 32(%rbx,%rax),      %ymm7   # ymm7= x1_im = x1[c*8+4]

    vfnmadd231pd %ymm4,%ymm2,%ymm6      # ymm6 = x1_re -= x0_re*f_re
    vfmadd231pd  %ymm4,%ymm3,%ymm7      # ymm7 = x1_im += x0_re*f_im
    vfnmadd231pd %ymm5,%ymm3,%ymm6      # ymm6 = x1_re -= x0_im*f_im
    vfnmadd231pd %ymm5,%ymm2,%ymm7      # ymm7 = x1_im -= x0_im*f_re

    vmulpd    %ymm1, %ymm6, %ymm6       # ymm6 = x1_re *= aa1InvSqrt
    vmulpd    %ymm1, %ymm7, %ymm7       # ymm7 = x1_im *= aa1InvSqrt

    vmovupd   %ymm4,  (%rbx)            # x0[c*8+0] = x0_re
    vmovupd   %ymm5, 32(%rbx)           # x0[c*8+4] = x0_im

    vmovupd   %ymm6,   (%rbx,%rax)      # x1[c*8+0] = x1_re
    vmovupd   %ymm7, 32(%rbx,%rax)      # x1[c*8+4] = x1_im

    add $64, %rbx                       # c += 1
  dec %edi
  jnz .top_rows_loop

  # store last diag and near diag
  lea (%rcx,%r9),   %rsi                # rsi  = x0 = &triang[xi]
  vblendpd $1, %xmm11,%xmm2, %xmm2      # ymm2 = aa0Sqrt,   f_re,0,0
  vblendpd $1, %xmm0, %xmm3, %xmm3      # ymm3 = aa0InvSqrt,f_im,0,0
  vmovupd  %xmm2,   (%rsi)              # x0[0:1] = [aa0Sqrt,    f_re]
  vmovupd  %xmm3, 32(%rsi)              # x0[4:5] = [aa0InvSqrt, f_im]
  vmovupd  %xmm10,  (%rsi,%rax)         # x1[0:1] = [x        aa1Sqrt]
  vmovsd   %xmm1, 40(%rsi,%rax)         # x1[5]   = aa1InvSqrt

  lea (%rcx,%r9, 4), %r10               # r10  = x0 = &triang[(rlen & 2)*4]
  lea (%rcx,%rax,2), %rcx               # rcx  = triang += xlenx*2

  # subtract outer product of two top rows from lower part of the matrix
  # process two output rows together
# RAX  - xlenx = (x1-x0)*sizeof(double)
# RCX  - triang
# RDX  - rlenx = rlen*sizeof(double) >= 4*sizeof(double)
# R10  - x0
# R11  - ltMask
# RBX, RSI, RDI, RBP, R8, R9 - free
  lea (,%edx,4),%edi                    # RDI  = yix == (y0-x0)*sizeof(double) = rlenx*4
  sub $16, %edx                         # RDX  = rlenx -= 2*sizeof(double)
  jz .caxpy2x2_outer_loop_done
  mov %edx,%esi                         # RSI  = clenx = clen*sizeof(double) = rlenx
  mov %edx,%r9d
  and $16, %r9d                         # R9   = xix = (clen & 2)*sizeof(double)
  .caxpy2x2_outer_loop:
    lea (%esi,%r9d),   %r8d             # R8   = ylenx = ((clen+2) & -4)*sizeof(double) = clenx+xix = ((y1-y0)/2)*sizeof(double)
    lea (%edi,%r8d, 2),%ebp             # RBP  = y1-x0 = (y0-x0)+(y1-y0) = (y0-x0)+ylenx*2
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

      vfnmadd231pd %ymm12,%ymm0, %ymm8  # ymm8 = y0_re -= x0_re*f00_re
      vfmadd231pd  %ymm12,%ymm1, %ymm9  # ymm9 = y0_im += x0_re*f00_im
      vfnmadd231pd %ymm12,%ymm2, %ymm10 # ymm10= y1_re -= x0_re*f01_re
      vfmadd231pd  %ymm12,%ymm3, %ymm11 # ymm11= y1_im += x0_re*f01_im

      vmovupd (%rbx,%rax),   %ymm12     # ymm12= x1_re = x1[c*8+0]
      vfnmadd231pd %ymm12,%ymm4, %ymm8  # ymm8 = y0_re -= x1_re*f10_re
      vfmadd231pd  %ymm12,%ymm5, %ymm9  # ymm9 = y0_im += x1_re*f10_im
      vfnmadd231pd %ymm12,%ymm6, %ymm10 # ymm10= y1_re -= x1_re*f11_re
      vfmadd231pd  %ymm12,%ymm7, %ymm11 # ymm11= y1_im += x1_re*f11_im

      vmovupd 32(%rbx),      %ymm12     # ymm12= x0_im = x0[c*8+4]
      vfnmadd231pd %ymm12,%ymm1, %ymm8  # ymm8 = y0_re -= x0_im*f00_im
      vfnmadd231pd %ymm12,%ymm0, %ymm9  # ymm9 = y0_im -= x0_im*f00_re
      vfnmadd231pd %ymm12,%ymm3, %ymm10 # ymm10= y1_re -= x0_im*f01_im
      vfnmadd231pd %ymm12,%ymm2, %ymm11 # ymm11= y1_im -= x0_im*f01_re

      vmovupd 32(%rbx,%rax), %ymm12     # ymm12= x1_im = x1[c*8+4]
      vfnmadd231pd %ymm12,%ymm5, %ymm8  # ymm8 = y0_re -= x1_im*f10_im
      vfnmadd231pd %ymm12,%ymm4, %ymm9  # ymm9 = y0_im -= x1_im*f10_re
      vfnmadd231pd %ymm12,%ymm7, %ymm10 # ymm10= y1_re -= x1_im*f11_im
      vfnmadd231pd %ymm12,%ymm6, %ymm11 # ymm11= y1_im -= x1_im*f11_re

      vmovupd %ymm8,   (%rbx,%rdi)      # y0[c*8+0] = y0_re
      vmovupd %ymm9, 32(%rbx,%rdi)      # y0[c*8+4] = y0_im
      vmovupd %ymm10,  (%rbx,%rbp)      # y1[c*8+0] = y1_re
      vmovupd %ymm11,32(%rbx,%rbp)      # y1[c*8+4] = y1_im
      add $64, %rbx                     # c += 1
    sub $32,   %r8d                     # ylenx -= 4*sizeof(double)
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
# R11 - ltMask
# RAX, RBX, RDI, RSI, RBP, R8, R9, R10 - free
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
  vfnmadd231sd %xmm2, %xmm2, %xmm1      # ymm1 = aa1 -= f_re*f_re
  vfnmadd231sd %xmm3, %xmm3, %xmm1      # ymm1 = aa1 -= f_im*f_im
  vunpcklpd    %xmm1, %xmm0, %xmm0      # ymm0 = aa0, aa1, 0, 0

  # check that we are positive definite
  vmovddup .LDIAG_MIN(%rip), %xmm1
  vcmpltpd  %xmm1, %xmm0, %xmm1         # xmm1 = ltmsk =[aa0<DIAG_MIN aa1<DIAG_MIN]
  vmovddup  .LDIAG_SUBST(%rip), %xmm4
  vblendvpd %xmm1, %xmm4, %xmm0,%xmm0   # xmm0 = aa[0:1] < DIAG_MIN ? DIAG_SUBST : aa[0:1]
  vmovmskpd %xmm1, %esi
  or        %esi,  %r11d
  vsqrtpd   %xmm0, %xmm11               # ymm11= sqrt(aa0), sqrt(aa1), 0, 0
  vmovddup  .LONE(%rip), %xmm1
  vdivpd    %xmm11,%xmm1,%xmm1          # xmm1 = [1/sqrt(aa0), 1/sqrt(aa1)]

  vmulsd    %xmm1, %xmm2,%xmm2          # xmm2  = f_re*=aa0InvSqrt 0 0 0
  vmulsd    %xmm1, %xmm3,%xmm3          # xmm3  = f_im*=aa0InvSqrt 0 0 0

  vbroadcastsd %xmm1, %ymm0             # ymm0 = 1/sqrt(aa0), 1/sqrt(aa0) = aa0InvSqrt x4
  vunpckhpd %xmm1, %xmm1, %xmm1         # ymm1 = 1/sqrt(aa1), 1/sqrt(aa1) 0,0
  vmulsd    %xmm11,%xmm1, %xmm1         # ymm1 = aa1InvSqrt = 1/sqrt(aa1)*aa0Sqrt, 0,0,0
  vmulpd    %xmm0, %xmm11,%xmm10        # ymm10= x aa1Sqrt=sqrt(aa1)*aa0InvSqrt, 0,0
  vbroadcastsd %xmm1, %ymm1             # ymm1 = aa1InvSqrt x4
  vbroadcastsd %xmm2, %ymm2             # ymm2 = f_re       x4
  vbroadcastsd %xmm3, %ymm3             # ymm3 = f_im       x4

  cmp       $16, %edx
  jg       .even_loop                   # x1 was not a last row

  # store last diag and near diag
  vblendpd $1, %xmm11,%xmm2, %xmm2      # ymm2 = aa0Sqrt,   f_re,0,0
  vblendpd $1, %xmm0, %xmm3, %xmm3      # ymm3 = aa0InvSqrt,f_im,0,0
  vmovupd  %xmm2,   (%rbx)              # x0[0:1] = [aa0Sqrt,    f_re]
  vmovupd  %xmm3, 32(%rbx)              # x0[4:5] = [aa0InvSqrt, f_im]
  vmovupd  %xmm10,  (%rbx,%rax)         # x1[0:1] = [x        aa1Sqrt]
  vmovsd   %xmm1, 40(%rbx,%rax)         # x1[5]   = aa1InvSqrt

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
  addq  $16*8+8, %rsp
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
