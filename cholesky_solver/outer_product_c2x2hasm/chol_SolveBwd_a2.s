  .text
  .def   @feat.00;
  .scl  3;
  .type 0;
  .endef
  .globl  @feat.00
.set @feat.00, 0
  .file "chol_SolveBwd_a1.s"
  .def   chol_SolveBwd_a;
  .scl  2;
  .type 32;
  .endef

  .text
  .globl  chol_SolveBwd_a
  .p2align  4, 0x90
# extern "C" void chol_SolveBwd(double *x, unsigned N, const double* triang)
chol_SolveBwd_a:
# RCX - x
# RDX - N
# R8  - triang
.seh_proc chol_SolveBwd_a
  pushq %rbx
  .seh_pushreg %rbx
  subq  $16*6, %rsp
  .seh_stackalloc 16*6
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

  cmp  $1,     %edx
  jbe  .short_N
  # N > 1
  lea  2(%edx),         %eax    # eax  = N+2
  imul %eax,            %eax    # eax  = (N+2)*(N+2)
  shr  $3,              %eax    # eax  = trlenq = (N+2)*(N+2)/8 = number of quad-elements in triang array
  shl  $6,              %eax    # eax  = trlenq*4*2*sizeof(double)
  lea  -128(%rax,%r8),  %r8     # r8   = triang = &triang[(N+1)*(N+1)/2 - 16]; // points line in triang that corresponds to rlen=2
  lea -1(%edx),         %eax    # eax  = N-1
  and  $-4,             %eax    # eax  = (N-1) & -4
  shl  $4,              %eax    # eax  = ((N-1) & -4)*2*sizeof(double)
  add  %rax,            %rcx    # rcx  = x += ((N-1) & -4)*2;                  // points to last quad of x
  shl  $3,              %edx    # edx  = Nx = N*sizeof(double)
  mov  %edx,            %r11d   # r11  = Nx = N*sizeof(double)

  mov  $16,             %edx    # edx  = rlenx = rlen*sizeof(double)=2*sizeof(double)
  mov  %edx,            %eax    # eax  = xix = 2*sizeof(double)
  vmovupd 16(%rcx),     %xmm1   # ymm1 = acc_re = acc0_re acc1_re 0 0
  vmovupd 48(%rcx),     %xmm2   # ymm2 = acc_im = acc0_im acc1_im 0 0
  vunpcklpd %xmm2,%xmm1,%xmm0   # ymm0 = acc0 = acc0_re acc0_im 0 0
  vunpckhpd %xmm2,%xmm1,%xmm1   # ymm1 = acc1 = acc1_re acc1_im 0 0
  sub  %rcx,            %r8     # r8   = y0-x
  lea 64(%r8),          %r9     # r9   = y1-x
  jmp       .rloop_entry

  # process two rows per iteration
  .rloop:
    # RAX  = xix
    # RCX  = x
    # RDX  = rlenx=rlen*sizeof(double) > 2*sizeof(double)
    # R8   = y0-x
    # R11  = Nx   = N*sizeof(double)
    # RBX, R9, R10 - free

    lea     (,%eax,4),         %ebx
    sub       %rbx,            %rcx    # rcx  = x -= xi*4
    lea     (,%edx,4),         %ebx
    sub       %rbx,            %r8     # r8   = y0-x -= rlen*4
    lea      (%eax,%edx),      %ebx    # rbx  = ylenx = rlenx+xix = (y1-y0)*sizeof(double)/2
    lea      (%r8, %rbx,2),    %r9     # r9   = y1-x

    vmovsd   (%rcx,%rax),      %xmm0   # ymm0 = acc0_rere = x[xix+0+0] 0 0 0
    vmovsd 32(%rcx,%rax),      %xmm1   # ymm1 = acc0_imre = x[xix+0+4] 0 0 0
    vmovsd  8(%rcx,%rax),      %xmm2   # ymm2 = acc1_rere = x[xix+1+0] 0 0 0
    vmovsd 40(%rcx,%rax),      %xmm3   # ymm3 = acc1_imre = x[xix+1+4] 0 0 0
    vxorpd     %ymm4, %ymm4,   %ymm4   # ymm4 = acc0_reim = zeros
    vxorpd     %ymm5, %ymm5,   %ymm5   # ymm5 = acc0_imim = zeros
    vxorpd     %ymm6, %ymm6,   %ymm6   # ymm6 = acc1_reim = zeros
    vxorpd     %ymm7, %ymm7,   %ymm7   # ymm7 = acc1_reim = zeros

    vmovupd    %xmm4,      (%rcx,%rax) # x[xix+0+0:1] = 0
    vmovupd    %xmm4,    32(%rcx,%rax) # x[xix+4+0:1] = 0

    lea       (%rcx,%rax,4),   %rbx    # rbx  = &x[c] = &x[xi]
    mov       %edx,            %r10d   # r10d = rlenx
    shr       $5,              %edx    # edx  = cnt = rlen/4
    .cloop:
      # acc0 -= x[c] * conj(y0[c]);
      # acc1 -= x[c] * conj(y1[c]);
      vmovupd   (%rbx, %r8),   %ymm10 # ymm10 = y0_re = y0[c+0:3]
      vmovupd   (%rbx),        %ymm8  # ymm8  = x_re  = x[c+0:3]
      vmovupd 32(%rbx),        %ymm9  # ymm9  = x_im  = x[c+4:7]

      vmulpd     %ymm8, %ymm10,%ymm11 # ymm11 = x_re*y0_re
      vsubpd     %ymm11,%ymm0, %ymm0  # ymm0  = acc0_rere -= x_re*y0_re
      vmulpd     %ymm9, %ymm10,%ymm11 # ymm11 = x_im*y0_re
      vsubpd     %ymm11,%ymm1, %ymm1  # ymm1  = acc0_imre -= x_im*y0_re

      vmovupd 32(%rbx, %r8),   %ymm10 # ymm10 = y0_im = y0[c+4:7]
      vmulpd     %ymm8, %ymm10,%ymm11 # ymm11 = x_re*y0_im
      vaddpd     %ymm11,%ymm4, %ymm4  # ymm4  = acc0_reim += x_re*y0_im
      vmulpd     %ymm9, %ymm10,%ymm11 # ymm11 = x_im*y0_im
      vaddpd     %ymm11,%ymm5, %ymm5  # ymm5  = acc0_imim += x_im*y0_im

      vmovupd   (%rbx, %r9),   %ymm10 # ymm10 = y1_re = y1[c+0:3]
      vmulpd     %ymm8, %ymm10,%ymm11 # ymm11 = x_re*y1_re
      vsubpd     %ymm11,%ymm2, %ymm2  # ymm2  = acc1_rere -= x_re*y1_re
      vmulpd     %ymm9, %ymm10,%ymm11 # ymm11 = x_im*y1_re
      vsubpd     %ymm11,%ymm3, %ymm3  # ymm3  = acc1_imre -= x_im*y1_re

      vmovupd 32(%rbx, %r9),   %ymm10 # ymm10 = y1_im = y1[c+4:7]
      vmulpd     %ymm8, %ymm10,%ymm11 # ymm11 = x_re*y1_im
      vaddpd     %ymm11,%ymm6, %ymm6  # ymm6  = acc1_reim += x_re*y1_im
      vmulpd     %ymm9, %ymm10,%ymm11 # ymm11 = x_im*y1_im
      vaddpd     %ymm11,%ymm7, %ymm7  # ymm7  = acc1_imim += x_im*y1_im

      add      $64,   %rbx            # rbx  = &x[c] += 8
    dec               %edx            # edx  = cnt -= 1
    jnz .cloop
    # reduce accumulators and pack re/im
    vsubpd     %ymm5, %ymm0, %ymm0         # ymm0 = acc0_re = acc0_rere - acc0_imim
    vaddpd     %ymm4, %ymm1, %ymm1         # ymm1 = acc0_im = acc0_imre + acc0_reim
    vsubpd     %ymm7, %ymm2, %ymm2         # ymm2 = acc1_re = acc1_rere - acc1_imim
    vaddpd     %ymm6, %ymm3, %ymm3         # ymm3 = acc1_im = acc1_imre + acc1_reim

    vhaddpd    %ymm1, %ymm0, %ymm0         # ymm0 = acc0a_re acc0a_im acc0b_re acc0b_im
    vhaddpd    %ymm3, %ymm2, %ymm1         # ymm1 = acc1a_re acc1a_im acc1b_re acc1b_im
    vperm2f128 $1,    %ymm0, %ymm0, %ymm2  # ymm2 = acc0b_re acc0b_im x x
    vaddpd     %xmm2, %xmm0, %xmm0         # ymm0 = acc0_re  acc0_im  0 0
    vperm2f128 $1,    %ymm1, %ymm1, %ymm2  # ymm2 = acc1b_re acc1b_im x x
    vaddpd     %xmm2, %xmm1, %xmm1         # ymm1 = acc1_re  acc1_im  0 0

    mov        %r10d, %edx                 # rdx  = rlenx

    .rloop_entry:
    # RAX  = xix
    # RCX  = x
    # RDX  = rlenx=rlen*sizeof(double)
    # R8   = y0-x
    # R9   = y1-x
    # R10  = rlenx
    # R11  = Nx   = N*sizeof(double)
    # YMM0 = acc0 = acc0_re acc0_im 0 0
    # YMM1 = acc1 = acc1_re acc1_im 0 0
    # RBX  - free
    lea           (%rcx,%rax),      %rbx  # RBX  = &x[xi]
    vmovddup    40(%rbx,%r9),       %xmm2 # ymm2 = aa1InvSqrt x2 0 0 = y1[xi+1+4]
    vmulpd         %xmm2, %xmm1,    %xmm1 # acc1 *= aa1InvSqrt
    # acc0 -= acc1*conj(y0[1]);
    vmovddup     8(%rbx,%r8),       %xmm2 # ymm2 = y_re       x2 0 0 = y0[xi+1+0]
    vmovddup    40(%rbx,%r8),       %xmm3 # ymm3 = y_im       x2 0 0 = y0[xi+1+4]
    vmovddup    32(%rbx,%r8),       %xmm4 # ymm4 = aa0InvSqrt x2 0 0 = y0[xi+0+4]
    vmulpd         %xmm1, %xmm2,    %xmm2 # ymm2 = acc1_re*y_re  acc1_im*y_re 0 0
    vsubpd         %xmm2, %xmm0,    %xmm0 # ymm0 = acc0_re -= acc1_re*y_re acc0_im -= acc1_im*y_re
    vpermilpd $1,  %xmm1,           %xmm2 # ymm2 = acc1_im  acc1_re  0 0
    vmulpd         %xmm3, %xmm2,    %xmm2 # ymm2 = acc1_im*y_im  acc1_re*y0_im 0 0
    vaddsubpd      %xmm2, %xmm0,    %xmm0 # ymm0 = acc0_re -= acc1_im*y_im acc0_im -= acc1_re*y_im
    vmulpd         %xmm4, %xmm0,    %xmm0 # ymm0 = acc0 *= aa0InvSqrt

    vunpcklpd      %xmm1, %xmm0,    %xmm2 # ymm2 = acc_re = acc0_re acc1_re 0 0
    vunpckhpd      %xmm1, %xmm0,    %xmm0 # ymm0 = acc_im = acc0_im acc1_im 0 0
    vmovupd        %xmm2,          (%rbx) # x[xi+0:1] = [acc0_re acc1_re]
    vmovupd        %xmm0,        32(%rbx) # x[xi+4:5] = [acc0_im acc1_im]

    add            $16,             %edx  # rlenx += 2*sizeof(double)
    xor            $16,             %eax  # xix   ^= 2*sizeof(double)
  cmp              %r11d,           %edx
  jbe .rloop                              # rlenx <= Nx

  je  .done
    # N % 2 == 1
    # handle the first row of matrix with odd number of elements
    # RAX  = xix
    # RCX  = x
    # RDX  = rlenx=rlen*sizeof(double) == (N+1)*sizeof(double)
    # R8   = y0-x
    # RBX, R9, R10, R11 - free

    lea     (,%eax,4),         %ebx
    sub       %rbx,            %rcx    # rcx  = x -= xi*4
    lea     (,%edx,4),         %ebx
    sub       %rbx,            %r8     # r8   = y0-x -= rlen*4
    lea      (%eax,%edx),      %ebx    # rbx  = ylenx = rlenx+xix = (y1-y0)*sizeof(double)/2
    lea      (%r8, %rbx,2),    %r9     # r9   = y1-x

    vmovsd  8(%rcx,%rax),      %xmm2   # ymm2 = acc1_rere = x[xix+1+0] 0 0 0
    vmovsd 40(%rcx,%rax),      %xmm3   # ymm3 = acc1_imre = x[xix+1+4] 0 0 0
    vxorpd     %ymm6, %ymm6,   %ymm6   # ymm6 = acc1_reim = zeros
    vxorpd     %ymm7, %ymm7,   %ymm7   # ymm7 = acc1_reim = zeros

    vmovupd    %xmm6,      (%rcx,%rax) # x[xix+0+0:1] = 0
    vmovupd    %xmm6,    32(%rcx,%rax) # x[xix+4+0:1] = 0

    lea       (%rcx,%rax,4),   %rbx    # rbx  = &x[c] = &x[xi]
    shr       $5,              %edx    # edx  = cnt = rlen/4
    .last_cloop:
      # acc1 -= x[c] * conj(y1[c]);
      vmovupd   (%rbx, %r9),   %ymm0   # ymm0  = y1_re = y1[c+0:3]
      vmovupd   (%rbx),        %ymm4   # ymm4  = x_re  = x[c+0:3]
      vmovupd 32(%rbx),        %ymm5   # ymm5  = x_im  = x[c+4:7]

      vmulpd     %ymm4, %ymm0, %ymm1   # ymm1  = x_re*y1_re
      vsubpd     %ymm1, %ymm2, %ymm2   # ymm2  = acc1_rere -= x_re*y1_re
      vmulpd     %ymm5, %ymm0, %ymm1   # ymm1  = x_im*y1_re
      vsubpd     %ymm1, %ymm3, %ymm3   # ymm3  = acc1_imre -= x_im*y1_re

      vmovupd 32(%rbx, %r9),   %ymm0   # ymm0  = y1_im = y1[c+4:7]
      vmulpd     %ymm4, %ymm0, %ymm1   # ymm1  = x_re*y1_im
      vaddpd     %ymm1, %ymm6, %ymm6   # ymm6  = acc1_reim += x_re*y1_im
      vmulpd     %ymm5, %ymm0, %ymm1   # ymm1  = x_im*y1_im
      vaddpd     %ymm1, %ymm7, %ymm7   # ymm7  = acc1_imim += x_im*y1_im

      add      $64,   %rbx             # rbx  = &x[c] += 8
    dec               %edx             # edx  = cnt -= 1
    jnz .last_cloop
    # reduce accumulators and pack re/im
    vsubpd     %ymm7, %ymm2, %ymm2        # ymm2 = acc1_re = acc1_rere - acc1_imim
    vaddpd     %ymm6, %ymm3, %ymm3        # ymm3 = acc1_im = acc1_imre + acc1_reim

    vhaddpd    %ymm3, %ymm2, %ymm1        # ymm1 = acc1a_re acc1a_im acc1b_re acc1b_im
    vperm2f128 $1,    %ymm1, %ymm1, %ymm2 # ymm2 = acc1b_re acc1b_im x x
    vaddpd     %xmm2, %xmm1, %xmm1        # ymm1 = acc1_re  acc1_im  0 0

    lea           (%rcx,%rax),      %rbx  # RBX  = &x[xi]
    vmovddup    40(%rbx,%r9),       %xmm2 # ymm2 = aa1InvSqrt x2 0 0 = y1[xi+1+4]
    vmulpd         %xmm2, %xmm1,    %xmm0 # ymm0 = acc1_re  acc1_im  0 0 *= aa1InvSqrt

    vunpckhpd      %xmm0, %xmm0,    %xmm1 # ymm0 = acc1_im acc1_im 0 0
    vmovsd         %xmm0,         8(%rbx) # x[xi+1] = acc1_re
    vmovsd         %xmm1,        40(%rbx) # x[xi+5] = acc1_im

  .done:
  vmovaps 0*16(%rsp), %xmm6
  vmovaps 1*16(%rsp), %xmm7
  vmovaps 2*16(%rsp), %xmm8
  vmovaps 3*16(%rsp), %xmm9
  vmovaps 4*16(%rsp), %xmm10
  vmovaps 5*16(%rsp), %xmm11
  addq  $16*6, %rsp
  pop     %rbx
  vzeroupper
  retq

  .short_N:  # N <= 1
  test        %edx,   %edx
  jz   .done
  # N ==1
  vmovsd  7*8(%r8),  %xmm0        # ymm0 = aaInvSqrt
  vmulsd  3*8(%rcx), %xmm0, %xmm1 # ymm1 = x0_re = x[3+0]*aaInvSqrt
  vmulsd  7*8(%rcx), %xmm0, %xmm0 # ymm0 = x0_im = x[3+4]*aaInvSqrt
  vmovsd  %xmm1,        3*8(%rcx) # ymm1 = x[3+0] = x0_re
  vmovsd  %xmm0,        7*8(%rcx) # ymm0 = x[3+4] = x0_im
  jmp  .done


.seh_handlerdata
  .text
  .seh_endproc
                                        # -- End function
  .addrsig
