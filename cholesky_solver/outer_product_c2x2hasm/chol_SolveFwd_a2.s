	.text
	.def	 @feat.00;
	.scl	3;
	.type	0;
	.endef
	.globl	@feat.00
.set @feat.00, 0
	.file	"chol_SolveFwd_a2.s"
	.def	 chol_SolveFwd_a;
	.scl	2;
	.type	32;
	.endef
	.globl	chol_SolveFwd_a         # -- Begin function chol_SolveFwd
	.p2align	4, 0x90
# extern "C" void chol_SolveFwd(double *x, unsigned N, const double* triang)
chol_SolveFwd_a:
  .seh_proc chol_SolveFwd_a
	subq	$24, %rsp
	.seh_stackalloc 24
	vmovapd	%xmm6, 0(%rsp)          # 16-byte Spill
	.seh_savexmm %xmm6, 0
	.seh_endprologue
# RCX - x
# RDX - N
# R8  - triang

  ## x = R' \ conj(B);
  ## R' is lower triangle - solve by forward propagation (caxpy-like)
  ## x = B;
  ## for r=1:N
  ##   x(r) = x(r)/R(r,r);
  ##   x(r+1:N) -= R(r, r+1:N).'*x(r);
  ## end

  shl $3,  %edx                         # edx = rlenx = rlen*sizeof(double) = N*sizeof(double)
  jz  .done
  test $8, %edx
  jz .even_loop_entry
    # special handling for the first row of matrix with odd number of elements
    add        $8,                 %edx  # edx = rlenx += 1*sizeof(double) = (N+1)*sizeof(double)
    mov        %edx,               %eax
    and        $16,                %eax  # RAX  = xix = (rlen & 2) * sizeof(double)
    vmovsd   40(%r8, %rax), %xmm2        # ymm2 = aaInvSqrt,0,0,0 - imag() of diag element contains inverse of real() part
    vmulsd    8(%rcx,%rax), %xmm2, %xmm0 # ymm0 = x0_re = x[xix+1+0] * aaInvSqrt
    vmulsd   40(%rcx,%rax), %xmm2, %xmm1 # ymm1 = x0_im = x[xix+1+4] * aaInvSqrt


    sub $16,   %edx                      # edx  = rlenx -= 2*sizeof(double) = (N-1)*sizeof(double)
    jz .odd_done                         # N==1, nothing more to do

    mov %rcx,  %r9                       # r9   = x
    sub %rcx,  %r8                       # r8   = (triang-x)*sizeof(double)
    mov %edx,  %r10d                     # r10  = rlenx

    vbroadcastsd %xmm0,            %ymm0 # ymm0 =  x0_re   x4
    vbroadcastsd %xmm1,            %ymm1 # ymm1 =  x0_im   x4

    lea        (%rcx, %rax,4),     %rcx  # rcx = &x[c] = &x[xi*4]
    .odd_loop:
      vmovupd   (%rcx,%r8),        %ymm4 # ymm4 = y_re = y[c+0]
      vmovupd   (%rcx),            %ymm2 # ymm2 = x_re = x[c+0]
      vmovupd 32(%rcx),            %ymm3 # ymm3 = x_im = x[c+4]
      vmovupd 32(%rcx,%r8),        %ymm5 # ymm5 = y_im = y[c+4]

      vfnmadd231pd %ymm4, %ymm0,   %ymm2 # ymm2 = x_re -= y_re*x0_re
      vfnmadd231pd %ymm4, %ymm1,   %ymm3 # ymm2 = x_im -= y_re*x0_im

      vfmadd231pd  %ymm5, %ymm1,   %ymm2 # ymm2 = x_re += y_im*x0_im
      vfnmadd231pd %ymm5, %ymm0,   %ymm3 # ymm2 = x_im -= y_im*x0_re

      vmovupd    %ymm2,   (%rcx)         # x[c+0] = x_re
      vmovupd    %ymm3, 32(%rcx)         # x[c+4] = x_im
      add        $64, %rcx               # rcx  = &x[c] += 8
    sub          $32, %edx               # edx  = cnt -= 4*sizeof(double)
    jg  .odd_loop

    add %rcx, %r8                        # r8   = triang += ((N+3)/4)*4*2
    mov %r9,  %rcx                       # rcx  = x
    mov %r10d,%edx                       # rdx  = rlenx
    # store x[0]
    vmovsd  %xmm0,  8(%rcx,%rax)         # x[xix+1+0] = x0_re
    vmovsd  %xmm1, 40(%rcx,%rax)         # x[xix+1+4] = x0_im
    lea        (%rcx, %rax,4),     %rcx  # rcx = x += xi*4
  jmp .even_loop_entry

  .odd_done: # finish case of N==1
  vmovsd  %xmm0,  8(%rcx,%rax)           # x[xix+1+0] = x0_re
  vmovsd  %xmm1, 40(%rcx,%rax)           # x[xix+1+4] = x0_im
  jmp .done

  .even_loop:
    # process two rows per iteration
    lea          (%r10,%rax,4),    %rcx  # RCX  = &x[xi*4]
    vperm2f128 $0,%ymm0,%ymm0,     %ymm0 # ymm0 = x0_re                    x4
    vperm2f128 $0,%ymm1,%ymm1,     %ymm1 # ymm1 = x0_im                    x4
    vperm2f128 $0,%ymm2,%ymm2,     %ymm2 # ymm2 = x1_re                    x4
    vperm2f128 $0,%ymm3,%ymm3,     %ymm3 # ymm3 = x1_im                    x4
    mov           %edx,            %r11d # R11  = rlen
    .inner_loop:
      vmovupd   (%rcx,%r8),        %ymm6 # ymm6 = y0_re = y0[c+0]
      vmovupd   (%rcx),            %ymm4 # ymm4 = x_re  = x[c+0]
      vmovupd 32(%rcx),            %ymm5 # ymm5 = x_im  = x[c+4]

      vfnmadd231pd %ymm6, %ymm0,   %ymm4 # ymm4 = x_re -= y0_re*x0_re
      vfnmadd231pd %ymm6, %ymm1,   %ymm5 # ymm5 = x_im -= y0_re*x0_im

      vmovupd 32(%rcx,%r8),        %ymm6 # ymm5 = y0_im = y[c+4]
      vfmadd231pd  %ymm6, %ymm1,   %ymm4 # ymm4 = x_re += y0_im*x0_im
      vfnmadd231pd %ymm6, %ymm0,   %ymm5 # ymm5 = x_im -= y0_im*x0_re

      vmovupd   (%rcx,%r9),        %ymm6 # ymm6 = y1_re = y1[c+0]
      vfnmadd231pd %ymm6, %ymm2,   %ymm4 # ymm4 = x_re -= y1_re*x1_re
      vfnmadd231pd %ymm6, %ymm3,   %ymm5 # ymm5 = x_im -= y1_re*x1_im

      vmovupd 32(%rcx,%r9),        %ymm6 # ymm5 = y1_im = y1[c+4]
      vfmadd231pd  %ymm6, %ymm3,   %ymm4 # ymm4 = x_re += y1_im*x1_im
      vfnmadd231pd %ymm6, %ymm2,   %ymm5 # ymm5 = x_im -= y1_im*x1_re

      vmovupd    %ymm4,   (%rcx)         # x[c+0] = x_re
      vmovupd    %ymm5, 32(%rcx)         # x[c+0] = x_im
      add        $64, %rcx               # rcx  = &x[c] += 8
    sub          $32, %edx               # edx  = cnt -= 4*sizeof(double)
    jg  .inner_loop

    lea        (%rcx,%r9), %r8           # r8   = triang = &y1[((rlen+2)/4)*4]
    mov         %r11d,     %edx          # rdx  = rlenx
    # store x[0:1]
    vblendpd $1, %xmm0, %xmm2,     %xmm0 # ymm0 = x0_re,x1_re,0,0
    vblendpd $1, %xmm1, %xmm3,     %xmm1 # ymm1 = x0_im,x1_im,0,0
    vmovupd      %xmm0,   (%r10,%rax)    # pr[0:1] = [x0_re x1_re]
    vmovupd      %xmm1, 32(%r10,%rax)    # pr[4:5] = [x0_im x1_im]
    lea         (%r10, %rax,4),    %rcx  # rcx = x += xi*4

    .even_loop_entry:
    mov        %edx,               %eax
    mov        %rcx,               %r10  # R10  = x
    and        $16,                %eax  # RAX  = xix = (rlen & 2) * sizeof(double)
    lea       (%edx,%eax),         %r9d  # R9   = ylenx = rlenx+xix = (rlen+2)/2*2*sizeof(double) = (y1-y0)*sizeof(double)/2
    sub        %rcx,               %r8   # R8   = (y0-x) * sizeof(double)
    add        %rax,               %rcx  # RCX  = &x[xi]
    lea       (%r8, %r9,2),        %r9   # R9   = (y1-x) * sizeof(double) = ((y0-x)+(y1-y0))*sizeof(double)

    # process 2 top results
    vmovupd       (%rcx),          %xmm5 # ymm5 = x0_re,x1_re,0,0 = x[xi+0:1]
    vmovupd     32(%rcx),          %xmm6 # ymm6 = x0_im,x1_im,0,0 = x[xi+4:5]
    vmovupd     32(%rcx,%r8),      %xmm1 # ymm1 = aa0InvSqrt,f_im,0,0 = y0[xi+4:5]
    vmovddup     8(%rcx,%r8),      %xmm0 # ymm0 = f_re       x2 0 0 = y1[xi+5]
    vunpcklpd      %xmm6, %xmm5,   %xmm4 # ymm4 = x0_re,x0_im,0,0
    vunpckhpd      %xmm6, %xmm5,   %xmm5 # ymm5 = x1_re,x1_im,0,0
    vmovddup       %xmm1,          %xmm2 # ymm2 = aa0InvSqrt x2 0 0
    vmulpd         %xmm2, %xmm4,   %xmm4 # ymm4 = x0_re,x0_im,0,0 *= aa0InvSqrt
    vunpckhpd      %xmm1, %xmm1,   %xmm1 # ymm1 = f_im       x2, 0,0
    vmovddup    40(%rcx,%r9),      %xmm2 # ymm2 = aa1InvSqrt x2 0 0 = y0[xi+5]
    # r1 -= r0*f
    vmulpd         %xmm4, %xmm0,   %xmm6 # ymm6 = x0_re*f_re,x0_im*f_re,0,0
    vpermilpd $1,  %xmm4,          %xmm3 # ymm3 = x0_im,x0_re,0,0
    vmulpd         %xmm3, %xmm1,   %xmm3 # ymm3 = x0_im*f_im,x0_re*f_im,0,0
    vaddsubpd      %xmm3, %xmm6,   %xmm6 # ymm6 = x0*f = x0_re*f_re-x0_im*f_im,x0_im*f_re+x0_re*f_im,0,0
    vsubpd         %xmm6, %xmm5,   %xmm5 # ymm5 = x1_re,x1_im,0,0 -= r0*f
    # r1 *= aa1InvSqrt
    vmulpd         %xmm2, %xmm5,   %xmm5 # ymm5 = r1_re,r1_im,0,0 *= aa1InvSqrt

    vmovddup       %xmm4,          %xmm0 # ymm0 = x0_re       x2, 0,0
    vunpckhpd      %xmm4, %xmm4,   %xmm1 # ymm1 = x0_im       x2, 0,0
    vmovddup       %xmm5,          %xmm2 # ymm2 = x1_re       x2, 0,0
    vunpckhpd      %xmm5, %xmm5,   %xmm3 # ymm3 = x1_im       x2, 0,0
  sub  $16, %edx
  jnz  .even_loop

  # store last pair of results
  vblendpd $1, %xmm0, %xmm2,       %xmm0 # ymm0 = x0_re,x1_re,0,0
  vblendpd $1, %xmm1, %xmm3,       %xmm1 # ymm1 = x0_im,x1_im,0,0
  vmovupd      %xmm0,   (%rcx)           # x[xi+0:1] = [x0_re x1_re]
  vmovupd      %xmm1, 32(%rcx)           # x[xi+4:5] = [x0_im x1_im]

  .done:
  vmovaps	 0(%rsp), %xmm6
  add     $24, %rsp
  vzeroupper
  retq
.seh_handlerdata
.text
.seh_endproc
.addrsig

