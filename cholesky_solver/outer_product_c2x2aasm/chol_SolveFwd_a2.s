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
# extern "C" void chol_SolveFwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
chol_SolveFwd_a:
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

  test $1, %edx
  jz .start_even
    # special handling for the first row of matrix with odd number of elements
    vmovddup 24(%r8),%xmm0              # imag() of diag element contains inverse of real() part
    vmulpd   (%rcx), %xmm0, %xmm0       # xmm0 = xr0 = x[0] * y[1].imag()
    vmovupd   %xmm0, (%rcx)             # x[0] = xr0

    shr $1, %edx                        # edx = rhlen = N/2
    jz  .done

    vunpckhpd %xmm0, %xmm0, %xmm1       # ymm3 =  xr1.im xr1.im 0 0
    vmovddup  %xmm0, %xmm0              # ymm2 =  xr1.im xr1.im 0 0
    vxorpd    %ymm4, %ymm4, %ymm4       # ymm4 = 0 0 0 0
    vaddsubpd %xmm1, %xmm4, %xmm1       # ymm1 = -xr0.im xr0.im 0 0
    vperm2f128 $0,%ymm0,%ymm0,%ymm0     # ymm0 =  xr0.re                    x4
    vperm2f128 $0,%ymm1,%ymm1,%ymm1     # ymm1 = -xr0.im xr0.im -xr0.im xr0.im

    add $32, %r8                        # r8  = triang += 2
    add $16, %rcx                       # rcx = x += 1
    mov %edx,%r9d                       # r9d = cnt
    mov %rcx,%rax
    sub %r8, %rax                       # rax = (x-tring)*sizeof(DComplex)
    .odd_loop:
      vmovupd (%r8),     %ymm2          # ymm2 = vy0 = y0[c] (re im re im)
      vmovupd (%r8,%rax),%ymm3          # ymm3 = vx  = x[c]  (re im re im)
      vfnmadd231pd %ymm2,%ymm0,%ymm3    # ymm3 = vx.re-vy0.re*xr0.re vx.im-vy0.im*vr0.re
      vpermilpd $5,%ymm2,%ymm2          # ymm2 = vy0         { im re im re }
      vfnmadd231pd %ymm2,%ymm1,%ymm3    # ymm3 = vx.re+vy0.im*xr0.im vx.re-vy0.im*vr0.im
      vmovupd %ymm3, (%r8,%rax)         # x[c] = vx == x[c] - y0[c]*xr0 - y1[c]*xr1
      add $32, %r8                      # triang += 2
    dec %r9d
    jnz .odd_loop
  jmp .even_loop_entry

.start_even:
  shr $1, %edx                          # edx = rhlen = N/2
  jnz .even_loop_entry
  jmp .done

  .even_loop:
    # process two rows per iteration
    add $32, %r8                        # triang += 2
    add $32, %rcx                       # x      += 2
    vunpckhpd %xmm2, %xmm2, %xmm3       # ymm3 =  xr1.im xr1.im 0 0
    vmovddup  %xmm2, %xmm2              # ymm2 =  xr1.im xr1.im 0 0
    vxorpd    %ymm4, %ymm4, %ymm4       # ymm4 = 0 0 0 0
    vaddsubpd %xmm1, %xmm4, %xmm1       # ymm1 = -xr0.im xr0.im 0 0
    vaddsubpd %xmm3, %xmm4, %xmm3       # ymm3 = -xr1.im xr1.im 0 0
    vperm2f128 $0,%ymm0,%ymm0,%ymm0     # ymm0 =  xr0.re                    x4
    vperm2f128 $0,%ymm1,%ymm1,%ymm1     # ymm1 = -xr0.im xr0.im -xr0.im xr0.im
    vperm2f128 $0,%ymm2,%ymm2,%ymm2     # ymm2 =  xr1.re                    x4
    vperm2f128 $0,%ymm3,%ymm3,%ymm3     # ymm3 = -xr1.im xr1.im -xr1.im xr1.im
    mov %rcx, %r10                      # save x
    sub %r8,  %rcx                      # rcx = (x-tring)*sizeof(DComplex)
    mov %edx, %r9d                      # r9d = cnt
    .inner_loop:
      vmovupd (%r8),       %ymm4        # ymm4 = vy0 = y0[c] (re im re im)
      vmovupd (%r8,%rcx),  %ymm5        # ymm5 = vx  = x[c]  (re im re im)
      vfnmadd231pd %ymm4,%ymm0,%ymm5    # ymm5 = vx.re-vy0.re*xr0.re vx.im-vy0.im*vr0.re
      vpermilpd $5,%ymm4,%ymm4          # ymm4 = vy0         { im re im re }
      vfnmadd231pd %ymm4,%ymm1,%ymm5    # ymm5 = vx.re+vy0.im*xr0.im vx.re-vy0.im*vr0.im

      vmovupd (%r8,%rax),  %ymm4        # ymm5 = vy1 = y1[c] (re im re im)
      vfnmadd231pd %ymm4,%ymm2,%ymm5    # ymm5 = vx.re-vy1.re*xr1.re vx.im-vy1.im*vr1.re
      vpermilpd $5,%ymm4,%ymm4          # ymm4 = vy0         { im re im re }
      vfnmadd231pd %ymm4,%ymm3,%ymm5    # ymm5 = vx.re+vy1.im*xr1.im vx.re-vy1.im*vr1.im

      vmovupd %ymm5, (%r8,%rcx)         # x[c] = vx == x[c] - y0[c]*xr0 - y1[c]*xr1
      add $32, %r8                      # triang += 2
    dec %r9d
    jnz .inner_loop
    mov %r10, %rcx
    add %rax, %r8                       # triang += (y1-y0)

    .even_loop_entry:
    vmovddup 8(%r8), %xmm0              # imag() of diag element contains inverse of real() part
    vmulpd   (%rcx), %xmm0, %xmm0       # xmm0 = xr0 = x[0] * y0[0].imag()
    mov      %edx, %eax
    shl      $5,   %eax                 # eax = (y1-y0)*sizeof(DComplex) = rhlen*2*sizeof(DComplex)
    # xr1 = (x[1] - y0[1]*xr0) * y1[1].imag(); // odd row is padded
    vmovupd  16(%r8),     %xmm2         # xmm2 = y0[1].re y0[1].im
    vmovupd  16(%rcx),    %xmm3         # xmm3 = x[1]
    vmovddup 24(%r8,%rax),%xmm5         # xmm5 = y1[1].im y1[1].im
    vunpckhpd      %xmm0, %xmm0, %xmm1  # xmm1 = xr0.im xr0.im
    vmovupd        %xmm0, (%rcx)        # x[0] = xr0
    vmovddup       %xmm0, %xmm0         # xmm0 = xr0.re xr0.re
    vpermilpd $1,  %xmm2, %xmm4         # xmm4 = y0[1].im y0[1].re
    vmulpd         %xmm4, %xmm1, %xmm4  # xmm2 = y0[1].im*xr0.im y0[1].re*xr0.im
    vfmaddsub231pd %xmm2, %xmm0, %xmm4  # xmm4 = y0[1]*xr0 = y0[1].re*xr0.re-y0[1].im*xr0.im y0[1].im*xr0.re+y0[1].re*xr0.im
    vsubpd         %xmm4, %xmm3, %xmm2  # xmm2 = x[1] - y0[1]*xr0
    vmulpd         %xmm2, %xmm5, %xmm2  # xmm2 = xr1 = (x[1] - y0[1]*xr0) * y1[1].imag()
    vmovupd        %xmm2, 16(%rcx)      # x[1] = xr1
  dec  %edx
  jnz  .even_loop

  .done:
  vzeroupper
  retq
.text
.addrsig

