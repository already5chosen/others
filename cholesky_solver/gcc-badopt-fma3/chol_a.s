	.file	"chol.cpp"
	.text
	.p2align 4
	.def	_ZN12_GLOBAL__N_125chol_FactorizeAndSolveFwdEPSt7complexIdEjS2_b;	.scl	3;	.type	32;	.endef
	.seh_proc	_ZN12_GLOBAL__N_125chol_FactorizeAndSolveFwdEPSt7complexIdEjS2_b
_ZN12_GLOBAL__N_125chol_FactorizeAndSolveFwdEPSt7complexIdEjS2_b:
.LFB7707:
	pushq	%r15
	.seh_pushreg	%r15
	pushq	%r14
	.seh_pushreg	%r14
	pushq	%r13
	.seh_pushreg	%r13
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbp
	.seh_pushreg	%rbp
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$184, %rsp
	.seh_stackalloc	184
	vmovaps	%xmm6, 16(%rsp)
	.seh_savexmm	%xmm6, 16
	vmovaps	%xmm7, 32(%rsp)
	.seh_savexmm	%xmm7, 32
	vmovaps	%xmm8, 48(%rsp)
	.seh_savexmm	%xmm8, 48
	vmovaps	%xmm9, 64(%rsp)
	.seh_savexmm	%xmm9, 64
	vmovaps	%xmm10, 80(%rsp)
	.seh_savexmm	%xmm10, 80
	vmovaps	%xmm11, 96(%rsp)
	.seh_savexmm	%xmm11, 96
	vmovaps	%xmm12, 112(%rsp)
	.seh_savexmm	%xmm12, 112
	vmovaps	%xmm13, 128(%rsp)
	.seh_savexmm	%xmm13, 128
	vmovaps	%xmm14, 144(%rsp)
	.seh_savexmm	%xmm14, 144
	vmovaps	%xmm15, 160(%rsp)
	.seh_savexmm	%xmm15, 160
	.seh_endprologue
	movl	%edx, %edi
	movb	%r9b, 14(%rsp)
	vmovsd	16(%rcx), %xmm4
	movq	%rcx, %r11
	movq	%r8, %rbx
	movl	%r9d, %eax
	shrl	%edi
	testb	$1, %dl
	je	.L27
	vmovsd	.LC2(%rip), %xmm0
	vcomisd	%xmm4, %xmm0
	ja	.L28
	vsqrtsd	%xmm4, %xmm4, %xmm4
	vmovsd	.LC3(%rip), %xmm0
	movb	$1, 15(%rsp)
	vdivsd	%xmm4, %xmm0, %xmm0
.L3:
	vunpcklpd	%xmm0, %xmm4, %xmm4
	vmovupd	%xmm4, 16(%r11)
	cmpl	$1, %edx
	je	.L1
	leaq	32(%r11), %r8
	testb	%al, %al
	jne	.L5
	testl	$-2, %edx
	jle	.L8
	movl	%edi, %ecx
	salq	$5, %rcx
	vbroadcastsd	%xmm0, %ymm0
	movq	%r11, %rax
	addq	%r11, %rcx
.L9:
	vmulpd	32(%rax), %ymm0, %ymm1
	addq	$32, %rax
	vmovupd	%ymm1, (%rax)
	cmpq	%rax, %rcx
	jne	.L9
.L8:
	movq	%rbx, %rsi
	leal	-1(%rdi), %ebp
.L7:
	leal	1(%rdx), %r10d
	salq	$4, %r10
	movl	%ebp, %eax
	addq	%r11, %r10
	salq	$5, %rax
	movq	%r10, %rbx
	movl	%edi, %r9d
	leaq	64(%r11,%rax), %r11
	vxorpd	%xmm7, %xmm7, %xmm7
.L12:
	vbroadcastsd	8(%r8), %ymm4
	vbroadcastsd	24(%r8), %ymm3
	vaddsubpd	%ymm4, %ymm7, %ymm4
	vaddsubpd	%ymm3, %ymm7, %ymm3
	leal	(%r9,%r9), %ecx
	vbroadcastsd	(%r8), %ymm6
	vbroadcastsd	16(%r8), %ymm5
	salq	$4, %rcx
	movq	%r8, %rdx
	movq	%rbx, %rax
	.p2align 4,,10
	.p2align 3
.L11:
	vmovupd	(%rdx), %ymm0
	addq	$32, %rdx
	vmovapd	%ymm0, %ymm1
	vmovapd	%ymm0, %ymm2
	vfnmadd213pd	(%rax), %ymm6, %ymm1
	vfnmadd213pd	(%rax,%rcx), %ymm5, %ymm2
	vpermilpd	$5, %ymm0, %ymm0
	addq	$32, %rax
	vfmadd231pd	%ymm4, %ymm0, %ymm1
	vfmadd132pd	%ymm3, %ymm2, %ymm0
	vmovupd	%ymm1, -32(%rax)
	vmovupd	%ymm0, -32(%rax,%rcx)
	cmpq	%r11, %rdx
	jne	.L11
	leal	0(,%r9,4), %eax
	salq	$4, %rax
	addq	$32, %r8
	addq	%rax, %rbx
	decl	%r9d
	jne	.L12
	vmovsd	16(%r10), %xmm4
	movq	%rsi, %rbx
	movq	%r10, %r11
	jmp	.L2
.L27:
	movb	$1, 15(%rsp)
	leal	-1(%rdi), %ebp
.L2:
	vmovapd	.LC4(%rip), %xmm15
	vmovapd	.LC5(%rip), %xmm14
	vmovapd	.LC6(%rip), %xmm13
	leal	(%rbp,%rbp), %esi
	leal	-2(%rdi), %r12d
	xorl	%r13d, %r13d
	vxorpd	%xmm12, %xmm12, %xmm12
	vmovsd	%xmm4, %xmm4, %xmm5
	.p2align 4,,10
	.p2align 3
.L23:
	vmovsd	24(%r11), %xmm2
	leal	2(%rsi), %eax
	vmulsd	%xmm2, %xmm2, %xmm0
	vmovsd	(%r11), %xmm1
	salq	$4, %rax
	leaq	(%r11,%rax), %rdx
	vfmadd231sd	%xmm5, %xmm5, %xmm0
	vfmsub231sd	16(%rdx), %xmm1, %xmm0
	vunpcklpd	%xmm0, %xmm1, %xmm0
	vcmpltpd	%xmm15, %xmm0, %xmm1
	vblendvpd	%xmm1, %xmm14, %xmm0, %xmm0
	vsqrtpd	%xmm0, %xmm0
	vdivpd	%xmm0, %xmm13, %xmm6
	vmovsd	%xmm0, (%r11)
	vmovmskpd	%xmm1, %ecx
	orl	%ecx, %r13d
	cmpb	$0, 14(%rsp)
	vunpckhpd	%xmm6, %xmm6, %xmm4
	vmulsd	%xmm0, %xmm4, %xmm4
	vunpckhpd	%xmm0, %xmm0, %xmm0
	vmulsd	%xmm6, %xmm0, %xmm0
	vmulsd	%xmm5, %xmm6, %xmm1
	vmulsd	%xmm2, %xmm6, %xmm2
	vmovsd	%xmm6, 8(%r11)
	vmovsd	%xmm6, %xmm6, %xmm3
	vunpcklpd	%xmm4, %xmm0, %xmm0
	vmovupd	%xmm0, 16(%rdx)
	vmovsd	%xmm1, 16(%r11)
	vmovsd	%xmm2, 24(%r11)
	je	.L13
	vmulsd	(%rbx), %xmm6, %xmm8
	vmulsd	8(%rbx), %xmm6, %xmm6
	vmovsd	%xmm2, %xmm2, %xmm7
	vmovsd	24(%rbx), %xmm0
	vmulsd	%xmm8, %xmm2, %xmm5
	vfmadd213sd	16(%rbx), %xmm6, %xmm7
	vmovsd	%xmm8, (%rbx)
	vmovsd	%xmm6, 8(%rbx)
	vfmadd231sd	%xmm6, %xmm1, %xmm5
	vfnmadd231sd	%xmm8, %xmm1, %xmm7
	vsubsd	%xmm5, %xmm0, %xmm5
	vmulsd	%xmm7, %xmm4, %xmm7
	vmulsd	%xmm4, %xmm5, %xmm5
	vmovsd	%xmm7, 16(%rbx)
	vmovsd	%xmm5, 24(%rbx)
	cmpl	$1, %edi
	jbe	.L14
	vbroadcastsd	%xmm2, %ymm2
	vbroadcastsd	%xmm6, %ymm6
	vbroadcastsd	%xmm5, %ymm5
	movl	%r12d, %ecx
	vaddsubpd	%ymm2, %ymm12, %ymm2
	vaddsubpd	%ymm6, %ymm12, %ymm6
	vaddsubpd	%ymm5, %ymm12, %ymm5
	leaq	2(%rcx), %r10
	leaq	32(%r11), %r8
	leaq	32(%rdx), %r9
	vbroadcastsd	%xmm1, %ymm1
	vbroadcastsd	%xmm8, %ymm8
	vbroadcastsd	%xmm7, %ymm7
	vbroadcastsd	%xmm3, %ymm3
	vbroadcastsd	%xmm4, %ymm4
	leaq	32(%rbx), %r14
	salq	$5, %r10
	movl	$32, %eax
	.p2align 4,,10
	.p2align 3
.L20:
	vmulpd	(%r11,%rax), %ymm3, %ymm9
	vmovapd	%ymm9, %ymm0
	vfnmadd213pd	(%rdx,%rax), %ymm1, %ymm0
	vmovapd	%ymm9, %ymm10
	vfnmadd213pd	(%rbx,%rax), %ymm8, %ymm10
	vmovupd	%ymm9, (%r11,%rax)
	vpermilpd	$5, %ymm9, %ymm9
	vfmadd231pd	%ymm2, %ymm9, %ymm0
	vfnmadd132pd	%ymm6, %ymm10, %ymm9
	vmulpd	%ymm4, %ymm0, %ymm0
	vfnmadd231pd	%ymm7, %ymm0, %ymm9
	vmovupd	%ymm0, (%rdx,%rax)
	vpermilpd	$5, %ymm0, %ymm0
	vfnmadd132pd	%ymm5, %ymm9, %ymm0
	vmovupd	%ymm0, (%rbx,%rax)
	addq	$32, %rax
	cmpq	%rax, %r10
	jne	.L20
	movq	%r14, %rbx
.L19:
	leal	0(,%rdi,4), %eax
	salq	$4, %rax
	addq	%rax, %r11
	rdtsc
	salq	$32, %rdx
	orq	%rdx, %rax
	subq	%rax, gl_ticks(%rip)
	movl	%esi, %r15d
	leaq	1(%rcx), %r10
	salq	$4, %r15
	salq	$5, %r10
	addq	%r8, %r15
	movl	%ebp, %r14d
	movq	%r11, %rdx
	.p2align 4,,10
	.p2align 3
.L22:
	vbroadcastsd	8(%r8), %ymm5
	vbroadcastsd	8(%r9), %ymm4
	vbroadcastsd	24(%r8), %ymm3
	vbroadcastsd	24(%r9), %ymm2
	movq	%r15, %rcx
	vaddsubpd	%ymm5, %ymm12, %ymm5
	vaddsubpd	%ymm4, %ymm12, %ymm4
	vaddsubpd	%ymm3, %ymm12, %ymm3
	vaddsubpd	%ymm2, %ymm12, %ymm2
	subq	%r8, %rcx
	vbroadcastsd	(%r8), %ymm9
	vbroadcastsd	(%r9), %ymm8
	vbroadcastsd	16(%r8), %ymm7
	vbroadcastsd	16(%r9), %ymm6
	xorl	%eax, %eax
	addq	%rdx, %rcx
	.p2align 4,,10
	.p2align 3
.L21:
	vmovupd	(%r8,%rax), %ymm1
	vmovupd	(%r9,%rax), %ymm0
	vmovapd	%ymm1, %ymm10
	vmovapd	%ymm1, %ymm11
	vfnmadd213pd	(%rdx,%rax), %ymm9, %ymm10
	vfnmadd213pd	(%rcx,%rax), %ymm7, %ymm11
	vpermilpd	$5, %ymm1, %ymm1
	vfnmadd231pd	%ymm8, %ymm0, %ymm10
	vfnmadd231pd	%ymm6, %ymm0, %ymm11
	vpermilpd	$5, %ymm0, %ymm0
	vfmadd231pd	%ymm5, %ymm1, %ymm10
	vfmadd132pd	%ymm3, %ymm11, %ymm1
	vfmadd231pd	%ymm4, %ymm0, %ymm10
	vfmadd132pd	%ymm2, %ymm1, %ymm0
	vmovupd	%ymm10, (%rdx,%rax)
	vmovupd	%ymm0, (%rcx,%rax)
	addq	$32, %rax
	cmpq	%rax, %r10
	jne	.L21
	leal	0(,%r14,4), %eax
	salq	$4, %rax
	addq	$32, %r8
	addq	$32, %r9
	addq	%rax, %rdx
	subq	$32, %r10
	decl	%r14d
	jne	.L22
	rdtsc
	salq	$32, %rdx
	orq	%rdx, %rax
	addq	%rax, gl_ticks(%rip)
	vmovsd	16(%r11), %xmm5
	subl	$2, %esi
	decl	%ebp
	decl	%edi
	decl	%r12d
	jmp	.L23
.L13:
	cmpl	$1, %edi
	jbe	.L14
	vmulsd	%xmm2, %xmm3, %xmm2
	leaq	32(%rdx), %r9
	leaq	32(%r11), %r8
	leaq	32(%r11,%rax), %rdx
	leaq	64(%r11,%rax), %rax
	vmulsd	%xmm1, %xmm3, %xmm1
	vxorpd	.LC7(%rip), %xmm2, %xmm0
	cmpq	%rax, %r8
	jnb	.L24
	leaq	64(%r11), %rax
	cmpq	%rax, %rdx
	jb	.L47
.L24:
	movl	%esi, %ecx
	vunpcklpd	%xmm1, %xmm0, %xmm5
	shrl	%ecx
	vunpcklpd	%xmm0, %xmm1, %xmm1
	vinsertf128	$1, %xmm5, %ymm5, %ymm5
	vinsertf128	$1, %xmm1, %ymm1, %ymm1
	vbroadcastsd	%xmm4, %ymm4
	vbroadcastsd	%xmm3, %ymm3
	salq	$5, %rcx
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L18:
	vmovupd	32(%r11,%rax), %ymm0
	vmovupd	(%rdx,%rax), %ymm2
	vmulpd	%ymm0, %ymm3, %ymm6
	vpermilpd	$0, %ymm0, %ymm7
	vpermilpd	$15, %ymm0, %ymm0
	vmulpd	%ymm5, %ymm0, %ymm0
	vmovupd	%ymm6, 32(%r11,%rax)
	vmovapd	%ymm1, %ymm6
	vfnmadd132pd	%ymm7, %ymm0, %ymm6
	vfmadd231pd	%ymm7, %ymm1, %ymm0
	vshufpd	$10, %ymm0, %ymm6, %ymm6
	vaddpd	%ymm2, %ymm6, %ymm0
	vsubpd	%ymm6, %ymm2, %ymm2
	vshufpd	$10, %ymm2, %ymm0, %ymm0
	vmulpd	%ymm4, %ymm0, %ymm0
	vmovupd	%ymm0, (%rdx,%rax)
	addq	$32, %rax
	cmpq	%rax, %rcx
	jne	.L18
	movl	%r12d, %ecx
	jmp	.L19
.L14:
	movzbl	15(%rsp), %esi
	testl	%r13d, %r13d
	movl	$0, %eax
	cmovne	%eax, %esi
	movb	%sil, 15(%rsp)
	vzeroupper
.L1:
	movzbl	15(%rsp), %eax
	vmovaps	16(%rsp), %xmm6
	vmovaps	32(%rsp), %xmm7
	vmovaps	48(%rsp), %xmm8
	vmovaps	64(%rsp), %xmm9
	vmovaps	80(%rsp), %xmm10
	vmovaps	96(%rsp), %xmm11
	vmovaps	112(%rsp), %xmm12
	vmovaps	128(%rsp), %xmm13
	vmovaps	144(%rsp), %xmm14
	vmovaps	160(%rsp), %xmm15
	addq	$184, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L47:
	leal	-1(%rsi), %ecx
	salq	$4, %rcx
	vunpcklpd	%xmm0, %xmm1, %xmm3
	movq	%r8, %rax
	movq	%r9, %rdx
	leaq	48(%r11,%rcx), %rcx
	vmovddup	%xmm4, %xmm4
	vpermilpd	$0, %xmm6, %xmm6
	vunpcklpd	%xmm1, %xmm0, %xmm0
	.p2align 4,,10
	.p2align 3
.L17:
	vmovupd	(%rax), %xmm2
	vmovupd	(%rdx), %xmm1
	vpermilpd	$0, %xmm2, %xmm5
	vmulpd	%xmm2, %xmm6, %xmm7
	vpermilpd	$3, %xmm2, %xmm2
	vmulpd	%xmm0, %xmm2, %xmm2
	addq	$16, %rax
	addq	$16, %rdx
	vmovupd	%xmm7, -16(%rax)
	vmovapd	%xmm3, %xmm7
	vfnmadd132pd	%xmm5, %xmm2, %xmm7
	vfmadd231pd	%xmm5, %xmm3, %xmm2
	vmovsd	%xmm7, %xmm2, %xmm2
	vaddpd	%xmm2, %xmm1, %xmm5
	vsubpd	%xmm2, %xmm1, %xmm1
	vmovsd	%xmm5, %xmm1, %xmm1
	vmulpd	%xmm4, %xmm1, %xmm1
	vmovupd	%xmm1, -16(%rdx)
	cmpq	%rax, %rcx
	jne	.L17
	movl	%r12d, %ecx
	jmp	.L19
.L28:
	movb	$0, 15(%rsp)
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	.LC1(%rip), %xmm4
	jmp	.L3
.L5:
	vmulsd	8(%rbx), %xmm0, %xmm1
	vmulsd	(%rbx), %xmm0, %xmm3
	leal	-1(%rdi), %ecx
	vxorpd	%xmm2, %xmm2, %xmm2
	movq	%rcx, %rbp
	vmovsd	%xmm1, 8(%rbx)
	incq	%rcx
	vbroadcastsd	%xmm1, %ymm1
	vaddsubpd	%ymm1, %ymm2, %ymm2
	vmovsd	%xmm3, (%rbx)
	leaq	16(%rbx), %rsi
	vbroadcastsd	%xmm3, %ymm3
	vbroadcastsd	%xmm0, %ymm1
	salq	$5, %rcx
	xorl	%eax, %eax
.L10:
	vmulpd	32(%r11,%rax), %ymm1, %ymm0
	vmovapd	%ymm0, %ymm4
	vfnmadd213pd	16(%rbx,%rax), %ymm3, %ymm4
	vmovupd	%ymm0, 32(%r11,%rax)
	vpermilpd	$5, %ymm0, %ymm0
	vfnmadd132pd	%ymm2, %ymm4, %ymm0
	vmovupd	%ymm0, 16(%rbx,%rax)
	addq	$32, %rax
	cmpq	%rax, %rcx
	jne	.L10
	jmp	.L7
	.seh_endproc
	.p2align 4
	.def	_ZN12_GLOBAL__N_113chol_SolveBwdEPSt7complexIdEjPKS1_;	.scl	3;	.type	32;	.endef
	.seh_proc	_ZN12_GLOBAL__N_113chol_SolveBwdEPSt7complexIdEjPKS1_
_ZN12_GLOBAL__N_113chol_SolveBwdEPSt7complexIdEjPKS1_:
.LFB7711:
	pushq	%r13
	.seh_pushreg	%r13
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbp
	.seh_pushreg	%rbp
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$24, %rsp
	.seh_stackalloc	24
	vmovaps	%xmm6, (%rsp)
	.seh_savexmm	%xmm6, 0
	.seh_endprologue
	movl	%edx, %esi
	leal	1(%rsi), %eax
	imull	%eax, %eax
	movq	%rsi, %r10
	movl	%r10d, %edi
	shrl	%eax
	salq	$4, %rax
	movl	%r10d, %r9d
	salq	$4, %rsi
	addq	%rax, %r8
	andl	$1, %edi
	shrl	%r9d
	movq	%rcx, %rbx
	leaq	(%rcx,%rsi), %r11
	je	.L49
	movq	%r11, %rdx
	movl	$1, %ecx
	.p2align 4,,10
	.p2align 3
.L53:
	leal	0(,%rcx,4), %eax
	salq	$4, %rax
	leal	(%rcx,%rcx), %ebp
	movq	%rdx, %r12
	subq	%rax, %r8
	salq	$4, %rbp
	subq	$32, %rdx
	addq	%r8, %rbp
	vmovsd	(%rdx), %xmm4
	vmovsd	8(%rdx), %xmm2
	vmovsd	16(%rdx), %xmm3
	vmovsd	24(%rdx), %xmm1
	cmpl	$1, %ecx
	je	.L50
	movl	%ecx, %r13d
	vmovq	%xmm4, %xmm4
	vmovq	%xmm3, %xmm3
	vmovq	%xmm2, %xmm2
	vmovq	%xmm1, %xmm1
	salq	$5, %r13
	movl	$32, %eax
	.p2align 4,,10
	.p2align 3
.L51:
	vmovupd	-32(%r12,%rax), %ymm0
	vmovupd	(%r8,%rax), %ymm6
	vmovupd	0(%rbp,%rax), %ymm5
	vfnmadd231pd	%ymm0, %ymm6, %ymm4
	vfnmadd231pd	%ymm0, %ymm5, %ymm3
	addq	$32, %rax
	vpermilpd	$5, %ymm0, %ymm0
	vfnmadd231pd	%ymm0, %ymm6, %ymm2
	vfnmadd231pd	%ymm0, %ymm5, %ymm1
	cmpq	%r13, %rax
	jne	.L51
	vhaddpd	%ymm3, %ymm4, %ymm3
	vhsubpd	%ymm1, %ymm2, %ymm1
	vunpcklpd	%ymm1, %ymm3, %ymm0
	vunpckhpd	%ymm1, %ymm3, %ymm1
	vperm2f128	$33, %ymm1, %ymm0, %ymm2
	vblendpd	$12, %ymm1, %ymm0, %ymm1
	vaddpd	%ymm2, %ymm1, %ymm1
	vmovsd	%xmm1, %xmm1, %xmm4
	vunpckhpd	%xmm1, %xmm1, %xmm2
	vextractf128	$0x1, %ymm1, %xmm1
	vmovsd	%xmm1, %xmm1, %xmm3
	vunpckhpd	%xmm1, %xmm1, %xmm1
.L50:
	vmovsd	24(%rbp), %xmm0
	vmovsd	24(%r8), %xmm5
	vmulsd	%xmm0, %xmm1, %xmm1
	vmulsd	%xmm0, %xmm3, %xmm3
	vmovsd	16(%r8), %xmm0
	leal	1(%rcx), %eax
	vmulsd	%xmm1, %xmm0, %xmm6
	vfnmadd231sd	%xmm1, %xmm5, %xmm4
	vunpcklpd	%xmm1, %xmm3, %xmm1
	vfnmadd231sd	%xmm3, %xmm5, %xmm6
	vfnmadd132sd	%xmm3, %xmm4, %xmm0
	vmovsd	8(%r8), %xmm4
	vmovupd	%xmm1, 16(%rdx)
	vsubsd	%xmm6, %xmm2, %xmm2
	vmulsd	%xmm0, %xmm4, %xmm0
	vmulsd	%xmm4, %xmm2, %xmm2
	vmovsd	%xmm0, (%rdx)
	vmovsd	%xmm2, 8(%rdx)
	cmpl	%ecx, %r9d
	je	.L52
	movl	%eax, %ecx
	jmp	.L53
	.p2align 4,,10
	.p2align 3
.L52:
	movl	%r9d, %eax
	negq	%rax
	salq	$5, %rax
	addq	%rax, %r11
	testl	%edi, %edi
	je	.L68
	andl	$-2, %r10d
	salq	$4, %r10
	leal	(%r9,%r9), %ecx
	movl	$1, %edx
	subq	%r10, %r8
	testl	%ecx, %ecx
	cmovle	%edx, %ecx
	addq	%rsi, %rax
	movl	%ecx, %edx
	shrl	%edx
	addq	%rax, %rbx
	salq	$5, %rdx
	xorl	%eax, %eax
	vxorpd	%xmm4, %xmm4, %xmm4
	.p2align 4,,10
	.p2align 3
.L56:
	vmovupd	(%r8,%rax), %ymm0
	vmovupd	(%rbx,%rax), %ymm2
	vpermilpd	$0, %ymm0, %ymm3
	vpermilpd	$15, %ymm0, %ymm0
	vmulpd	%ymm2, %ymm0, %ymm0
	vpermilpd	$5, %ymm2, %ymm1
	vmovapd	%ymm1, %ymm2
	addq	$32, %rax
	vfmsub132pd	%ymm3, %ymm0, %ymm2
	vfmadd231pd	%ymm3, %ymm1, %ymm0
	vshufpd	$10, %ymm0, %ymm2, %ymm0
	vaddpd	%ymm0, %ymm4, %ymm4
	cmpq	%rax, %rdx
	jne	.L56
	vextractf128	$0x1, %ymm4, %xmm0
	vaddpd	%xmm4, %xmm0, %xmm0
	movl	%ecx, %eax
	andl	$-2, %eax
	andl	$1, %ecx
	vmovsd	%xmm0, %xmm0, %xmm1
	vunpckhpd	%xmm0, %xmm0, %xmm0
	je	.L66
	salq	$4, %rax
	leaq	(%r8,%rax), %rdx
	addq	%r11, %rax
	vmovsd	(%rdx), %xmm5
	vmovsd	(%rax), %xmm2
	vmovsd	8(%rax), %xmm4
	vfmadd231sd	%xmm5, %xmm2, %xmm0
	vfmadd231sd	%xmm5, %xmm4, %xmm1
	vmovsd	8(%rdx), %xmm3
	vfmadd231sd	%xmm4, %xmm3, %xmm0
	vfnmadd231sd	%xmm2, %xmm3, %xmm1
	vzeroupper
.L59:
	vmovupd	-16(%r11), %xmm4
	vunpcklpd	%xmm1, %xmm0, %xmm0
	vsubpd	%xmm0, %xmm4, %xmm0
	vmovddup	-8(%r8), %xmm1
	vmulpd	%xmm1, %xmm0, %xmm0
	vmovupd	%xmm0, -16(%r11)
	jmp	.L67
	.p2align 4,,10
	.p2align 3
.L68:
	vzeroupper
.L67:
	vmovaps	(%rsp), %xmm6
	addq	$24, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	ret
.L66:
	vzeroupper
	jmp	.L59
.L49:
	testl	%edi, %edi
	je	.L67
	vxorpd	%xmm1, %xmm1, %xmm1
	vmovsd	%xmm1, %xmm1, %xmm0
	jmp	.L59
	.seh_endproc
	.p2align 4
	.globl	_Z22chol_getWorkBufferSizei
	.def	_Z22chol_getWorkBufferSizei;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z22chol_getWorkBufferSizei
_Z22chol_getWorkBufferSizei:
.LFB7712:
	.seh_endprologue
	leal	1(%rcx), %eax
	addl	$3, %ecx
	imull	%eax, %ecx
	movl	%ecx, %eax
	shrl	$31, %eax
	addl	%ecx, %eax
	sarl	%eax
	sall	$4, %eax
	ret
	.seh_endproc
	.p2align 4
	.globl	_Z4cholPSt7complexIdEPKS0_iPvi
	.def	_Z4cholPSt7complexIdEPKS0_iPvi;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z4cholPSt7complexIdEPKS0_iPvi
_Z4cholPSt7complexIdEPKS0_iPvi:
.LFB7713:
	pushq	%r15
	.seh_pushreg	%r15
	pushq	%r14
	.seh_pushreg	%r14
	pushq	%r13
	.seh_pushreg	%r13
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbp
	.seh_pushreg	%rbp
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$72, %rsp
	.seh_stackalloc	72
	vmovaps	%xmm6, 48(%rsp)
	.seh_savexmm	%xmm6, 48
	.seh_endprologue
	leal	-1(%r8), %r14d
	movl	%r8d, %edi
	movl	%r8d, %r13d
	movl	176(%rsp), %eax
	movq	%rcx, %rsi
	movq	%rdx, %r10
	xorl	%r8d, %r8d
	cmpl	$16382, %r14d
	jbe	.L126
.L70:
	vmovaps	48(%rsp), %xmm6
	movl	%r8d, %eax
	addq	$72, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
	.p2align 4,,10
	.p2align 3
.L126:
	movl	%edi, %ebx
	andl	$1, %ebx
	movslq	%edi, %rdx
	leal	1(%rdi), %r12d
	addq	%rdx, %rbx
	movl	%edi, %r15d
	salq	$4, %rbx
	movl	%edi, %ebp
	movq	%r12, %rcx
	salq	$4, %r12
	movq	%r12, 32(%rsp)
	andl	$1, %r15d
	addq	%r9, %rbx
	salq	$4, %rbp
	cmpl	$80, %eax
	je	.L72
	cmpl	$81, %eax
	je	.L73
	movl	%r15d, %r9d
	salq	$4, %r9
	addq	%rbx, %r9
	cmpl	$67, %eax
	je	.L127
.L96:
	movq	%r9, %rcx
	movq	%r10, %rdx
	movl	%edi, %eax
	.p2align 4,,10
	.p2align 3
.L95:
	vmovsd	(%rdx), %xmm0
	movl	%eax, %r8d
	vmovsd	%xmm0, (%rcx)
	vmovsd	8(%rdx), %xmm0
	andl	$-2, %r8d
	salq	$4, %r8
	decl	%eax
	vmovsd	%xmm0, 8(%rcx)
	addq	$16, %rdx
	addq	%r8, %rcx
	cmpl	%eax, %r14d
	jne	.L95
	addq	%rbp, %r10
	addq	$16, %r9
	subl	$1, %r14d
	jnb	.L96
	leal	(%rdi,%rdi), %eax
	movq	%rbx, %rcx
	testl	%r15d, %r15d
	jne	.L97
	leaq	(%rbx,%rbp), %rcx
	subl	$2, %eax
.L97:
	cmpl	$3, %eax
	jbe	.L80
	vxorpd	%xmm0, %xmm0, %xmm0
	.p2align 4,,10
	.p2align 3
.L98:
	movl	%eax, %edx
	salq	$4, %rdx
	subl	$4, %eax
	vmovupd	%xmm0, (%rcx)
	addq	%rdx, %rcx
	cmpl	$3, %eax
	ja	.L98
.L80:
	xorl	%r8d, %r8d
	xorl	%r9d, %r9d
	movl	%edi, %edx
	movq	%rbx, %rcx
	call	_ZN12_GLOBAL__N_125chol_FactorizeAndSolveFwdEPSt7complexIdEjS2_b
	movl	%eax, %r8d
	cmpl	$1, %edi
	jne	.L128
	vmovsd	16(%rbx), %xmm0
	movq	$0x000000000, 8(%rsi)
	vmovsd	%xmm0, (%rsi)
	jmp	.L70
.L127:
	movl	%edi, %r11d
	shrl	%r11d
	movq	%rbx, %r8
	testl	%r15d, %r15d
	je	.L75
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovupd	%xmm0, (%rbx)
	vmovsd	(%r10), %xmm0
	leal	(%r11,%r11), %r8d
	vmovsd	%xmm0, 16(%rbx)
	vmovsd	8(%r10), %xmm0
	movslq	%r8d, %rdx
	vmovsd	%xmm0, 24(%rbx)
	salq	$4, %rdx
	xorl	%eax, %eax
	testl	%r8d, %r8d
	je	.L77
	.p2align 4,,10
	.p2align 3
.L78:
	vmovsd	16(%r10,%rax), %xmm0
	vmovsd	%xmm0, 32(%rbx,%rax)
	vmovsd	24(%r10,%rax), %xmm0
	vmovsd	%xmm0, 40(%rbx,%rax)
	addq	$16, %rax
	cmpq	%rax, %rdx
	jne	.L78
.L77:
	movq	32(%rsp), %rax
	leaq	(%rbx,%rax), %r8
	addq	%rax, %r10
.L75:
	vmovsd	(%r10), %xmm0
	leal	(%r11,%r11), %edx
	vmovsd	%xmm0, (%r8)
	vmovsd	8(%r10), %xmm0
	salq	$4, %rdx
	vmovsd	%xmm0, 8(%r8)
	vmovsd	16(%r10), %xmm0
	addq	%r8, %rdx
	vmovsd	%xmm0, 16(%r8)
	vmovsd	24(%r10), %xmm0
	leaq	(%r10,%rbp), %r9
	vmovsd	%xmm0, 24(%r8)
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovupd	%xmm0, (%rdx)
	vmovsd	16(%r9), %xmm0
	vmovsd	%xmm0, 16(%rdx)
	vmovsd	24(%r9), %xmm0
	vmovsd	%xmm0, 24(%rdx)
	cmpl	$1, %r11d
	je	.L80
	leal	(%rcx,%rcx), %r15d
	salq	$4, %r15
	leaq	(%r15,%rbp), %rax
	addq	%r10, %rax
	subq	%r9, %rax
	leal	0(,%r11,4), %r14d
	movq	%r10, %rcx
	movq	%rax, %r12
	vxorpd	%xmm1, %xmm1, %xmm1
	.p2align 4,,10
	.p2align 3
.L83:
	decl	%r11d
	leal	(%r11,%r11), %eax
	testl	%eax, %eax
	jle	.L84
	leal	2(%rax), %r10d
	salq	$4, %r10
	movl	$32, %eax
	.p2align 4,,10
	.p2align 3
.L85:
	vmovsd	(%rcx,%rax), %xmm0
	vmovsd	%xmm0, (%r8,%rax)
	vmovsd	8(%rcx,%rax), %xmm0
	vmovsd	%xmm0, 8(%r8,%rax)
	vmovsd	(%r9,%rax), %xmm0
	vmovsd	%xmm0, (%rdx,%rax)
	vmovsd	8(%r9,%rax), %xmm0
	vmovsd	%xmm0, 8(%rdx,%rax)
	addq	$16, %rax
	cmpq	%rax, %r10
	jne	.L85
.L84:
	addq	%r15, %rcx
	movl	%r14d, %eax
	vmovsd	(%rcx), %xmm0
	salq	$4, %rax
	addq	%rax, %r8
	vmovsd	%xmm0, (%r8)
	vmovsd	8(%rcx), %xmm0
	leal	(%r11,%r11), %edx
	vmovsd	%xmm0, 8(%r8)
	vmovsd	16(%rcx), %xmm0
	salq	$4, %rdx
	vmovsd	%xmm0, 16(%r8)
	vmovsd	24(%rcx), %xmm0
	addq	%r8, %rdx
	vmovsd	%xmm0, 24(%r8)
	vmovupd	%xmm1, (%rdx)
	vmovsd	16(%rcx,%rbp), %xmm0
	addq	%r12, %r9
	vmovsd	%xmm0, 16(%rdx)
	vmovsd	24(%rcx,%rbp), %xmm0
	subl	$4, %r14d
	vmovsd	%xmm0, 24(%rdx)
	cmpl	$1, %r11d
	jne	.L83
	jmp	.L80
.L73:
	movl	%r15d, %eax
	salq	$4, %rax
	movl	$16, %r9d
	leaq	-16(%rbx,%rax), %r11
	.p2align 4,,10
	.p2align 3
.L91:
	leaq	(%r11,%r9), %rcx
	movq	%r10, %rdx
	movl	%edi, %eax
	.p2align 4,,10
	.p2align 3
.L89:
	vmovsd	(%rdx), %xmm0
	movl	%eax, %r8d
	vmovsd	%xmm0, (%rcx)
	vmovsd	8(%rdx), %xmm0
	andl	$-2, %r8d
	salq	$4, %r8
	decl	%eax
	vmovsd	%xmm0, 8(%rcx)
	addq	$16, %rdx
	addq	%r8, %rcx
	cmpl	%eax, %r14d
	jne	.L89
	addq	%r9, %r10
	addq	$16, %r9
	subl	$1, %r14d
	jnb	.L91
	leal	(%rdi,%rdi), %eax
	movq	%rbx, %rcx
	testl	%r15d, %r15d
	jne	.L93
	leaq	(%rbx,%rbp), %rcx
	subl	$2, %eax
.L93:
	cmpl	$3, %eax
	jbe	.L80
	vxorpd	%xmm0, %xmm0, %xmm0
	.p2align 4,,10
	.p2align 3
.L94:
	movl	%eax, %edx
	salq	$4, %rdx
	subl	$4, %eax
	vmovupd	%xmm0, (%rcx)
	addq	%rdx, %rcx
	cmpl	$3, %eax
	ja	.L94
	jmp	.L80
.L72:
	movq	%rbx, %rax
	testl	%r15d, %r15d
	je	.L129
.L86:
	shrl	%r14d
	je	.L87
	leal	-1(%r14), %edx
	leal	1(,%r14,4), %r15d
	salq	$2, %rdx
	movq	%r15, %r12
	subq	%rdx, %r15
	movq	%rbx, 40(%rsp)
	salq	$4, %r12
	salq	$4, %r15
	vxorpd	%xmm6, %xmm6, %xmm6
	movq	%r10, %rbx
	.p2align 4,,10
	.p2align 3
.L88:
	movq	%r12, %r8
	movq	%rbx, %rdx
	leaq	16(%rax), %rcx
	vmovupd	%xmm6, (%rax)
	movq	%r12, %r14
	call	memcpy
	addq	%r12, %rax
	addq	%r12, %rbx
	subq	$64, %r12
	cmpq	%r14, %r15
	jne	.L88
	movq	%rbx, %r10
	movq	40(%rsp), %rbx
.L87:
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovupd	%xmm0, (%rax)
	vmovsd	(%r10), %xmm0
	vmovsd	%xmm0, 16(%rax)
	vmovsd	8(%r10), %xmm0
	vmovsd	%xmm0, 24(%rax)
	jmp	.L80
.L128:
	movq	32(%rsp), %r9
	subl	$2, %edi
	.p2align 4,,10
	.p2align 3
.L101:
	movl	%r13d, %ecx
	andl	$1, %ecx
	salq	$4, %rcx
	addq	%rbx, %rcx
	vmovsd	(%rcx), %xmm0
	leaq	16(%rcx), %rax
	movq	$0x000000000, 8(%rsi)
	vmovsd	%xmm0, (%rsi)
	leaq	(%rsi,%rbp), %rdx
	cmpl	$1, %r13d
	je	.L70
	movl	%edi, %ebx
	addq	$2, %rbx
	salq	$4, %rbx
	addq	%rcx, %rbx
	.p2align 4,,10
	.p2align 3
.L100:
	vmovsd	(%rax), %xmm0
	addq	$16, %rax
	vmovsd	%xmm0, (%rdx)
	vmovsd	-8(%rax), %xmm0
	vmovsd	%xmm0, 8(%rdx)
	addq	%rbp, %rdx
	cmpq	%rbx, %rax
	jne	.L100
	addq	%r9, %rsi
	decl	%r13d
	decl	%edi
	jmp	.L101
.L129:
	movq	%r10, %rdx
	movq	%rbp, %r8
	movq	%rbx, %rcx
	movq	%r10, 152(%rsp)
	call	memcpy
	movq	152(%rsp), %r10
	leaq	(%rbx,%rbp), %rax
	addq	%rbp, %r10
	jmp	.L86
	.seh_endproc
	.p2align 4
	.globl	_Z11chol_solverPSt7complexIdEPKS0_S3_iPvi
	.def	_Z11chol_solverPSt7complexIdEPKS0_S3_iPvi;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z11chol_solverPSt7complexIdEPKS0_S3_iPvi
_Z11chol_solverPSt7complexIdEPKS0_S3_iPvi:
.LFB7714:
	pushq	%r15
	.seh_pushreg	%r15
	pushq	%r14
	.seh_pushreg	%r14
	pushq	%r13
	.seh_pushreg	%r13
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbp
	.seh_pushreg	%rbp
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$72, %rsp
	.seh_stackalloc	72
	vmovaps	%xmm6, 48(%rsp)
	.seh_savexmm	%xmm6, 48
	.seh_endprologue
	leal	-1(%r9), %edi
	movq	%rcx, 144(%rsp)
	movl	184(%rsp), %eax
	movq	%rdx, %r14
	movq	%r8, %r13
	movl	%r9d, %r15d
	cmpl	$16382, %edi
	jbe	.L131
.L160:
	xorl	%r13d, %r13d
.L130:
	vmovaps	48(%rsp), %xmm6
	movl	%r13d, %eax
	addq	$72, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
	.p2align 4,,10
	.p2align 3
.L131:
	movl	%r15d, %r12d
	andl	$1, %r12d
	salq	$4, %r12
	movslq	%r15d, %rsi
	addq	176(%rsp), %r12
	salq	$4, %rsi
	andl	$1, %r9d
	leaq	(%r12,%rsi), %rbx
	cmpl	$80, %eax
	je	.L133
	cmpl	$81, %eax
	je	.L134
	cmpl	$67, %eax
	je	.L180
	movl	%r9d, %r11d
	movl	%r15d, %ebp
	salq	$4, %r11
	salq	$4, %rbp
	addq	%rbx, %r11
	.p2align 4,,10
	.p2align 3
.L157:
	movq	%r11, %rcx
	movq	%r14, %rdx
	movl	%r15d, %eax
	.p2align 4,,10
	.p2align 3
.L156:
	vmovsd	(%rdx), %xmm0
	movl	%eax, %r8d
	vmovsd	%xmm0, (%rcx)
	vmovsd	8(%rdx), %xmm0
	andl	$-2, %r8d
	salq	$4, %r8
	decl	%eax
	vmovsd	%xmm0, 8(%rcx)
	addq	$16, %rdx
	addq	%r8, %rcx
	cmpl	%edi, %eax
	jne	.L156
	addq	%rbp, %r14
	addq	$16, %r11
	leal	-1(%rax), %edi
	testl	%eax, %eax
	jne	.L157
	leal	(%r15,%r15), %eax
	movq	%rbx, %rcx
	testl	%r9d, %r9d
	jne	.L158
	leaq	(%rbx,%rbp), %rcx
	subl	$2, %eax
.L158:
	cmpl	$3, %eax
	jbe	.L141
	vxorpd	%xmm0, %xmm0, %xmm0
	.p2align 4,,10
	.p2align 3
.L159:
	movl	%eax, %edx
	salq	$4, %rdx
	subl	$4, %eax
	vmovupd	%xmm0, (%rcx)
	addq	%rdx, %rcx
	cmpl	$3, %eax
	ja	.L159
.L141:
	movq	%r13, %rdx
	movq	%rsi, %r8
	movq	%r12, %rcx
	call	memcpy
	movl	$1, %r9d
	movq	%r12, %r8
	movl	%r15d, %edx
	movq	%rbx, %rcx
	call	_ZN12_GLOBAL__N_125chol_FactorizeAndSolveFwdEPSt7complexIdEjS2_b
	movl	%eax, %r13d
	testb	%al, %al
	je	.L160
	movq	%rbx, %r8
	movl	%r15d, %edx
	movq	%r12, %rcx
	call	_ZN12_GLOBAL__N_113chol_SolveBwdEPSt7complexIdEjPKS1_
	movq	144(%rsp), %rcx
	movq	%rsi, %r8
	movq	%r12, %rdx
	call	memcpy
	jmp	.L130
	.p2align 4,,10
	.p2align 3
.L180:
	movl	%r15d, %r11d
	shrl	%r11d
	movq	%rbx, %r8
	testl	%r9d, %r9d
	je	.L136
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovupd	%xmm0, (%rbx)
	vmovsd	(%rdx), %xmm0
	leal	(%r11,%r11), %r8d
	vmovsd	%xmm0, 16(%rbx)
	vmovsd	8(%rdx), %xmm0
	movslq	%r8d, %rcx
	vmovsd	%xmm0, 24(%rbx)
	salq	$4, %rcx
	xorl	%eax, %eax
	testl	%r8d, %r8d
	je	.L138
	.p2align 4,,10
	.p2align 3
.L139:
	vmovsd	16(%r14,%rax), %xmm0
	vmovsd	%xmm0, 32(%rbx,%rax)
	vmovsd	24(%r14,%rax), %xmm0
	vmovsd	%xmm0, 40(%rbx,%rax)
	addq	$16, %rax
	cmpq	%rax, %rcx
	jne	.L139
.L138:
	leal	1(%r15), %eax
	salq	$4, %rax
	leaq	(%rbx,%rax), %r8
	addq	%rax, %r14
.L136:
	vmovsd	(%r14), %xmm0
	leal	(%r11,%r11), %ecx
	vmovsd	%xmm0, (%r8)
	vmovsd	8(%r14), %xmm0
	salq	$4, %rcx
	vmovsd	%xmm0, 8(%r8)
	vmovsd	16(%r14), %xmm0
	movl	%r15d, %ebp
	vmovsd	%xmm0, 16(%r8)
	vmovsd	24(%r14), %xmm0
	addq	%r8, %rcx
	vmovsd	%xmm0, 24(%r8)
	salq	$4, %rbp
	vxorpd	%xmm0, %xmm0, %xmm0
	leaq	(%r14,%rbp), %r9
	vmovupd	%xmm0, (%rcx)
	vmovsd	16(%r9), %xmm0
	vmovsd	%xmm0, 16(%rcx)
	vmovsd	24(%r9), %xmm0
	vmovsd	%xmm0, 24(%rcx)
	cmpl	$1, %r11d
	je	.L141
	leal	2(%r15,%r15), %eax
	movq	%rax, %rdx
	salq	$4, %rdx
	leaq	(%rdx,%rbp), %rax
	addq	%r14, %rax
	subq	%r9, %rax
	movq	%rsi, 40(%rsp)
	movq	%r14, %r10
	leal	0(,%r11,4), %edi
	movq	%rax, %r14
	vxorpd	%xmm1, %xmm1, %xmm1
	movq	%rdx, %rsi
	.p2align 4,,10
	.p2align 3
.L144:
	decl	%r11d
	leal	(%r11,%r11), %edx
	testl	%edx, %edx
	jle	.L145
	addl	$2, %edx
	salq	$4, %rdx
	movl	$32, %eax
	.p2align 4,,10
	.p2align 3
.L146:
	vmovsd	(%r10,%rax), %xmm0
	vmovsd	%xmm0, (%r8,%rax)
	vmovsd	8(%r10,%rax), %xmm0
	vmovsd	%xmm0, 8(%r8,%rax)
	vmovsd	(%r9,%rax), %xmm0
	vmovsd	%xmm0, (%rcx,%rax)
	vmovsd	8(%r9,%rax), %xmm0
	vmovsd	%xmm0, 8(%rcx,%rax)
	addq	$16, %rax
	cmpq	%rax, %rdx
	jne	.L146
.L145:
	addq	%rsi, %r10
	movl	%edi, %eax
	vmovsd	(%r10), %xmm0
	salq	$4, %rax
	addq	%rax, %r8
	vmovsd	%xmm0, (%r8)
	vmovsd	8(%r10), %xmm0
	leal	(%r11,%r11), %ecx
	vmovsd	%xmm0, 8(%r8)
	vmovsd	16(%r10), %xmm0
	salq	$4, %rcx
	vmovsd	%xmm0, 16(%r8)
	vmovsd	24(%r10), %xmm0
	addq	%r8, %rcx
	vmovsd	%xmm0, 24(%r8)
	vmovupd	%xmm1, (%rcx)
	vmovsd	16(%r10,%rbp), %xmm0
	addq	%r14, %r9
	vmovsd	%xmm0, 16(%rcx)
	vmovsd	24(%r10,%rbp), %xmm0
	subl	$4, %edi
	vmovsd	%xmm0, 24(%rcx)
	cmpl	$1, %r11d
	jne	.L144
	movq	40(%rsp), %rsi
	jmp	.L141
	.p2align 4,,10
	.p2align 3
.L134:
	movl	%r9d, %eax
	salq	$4, %rax
	movl	$16, %r11d
	leaq	-16(%rbx,%rax), %rbp
	.p2align 4,,10
	.p2align 3
.L152:
	leaq	0(%rbp,%r11), %rcx
	movq	%r14, %rdx
	movl	%r15d, %eax
	.p2align 4,,10
	.p2align 3
.L150:
	vmovsd	(%rdx), %xmm0
	movl	%eax, %r8d
	vmovsd	%xmm0, (%rcx)
	vmovsd	8(%rdx), %xmm0
	andl	$-2, %r8d
	salq	$4, %r8
	decl	%eax
	vmovsd	%xmm0, 8(%rcx)
	addq	$16, %rdx
	addq	%r8, %rcx
	cmpl	%eax, %edi
	jne	.L150
	addq	%r11, %r14
	addq	$16, %r11
	subl	$1, %edi
	jnb	.L152
	leal	(%r15,%r15), %eax
	movq	%rbx, %rcx
	testl	%r9d, %r9d
	jne	.L154
	movl	%r15d, %ecx
	salq	$4, %rcx
	addq	%rbx, %rcx
	subl	$2, %eax
.L154:
	cmpl	$3, %eax
	jbe	.L141
	vxorpd	%xmm0, %xmm0, %xmm0
	.p2align 4,,10
	.p2align 3
.L155:
	movl	%eax, %edx
	salq	$4, %rdx
	subl	$4, %eax
	vmovupd	%xmm0, (%rcx)
	addq	%rdx, %rcx
	cmpl	$3, %eax
	ja	.L155
	jmp	.L141
	.p2align 4,,10
	.p2align 3
.L133:
	movq	%rbx, %rax
	testl	%r9d, %r9d
	je	.L181
.L147:
	shrl	%edi
	je	.L148
	leal	-1(%rdi), %ecx
	leal	1(,%rdi,4), %edx
	salq	$2, %rcx
	movq	%rdx, %rbp
	subq	%rcx, %rdx
	movq	%rdx, %rdi
	salq	$4, %rdi
	movq	%rsi, 40(%rsp)
	salq	$4, %rbp
	movq	%r14, %rsi
	vxorpd	%xmm6, %xmm6, %xmm6
	movq	%r13, %r14
	movq	%rbx, %r13
	movq	%rdi, %rbx
	.p2align 4,,10
	.p2align 3
.L149:
	movq	%rbp, %r8
	movq	%rsi, %rdx
	leaq	16(%rax), %rcx
	vmovupd	%xmm6, (%rax)
	movq	%rbp, %rdi
	call	memcpy
	addq	%rbp, %rax
	addq	%rbp, %rsi
	subq	$64, %rbp
	cmpq	%rdi, %rbx
	jne	.L149
	movq	%r13, %rbx
	movq	%r14, %r13
	movq	%rsi, %r14
	movq	40(%rsp), %rsi
.L148:
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovupd	%xmm0, (%rax)
	vmovsd	(%r14), %xmm0
	vmovsd	%xmm0, 16(%rax)
	vmovsd	8(%r14), %xmm0
	vmovsd	%xmm0, 24(%rax)
	jmp	.L141
.L181:
	movl	%r15d, %ebp
	salq	$4, %rbp
	movq	%rbp, %r8
	movq	%rbx, %rcx
	call	memcpy
	addq	%rbp, %r14
	leaq	(%rbx,%rbp), %rax
	jmp	.L147
	.seh_endproc
	.p2align 4
	.globl	_Z16chol_trif_solverPSt7complexIdEPKS0_iPv
	.def	_Z16chol_trif_solverPSt7complexIdEPKS0_iPv;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z16chol_trif_solverPSt7complexIdEPKS0_iPv
_Z16chol_trif_solverPSt7complexIdEPKS0_iPv:
.LFB7715:
	pushq	%r14
	.seh_pushreg	%r14
	pushq	%r13
	.seh_pushreg	%r13
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbp
	.seh_pushreg	%rbp
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$64, %rsp
	.seh_stackalloc	64
	vmovaps	%xmm6, 32(%rsp)
	.seh_savexmm	%xmm6, 32
	vmovaps	%xmm7, 48(%rsp)
	.seh_savexmm	%xmm7, 48
	.seh_endprologue
	leal	-1(%r8), %eax
	movq	%rcx, %r12
	movl	%r8d, %r13d
	movq	%r9, %rsi
	cmpl	$16382, %eax
	jbe	.L204
	vmovaps	32(%rsp), %xmm6
	vmovaps	48(%rsp), %xmm7
	addq	$64, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	ret
	.p2align 4,,10
	.p2align 3
.L204:
	movl	%r8d, %ebp
	andl	$1, %ebp
	movslq	%r8d, %rbx
	salq	$4, %rbp
	leaq	(%r9,%rbp), %r14
	movl	%r8d, %edi
	salq	$4, %rbx
	movq	%rbx, %r8
	movq	%r14, %rcx
	andl	$1, %edi
	call	memcpy
	leaq	(%r14,%rbx), %r8
	testl	%edi, %edi
	jne	.L184
	movl	%r13d, %r11d
	shrl	%r11d
	movq	%r8, %r9
	movq	%r14, %r10
.L185:
	vmovsd	8(%r9), %xmm5
	vmovsd	16(%r9), %xmm2
	vmulsd	(%r10), %xmm5, %xmm6
	vmulsd	8(%r10), %xmm5, %xmm5
	vmovsd	24(%r9), %xmm1
	vmovsd	24(%r10), %xmm0
	leal	(%r11,%r11), %ecx
	salq	$4, %rcx
	vmulsd	%xmm2, %xmm5, %xmm3
	addq	%r9, %rcx
	vmovsd	24(%rcx), %xmm4
	vmovsd	%xmm6, (%r10)
	vmovsd	%xmm5, 8(%r10)
	vfmadd231sd	%xmm1, %xmm6, %xmm3
	vfmadd213sd	16(%r10), %xmm5, %xmm1
	vsubsd	%xmm3, %xmm0, %xmm3
	vfnmadd231sd	%xmm2, %xmm6, %xmm1
	vmulsd	%xmm4, %xmm3, %xmm3
	vmulsd	%xmm1, %xmm4, %xmm1
	vmovsd	%xmm3, 24(%r10)
	vmovsd	%xmm1, 16(%r10)
	cmpl	$1, %r11d
	je	.L202
	leal	0(,%r11,4), %esi
	vxorpd	%xmm4, %xmm4, %xmm4
	.p2align 4,,10
	.p2align 3
.L193:
	movl	%r11d, %eax
	vbroadcastsd	%xmm5, %ymm5
	vbroadcastsd	%xmm3, %ymm3
	leal	-2(%rax), %edi
	vaddsubpd	%ymm5, %ymm4, %ymm5
	vaddsubpd	%ymm3, %ymm4, %ymm3
	addq	$2, %rdi
	movq	%r10, %rdx
	vbroadcastsd	%xmm6, %ymm6
	addq	$32, %r10
	vbroadcastsd	%xmm1, %ymm7
	decl	%r11d
	salq	$5, %rdi
	movl	$32, %eax
	.p2align 4,,10
	.p2align 3
.L192:
	vmovupd	(%r9,%rax), %ymm1
	vmovupd	(%rcx,%rax), %ymm0
	vmovapd	%ymm1, %ymm2
	vfnmadd213pd	(%rdx,%rax), %ymm6, %ymm2
	vpermilpd	$5, %ymm1, %ymm1
	vfnmadd231pd	%ymm7, %ymm0, %ymm2
	vpermilpd	$5, %ymm0, %ymm0
	vfnmadd132pd	%ymm5, %ymm2, %ymm1
	vfnmadd132pd	%ymm3, %ymm1, %ymm0
	vmovupd	%ymm0, (%rdx,%rax)
	addq	$32, %rax
	cmpq	%rdi, %rax
	jne	.L192
	movl	%esi, %eax
	salq	$4, %rax
	addq	%rax, %r9
	vmovsd	8(%r9), %xmm5
	vmovsd	16(%r9), %xmm2
	vmulsd	(%r10), %xmm5, %xmm6
	vmulsd	8(%r10), %xmm5, %xmm5
	vmovsd	24(%r9), %xmm1
	vmovsd	24(%r10), %xmm0
	leal	(%r11,%r11), %ecx
	salq	$4, %rcx
	vmulsd	%xmm2, %xmm5, %xmm3
	addq	%r9, %rcx
	vmovsd	24(%rcx), %xmm7
	subl	$4, %esi
	vmovsd	%xmm6, (%r10)
	vfmadd231sd	%xmm1, %xmm6, %xmm3
	vfmadd213sd	16(%r10), %xmm5, %xmm1
	vmovsd	%xmm5, 8(%r10)
	vsubsd	%xmm3, %xmm0, %xmm3
	vfnmadd231sd	%xmm2, %xmm6, %xmm1
	vmulsd	%xmm7, %xmm3, %xmm3
	vmulsd	%xmm1, %xmm7, %xmm1
	vmovsd	%xmm3, 24(%r10)
	vmovsd	%xmm1, 16(%r10)
	cmpl	$1, %r11d
	jne	.L193
.L202:
	vzeroupper
.L191:
	movl	%r13d, %edx
	movq	%r14, %rcx
	call	_ZN12_GLOBAL__N_113chol_SolveBwdEPSt7complexIdEjPKS1_
	nop
	vmovaps	32(%rsp), %xmm6
	vmovaps	48(%rsp), %xmm7
	movq	%rbx, %r8
	movq	%r14, %rdx
	movq	%r12, %rcx
	addq	$64, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	jmp	memcpy
.L184:
	vmovsd	24(%r8), %xmm0
	vmulsd	(%r14), %xmm0, %xmm4
	vmulsd	8(%r14), %xmm0, %xmm5
	vmovsd	%xmm4, (%r14)
	vmovsd	%xmm5, 8(%r14)
	cmpl	$1, %r13d
	je	.L191
	leaq	16(%rsi,%rbp), %rdx
	leaq	48(%r8), %r9
	movq	%rdx, %rcx
	subq	%r9, %rcx
	movl	%r13d, %r11d
	movl	%r13d, %eax
	addq	$8, %rcx
	leaq	16(%r14), %r10
	shrl	%r11d
	andl	$-2, %eax
	cmpq	$16, %rcx
	jbe	.L187
	shrl	%eax
	movl	%eax, %ecx
	vbroadcastsd	%xmm5, %ymm5
	vbroadcastsd	%xmm4, %ymm4
	salq	$5, %rcx
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L189:
	vmovupd	32(%r8,%rax), %ymm0
	vmovapd	%ymm4, %ymm2
	vpermilpd	$5, %ymm0, %ymm3
	vmulpd	%ymm5, %ymm3, %ymm3
	vmovupd	(%rdx,%rax), %ymm1
	vfnmadd132pd	%ymm0, %ymm3, %ymm2
	vfmadd132pd	%ymm4, %ymm3, %ymm0
	vshufpd	$10, %ymm0, %ymm2, %ymm0
	vaddpd	%ymm1, %ymm0, %ymm2
	vsubpd	%ymm0, %ymm1, %ymm1
	vshufpd	$10, %ymm1, %ymm2, %ymm1
	vmovupd	%ymm1, (%rdx,%rax)
	addq	$32, %rax
	cmpq	%rax, %rcx
	jne	.L189
.L188:
	leal	1(%r13), %r9d
	salq	$4, %r9
	addq	%r8, %r9
	jmp	.L185
.L187:
	decl	%eax
	salq	$4, %rax
	leaq	32(%r8), %rcx
	addq	%rax, %r9
	vmovddup	%xmm4, %xmm0
	vmovddup	%xmm5, %xmm5
.L190:
	vmovupd	(%rcx), %xmm1
	vmovupd	(%rdx), %xmm2
	vpermilpd	$1, %xmm1, %xmm3
	vmulpd	%xmm5, %xmm3, %xmm3
	vmovapd	%xmm1, %xmm4
	addq	$16, %rcx
	addq	$16, %rdx
	vfnmadd132pd	%xmm0, %xmm3, %xmm4
	vfmadd132pd	%xmm0, %xmm3, %xmm1
	vmovsd	%xmm4, %xmm1, %xmm1
	vaddpd	%xmm1, %xmm2, %xmm3
	vsubpd	%xmm1, %xmm2, %xmm1
	vmovsd	%xmm3, %xmm1, %xmm4
	vmovupd	%xmm4, -16(%rdx)
	cmpq	%r9, %rcx
	jne	.L190
	jmp	.L188
	.seh_endproc
	.globl	gl_ticks
	.bss
	.align 8
gl_ticks:
	.space 8
	.section .rdata,"dr"
	.align 8
.LC0:
	.long	210911779
	.long	1002937505
	.align 8
.LC1:
	.long	2025163840
	.long	1142271773
	.align 8
.LC2:
	.long	0
	.long	940572672
	.align 8
.LC3:
	.long	0
	.long	1072693248
	.align 16
.LC4:
	.long	0
	.long	940572672
	.long	0
	.long	940572672
	.align 16
.LC5:
	.long	-238855003
	.long	1211982633
	.long	-238855003
	.long	1211982633
	.align 16
.LC6:
	.long	0
	.long	1072693248
	.long	0
	.long	1072693248
	.align 16
.LC7:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.ident	"GCC: (Rev1, Built by MSYS2 project) 10.2.0"
	.def	memcpy;	.scl	2;	.type	32;	.endef
