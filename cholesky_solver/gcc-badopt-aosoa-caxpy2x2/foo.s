	.file	"foo.c"
	.text
	.p2align 4
	.globl	foo
	.def	foo;	.scl	2;	.type	32;	.endef
	.seh_proc	foo
foo:
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
	subq	$1192, %rsp
	.seh_stackalloc	1192
	vmovaps	%xmm6, 1024(%rsp)
	.seh_savexmm	%xmm6, 1024
	vmovaps	%xmm7, 1040(%rsp)
	.seh_savexmm	%xmm7, 1040
	vmovaps	%xmm8, 1056(%rsp)
	.seh_savexmm	%xmm8, 1056
	vmovaps	%xmm9, 1072(%rsp)
	.seh_savexmm	%xmm9, 1072
	vmovaps	%xmm10, 1088(%rsp)
	.seh_savexmm	%xmm10, 1088
	vmovaps	%xmm11, 1104(%rsp)
	.seh_savexmm	%xmm11, 1104
	vmovaps	%xmm12, 1120(%rsp)
	.seh_savexmm	%xmm12, 1120
	vmovaps	%xmm13, 1136(%rsp)
	.seh_savexmm	%xmm13, 1136
	vmovaps	%xmm14, 1152(%rsp)
	.seh_savexmm	%xmm14, 1152
	vmovaps	%xmm15, 1168(%rsp)
	.seh_savexmm	%xmm15, 1168
	.seh_endprologue
	movq	%rdx, %r10
	movl	%r9d, %r11d
	movl	%r9d, %edx
	andl	$2, %edx
	andl	$2, %r11d
	leaq	0(,%rdx,8), %rax
	addl	%r11d, %r9d
	vmovsd	(%r10,%rdx,8), %xmm7
	vmovsd	(%r8,%rdx,8), %xmm6
	vmovsd	8(%r10,%rax), %xmm5
	vmovsd	8(%r8,%rax), %xmm4
	vmovsd	32(%r10,%rax), %xmm1
	vmovsd	32(%r8,%rax), %xmm3
	vmovsd	40(%r10,%rax), %xmm2
	vmovsd	40(%r8,%rax), %xmm0
	addl	%r9d, %r9d
	movslq	%r9d, %r12
	vmovsd	%xmm7, 280(%rsp)
	vmovsd	%xmm6, 992(%rsp)
	vmovsd	%xmm5, 160(%rsp)
	vmovsd	%xmm4, 272(%rsp)
	vmovsd	%xmm1, 1000(%rsp)
	vmovsd	%xmm3, 1008(%rsp)
	vmovsd	%xmm2, 1016(%rsp)
	vmovsd	%xmm0, 168(%rsp)
	leaq	(%rcx,%r12,8), %rdx
	testl	%r9d, %r9d
	jle	.L22
	leaq	256(,%r12,8), %rax
	decl	%r9d
	testq	%rax, %rax
	leaq	256(%rcx), %rax
	setle	%r11b
	cmpq	%rax, %rdx
	setnb	%al
	orb	%al, %r11b
	je	.L3
	cmpl	$7, %r9d
	jbe	.L3
	movl	%r9d, %r13d
	shrl	$3, %r13d
	leal	1(%r13), %ebp
	cmpl	$23, %r9d
	jbe	.L11
	vbroadcastsd	%xmm6, %ymm6
	vbroadcastsd	%xmm7, %ymm7
	vmovupd	%ymm6, 320(%rsp)
	movl	%ebp, %edi
	vbroadcastsd	%xmm3, %ymm6
	vmovupd	%ymm7, 176(%rsp)
	vmovupd	%ymm6, 64(%rsp)
	vbroadcastsd	%xmm1, %ymm7
	vbroadcastsd	%xmm4, %ymm6
	shrl	$2, %edi
	vmovupd	%ymm7, (%rsp)
	vbroadcastsd	%xmm5, %ymm5
	vbroadcastsd	%xmm2, %ymm7
	vmovupd	%ymm6, 208(%rsp)
	salq	$8, %rdi
	vbroadcastsd	%xmm0, %ymm6
	vmovupd	%ymm5, 32(%rsp)
	vmovupd	%ymm7, 128(%rsp)
	vmovupd	%ymm6, 96(%rsp)
	movq	%r10, %r9
	movq	%rcx, %rax
	movq	%rdx, %rsi
	movq	%r8, %rbx
	addq	%r10, %rdi
	movq	%rdx, %r11
	.p2align 4,,10
	.p2align 3
.L5:
	vmovupd	(%r9), %ymm7
	vmovupd	64(%r9), %ymm6
	vunpcklpd	32(%r9), %ymm7, %ymm2
	vunpckhpd	32(%r9), %ymm7, %ymm0
	vmovupd	64(%r9), %ymm7
	vmovupd	192(%r9), %ymm4
	vunpckhpd	96(%r9), %ymm7, %ymm5
	vmovupd	128(%r9), %ymm7
	vunpcklpd	96(%r9), %ymm6, %ymm6
	vunpcklpd	160(%r9), %ymm7, %ymm3
	vunpckhpd	160(%r9), %ymm7, %ymm1
	vmovupd	192(%r9), %ymm7
	vunpcklpd	224(%r9), %ymm4, %ymm4
	vunpckhpd	224(%r9), %ymm7, %ymm8
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm8, %ymm8
	vpermpd	$216, %ymm2, %ymm2
	vpermpd	$216, %ymm0, %ymm0
	vpermpd	$216, %ymm3, %ymm3
	vpermpd	$216, %ymm1, %ymm1
	vunpcklpd	%ymm6, %ymm2, %ymm7
	vunpckhpd	%ymm6, %ymm2, %ymm2
	vunpcklpd	%ymm4, %ymm3, %ymm6
	vunpckhpd	%ymm4, %ymm3, %ymm3
	vunpcklpd	%ymm5, %ymm0, %ymm4
	vunpckhpd	%ymm5, %ymm0, %ymm0
	vunpcklpd	%ymm8, %ymm1, %ymm5
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm3, %ymm3
	vunpcklpd	%ymm5, %ymm4, %ymm11
	vpermpd	$216, %ymm2, %ymm2
	vunpckhpd	%ymm5, %ymm4, %ymm4
	vunpckhpd	%ymm8, %ymm1, %ymm1
	vpermpd	$216, %ymm0, %ymm0
	vpermpd	$216, %ymm4, %ymm8
	vpermpd	$216, %ymm1, %ymm1
	vunpcklpd	%ymm3, %ymm2, %ymm4
	vunpckhpd	%ymm3, %ymm2, %ymm2
	vpermpd	$216, %ymm2, %ymm5
	vunpcklpd	%ymm1, %ymm0, %ymm2
	vpermpd	$216, %ymm4, %ymm10
	vpermpd	$216, %ymm2, %ymm4
	vmovupd	64(%rax), %ymm2
	vmovupd	(%rax), %ymm3
	vmovupd	%ymm4, 448(%rsp)
	vunpckhpd	96(%rax), %ymm2, %ymm4
	vmovupd	128(%rax), %ymm2
	vpermpd	$216, %ymm6, %ymm6
	vunpckhpd	%ymm1, %ymm0, %ymm1
	vpermpd	$216, %ymm7, %ymm7
	vunpcklpd	32(%rax), %ymm3, %ymm9
	vunpckhpd	32(%rax), %ymm3, %ymm14
	vunpckhpd	160(%rax), %ymm2, %ymm0
	vmovupd	64(%rax), %ymm3
	vunpcklpd	%ymm6, %ymm7, %ymm12
	vunpckhpd	%ymm6, %ymm7, %ymm7
	vpermpd	$216, %ymm1, %ymm6
	vunpcklpd	160(%rax), %ymm2, %ymm1
	vmovupd	192(%rax), %ymm2
	vunpcklpd	96(%rax), %ymm3, %ymm3
	vmovupd	%ymm5, 416(%rsp)
	vunpcklpd	224(%rax), %ymm2, %ymm5
	vunpckhpd	224(%rax), %ymm2, %ymm2
	vpermpd	$216, %ymm3, %ymm3
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm9, %ymm9
	vpermpd	$216, %ymm1, %ymm1
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm0, %ymm0
	vmovupd	%ymm10, 384(%rsp)
	vpermpd	$216, %ymm14, %ymm14
	vunpcklpd	%ymm3, %ymm9, %ymm10
	vpermpd	$216, %ymm2, %ymm2
	vunpckhpd	%ymm3, %ymm9, %ymm9
	vunpcklpd	%ymm5, %ymm1, %ymm3
	vpermpd	$216, %ymm3, %ymm3
	vmovupd	%ymm8, 288(%rsp)
	vpermpd	$216, %ymm10, %ymm10
	vunpcklpd	%ymm4, %ymm14, %ymm8
	vunpckhpd	%ymm4, %ymm14, %ymm14
	vunpcklpd	%ymm2, %ymm0, %ymm4
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm8, %ymm8
	vunpckhpd	%ymm2, %ymm0, %ymm2
	vunpcklpd	%ymm3, %ymm10, %ymm0
	vpermpd	$216, %ymm0, %ymm13
	vunpcklpd	%ymm4, %ymm8, %ymm0
	vunpckhpd	%ymm4, %ymm8, %ymm8
	vpermpd	$216, %ymm2, %ymm2
	vunpckhpd	%ymm3, %ymm10, %ymm10
	vpermpd	$216, %ymm14, %ymm14
	vpermpd	$216, %ymm0, %ymm3
	vpermpd	$216, %ymm8, %ymm0
	vmovupd	%ymm6, 480(%rsp)
	vunpckhpd	%ymm5, %ymm1, %ymm1
	vmovupd	%ymm3, 512(%rsp)
	vmovupd	(%rsi), %ymm3
	vmovupd	%ymm0, 544(%rsp)
	vunpcklpd	%ymm2, %ymm14, %ymm0
	vpermpd	$216, %ymm1, %ymm1
	vpermpd	$216, %ymm0, %ymm4
	vpermpd	$216, %ymm9, %ymm9
	vunpcklpd	%ymm1, %ymm9, %ymm6
	vmovupd	%ymm4, 640(%rsp)
	vunpckhpd	%ymm1, %ymm9, %ymm9
	vunpcklpd	32(%rsi), %ymm3, %ymm4
	vunpckhpd	32(%rsi), %ymm3, %ymm1
	vmovupd	64(%rsi), %ymm3
	vunpckhpd	%ymm2, %ymm14, %ymm14
	vunpcklpd	96(%rsi), %ymm3, %ymm8
	vunpckhpd	96(%rsi), %ymm3, %ymm5
	vmovupd	128(%rsi), %ymm3
	vpermpd	$216, %ymm14, %ymm2
	vunpckhpd	160(%rsi), %ymm3, %ymm0
	vmovupd	%ymm2, 672(%rsp)
	vunpcklpd	160(%rsi), %ymm3, %ymm2
	vmovupd	192(%rsi), %ymm3
	vmovupd	192(%rsi), %ymm14
	vunpcklpd	224(%rsi), %ymm3, %ymm3
	vpermpd	$216, %ymm9, %ymm9
	vmovupd	%ymm9, 608(%rsp)
	vunpckhpd	224(%rsi), %ymm14, %ymm9
	vpermpd	$216, %ymm8, %ymm8
	vpermpd	$216, %ymm3, %ymm3
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm2, %ymm2
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm9, %ymm9
	vmovupd	%ymm6, 576(%rsp)
	vpermpd	$216, %ymm1, %ymm1
	vunpcklpd	%ymm8, %ymm4, %ymm6
	vpermpd	$216, %ymm0, %ymm0
	vunpckhpd	%ymm8, %ymm4, %ymm4
	vunpcklpd	%ymm3, %ymm2, %ymm8
	vpermpd	$216, %ymm8, %ymm8
	vpermpd	$216, %ymm6, %ymm6
	vunpckhpd	%ymm3, %ymm2, %ymm2
	vunpcklpd	%ymm5, %ymm1, %ymm3
	vunpckhpd	%ymm5, %ymm1, %ymm1
	vunpcklpd	%ymm9, %ymm0, %ymm5
	vpermpd	$216, %ymm2, %ymm2
	vpermpd	$216, %ymm5, %ymm5
	vunpcklpd	%ymm8, %ymm6, %ymm14
	vpermpd	$216, %ymm4, %ymm4
	vunpckhpd	%ymm8, %ymm6, %ymm6
	vpermpd	$216, %ymm3, %ymm3
	vunpckhpd	%ymm9, %ymm0, %ymm0
	vpermpd	$216, %ymm6, %ymm9
	vunpcklpd	%ymm5, %ymm3, %ymm6
	vunpckhpd	%ymm5, %ymm3, %ymm3
	vunpcklpd	%ymm2, %ymm4, %ymm5
	vunpckhpd	%ymm2, %ymm4, %ymm4
	vpermpd	$216, %ymm0, %ymm0
	vpermpd	$216, %ymm4, %ymm2
	vpermpd	$216, %ymm1, %ymm1
	vmovupd	%ymm2, 832(%rsp)
	vunpcklpd	%ymm0, %ymm1, %ymm2
	vunpckhpd	%ymm0, %ymm1, %ymm1
	vpermpd	$216, %ymm1, %ymm0
	vmovupd	%ymm0, 896(%rsp)
	vmovupd	(%rbx), %ymm0
	vpermpd	$216, %ymm2, %ymm4
	vunpckhpd	32(%rbx), %ymm0, %ymm1
	vunpcklpd	32(%rbx), %ymm0, %ymm2
	vmovupd	64(%rbx), %ymm0
	vpermpd	$216, %ymm5, %ymm5
	vmovupd	%ymm5, 800(%rsp)
	vmovupd	%ymm4, 864(%rsp)
	vunpcklpd	96(%rbx), %ymm0, %ymm5
	vunpckhpd	96(%rbx), %ymm0, %ymm4
	vmovupd	128(%rbx), %ymm0
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm3, %ymm3
	vmovupd	%ymm9, 704(%rsp)
	vmovupd	%ymm6, 736(%rsp)
	vmovupd	%ymm3, 768(%rsp)
	vunpcklpd	160(%rbx), %ymm0, %ymm3
	vmovupd	192(%rbx), %ymm8
	vunpckhpd	160(%rbx), %ymm0, %ymm0
	vunpcklpd	224(%rbx), %ymm8, %ymm6
	vunpckhpd	224(%rbx), %ymm8, %ymm9
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm9, %ymm9
	vpermpd	$216, %ymm2, %ymm2
	vpermpd	$216, %ymm1, %ymm1
	vpermpd	$216, %ymm3, %ymm3
	vpermpd	$216, %ymm0, %ymm0
	vunpcklpd	%ymm5, %ymm2, %ymm8
	vunpckhpd	%ymm5, %ymm2, %ymm2
	vunpcklpd	%ymm6, %ymm3, %ymm5
	vunpckhpd	%ymm6, %ymm3, %ymm3
	vunpcklpd	%ymm4, %ymm1, %ymm6
	vunpckhpd	%ymm4, %ymm1, %ymm1
	vunpcklpd	%ymm9, %ymm0, %ymm4
	vunpckhpd	%ymm9, %ymm0, %ymm0
	vpermpd	$216, %ymm5, %ymm5
	vpermpd	$216, %ymm3, %ymm3
	vpermpd	$216, %ymm4, %ymm4
	vpermpd	$216, %ymm0, %ymm0
	vpermpd	$216, %ymm8, %ymm8
	vpermpd	$216, %ymm2, %ymm2
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm1, %ymm1
	vunpcklpd	%ymm5, %ymm8, %ymm9
	vunpckhpd	%ymm5, %ymm8, %ymm8
	vunpcklpd	%ymm4, %ymm6, %ymm5
	vunpckhpd	%ymm4, %ymm6, %ymm6
	vunpcklpd	%ymm3, %ymm2, %ymm4
	vunpckhpd	%ymm3, %ymm2, %ymm2
	vunpcklpd	%ymm0, %ymm1, %ymm3
	vunpckhpd	%ymm0, %ymm1, %ymm1
	vpermpd	$216, %ymm9, %ymm9
	vpermpd	$216, %ymm8, %ymm8
	vpermpd	$216, %ymm1, %ymm0
	vpermpd	$216, %ymm10, %ymm15
	vmovupd	%ymm0, 240(%rsp)
	vmulpd	320(%rsp), %ymm9, %ymm10
	vmulpd	64(%rsp), %ymm8, %ymm0
	vmovupd	(%rsp), %ymm1
	vpermpd	$216, %ymm12, %ymm12
	vpermpd	$216, %ymm7, %ymm7
	vfmadd231pd	176(%rsp), %ymm12, %ymm10
	vfmadd231pd	%ymm1, %ymm7, %ymm0
	vpermpd	$216, %ymm14, %ymm14
	vpermpd	$216, %ymm11, %ymm11
	vpermpd	$216, %ymm6, %ymm6
	vpermpd	$216, %ymm5, %ymm5
	vaddpd	%ymm10, %ymm0, %ymm0
	vmulpd	64(%rsp), %ymm9, %ymm10
	vpermpd	$216, %ymm2, %ymm2
	vsubpd	%ymm0, %ymm13, %ymm0
	vmulpd	320(%rsp), %ymm8, %ymm13
	vpermpd	$216, %ymm4, %ymm4
	vfmadd231pd	%ymm1, %ymm12, %ymm10
	vmovupd	%ymm0, 352(%rsp)
	vmovupd	208(%rsp), %ymm0
	vfmadd231pd	176(%rsp), %ymm7, %ymm13
	vpermpd	$216, %ymm3, %ymm3
	addq	$256, %r9
	addq	$256, %rax
	addq	$256, %rsi
	vsubpd	%ymm13, %ymm10, %ymm10
	vmulpd	%ymm0, %ymm9, %ymm13
	vmulpd	96(%rsp), %ymm9, %ymm9
	vaddpd	%ymm15, %ymm10, %ymm10
	vmulpd	96(%rsp), %ymm8, %ymm15
	vmulpd	%ymm0, %ymm8, %ymm8
	vmovupd	%ymm10, 928(%rsp)
	vmovupd	128(%rsp), %ymm10
	vfmadd231pd	32(%rsp), %ymm12, %ymm13
	vfmadd231pd	%ymm10, %ymm12, %ymm9
	vmovupd	32(%rsp), %ymm12
	vfmadd231pd	%ymm10, %ymm7, %ymm15
	vfmadd231pd	%ymm12, %ymm7, %ymm8
	vmovupd	(%rsp), %ymm7
	addq	$256, %rbx
	addq	$256, %r11
	vaddpd	%ymm15, %ymm13, %ymm13
	vsubpd	%ymm8, %ymm9, %ymm9
	vmovapd	%ymm10, %ymm15
	vmovupd	288(%rsp), %ymm10
	vsubpd	%ymm13, %ymm14, %ymm1
	vaddpd	704(%rsp), %ymm9, %ymm13
	vmulpd	%ymm15, %ymm10, %ymm9
	vmulpd	%ymm10, %ymm7, %ymm7
	vmovupd	176(%rsp), %ymm14
	vmovupd	320(%rsp), %ymm15
	vmovupd	%ymm1, 960(%rsp)
	vfmadd231pd	%ymm12, %ymm11, %ymm9
	vmovupd	64(%rsp), %ymm12
	vfmadd231pd	%ymm14, %ymm11, %ymm7
	vmulpd	%ymm12, %ymm6, %ymm8
	vmovupd	512(%rsp), %ymm1
	vfmadd231pd	%ymm15, %ymm5, %ymm8
	vaddpd	%ymm8, %ymm7, %ymm7
	vmulpd	%ymm12, %ymm5, %ymm8
	vmulpd	%ymm15, %ymm6, %ymm12
	vsubpd	%ymm7, %ymm1, %ymm7
	vfmadd231pd	(%rsp), %ymm11, %ymm8
	vfmadd231pd	%ymm14, %ymm10, %ymm12
	vsubpd	%ymm12, %ymm8, %ymm8
	vmulpd	96(%rsp), %ymm6, %ymm12
	vmulpd	%ymm0, %ymm6, %ymm6
	vaddpd	544(%rsp), %ymm8, %ymm8
	vmovupd	736(%rsp), %ymm1
	vmovupd	288(%rsp), %ymm10
	vfmadd231pd	%ymm0, %ymm5, %ymm12
	vmulpd	96(%rsp), %ymm5, %ymm5
	vmovupd	416(%rsp), %ymm0
	vaddpd	%ymm12, %ymm9, %ymm9
	vmovupd	32(%rsp), %ymm12
	vsubpd	%ymm9, %ymm1, %ymm9
	vmovupd	128(%rsp), %ymm1
	vfmadd231pd	%ymm12, %ymm10, %ymm6
	vfmadd231pd	%ymm1, %ymm11, %ymm5
	vmovupd	384(%rsp), %ymm10
	vmovupd	%ymm9, 512(%rsp)
	vsubpd	%ymm6, %ymm5, %ymm11
	vmulpd	%ymm1, %ymm0, %ymm5
	vmovupd	(%rsp), %ymm6
	vmovupd	576(%rsp), %ymm1
	vmulpd	%ymm0, %ymm6, %ymm9
	vaddpd	768(%rsp), %ymm11, %ymm11
	vfmadd231pd	%ymm12, %ymm10, %ymm5
	vmovupd	64(%rsp), %ymm12
	vmulpd	%ymm12, %ymm2, %ymm6
	vfmadd231pd	%ymm10, %ymm14, %ymm9
	vfmadd231pd	%ymm15, %ymm4, %ymm6
	vaddpd	%ymm9, %ymm6, %ymm6
	vmulpd	%ymm12, %ymm4, %ymm9
	vmulpd	%ymm15, %ymm2, %ymm12
	vsubpd	%ymm6, %ymm1, %ymm6
	vmovupd	800(%rsp), %ymm1
	vfmadd231pd	(%rsp), %ymm10, %ymm9
	vfmadd231pd	%ymm14, %ymm0, %ymm12
	vsubpd	%ymm12, %ymm9, %ymm9
	vmulpd	96(%rsp), %ymm2, %ymm12
	vmulpd	208(%rsp), %ymm2, %ymm2
	vaddpd	608(%rsp), %ymm9, %ymm9
	vfmadd231pd	208(%rsp), %ymm4, %ymm12
	vmulpd	96(%rsp), %ymm4, %ymm4
	vaddpd	%ymm12, %ymm5, %ymm5
	vfmadd231pd	128(%rsp), %ymm10, %ymm4
	vmovupd	480(%rsp), %ymm10
	vsubpd	%ymm5, %ymm1, %ymm5
	vmovapd	%ymm0, %ymm1
	vmovupd	32(%rsp), %ymm0
	vfmadd231pd	%ymm0, %ymm1, %ymm2
	vmovupd	448(%rsp), %ymm1
	vsubpd	%ymm2, %ymm4, %ymm4
	vmovupd	(%rsp), %ymm2
	vmulpd	%ymm10, %ymm2, %ymm12
	vmulpd	128(%rsp), %ymm10, %ymm2
	vaddpd	832(%rsp), %ymm4, %ymm4
	vfmadd231pd	%ymm1, %ymm14, %ymm12
	vfmadd231pd	%ymm0, %ymm1, %ymm2
	vmovupd	240(%rsp), %ymm0
	vmulpd	64(%rsp), %ymm0, %ymm14
	vfmadd231pd	%ymm15, %ymm3, %ymm14
	vmulpd	%ymm0, %ymm15, %ymm15
	vaddpd	%ymm14, %ymm12, %ymm12
	vmovupd	640(%rsp), %ymm14
	vfmadd231pd	176(%rsp), %ymm10, %ymm15
	vsubpd	%ymm12, %ymm14, %ymm12
	vmulpd	64(%rsp), %ymm3, %ymm14
	vfmadd231pd	(%rsp), %ymm1, %ymm14
	vsubpd	%ymm15, %ymm14, %ymm14
	vaddpd	672(%rsp), %ymm14, %ymm14
	vmulpd	96(%rsp), %ymm0, %ymm15
	vmovupd	208(%rsp), %ymm0
	vfmadd231pd	%ymm0, %ymm3, %ymm15
	vmulpd	96(%rsp), %ymm3, %ymm3
	vaddpd	%ymm15, %ymm2, %ymm2
	vmovupd	864(%rsp), %ymm15
	vfmadd231pd	128(%rsp), %ymm1, %ymm3
	vsubpd	%ymm2, %ymm15, %ymm2
	vmovupd	240(%rsp), %ymm15
	vmulpd	%ymm0, %ymm15, %ymm1
	vpermpd	$68, 352(%rsp), %ymm15
	vpermpd	$238, 352(%rsp), %ymm0
	vfmadd231pd	32(%rsp), %ymm10, %ymm1
	vmovupd	928(%rsp), %ymm10
	vsubpd	%ymm1, %ymm3, %ymm1
	vpermpd	$68, %ymm10, %ymm3
	vpermpd	$238, %ymm10, %ymm10
	vshufpd	$12, %ymm3, %ymm15, %ymm3
	vshufpd	$12, %ymm10, %ymm0, %ymm10
	vpermpd	$68, %ymm7, %ymm15
	vpermpd	$68, %ymm8, %ymm0
	vpermpd	$238, %ymm7, %ymm7
	vpermpd	$238, %ymm8, %ymm8
	vshufpd	$12, %ymm0, %ymm15, %ymm15
	vshufpd	$12, %ymm8, %ymm7, %ymm7
	vpermpd	$68, %ymm9, %ymm0
	vpermpd	$68, %ymm6, %ymm8
	vshufpd	$12, %ymm0, %ymm8, %ymm8
	vpermpd	$238, %ymm6, %ymm6
	vpermpd	$238, %ymm9, %ymm0
	vshufpd	$12, %ymm0, %ymm6, %ymm0
	vpermpd	$68, %ymm12, %ymm9
	vpermpd	$68, %ymm14, %ymm6
	vpermpd	$238, %ymm12, %ymm12
	vpermpd	$238, %ymm14, %ymm14
	vshufpd	$12, %ymm6, %ymm9, %ymm6
	vshufpd	$12, %ymm14, %ymm12, %ymm12
	vpermpd	$68, %ymm8, %ymm9
	vpermpd	$68, %ymm3, %ymm14
	vpermpd	$238, %ymm8, %ymm8
	vpermpd	$238, %ymm3, %ymm3
	vshufpd	$12, %ymm9, %ymm14, %ymm9
	vshufpd	$12, %ymm8, %ymm3, %ymm8
	vpermpd	$68, %ymm10, %ymm14
	vpermpd	$68, %ymm0, %ymm3
	vpermpd	$238, %ymm10, %ymm10
	vpermpd	$238, %ymm0, %ymm0
	vshufpd	$12, %ymm3, %ymm14, %ymm3
	vshufpd	$12, %ymm0, %ymm10, %ymm0
	vpermpd	$68, %ymm15, %ymm14
	vpermpd	$68, %ymm6, %ymm10
	vpermpd	$238, %ymm15, %ymm15
	vpermpd	$238, %ymm6, %ymm6
	vshufpd	$12, %ymm10, %ymm14, %ymm10
	vshufpd	$12, %ymm6, %ymm15, %ymm15
	vpermpd	$68, %ymm7, %ymm14
	vpermpd	$68, %ymm12, %ymm6
	vpermpd	$238, %ymm7, %ymm7
	vpermpd	$238, %ymm12, %ymm12
	vshufpd	$12, %ymm6, %ymm14, %ymm6
	vshufpd	$12, %ymm12, %ymm7, %ymm7
	vpermpd	$68, %ymm10, %ymm14
	vpermpd	$68, %ymm9, %ymm12
	vpermpd	$238, %ymm10, %ymm10
	vpermpd	$238, %ymm9, %ymm9
	vshufpd	$12, %ymm10, %ymm9, %ymm9
	vpermpd	$68, %ymm15, %ymm10
	vmovupd	%ymm9, -224(%rax)
	vpermpd	$238, %ymm15, %ymm15
	vpermpd	$68, %ymm8, %ymm9
	vpermpd	$238, %ymm8, %ymm8
	vshufpd	$12, %ymm10, %ymm9, %ymm9
	vshufpd	$12, %ymm15, %ymm8, %ymm8
	vmovupd	%ymm9, -192(%rax)
	vmovupd	%ymm8, -160(%rax)
	vpermpd	$68, %ymm6, %ymm9
	vpermpd	$68, %ymm3, %ymm8
	vpermpd	$238, %ymm6, %ymm6
	vpermpd	$238, %ymm3, %ymm3
	vshufpd	$12, %ymm6, %ymm3, %ymm3
	vpermpd	$68, %ymm7, %ymm6
	vmovupd	%ymm3, -96(%rax)
	vpermpd	$238, %ymm7, %ymm7
	vpermpd	$68, %ymm0, %ymm3
	vpermpd	$238, %ymm0, %ymm0
	vshufpd	$12, %ymm7, %ymm0, %ymm0
	vmovupd	960(%rsp), %ymm7
	vshufpd	$12, %ymm6, %ymm3, %ymm3
	vshufpd	$12, %ymm14, %ymm12, %ymm12
	vmovupd	%ymm3, -64(%rax)
	vpermpd	$238, %ymm7, %ymm14
	vpermpd	$68, %ymm7, %ymm3
	vmovupd	512(%rsp), %ymm7
	vmovupd	%ymm0, -32(%rax)
	vaddpd	896(%rsp), %ymm1, %ymm1
	vpermpd	$68, %ymm13, %ymm0
	vshufpd	$12, %ymm0, %ymm3, %ymm3
	vpermpd	$68, %ymm7, %ymm6
	vpermpd	$68, %ymm11, %ymm0
	vshufpd	$12, %ymm9, %ymm8, %ymm8
	vshufpd	$12, %ymm0, %ymm6, %ymm6
	vmovupd	%ymm8, -128(%rax)
	vpermpd	$238, %ymm7, %ymm9
	vpermpd	$68, %ymm4, %ymm0
	vpermpd	$68, %ymm5, %ymm8
	vpermpd	$238, %ymm11, %ymm11
	vshufpd	$12, %ymm11, %ymm9, %ymm11
	vshufpd	$12, %ymm0, %ymm8, %ymm8
	vpermpd	$68, %ymm2, %ymm9
	vpermpd	$68, %ymm1, %ymm0
	vshufpd	$12, %ymm0, %ymm9, %ymm9
	vpermpd	$238, %ymm5, %ymm5
	vpermpd	$68, %ymm8, %ymm0
	vpermpd	$68, %ymm3, %ymm7
	vpermpd	$238, %ymm13, %ymm13
	vpermpd	$238, %ymm4, %ymm4
	vshufpd	$12, %ymm4, %ymm5, %ymm4
	vshufpd	$12, %ymm0, %ymm7, %ymm7
	vpermpd	$238, %ymm2, %ymm2
	vpermpd	$68, %ymm4, %ymm0
	vshufpd	$12, %ymm13, %ymm14, %ymm13
	vpermpd	$238, %ymm4, %ymm4
	vpermpd	$68, %ymm13, %ymm10
	vpermpd	$238, %ymm1, %ymm1
	vpermpd	$238, %ymm13, %ymm13
	vshufpd	$12, %ymm1, %ymm2, %ymm1
	vshufpd	$12, %ymm0, %ymm10, %ymm10
	vpermpd	$68, %ymm9, %ymm2
	vshufpd	$12, %ymm4, %ymm13, %ymm0
	vpermpd	$68, %ymm6, %ymm4
	vshufpd	$12, %ymm2, %ymm4, %ymm4
	vpermpd	$238, %ymm8, %ymm8
	vpermpd	$68, %ymm1, %ymm2
	vpermpd	$68, %ymm11, %ymm5
	vpermpd	$238, %ymm3, %ymm3
	vshufpd	$12, %ymm8, %ymm3, %ymm3
	vshufpd	$12, %ymm2, %ymm5, %ymm5
	vpermpd	$68, %ymm4, %ymm8
	vpermpd	$68, %ymm7, %ymm2
	vpermpd	$238, %ymm4, %ymm4
	vpermpd	$238, %ymm6, %ymm6
	vpermpd	$238, %ymm9, %ymm9
	vpermpd	$238, %ymm7, %ymm7
	vmovupd	%ymm12, -256(%rax)
	vshufpd	$12, %ymm9, %ymm6, %ymm6
	vshufpd	$12, %ymm8, %ymm2, %ymm2
	vshufpd	$12, %ymm4, %ymm7, %ymm7
	vmovupd	%ymm2, -256(%r11)
	vpermpd	$68, %ymm6, %ymm4
	vpermpd	$68, %ymm3, %ymm2
	vpermpd	$238, %ymm6, %ymm6
	vpermpd	$238, %ymm3, %ymm3
	vshufpd	$12, %ymm4, %ymm2, %ymm2
	vshufpd	$12, %ymm6, %ymm3, %ymm3
	vmovupd	%ymm2, -192(%r11)
	vmovupd	%ymm3, -160(%r11)
	vpermpd	$68, %ymm10, %ymm2
	vpermpd	$68, %ymm5, %ymm3
	vpermpd	$238, %ymm11, %ymm11
	vpermpd	$238, %ymm1, %ymm1
	vshufpd	$12, %ymm3, %ymm2, %ymm2
	vshufpd	$12, %ymm1, %ymm11, %ymm1
	vmovupd	%ymm2, -128(%r11)
	vpermpd	$68, %ymm1, %ymm3
	vpermpd	$68, %ymm0, %ymm2
	vpermpd	$238, %ymm10, %ymm10
	vpermpd	$238, %ymm5, %ymm5
	vpermpd	$238, %ymm0, %ymm0
	vpermpd	$238, %ymm1, %ymm1
	vshufpd	$12, %ymm5, %ymm10, %ymm5
	vshufpd	$12, %ymm3, %ymm2, %ymm2
	vshufpd	$12, %ymm1, %ymm0, %ymm1
	vmovupd	%ymm7, -224(%r11)
	vmovupd	%ymm5, -96(%r11)
	vmovupd	%ymm2, -64(%r11)
	vmovupd	%ymm1, -32(%r11)
	cmpq	%r9, %rdi
	jne	.L5
	movl	%ebp, %eax
	andl	$-4, %eax
	leal	0(,%rax,8), %edi
	cmpl	%eax, %ebp
	je	.L20
	subl	%eax, %ebp
	cmpl	%eax, %r13d
	je	.L24
	vzeroupper
.L4:
	movq	%rax, %r11
	salq	$6, %r11
	leaq	(%r10,%r11), %rbx
	vmovupd	16(%rbx), %xmm1
	vmovupd	32(%rbx), %xmm7
	vmovupd	64(%rbx), %xmm3
	vmovupd	96(%rbx), %xmm4
	vmovupd	(%rbx), %xmm8
	vmovupd	48(%rbx), %xmm5
	vmovupd	80(%rbx), %xmm0
	vmovupd	112(%rbx), %xmm2
	vunpcklpd	%xmm1, %xmm8, %xmm6
	vunpcklpd	%xmm0, %xmm3, %xmm13
	vunpckhpd	%xmm1, %xmm8, %xmm8
	vunpckhpd	%xmm0, %xmm3, %xmm0
	vunpcklpd	%xmm5, %xmm7, %xmm1
	vunpcklpd	%xmm2, %xmm4, %xmm3
	vunpckhpd	%xmm5, %xmm7, %xmm5
	vunpckhpd	%xmm2, %xmm4, %xmm2
	leaq	(%rcx,%r11), %r9
	vunpcklpd	%xmm3, %xmm13, %xmm14
	vunpckhpd	%xmm3, %xmm13, %xmm13
	vunpcklpd	%xmm5, %xmm8, %xmm3
	vunpckhpd	%xmm5, %xmm8, %xmm8
	vunpcklpd	%xmm2, %xmm0, %xmm5
	vunpcklpd	%xmm1, %xmm6, %xmm4
	vunpckhpd	%xmm2, %xmm0, %xmm0
	vunpckhpd	%xmm1, %xmm6, %xmm1
	vmovupd	32(%r9), %xmm2
	vunpcklpd	%xmm5, %xmm3, %xmm6
	vmovapd	%xmm6, 32(%rsp)
	vunpckhpd	%xmm5, %xmm3, %xmm6
	vunpcklpd	%xmm13, %xmm1, %xmm5
	vmovhpd	48(%r9), %xmm2, %xmm3
	vmovupd	96(%r9), %xmm15
	vmovupd	48(%r9), %xmm2
	vmovapd	%xmm5, 176(%rsp)
	vunpckhpd	%xmm13, %xmm1, %xmm5
	vmovupd	(%r9), %xmm1
	vunpcklpd	%xmm14, %xmm4, %xmm9
	vunpcklpd	%xmm0, %xmm8, %xmm7
	vunpckhpd	%xmm14, %xmm4, %xmm14
	vmovapd	%xmm6, 64(%rsp)
	vunpckhpd	%xmm0, %xmm8, %xmm4
	vmovhpd	16(%r9), %xmm1, %xmm6
	vmovlpd	40(%r9), %xmm2, %xmm8
	vmovupd	16(%r9), %xmm1
	vmovupd	64(%r9), %xmm2
	vmovupd	80(%r9), %xmm0
	vmovapd	%xmm5, 96(%rsp)
	vmovhpd	112(%r9), %xmm15, %xmm5
	vmovupd	112(%r9), %xmm15
	vmovddup	1000(%rsp), %xmm11
	vmovlpd	8(%r9), %xmm1, %xmm1
	vmovlpd	72(%r9), %xmm0, %xmm0
	vmovhpd	80(%r9), %xmm2, %xmm2
	vmovapd	%xmm9, 128(%rsp)
	vmovapd	%xmm14, (%rsp)
	vmovapd	%xmm7, 208(%rsp)
	vmovapd	%xmm4, 288(%rsp)
	vmovlpd	104(%r9), %xmm15, %xmm4
	leaq	(%r12,%rax,8), %rax
	leaq	(%rcx,%rax,8), %rax
	vunpcklpd	%xmm3, %xmm6, %xmm7
	vunpckhpd	%xmm3, %xmm6, %xmm3
	vunpcklpd	%xmm5, %xmm2, %xmm6
	vunpckhpd	%xmm5, %xmm2, %xmm2
	vunpcklpd	%xmm8, %xmm1, %xmm5
	vunpckhpd	%xmm8, %xmm1, %xmm1
	vunpcklpd	%xmm4, %xmm0, %xmm8
	vunpckhpd	%xmm4, %xmm0, %xmm0
	vunpcklpd	%xmm2, %xmm3, %xmm4
	vunpckhpd	%xmm2, %xmm3, %xmm3
	vunpcklpd	%xmm0, %xmm1, %xmm2
	vunpckhpd	%xmm0, %xmm1, %xmm1
	vmovupd	(%rax), %xmm0
	vunpcklpd	%xmm6, %xmm7, %xmm15
	vunpckhpd	%xmm6, %xmm7, %xmm12
	vmovhpd	16(%rax), %xmm0, %xmm6
	vmovupd	16(%rax), %xmm0
	vunpcklpd	%xmm8, %xmm5, %xmm10
	vunpckhpd	%xmm8, %xmm5, %xmm7
	vmovlpd	8(%rax), %xmm0, %xmm8
	vmovupd	32(%rax), %xmm0
	vmovapd	%xmm1, 544(%rsp)
	vmovhpd	48(%rax), %xmm0, %xmm1
	vmovupd	48(%rax), %xmm0
	vmovapd	%xmm15, 320(%rsp)
	vmovlpd	40(%rax), %xmm0, %xmm5
	vmovupd	96(%rax), %xmm15
	vmovupd	64(%rax), %xmm0
	vmovapd	%xmm3, 480(%rsp)
	vmovhpd	80(%rax), %xmm0, %xmm13
	vmovhpd	112(%rax), %xmm15, %xmm3
	vmovupd	80(%rax), %xmm0
	vmovupd	112(%rax), %xmm15
	vmovlpd	72(%rax), %xmm0, %xmm0
	vmovapd	%xmm2, 512(%rsp)
	vmovlpd	104(%rax), %xmm15, %xmm2
	vunpcklpd	%xmm3, %xmm13, %xmm14
	vunpckhpd	%xmm3, %xmm13, %xmm13
	vunpcklpd	%xmm5, %xmm8, %xmm3
	vunpckhpd	%xmm5, %xmm8, %xmm8
	vunpcklpd	%xmm2, %xmm0, %xmm5
	vmovapd	%xmm12, 352(%rsp)
	vmovapd	%xmm4, 448(%rsp)
	vunpckhpd	%xmm5, %xmm3, %xmm12
	vunpcklpd	%xmm1, %xmm6, %xmm4
	vunpckhpd	%xmm1, %xmm6, %xmm1
	vunpckhpd	%xmm2, %xmm0, %xmm0
	vmovapd	%xmm12, 640(%rsp)
	vunpcklpd	%xmm13, %xmm1, %xmm12
	vunpckhpd	%xmm13, %xmm1, %xmm1
	addq	%r8, %r11
	vunpckhpd	%xmm14, %xmm4, %xmm15
	vunpcklpd	%xmm5, %xmm3, %xmm2
	vmovapd	%xmm1, 704(%rsp)
	vunpcklpd	%xmm0, %xmm8, %xmm1
	vunpckhpd	%xmm0, %xmm8, %xmm0
	vmovapd	%xmm10, 384(%rsp)
	vmovapd	%xmm7, 416(%rsp)
	vunpcklpd	%xmm14, %xmm4, %xmm10
	vmovapd	%xmm15, 576(%rsp)
	vmovapd	%xmm2, 608(%rsp)
	vmovapd	%xmm12, 672(%rsp)
	vmovapd	%xmm1, 736(%rsp)
	vmovapd	%xmm0, 768(%rsp)
	vmovupd	(%r11), %xmm9
	vmovupd	16(%r11), %xmm2
	vmovupd	48(%r11), %xmm0
	vmovupd	64(%r11), %xmm7
	vmovupd	80(%r11), %xmm4
	vmovupd	96(%r11), %xmm5
	vmovupd	32(%r11), %xmm1
	vmovupd	112(%r11), %xmm3
	vunpcklpd	%xmm2, %xmm9, %xmm6
	vunpckhpd	%xmm4, %xmm7, %xmm8
	vunpckhpd	%xmm2, %xmm9, %xmm9
	vunpcklpd	%xmm0, %xmm1, %xmm2
	vunpckhpd	%xmm0, %xmm1, %xmm1
	vunpcklpd	%xmm4, %xmm7, %xmm0
	vunpcklpd	%xmm3, %xmm5, %xmm4
	vunpckhpd	%xmm3, %xmm5, %xmm3
	vunpckhpd	%xmm3, %xmm8, %xmm7
	vunpcklpd	%xmm2, %xmm6, %xmm13
	vunpckhpd	%xmm2, %xmm6, %xmm5
	vunpcklpd	%xmm1, %xmm9, %xmm12
	vunpcklpd	%xmm4, %xmm0, %xmm2
	vunpckhpd	%xmm1, %xmm9, %xmm1
	vunpcklpd	%xmm7, %xmm1, %xmm15
	vunpckhpd	%xmm4, %xmm0, %xmm4
	vunpckhpd	%xmm7, %xmm1, %xmm1
	vunpcklpd	%xmm3, %xmm8, %xmm0
	vunpcklpd	%xmm2, %xmm13, %xmm3
	vunpckhpd	%xmm2, %xmm13, %xmm13
	vunpckhpd	%xmm4, %xmm5, %xmm6
	vunpcklpd	%xmm0, %xmm12, %xmm2
	vmovapd	%xmm1, 240(%rsp)
	vunpckhpd	%xmm0, %xmm12, %xmm12
	vmovddup	1008(%rsp), %xmm1
	vunpcklpd	%xmm4, %xmm5, %xmm0
	vmulpd	%xmm1, %xmm13, %xmm7
	vmovddup	992(%rsp), %xmm5
	vmulpd	%xmm5, %xmm3, %xmm8
	vmovapd	128(%rsp), %xmm9
	vmovapd	(%rsp), %xmm14
	vmovddup	280(%rsp), %xmm4
	vfmadd231pd	%xmm14, %xmm11, %xmm7
	vfmadd231pd	%xmm9, %xmm4, %xmm8
	vaddpd	%xmm8, %xmm7, %xmm7
	vmovapd	320(%rsp), %xmm8
	vsubpd	%xmm7, %xmm8, %xmm8
	vmulpd	%xmm1, %xmm3, %xmm7
	vmovapd	%xmm8, 320(%rsp)
	vmulpd	%xmm5, %xmm13, %xmm8
	vfmadd231pd	%xmm9, %xmm11, %xmm7
	vfmadd231pd	%xmm14, %xmm4, %xmm8
	vsubpd	%xmm8, %xmm7, %xmm7
	vmovddup	160(%rsp), %xmm8
	vfnmadd231pd	%xmm9, %xmm8, %xmm10
	vaddpd	352(%rsp), %xmm7, %xmm7
	vmovddup	1016(%rsp), %xmm9
	vmovapd	%xmm7, 352(%rsp)
	vmovddup	272(%rsp), %xmm7
	vmulpd	%xmm7, %xmm3, %xmm14
	vfmadd231pd	(%rsp), %xmm9, %xmm14
	vsubpd	%xmm14, %xmm10, %xmm14
	vmovddup	168(%rsp), %xmm10
	vmulpd	%xmm10, %xmm3, %xmm3
	vfnmadd231pd	%xmm10, %xmm13, %xmm14
	vfmadd231pd	128(%rsp), %xmm9, %xmm3
	vmovapd	%xmm14, 800(%rsp)
	vmovapd	(%rsp), %xmm14
	vfnmadd213pd	576(%rsp), %xmm8, %xmm14
	vaddpd	%xmm14, %xmm3, %xmm3
	vmulpd	64(%rsp), %xmm11, %xmm14
	vfnmadd132pd	%xmm7, %xmm3, %xmm13
	vmovapd	32(%rsp), %xmm3
	vfmadd231pd	%xmm3, %xmm4, %xmm14
	vmovapd	%xmm13, (%rsp)
	vmovapd	608(%rsp), %xmm13
	vfnmadd231pd	%xmm3, %xmm8, %xmm13
	vmulpd	%xmm1, %xmm12, %xmm3
	vfmadd231pd	%xmm5, %xmm2, %xmm3
	vaddpd	%xmm14, %xmm3, %xmm3
	vmovapd	384(%rsp), %xmm14
	vsubpd	%xmm3, %xmm14, %xmm3
	vmulpd	%xmm5, %xmm12, %xmm14
	vmovapd	%xmm3, 128(%rsp)
	vmulpd	%xmm1, %xmm2, %xmm3
	vfmadd231pd	64(%rsp), %xmm4, %xmm14
	vfmadd231pd	32(%rsp), %xmm11, %xmm3
	vsubpd	%xmm14, %xmm3, %xmm3
	vmulpd	%xmm7, %xmm2, %xmm14
	vmulpd	%xmm10, %xmm2, %xmm2
	vaddpd	416(%rsp), %xmm3, %xmm3
	vfmadd231pd	64(%rsp), %xmm9, %xmm14
	vfmadd231pd	32(%rsp), %xmm9, %xmm2
	vsubpd	%xmm14, %xmm13, %xmm13
	vmovapd	176(%rsp), %xmm14
	vfnmadd231pd	%xmm10, %xmm12, %xmm13
	vmovapd	%xmm13, 384(%rsp)
	vmovapd	64(%rsp), %xmm13
	vfnmadd213pd	640(%rsp), %xmm8, %xmm13
	vaddpd	%xmm13, %xmm2, %xmm2
	vfnmadd132pd	%xmm7, %xmm2, %xmm12
	vmulpd	96(%rsp), %xmm11, %xmm2
	vmovapd	%xmm12, 32(%rsp)
	vmovapd	672(%rsp), %xmm12
	vfmadd231pd	%xmm14, %xmm4, %xmm2
	vfnmadd231pd	%xmm14, %xmm8, %xmm12
	vmovapd	%xmm12, %xmm13
	vmulpd	%xmm1, %xmm6, %xmm12
	vfmadd231pd	%xmm5, %xmm0, %xmm12
	vaddpd	%xmm2, %xmm12, %xmm12
	vmovapd	448(%rsp), %xmm2
	vsubpd	%xmm12, %xmm2, %xmm12
	vmulpd	%xmm1, %xmm0, %xmm2
	vfmadd231pd	%xmm14, %xmm11, %xmm2
	vmulpd	%xmm5, %xmm6, %xmm14
	vfmadd231pd	96(%rsp), %xmm4, %xmm14
	vsubpd	%xmm14, %xmm2, %xmm2
	vmulpd	%xmm7, %xmm0, %xmm14
	vmulpd	%xmm10, %xmm0, %xmm0
	vaddpd	480(%rsp), %xmm2, %xmm2
	vfmadd231pd	96(%rsp), %xmm9, %xmm14
	vfmadd231pd	176(%rsp), %xmm9, %xmm0
	vsubpd	%xmm14, %xmm13, %xmm13
	vmulpd	288(%rsp), %xmm11, %xmm14
	vfnmadd231pd	%xmm10, %xmm6, %xmm13
	vmovapd	%xmm13, 64(%rsp)
	vmovapd	96(%rsp), %xmm13
	vfnmadd213pd	704(%rsp), %xmm8, %xmm13
	vaddpd	%xmm13, %xmm0, %xmm0
	vmovapd	736(%rsp), %xmm13
	vfnmadd231pd	%xmm7, %xmm6, %xmm0
	vmovapd	208(%rsp), %xmm6
	vfmadd231pd	%xmm6, %xmm4, %xmm14
	vfnmadd231pd	%xmm6, %xmm8, %xmm13
	vmulpd	240(%rsp), %xmm1, %xmm6
	vmulpd	%xmm1, %xmm15, %xmm1
	vfmadd231pd	%xmm5, %xmm15, %xmm6
	vmulpd	240(%rsp), %xmm5, %xmm5
	vaddpd	%xmm14, %xmm6, %xmm6
	vmovapd	512(%rsp), %xmm14
	vsubpd	%xmm6, %xmm14, %xmm6
	vmovapd	208(%rsp), %xmm14
	vfmadd231pd	%xmm14, %xmm11, %xmm1
	vmovapd	288(%rsp), %xmm11
	vfmadd231pd	%xmm11, %xmm4, %xmm5
	vmulpd	%xmm7, %xmm15, %xmm4
	vfnmadd213pd	768(%rsp), %xmm11, %xmm8
	vsubpd	%xmm5, %xmm1, %xmm1
	vfmadd231pd	%xmm11, %xmm9, %xmm4
	vaddpd	544(%rsp), %xmm1, %xmm1
	vsubpd	%xmm4, %xmm13, %xmm13
	vmovapd	240(%rsp), %xmm4
	vfnmadd231pd	%xmm4, %xmm10, %xmm13
	vmulpd	%xmm10, %xmm15, %xmm10
	vfmadd231pd	%xmm14, %xmm9, %xmm10
	vaddpd	%xmm8, %xmm10, %xmm10
	vfnmadd132pd	%xmm4, %xmm10, %xmm7
	vmovapd	320(%rsp), %xmm4
	vmovapd	352(%rsp), %xmm11
	vunpckhpd	%xmm11, %xmm4, %xmm9
	vunpcklpd	%xmm11, %xmm4, %xmm10
	vmovapd	128(%rsp), %xmm4
	vunpcklpd	%xmm1, %xmm6, %xmm11
	vunpcklpd	%xmm3, %xmm4, %xmm8
	vunpckhpd	%xmm3, %xmm4, %xmm3
	vunpcklpd	%xmm2, %xmm12, %xmm4
	vunpckhpd	%xmm2, %xmm12, %xmm2
	vunpckhpd	%xmm1, %xmm6, %xmm1
	vunpcklpd	%xmm4, %xmm10, %xmm14
	vunpcklpd	%xmm2, %xmm9, %xmm5
	vunpckhpd	%xmm2, %xmm9, %xmm2
	vunpcklpd	%xmm11, %xmm8, %xmm9
	vunpcklpd	%xmm1, %xmm3, %xmm12
	vunpckhpd	%xmm11, %xmm8, %xmm8
	vunpckhpd	%xmm1, %xmm3, %xmm3
	vunpckhpd	%xmm4, %xmm10, %xmm4
	vunpcklpd	%xmm9, %xmm14, %xmm1
	vmovupd	%xmm1, (%r9)
	vunpcklpd	%xmm8, %xmm4, %xmm1
	vmovupd	%xmm1, 32(%r9)
	vunpcklpd	%xmm12, %xmm5, %xmm1
	vunpckhpd	%xmm12, %xmm5, %xmm5
	vmovupd	%xmm5, 80(%r9)
	vmovapd	(%rsp), %xmm6
	vmovapd	800(%rsp), %xmm5
	vunpckhpd	%xmm8, %xmm4, %xmm4
	vmovupd	%xmm4, 48(%r9)
	vunpckhpd	%xmm9, %xmm14, %xmm9
	vunpcklpd	%xmm6, %xmm5, %xmm4
	vunpckhpd	%xmm6, %xmm5, %xmm14
	vmovapd	384(%rsp), %xmm5
	vmovapd	32(%rsp), %xmm6
	vmovupd	%xmm1, 64(%r9)
	vunpcklpd	%xmm3, %xmm2, %xmm1
	vunpckhpd	%xmm3, %xmm2, %xmm2
	vmovupd	%xmm2, 112(%r9)
	vunpcklpd	%xmm6, %xmm5, %xmm3
	vunpckhpd	%xmm6, %xmm5, %xmm2
	vmovapd	64(%rsp), %xmm5
	vmovupd	%xmm1, 96(%r9)
	vunpcklpd	%xmm0, %xmm5, %xmm1
	vunpckhpd	%xmm0, %xmm5, %xmm0
	vunpcklpd	%xmm7, %xmm13, %xmm5
	vunpcklpd	%xmm1, %xmm4, %xmm6
	vunpckhpd	%xmm7, %xmm13, %xmm7
	vunpckhpd	%xmm1, %xmm4, %xmm1
	vunpcklpd	%xmm5, %xmm3, %xmm4
	vmovupd	%xmm9, 16(%r9)
	vunpckhpd	%xmm5, %xmm3, %xmm3
	vunpcklpd	%xmm7, %xmm2, %xmm9
	vunpckhpd	%xmm7, %xmm2, %xmm7
	vunpcklpd	%xmm4, %xmm6, %xmm2
	vunpcklpd	%xmm0, %xmm14, %xmm8
	vmovupd	%xmm2, (%rax)
	vunpcklpd	%xmm3, %xmm1, %xmm2
	vunpckhpd	%xmm3, %xmm1, %xmm1
	vunpckhpd	%xmm0, %xmm14, %xmm0
	vmovupd	%xmm1, 48(%rax)
	vunpcklpd	%xmm9, %xmm8, %xmm1
	vmovupd	%xmm1, 64(%rax)
	vunpckhpd	%xmm4, %xmm6, %xmm4
	vunpcklpd	%xmm7, %xmm0, %xmm1
	vunpckhpd	%xmm9, %xmm8, %xmm8
	vunpckhpd	%xmm7, %xmm0, %xmm0
	vmovupd	%xmm4, 16(%rax)
	vmovupd	%xmm2, 32(%rax)
	vmovupd	%xmm8, 80(%rax)
	vmovupd	%xmm1, 96(%rax)
	vmovupd	%xmm0, 112(%rax)
	movl	%ebp, %eax
	andl	$-2, %eax
	leal	(%rdi,%rax,8), %edi
	cmpl	%eax, %ebp
	je	.L22
.L7:
	movslq	%edi, %rdi
	leaq	0(,%rdi,8), %rax
	vmovsd	992(%rsp), %xmm13
	vmovsd	1008(%rsp), %xmm12
	vmovsd	(%r8,%rdi,8), %xmm3
	vmovsd	32(%r8,%rax), %xmm1
	vmulsd	%xmm13, %xmm3, %xmm9
	vmulsd	%xmm12, %xmm1, %xmm8
	vmovsd	1000(%rsp), %xmm11
	vmovsd	280(%rsp), %xmm14
	vmovsd	32(%r10,%rax), %xmm4
	vmovsd	(%r10,%rdi,8), %xmm2
	vfmadd231sd	%xmm11, %xmm4, %xmm8
	vfmadd231sd	%xmm14, %xmm2, %xmm9
	leaq	(%rcx,%rax), %rsi
	vmovsd	(%rsi), %xmm7
	leaq	32(%rax), %r9
	leaq	(%rcx,%r9), %rbx
	vaddsd	%xmm9, %xmm8, %xmm8
	vmovsd	%xmm11, %xmm11, %xmm5
	vfmadd213sd	(%rbx), %xmm2, %xmm5
	vsubsd	%xmm8, %xmm7, %xmm8
	vmulsd	%xmm13, %xmm1, %xmm7
	vmovsd	272(%rsp), %xmm10
	vmovsd	1016(%rsp), %xmm15
	vmovsd	160(%rsp), %xmm9
	addq	%rdx, %r9
	vfmadd231sd	%xmm14, %xmm4, %xmm7
	vmovsd	%xmm15, %xmm15, %xmm0
	vfmadd213sd	(%r9), %xmm2, %xmm0
	leaq	(%rdx,%rax), %r11
	vmovsd	(%r11), %xmm6
	vsubsd	%xmm7, %xmm5, %xmm5
	vmulsd	%xmm10, %xmm3, %xmm7
	vmovsd	%xmm8, (%rsi)
	vfmadd231sd	%xmm12, %xmm3, %xmm5
	vfmadd132sd	%xmm9, %xmm7, %xmm2
	vmulsd	168(%rsp), %xmm1, %xmm7
	vmulsd	%xmm10, %xmm1, %xmm1
	vmovsd	%xmm5, (%rbx)
	vfmadd231sd	%xmm15, %xmm4, %xmm7
	vfmadd231sd	%xmm9, %xmm4, %xmm1
	vmovsd	40(%r10,%rax), %xmm4
	vmulsd	%xmm11, %xmm4, %xmm10
	vaddsd	%xmm2, %xmm7, %xmm2
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	40(%r8,%rax), %xmm1
	vsubsd	%xmm2, %xmm6, %xmm2
	vfmadd231sd	168(%rsp), %xmm3, %xmm0
	vmulsd	%xmm12, %xmm1, %xmm8
	vmovsd	%xmm2, (%r11)
	vmovsd	8(%r8,%rax), %xmm2
	leaq	8(%rax), %r11
	vmovsd	%xmm0, (%r9)
	vmovsd	8(%r10,%rax), %xmm0
	vfmadd231sd	%xmm13, %xmm2, %xmm8
	vfmadd231sd	%xmm14, %xmm0, %xmm10
	leaq	(%rcx,%r11), %rsi
	vmovsd	(%rsi), %xmm3
	vmulsd	%xmm13, %xmm1, %xmm5
	vmulsd	%xmm15, %xmm4, %xmm9
	vaddsd	%xmm10, %xmm8, %xmm8
	leaq	40(%rax), %r9
	vmovsd	272(%rsp), %xmm10
	vsubsd	%xmm8, %xmm3, %xmm8
	vmulsd	%xmm12, %xmm2, %xmm3
	vfmadd231sd	%xmm14, %xmm4, %xmm5
	leaq	(%rcx,%r9), %rbx
	addq	%rdx, %r11
	addq	%rdx, %r9
	vfmadd231sd	%xmm11, %xmm0, %xmm3
	vfmadd231sd	160(%rsp), %xmm0, %xmm9
	vmovsd	(%r11), %xmm6
	vmovsd	(%r9), %xmm7
	vsubsd	%xmm5, %xmm3, %xmm3
	vmulsd	168(%rsp), %xmm1, %xmm5
	vmulsd	%xmm10, %xmm1, %xmm1
	vaddsd	(%rbx), %xmm3, %xmm3
	vfmadd231sd	%xmm10, %xmm2, %xmm5
	vmulsd	168(%rsp), %xmm2, %xmm2
	vmovsd	%xmm8, (%rsi)
	vmovsd	%xmm3, (%rbx)
	vaddsd	%xmm9, %xmm5, %xmm5
	vfmadd132sd	%xmm15, %xmm2, %xmm0
	vmovsd	160(%rsp), %xmm2
	vsubsd	%xmm5, %xmm6, %xmm5
	vfmadd231sd	%xmm2, %xmm4, %xmm1
	vmovsd	48(%r10,%rax), %xmm4
	vmovsd	%xmm5, (%r11)
	vmulsd	%xmm15, %xmm4, %xmm9
	vmulsd	%xmm11, %xmm4, %xmm10
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	48(%r8,%rax), %xmm1
	leaq	16(%rax), %r11
	vaddsd	%xmm7, %xmm0, %xmm0
	vmulsd	%xmm12, %xmm1, %xmm8
	leaq	(%rcx,%r11), %rsi
	vmovsd	%xmm0, (%r9)
	vmovsd	16(%r10,%rax), %xmm0
	vmovsd	(%rsi), %xmm3
	vfmadd231sd	%xmm2, %xmm0, %xmm9
	vmovsd	16(%r8,%rax), %xmm2
	vfmadd231sd	%xmm14, %xmm0, %xmm10
	vfmadd231sd	%xmm13, %xmm2, %xmm8
	vmulsd	%xmm13, %xmm1, %xmm5
	leaq	48(%rax), %r9
	leaq	(%rcx,%r9), %rbx
	addq	%rdx, %r11
	vaddsd	%xmm10, %xmm8, %xmm8
	vfmadd231sd	%xmm14, %xmm4, %xmm5
	vmovsd	%xmm14, %xmm14, %xmm10
	vsubsd	%xmm8, %xmm3, %xmm8
	vmulsd	%xmm12, %xmm2, %xmm3
	vmovsd	272(%rsp), %xmm14
	addq	%rdx, %r9
	vmovsd	(%r11), %xmm6
	vmovsd	(%r9), %xmm7
	vfmadd231sd	%xmm11, %xmm0, %xmm3
	vsubsd	%xmm5, %xmm3, %xmm3
	vmulsd	168(%rsp), %xmm1, %xmm5
	vmulsd	%xmm14, %xmm1, %xmm1
	vaddsd	(%rbx), %xmm3, %xmm3
	vmovsd	%xmm8, (%rsi)
	vfmadd231sd	%xmm14, %xmm2, %xmm5
	vmulsd	168(%rsp), %xmm2, %xmm2
	vfmadd231sd	160(%rsp), %xmm4, %xmm1
	vmovsd	%xmm3, (%rbx)
	vmovsd	56(%r8,%rax), %xmm3
	vaddsd	%xmm9, %xmm5, %xmm5
	vfmadd132sd	%xmm15, %xmm2, %xmm0
	vmovsd	56(%r10,%rax), %xmm2
	vsubsd	%xmm5, %xmm6, %xmm5
	vmulsd	%xmm11, %xmm2, %xmm6
	vmovsd	24(%r8,%rax), %xmm4
	vsubsd	%xmm1, %xmm0, %xmm0
	vmovsd	%xmm5, (%r11)
	vmulsd	%xmm15, %xmm2, %xmm9
	vaddsd	%xmm7, %xmm0, %xmm0
	vmulsd	%xmm12, %xmm3, %xmm7
	leaq	56(%rax), %r11
	vmovsd	%xmm0, (%r9)
	vmovsd	24(%r10,%rax), %xmm0
	leaq	24(%rax), %r9
	vfmadd231sd	%xmm13, %xmm4, %xmm7
	vfmadd231sd	%xmm10, %xmm0, %xmm6
	leaq	(%rcx,%r9), %r10
	vmovsd	(%r10), %xmm1
	addq	%rdx, %r9
	addq	%r11, %rcx
	vaddsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm13, %xmm3, %xmm7
	addq	%r11, %rdx
	vsubsd	%xmm6, %xmm1, %xmm6
	vmulsd	%xmm12, %xmm4, %xmm1
	vfmadd231sd	160(%rsp), %xmm0, %xmm9
	vfmadd231sd	%xmm10, %xmm2, %xmm7
	vmovsd	(%r9), %xmm5
	vmovsd	(%rdx), %xmm8
	vfmadd231sd	%xmm11, %xmm0, %xmm1
	vmovsd	168(%rsp), %xmm11
	vsubsd	%xmm7, %xmm1, %xmm1
	vmulsd	%xmm11, %xmm3, %xmm7
	vmulsd	%xmm14, %xmm3, %xmm3
	vaddsd	(%rcx), %xmm1, %xmm1
	vfmadd231sd	%xmm14, %xmm4, %xmm7
	vmulsd	%xmm11, %xmm4, %xmm4
	vfmadd132sd	160(%rsp), %xmm3, %xmm2
	vmovsd	%xmm6, (%r10)
	vmovsd	%xmm1, (%rcx)
	vaddsd	%xmm9, %xmm7, %xmm7
	vfmadd132sd	%xmm15, %xmm4, %xmm0
	vsubsd	%xmm7, %xmm5, %xmm5
	vsubsd	%xmm2, %xmm0, %xmm0
	vmovsd	%xmm5, (%r9)
	vaddsd	%xmm8, %xmm0, %xmm0
	vmovsd	%xmm0, (%rdx)
.L22:
	vmovaps	1024(%rsp), %xmm6
	vmovaps	1040(%rsp), %xmm7
	vmovaps	1056(%rsp), %xmm8
	vmovaps	1072(%rsp), %xmm9
	vmovaps	1088(%rsp), %xmm10
	vmovaps	1104(%rsp), %xmm11
	vmovaps	1120(%rsp), %xmm12
	vmovaps	1136(%rsp), %xmm13
	vmovaps	1152(%rsp), %xmm14
	vmovaps	1168(%rsp), %xmm15
	addq	$1192, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	popq	%rbp
	popq	%r12
	popq	%r13
	ret
.L20:
	vzeroupper
	jmp	.L22
.L3:
	shrl	$3, %r9d
	movl	%r9d, %eax
	salq	$6, %rax
	vmovsd	272(%rsp), %xmm11
	vmovsd	1000(%rsp), %xmm12
	vmovsd	1008(%rsp), %xmm13
	vmovsd	1016(%rsp), %xmm14
	vmovsd	168(%rsp), %xmm15
	leaq	64(%r10,%rax), %rax
	.p2align 4,,10
	.p2align 3
.L9:
	vmovsd	(%r8), %xmm4
	vmovsd	32(%r8), %xmm3
	vmovsd	992(%rsp), %xmm10
	vmulsd	%xmm3, %xmm13, %xmm6
	vmulsd	%xmm10, %xmm4, %xmm7
	vmovsd	(%r10), %xmm0
	vmovsd	32(%r10), %xmm2
	vmovsd	280(%rsp), %xmm9
	vfmadd231sd	%xmm2, %xmm12, %xmm6
	vfmadd231sd	%xmm9, %xmm0, %xmm7
	vmovsd	(%rcx), %xmm1
	vmovsd	32(%rdx), %xmm8
	vmovsd	(%rdx), %xmm5
	addq	$64, %r10
	vaddsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm10, %xmm3, %xmm7
	vmovsd	160(%rsp), %xmm10
	vsubsd	%xmm6, %xmm1, %xmm6
	vmulsd	%xmm4, %xmm13, %xmm1
	addq	$64, %rcx
	vfmadd231sd	%xmm9, %xmm2, %xmm7
	vmulsd	%xmm4, %xmm11, %xmm9
	vmulsd	%xmm4, %xmm15, %xmm4
	vfmadd231sd	%xmm0, %xmm12, %xmm1
	vmovsd	%xmm6, -64(%rcx)
	addq	$64, %rdx
	vfmadd231sd	%xmm10, %xmm0, %xmm9
	vfmadd132sd	%xmm14, %xmm4, %xmm0
	vsubsd	%xmm7, %xmm1, %xmm1
	vmulsd	%xmm3, %xmm15, %xmm7
	vmulsd	%xmm3, %xmm11, %xmm3
	vaddsd	-32(%rcx), %xmm1, %xmm1
	vmovsd	%xmm10, %xmm10, %xmm4
	addq	$64, %r8
	vfmadd231sd	%xmm2, %xmm14, %xmm7
	vfmadd132sd	%xmm10, %xmm3, %xmm2
	vmovsd	-24(%r8), %xmm3
	vmovsd	%xmm1, -32(%rcx)
	vmulsd	%xmm3, %xmm13, %xmm6
	vaddsd	%xmm9, %xmm7, %xmm7
	vsubsd	%xmm2, %xmm0, %xmm0
	vmovsd	-24(%r10), %xmm2
	vsubsd	%xmm7, %xmm5, %xmm5
	vaddsd	%xmm8, %xmm0, %xmm0
	vmulsd	%xmm2, %xmm14, %xmm9
	vmulsd	%xmm2, %xmm12, %xmm10
	vmovsd	%xmm0, -32(%rdx)
	vmovsd	-56(%r10), %xmm0
	vmovsd	992(%rsp), %xmm7
	vfmadd231sd	%xmm4, %xmm0, %xmm9
	vmovsd	-56(%r8), %xmm4
	vfmadd231sd	280(%rsp), %xmm0, %xmm10
	vfmadd231sd	%xmm7, %xmm4, %xmm6
	vmovsd	%xmm5, -64(%rdx)
	vmovsd	-56(%rcx), %xmm1
	vmulsd	%xmm7, %xmm3, %xmm7
	vmovsd	-24(%rdx), %xmm8
	vaddsd	%xmm10, %xmm6, %xmm6
	vmovsd	-56(%rdx), %xmm5
	vsubsd	%xmm6, %xmm1, %xmm6
	vmulsd	%xmm4, %xmm13, %xmm1
	vfmadd231sd	280(%rsp), %xmm2, %xmm7
	vmovsd	%xmm6, -56(%rcx)
	vfmadd231sd	%xmm0, %xmm12, %xmm1
	vsubsd	%xmm7, %xmm1, %xmm1
	vmulsd	%xmm3, %xmm15, %xmm7
	vmulsd	%xmm3, %xmm11, %xmm3
	vaddsd	-24(%rcx), %xmm1, %xmm1
	vfmadd231sd	%xmm4, %xmm11, %xmm7
	vmulsd	%xmm4, %xmm15, %xmm4
	vmovsd	%xmm1, -24(%rcx)
	vaddsd	%xmm9, %xmm7, %xmm7
	vfmadd132sd	%xmm14, %xmm4, %xmm0
	vmovsd	160(%rsp), %xmm4
	vsubsd	%xmm7, %xmm5, %xmm5
	vfmadd132sd	%xmm4, %xmm3, %xmm2
	vmovsd	%xmm5, -56(%rdx)
	vsubsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm8, %xmm0, %xmm0
	vmovsd	%xmm0, -24(%rdx)
	vmovsd	-48(%r10), %xmm0
	vmovsd	-16(%r10), %xmm2
	vmovsd	-16(%r8), %xmm3
	vmulsd	%xmm2, %xmm14, %xmm9
	vmulsd	%xmm2, %xmm12, %xmm10
	vmulsd	%xmm3, %xmm13, %xmm6
	vmovsd	992(%rsp), %xmm7
	vmovsd	-48(%rcx), %xmm1
	vfmadd231sd	%xmm4, %xmm0, %xmm9
	vmovsd	-48(%r8), %xmm4
	vfmadd231sd	280(%rsp), %xmm0, %xmm10
	vfmadd231sd	%xmm7, %xmm4, %xmm6
	vmovsd	-16(%rdx), %xmm8
	vmovsd	-48(%rdx), %xmm5
	vaddsd	%xmm10, %xmm6, %xmm6
	vmovsd	%xmm7, %xmm7, %xmm10
	vmulsd	%xmm7, %xmm3, %xmm7
	vsubsd	%xmm6, %xmm1, %xmm6
	vmulsd	%xmm4, %xmm13, %xmm1
	vmovsd	%xmm6, -48(%rcx)
	vfmadd231sd	280(%rsp), %xmm2, %xmm7
	vfmadd231sd	%xmm0, %xmm12, %xmm1
	vsubsd	%xmm7, %xmm1, %xmm1
	vmulsd	%xmm3, %xmm15, %xmm7
	vmulsd	%xmm3, %xmm11, %xmm3
	vaddsd	-16(%rcx), %xmm1, %xmm1
	vfmadd231sd	%xmm4, %xmm11, %xmm7
	vmulsd	%xmm4, %xmm15, %xmm4
	vmovsd	%xmm1, -16(%rcx)
	vaddsd	%xmm9, %xmm7, %xmm7
	vfmadd132sd	%xmm14, %xmm4, %xmm0
	vmovsd	160(%rsp), %xmm4
	vsubsd	%xmm7, %xmm5, %xmm5
	vfmadd132sd	%xmm4, %xmm3, %xmm2
	vmovsd	-8(%r8), %xmm3
	vmovsd	%xmm5, -48(%rdx)
	vmulsd	%xmm13, %xmm3, %xmm9
	vmovsd	-40(%rdx), %xmm5
	vsubsd	%xmm2, %xmm0, %xmm0
	vmovsd	-8(%r10), %xmm2
	vmovsd	-8(%rdx), %xmm7
	vaddsd	%xmm8, %xmm0, %xmm0
	vmulsd	%xmm2, %xmm14, %xmm8
	vmulsd	%xmm2, %xmm12, %xmm6
	vmovsd	%xmm0, -16(%rdx)
	vmovsd	-40(%r10), %xmm0
	vmovsd	-40(%rcx), %xmm1
	vfmadd231sd	%xmm4, %xmm0, %xmm8
	vmovsd	-40(%r8), %xmm4
	vfmadd231sd	280(%rsp), %xmm0, %xmm6
	vfmadd231sd	%xmm10, %xmm4, %xmm9
	vaddsd	%xmm9, %xmm6, %xmm6
	vmulsd	%xmm10, %xmm3, %xmm9
	vsubsd	%xmm6, %xmm1, %xmm6
	vmulsd	%xmm13, %xmm4, %xmm1
	vfmadd231sd	280(%rsp), %xmm2, %xmm9
	vmovsd	%xmm6, -40(%rcx)
	vfmadd231sd	%xmm0, %xmm12, %xmm1
	vsubsd	%xmm9, %xmm1, %xmm1
	vmulsd	%xmm15, %xmm3, %xmm9
	vmulsd	%xmm11, %xmm3, %xmm3
	vaddsd	-8(%rcx), %xmm1, %xmm1
	vfmadd231sd	%xmm11, %xmm4, %xmm9
	vmulsd	%xmm15, %xmm4, %xmm4
	vfmadd132sd	160(%rsp), %xmm3, %xmm2
	vmovsd	%xmm1, -8(%rcx)
	vaddsd	%xmm9, %xmm8, %xmm8
	vfmadd132sd	%xmm14, %xmm4, %xmm0
	vsubsd	%xmm8, %xmm5, %xmm5
	vsubsd	%xmm2, %xmm0, %xmm0
	vmovsd	%xmm5, -40(%rdx)
	vaddsd	%xmm7, %xmm0, %xmm0
	vmovsd	%xmm0, -8(%rdx)
	cmpq	%rax, %r10
	jne	.L9
	jmp	.L22
.L11:
	xorl	%eax, %eax
	xorl	%edi, %edi
	jmp	.L4
.L24:
	vzeroupper
	jmp	.L7
	.seh_endproc
	.ident	"GCC: (Rev3, Built by MSYS2 project) 10.2.0"
