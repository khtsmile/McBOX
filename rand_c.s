# mark_description "Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 14.0.3.174 Build 2014042";
# mark_description "2";
# mark_description "-openmp -pg -save-temps -c";
	.file "rand_c.c"
	.text
..TXTST0:
# -- Begin  rang
# mark_begin;
       .align    16,0x90
	.globl rang
rang:
..B1.1:                         # Preds ..B1.0
..___tag_value_rang.1:                                          #8.15
        pushq     %rbp                                          #8.15
..___tag_value_rang.3:                                          #
        movq      %rsp, %rbp                                    #8.15
..___tag_value_rang.4:                                          #
        subq      $16, %rsp                                     #8.15
        lea       _gprof_pack0(%rip), %rdx                      #8.15
        movq      %r15, -16(%rbp)                               #8.15
        movq      %r14, -8(%rbp)                                #8.15
..___tag_value_rang.6:                                          #8.15
        call      mcount                                        #8.15
..___tag_value_rang.7:                                          #
                                # LOE rbx r12 r13
..B1.12:                        # Preds ..B1.1
        movl      $.2.3_2_kmpc_loc_struct_pack.12, %edi         #8.15
        call      __kmpc_global_thread_num                      #8.15
                                # LOE rbx r12 r13 eax
..B1.11:                        # Preds ..B1.12
        movq      .2.3_2_kmpv_seed_V$2_cache_0.27(%rip), %rdx   #8.15
        movl      %eax, %esi                                    #8.15
        testq     %rdx, %rdx                                    #8.15
        je        ..B1.3        # Prob 50%                      #8.15
                                # LOE rdx rbx r12 r13 esi
..B1.2:                         # Preds ..B1.11
        lea       (,%rsi,8), %eax                               #8.15
        cltq                                                    #8.15
        movq      (%rax,%rdx), %r15                             #8.15
        testq     %r15, %r15                                    #8.15
        jne       ..B1.5        # Prob 50%                      #8.15
                                # LOE rbx r12 r13 r15 esi
..B1.3:                         # Preds ..B1.2 ..B1.11
        movl      $.2.3_2_kmpc_loc_struct_pack.20, %edi         #8.15
        movl      $seed, %edx                                   #8.15
        movl      $.2.3_2_kmpv_seed_V$2_cache_0.27, %r8d        #8.15
        xorl      %eax, %eax                                    #8.15
        movq      .2.3_2_kmpv_seed_V$2_size_0.28(%rip), %rcx    #8.15
..___tag_value_rang.10:                                         #8.15
        call      __kmpc_threadprivate_cached                   #8.15
..___tag_value_rang.11:                                         #
                                # LOE rax rbx r12 r13
..B1.13:                        # Preds ..B1.3
        movq      %rax, %r15                                    #8.15
                                # LOE rbx r12 r13 r15
..B1.5:                         # Preds ..B1.13 ..B1.2
        xorl      %edi, %edi                                    #11.20
        call      time                                          #11.20
                                # LOE rax rbx r12 r13 r15
..B1.6:                         # Preds ..B1.5
        imull     $125346, %eax, %r14d                          #13.26
        addl      (%r15), %r14d                                 #13.26
        movl      %r14d, (%r15)                                 #13.2
        call      omp_get_thread_num                            #14.21
                                # LOE rbx r12 r13 r15 eax r14d
..B1.7:                         # Preds ..B1.6
        imull     $125346, %eax, %eax                           #14.21
        movq      %r15, %rdi                                    #20.21
        addl      %eax, %r14d                                   #14.21
        movl      %r14d, (%r15)                                 #14.2
..___tag_value_rang.12:                                         #20.21
        call      rand_r                                        #20.21
..___tag_value_rang.13:                                         #
                                # LOE rbx r12 r13 eax
..B1.8:                         # Preds ..B1.7
        cvtsi2sd  %eax, %xmm0                                   #20.21
        divsd     .L_2il0floatpacket.29(%rip), %xmm0            #20.37
        movq      -8(%rbp), %r14                                #21.9
..___tag_value_rang.14:                                         #
        movq      -16(%rbp), %r15                               #21.9
..___tag_value_rang.15:                                         #
        movq      %rbp, %rsp                                    #21.9
        popq      %rbp                                          #21.9
..___tag_value_rang.16:                                         #
        ret                                                     #21.9
        .align    16,0x90
..___tag_value_rang.17:                                         #
                                # LOE
# mark_end;
	.type	rang,@function
	.size	rang,.-rang
	.bss
	.align 8
.2.3_2_kmpv_seed_V$2_cache_0.27:
	.type	.2.3_2_kmpv_seed_V$2_cache_0.27,@object
	.size	.2.3_2_kmpv_seed_V$2_cache_0.27,8
	.space 8	# pad
	.data
	.align 8
.2.3_2_kmpv_seed_V$2_size_0.28:
	.long	0x00000004,0x00000000
	.align 4
.2.3_2_kmpc_loc_struct_pack.12:
	.long	0
	.long	2
	.long	0
	.long	0
	.quad	.2.3_2__kmpc_loc_pack.11
	.align 4
.2.3_2__kmpc_loc_pack.11:
	.byte	59
	.byte	117
	.byte	110
	.byte	107
	.byte	110
	.byte	111
	.byte	119
	.byte	110
	.byte	59
	.byte	114
	.byte	97
	.byte	110
	.byte	103
	.byte	59
	.byte	56
	.byte	59
	.byte	56
	.byte	59
	.byte	59
	.space 1, 0x00 	# pad
	.align 4
.2.3_2_kmpc_loc_struct_pack.20:
	.long	0
	.long	2
	.long	0
	.long	0
	.quad	.2.3_2__kmpc_loc_pack.19
	.align 4
.2.3_2__kmpc_loc_pack.19:
	.byte	59
	.byte	117
	.byte	110
	.byte	107
	.byte	110
	.byte	111
	.byte	119
	.byte	110
	.byte	59
	.byte	114
	.byte	97
	.byte	110
	.byte	103
	.byte	59
	.byte	56
	.byte	59
	.byte	50
	.byte	49
	.byte	59
	.byte	59
	.data
# -- End  rang
	.data
	.align 4
	.globl seed
seed:
	.long	12357425
	.type	seed,@object
	.size	seed,4
	.section .rodata, "a"
	.align 8
.L_2il0floatpacket.29:
	.long	0xffc00000,0x41dfffff
	.type	.L_2il0floatpacket.29,@object
	.size	.L_2il0floatpacket.29,8
	.align 4
_gprof_pack0:
	.long	0
	.type	_gprof_pack0,@object
	.size	_gprof_pack0,4
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
	.4byte 0x00000014
	.8byte 0x7801000100000000
	.8byte 0x0000019008070c10
	.4byte 0x00000000
	.4byte 0x00000044
	.4byte 0x0000001c
	.8byte ..___tag_value_rang.1
	.8byte ..___tag_value_rang.17-..___tag_value_rang.1
	.byte 0x04
	.4byte ..___tag_value_rang.3-..___tag_value_rang.1
	.2byte 0x100e
	.byte 0x04
	.4byte ..___tag_value_rang.4-..___tag_value_rang.3
	.4byte 0x8610060c
	.2byte 0x0402
	.4byte ..___tag_value_rang.7-..___tag_value_rang.4
	.4byte 0x048f038e
	.byte 0x04
	.4byte ..___tag_value_rang.14-..___tag_value_rang.7
	.2byte 0x04ce
	.4byte ..___tag_value_rang.15-..___tag_value_rang.14
	.2byte 0x04cf
	.4byte ..___tag_value_rang.16-..___tag_value_rang.15
	.4byte 0x000000c6
	.byte 0x00
# End
