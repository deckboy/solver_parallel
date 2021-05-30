version: 1
creator: callgrind-3.8.1
pid: 26349
cmd:  ./../FTSim.exe
part: 1


desc: I1 cache: 
desc: D1 cache: 
desc: LL cache: 

desc: Timerange: Basic block 0 - 29514914869
desc: Trigger: Program termination

positions: line
events: Ir
summary: 252069463065


ob=(8) /usr/lib64/valgrind/vgpreload_core-amd64-linux.so
fl=(10) ???
fn=(334) 0x00000000000006d0
0 7
0 4

fn=(328) 0x0000000000000530
0 4
0 2

fn=(1104) 0x00000000000004b0
0 10
cob=(3) /lib64/libc-2.12.so
cfi=(3) ???
cfn=(1116) __cxa_finalize
calls=1 0 
0 68
cob=(1) /lib64/ld-2.12.so
cfi=(1) ???
cfn=(186) _dl_runtime_resolve
calls=1 0 
0 1026
0 8
0 5

fn=(322) 0x0000000000000490
0 4
0 2

ob=(10) /usr/lib64/valgrind/callgrind-amd64-linux
fl=(12) /builddir/build/BUILD/valgrind-3.8.1/coregrind/m_trampoline.S
fn=(426) vgPlain_amd64_linux_REDIR_FOR_vgettimeofday
176 213
+1 213
+1 213

ob=(2) ???
fl=(2) ???
fn=(256) 0x0000003001400cc0
0 2
cob=(5) /lib64/libdl-2.12.so
cfi=(5) ???
cfn=(258) 0x0000003001400de0
calls=1 0 
0 6
0 1
cob=(5)
cfi=(5)
cfn=(264) 0x0000003001400e80
calls=1 0 
0 6
0 1
cob=(5)
cfi=(5)
cfn=(270) 0x0000003001401960
calls=1 0 
0 11
0 2

fn=(358) 0x0000000000402478
0 2
cob=(9) /general/home/crriedel/csce657_project/TOUGH_HCGas_xfer/FTSim.exe
cfi=(11) ???
cfn=(360) 0x0000000000402bec
calls=1 0 
0 6
0 2

fn=(1146) 0x00000030010111f0
0 2
cob=(4) /lib64/libpthread-2.12.so
cfi=(4) ???
cfn=(1148) 0x0000003001005790
calls=1 0 
0 1463
0 2

fn=(174) 0x0000003001005300
0 2
cob=(4)
cfi=(4)
cfn=(176) __pthread_initialize_minimal
calls=1 0 
0 9448
0 1
cob=(4)
cfi=(4)
cfn=(224) 0x0000003001005810
calls=1 0 
0 6
0 1
cob=(4)
cfi=(4)
cfn=(230) 0x0000003001011160
calls=1 0 
0 11
0 2

fn=(1184) 0x0000000004a2a7e8
0 3

fn=(300) 0x0000003000c03d90
0 2
cob=(7) /lib64/libm-2.12.so
cfi=(9) ???
cfn=(302) 0x0000003000c03e70
calls=1 0 
0 6
0 1
cob=(7)
cfi=(9)
cfn=(308) 0x0000003000c03f10
calls=1 0 
0 6
0 1
cob=(7)
cfi=(9)
cfn=(314) 0x0000003000c43f80
calls=1 0 
0 11
0 2

fn=(278) 0x0000000004a1a748
0 2
cob=(6) /general/software/x86_64/easybuild/software/GCCcore/5.4.0/lib64/libgcc_s.so.1
cfi=(6) ???
cfn=(280) 0x0000000000002f1c
calls=1 0 
0 6
0 2

fn=(320) 0x0000000004801448
0 2
cob=(8)
cfi=(10)
cfn=(322)
calls=1 0 
0 6
0 1
cob=(8)
cfi=(10)
cfn=(328)
calls=1 0 
0 6
0 1
cob=(8)
cfi=(10)
cfn=(334)
calls=1 0 
0 11
0 2

fn=(1102) 0x0000000004801708
0 2
cob=(8)
cfi=(10)
cfn=(1104)
calls=1 0 
0 1117
0 2

fn=(1186) 0x0000003001401998
0 2
cob=(5)
cfi=(5)
cfn=(1188) 0x0000003001400e00
calls=1 0 
0 1476
0 2

fn=(1100) 0x0000000000523c3c
0 3

fn=(1126) 0x0000003000c43fb8
0 2
cob=(7)
cfi=(9)
cfn=(1128) 0x0000003000c03e90
calls=1 0 
0 1463
0 2

ob=(6)
fl=(8) /tmp/easybuild/GCCcore/5.4.0/dummy-/gcc-5.4.0/stage3_obj/gcc/include/cpuid.h
fn=(288) __get_cpuid_output
219 2
fi=(7) /tmp/easybuild/GCCcore/5.4.0/dummy-/gcc-5.4.0/stage3_obj/x86_64-unknown-linux-gnu/libgcc/../../../libgcc/config/i386/cpuinfo.c
336 4
fe=(8)
219 2
fi=(7)
336 2
fe=(8)
219 8
+20 4
+3 8
+1 2
-1 4
fi=(7)
+96 4
fe=(8)

fl=(6)
fn=(1178) 0x0000000000002f40
0 7
0 2

fn=(280)
0 4
0 2

fn=(292) 0x0000000000003010
0 11
0 6

fn=(1164) 0x0000000000002fd0
0 8
cob=(3)
cfi=(3)
cfn=(1116)
calls=1 0 
0 68
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1372
0 1
cfn=(1178)
calls=1 0 
0 9
0 3

fl=(7)
fn=(290) get_available_features
267 1
-1 1
+1 2
+5 4
+2 4
+2 4
+2 4
+2 4
+2 4
+2 4
+2 4
+2 4
+2 4
+3 2
+3 3
+2 4
+2 4
+2 4
+2 4
+6 2
+2 2
+2 2
+3 4
+2 4
+2 4
+3 2
+1 2

fn=(286) __cpu_indicator_init
349 6
+10 1
-10 1
+9 4
jcnd=1/1 +4 
* 
+71 9

ob=(5)
fl=(5)
fn=(264)
0 4
0 2

fn=(1202) fini
0 4
0 2

fn=(258)
0 4
0 2

fn=(276) init
0 6
cob=(4)
cfi=(4)
cfn=(554) pthread_key_create
calls=1 0 
0 20
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1233
0 2
0 5

fn=(540) _dlerror_run
0 121
cob=(4)
cfi=(4)
cfn=(548) pthread_once
calls=11 0 
0 1332
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1100
0 33
0 33
0 77
cob=(1)
cfi=(1)
cfn=(70) _dl_catch_error
calls=11 0 
0 11639
0 187
cob=(4)
cfi=(4)
cfn=(560) pthread_getspecific
calls=11 0 
0 192
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1221
0 33
0 3
cob=(3)
cfi=(3)
cfn=(566) calloc
calls=1 0 
0 248
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1172
0 5

ob=(6)
fl=(7)
fn=(286)
362 11
cfi=(8)
cfn=(288)
calls=1 219 
* 20
* 3
+7 1
-1 1
+3 2
+6 6
cfi=(8)
cfn=(288)
calls=1 219 
* 20
* 2
+6 2
+1 1
+2 1
-3 1
+1 1
+2 1
+1 1
-4 1
+1 1
+2 1
+3 1
-2 1
+2 1
jcnd=1/1 +3 
* 
+3 1
-7 1
+7 1

ob=(5)
fl=(5)
fn=(540)
0 3
cob=(4)
cfi=(4)
cfn=(572) pthread_setspecific
calls=1 0 
0 31
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1139
0 1

fn=(574) dlsym_doit
0 66
cob=(3)
cfi=(3)
cfn=(580) _dl_sym
calls=11 0 
0 9760
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1142
0 33

fn=(538) dlsym
0 154
cob=(4)
cfi=(4)
cfn=(458) pthread_mutex_lock
calls=11 0 
0 253
0 33
cfn=(540)
calls=11 0 
0 18570
0 44
cob=(4)
cfi=(4)
cfn=(460) pthread_mutex_unlock
calls=11 0 
0 286
0 55

ob=(6)
fl=(7)
fn=(286)
397 2
jcnd=1/1 +1 
* 
186 7
jump=1 +33 
* 
403 3
+25 1
-25 1
cfn=(290)
calls=1 267 
* 87

ob=(5)
fl=(5)
fn=(1188)
0 10
cob=(3)
cfi=(3)
cfn=(1116)
calls=1 0 
0 68
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1372
0 8

ob=(6)
fl=(7)
fn=(286)
404 3
jump=1 +25 
* 

ob=(5)
fl=(5)
fn=(1188)
0 4
cfn=(1202)
calls=1 0 
0 6
0 8

ob=(6)
fl=(7)
fn=(286)
398 1
178 2
jcnd=1/1 +8 
* 
+41 2
+1 2
jump=1 403 
* 

ob=(5)
fl=(5)
fn=(270)
0 7
0 4

ob=(4)
fl=(4)
fn=(560)
0 160
0 116

fn=(1148)
0 10
cob=(3)
cfi=(3)
cfn=(1116)
calls=1 0 
0 68
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1372
0 8
0 5

fn=(760) write
0 7575476
cfn=(762) __write_nocancel
calls=3787738 0 
0 18938690

fn=(802) __read_nocancel
0 16505

fn=(896) pthread_mutex_destroy
0 31349007

fn=(178) __libc_sigaction
0 558
0 94

fn=(496) sigaction
0 28
cfn=(178)
calls=7 0 
0 524

fn=(548)
0 22
0 33
0 3
cob=(5)
cfi=(5)
cfn=(276)
calls=1 0 
0 1263
0 11

fn=(754) __errno_location
0 11373105

fn=(460)
0 65501352
0 57313748
0 78

fn=(690) open
0 20
cfn=(692) __open_nocancel
calls=10 0 
0 50

fn=(1070) __close_nocancel
0 50

fn=(176)
0 74
0 16
0 2
cfn=(178)
calls=1 0 
0 64
0 7
cfn=(178)
calls=1 0 
0 64
0 11
cob=(1)
cfi=(1)
cfn=(190) _dl_get_tls_static_info
calls=1 0 
0 5
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 2170
0 4
0 11
cob=(3)
cfi=(3)
cfn=(196) __secure_getenv
calls=1 0 
0 691
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1352
0 6
cob=(3)
cfi=(3)
cfn=(206) getrlimit
calls=1 0 
0 5
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1341
0 2
0 14
cob=(3)
cfi=(3)
cfn=(212) __libc_dl_error_tsd
calls=1 0 
0 3
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1589
0 2
cob=(1)
cfi=(1)
cfn=(72) _dl_initial_error_catch_tsd
calls=1 0 
0 2
0 12
0 10
cob=(3)
cfi=(3)
cfn=(218) __libc_pthread_init
calls=1 0 
0 453
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1523
0 15

fn=(230)
0 7
0 4

fn=(732) lseek
0 408
cfn=(734) __lseek_nocancel
calls=204 0 
0 1020

fn=(554)
0 8
0 13
0 26

fn=(762)
0 18938690

fn=(224)
0 4
0 2

fn=(656) pthread_mutex_init
0 31349061
0 38315519
0 31349061
0 10449687

fn=(800) read
0 6602
cfn=(802)
calls=3301 0 
0 16505

fn=(458)
0 81876690
0 73689073

ob=(7)
fl=(9)
fn=(302)
0 4
0 2

fn=(1128)
0 6

ob=(4)
fl=(4)
fn=(572)
0 62

fn=(692)
0 40

ob=(7)
fl=(9)
fn=(1128)
0 4
cob=(3)
cfi=(3)
cfn=(1116)
calls=1 0 
0 68
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1372

ob=(4)
fl=(4)
fn=(692)
0 10

ob=(7)
fl=(9)
fn=(1128)
0 8
0 5

fn=(314)
0 7
0 4

fn=(308)
0 4

ob=(4)
fl=(4)
fn=(666) pthread_self
0 16375216

fn=(734)
0 408

ob=(7)
fl=(9)
fn=(308)
0 2

ob=(4)
fl=(4)
fn=(734)
0 612

fn=(1068) close
0 20
cfn=(1070)
calls=10 0 
0 50

ob=(1)
fl=(1)
fn=(12) sbrk
0 12
0 5

fn=(28) malloc
0 159
cfn=(32) __libc_memalign
calls=53 0 
0 1431

fn=(40) rtld_lock_default_unlock_recursive
0 14

fn=(62) _dl_debug_initialize
0 48
0 32
0 6

fn=(98) close
0 35

fn=(132) _dl_check_all_versions
0 14
0 14
0 40
0 39
cfn=(134) _dl_check_map_versions
calls=8 0 
0 12478
0 32

fn=(138) init_tls
0 6
0 8
cfn=(24) calloc
calls=1 0 
0 37
0 9
0 14
0 7
0 21
0 1
cfn=(140) _dl_determine_tlsoffset
calls=1 0 
0 55
0 1
cfn=(142) _dl_allocate_tls_storage
calls=1 0 
0 420
0 14

fn=(142)
0 7
cfn=(32)
calls=1 0 
0 50
0 18
0 288
0 10
cfn=(24)
calls=1 0 
0 37
0 10

fn=(148) do_lookup_x
0 7164
0 7774
0 30
0 17119
0 11850
0 4718
0 1980
0 1308
0 1524
0 2349
0 640
0 1752
0 1665
cfn=(80) _dl_name_match_p
calls=333 0 
0 15997
0 1665
0 868
0 2170
0 3000
cfn=(150) check_match.12447
calls=375 0 
0 10947
0 1875
0 2562
0 434
0 832
0 1656
0 1407
0 358
0 930
0 1005
0 2685
cfn=(150)
calls=179 0 
0 33740
0 1074
0 4
0 36
0 96
0 708
0 5735
0 1218
0 4467
0 14890
0 522
0 22

fn=(172) _dl_init
0 16
0 2
0 16
0 56
0 42
0 18
0 50
cob=(2)
cfi=(2)
cfn=(320)
calls=1 0 
0 29
cob=(2)
cfi=(2)
cfn=(300)
calls=1 0 
0 29
cob=(2)
cfi=(2)
cfn=(278)
calls=1 0 
0 10
cob=(2)
cfi=(2)
cfn=(256)
calls=1 0 
0 29
cob=(3)
cfi=(3)
cfn=(236) _init
calls=1 0 
0 2152
0 15
0 30
cob=(6)
cfi=(7)
cfn=(286)
calls=1 349 
0 215
cob=(5)
cfi=(5)
cfn=(276)
calls=1 0 
0 3
0 5
cob=(6)
cfi=(6)
cfn=(292)
calls=1 0 
0 17
0 15
0 4
0 10
0 4
0 23
cob=(2)
cfi=(2)
cfn=(174)
calls=1 0 
0 9471
0 4
0 2

fn=(10) strlen
0 184
0 116
0 489
0 315
0 923
0 615
0 476
0 525
0 376
0 445
0 352
0 20
0 120
0 114
0 72
0 54
0 60
0 30
0 21
0 93

fn=(24)
0 144
cfn=(28)
calls=16 0 
0 529

fn=(84) open_verify
0 1955
0 228
0 575
cfn=(86) open
calls=115 0 
0 1229
0 345
0 48
cfn=(88) read
calls=6 0 
0 30
0 66
0 46
0 6
0 42
0 42
0 156
0 256
0 72
0 65
0 50
0 24
0 20
0 84
0 9
0 80
0 1037
0 12
cfn=(120) bcmp
calls=2 0 
0 178
0 8

fn=(86)
0 464
0 770

fn=(96) memset
0 63
0 126
0 108
0 139
0 4
0 81
0 912
0 21

fn=(106) open_path
0 661
cfn=(52) mempcpy
calls=68 0 
0 4269
0 612
0 112
0 5
0 265
0 872
0 612
0 136
0 660
cfn=(52)
calls=110 0 
0 1946
0 440
cfn=(52)
calls=110 0 
0 2634
0 1210
cfn=(84)
calls=110 0 
0 5297
0 440
0 324
0 53
0 871
0 126
0 1008
cfn=(108) _xstat
calls=56 0 
0 686
0 112
0 56
0 210
0 42
0 3
cfn=(28)
calls=1 0 
0 28
0 9
cfn=(34) memcpy
calls=1 0 
0 64
0 67
0 42

fn=(140)
0 13
0 23
0 19

fn=(146) _dl_lookup_symbol_x
0 5970
0 21528
0 1791
0 2355
0 91
0 2985
cfn=(148)
calls=199 0 
0 170776
0 398
0 597
0 558
0 2374
0 2799
0 39
0 66
0 55
0 2

fn=(162) _dl_allocate_tls_init
0 32
0 15
cfn=(52)
calls=1 0 
0 23
0 4
cfn=(96)
calls=1 0 
0 27
0 10
0 9
0 10

fn=(4) _dl_setup_hash
0 24
0 148

fn=(34)
0 93
0 60
0 36
0 48
0 70
0 40
0 18
0 54
0 84
0 32
0 225
0 9
0 19
0 75
0 28
0 186
0 4

fn=(36) _dl_add_to_namespace_list
0 49
cfn=(38) rtld_lock_default_lock_recursive
calls=7 0 
0 14
0 42
0 30
0 84
0 84
cfn=(40)
calls=6 0 
0 12
0 12
cfn=(40)
calls=1 0 
0 2
0 2

fn=(52)
0 622
0 622
0 612
0 418
0 1074
0 264
0 462
0 468
0 1099
0 308
0 1540
0 108
0 270
0 1656
0 636
0 105
0 98

fn=(54) index
0 420
0 990
0 40
0 4
0 39

fn=(70)
0 216
cob=(3)
cfi=(3)
cfn=(212)
calls=11 0 
0 33
cfn=(72)
calls=13 0 
0 26
0 144
cfn=(74) __sigsetjmp
calls=24 0 
0 456
0 192
cob=(5)
cfi=(5)
cfn=(574)
calls=11 0 
0 11001
cfn=(104) openaux
calls=12 0 
0 46427
cfn=(76) map_doit
calls=1 0 
0 2152
0 312

fn=(72)
0 30

fn=(92) _fxstat
0 70

fn=(108)
0 504
0 182

fn=(160) _dl_add_to_slotinfo
0 9
0 3
0 12

fn=(0) 0x0000003000000b00
0 2
cfn=(2) _dl_start
calls=1 0 
0 276187
0 14
cfn=(172)
calls=1 0 
0 12267
0 3
cob=(9)
cfi=(11)
cfn=(342) 0x0000000000402bc0
calls=1 0 
0 252069174592

fn=(2)
0 31
0 82
0 44
0 32
0 35
0 3
0 71
cfn=(4)
calls=1 0 
0 24
0 17
cfn=(6) _dl_sysdep_start
calls=1 0 
0 275792
0 26
0 12
0 8
0 8
0 2

fn=(8) brk
0 11

fn=(20) _dl_new_object
0 105
cfn=(10)
calls=7 0 
0 328
0 28
0 28
cfn=(24)
calls=7 0 
0 315
0 91
cfn=(34)
calls=7 0 
0 193
0 84
0 28
0 9
0 60
0 77
0 36
0 3
0 54
0 6
0 42
cfn=(10)
calls=6 0 
0 329
0 24
0 139
0 4
0 14
0 6
0 4
0 12
cfn=(28)
calls=6 0 
0 168
0 48
cfn=(52)
calls=6 0 
0 227
0 6
0 364
0 18
0 24

fn=(50) _dl_important_hwcaps
0 18
0 3
0 2
0 37
cfn=(28)
calls=1 0 
0 28
0 13
0 13
0 18
0 25
0 32
0 8
0 6
0 8
0 6
0 9
0 14
0 29
0 5
cfn=(52)
calls=1 0 
0 22
0 5
cfn=(52)
calls=1 0 
0 22
0 3

fn=(58) strsep
0 60
0 7598
0 14820
0 5928
0 104

fn=(60) mmap
0 120

fn=(76)
0 9
cfn=(78) _dl_map_object
calls=1 0 
0 2140
0 3

fn=(94) mprotect
0 50

fn=(112) _dl_sysdep_read_whole_file
0 10
cfn=(86)
calls=1 0 
0 5
0 8
cfn=(92)
calls=1 0 
0 10
0 6
0 2
cfn=(98)
calls=1 0 
0 5
0 14
cfn=(60)
calls=1 0 
0 6
0 2

fn=(120)
0 16
0 16
0 126
0 4
0 16

fn=(128) _dl_receive_error
0 11
cfn=(72)
calls=1 0 
0 2
0 7
cfn=(130) version_check_doit
calls=1 0 
0 12627
0 10

fn=(130)
0 6
cfn=(132)
calls=1 0 
0 12617
0 2
0 2

fn=(1086) _dl_fini
0 20
0 2
cob=(4)
cfi=(4)
cfn=(458)
calls=1 0 
0 23
0 15
0 4
0 88
0 5
cfn=(1088) _dl_sort_fini
calls=1 0 
0 1090
0 2
cob=(4)
cfi=(4)
cfn=(460)
calls=1 0 
0 26
0 86
0 20
0 45
cob=(3)
cfi=(3)
cfn=(1210) __libc_fini
calls=1 0 
0 6
cob=(6)
cfi=(6)
cfn=(1164)
calls=1 0 
0 1461
cob=(9)
cfi=(11)
cfn=(1090) 0x0000000000402c90
calls=1 0 
0 16
0 33
0 18
cob=(2)
cfi=(2)
cfn=(1186)
calls=1 0 
0 1480
cob=(2)
cfi=(2)
cfn=(1184)
calls=1 0 
0 3
cob=(2)
cfi=(2)
cfn=(1146)
calls=1 0 
0 1467
cob=(2)
cfi=(2)
cfn=(1126)
calls=1 0 
0 1467
cob=(2)
cfi=(2)
cfn=(1102)
calls=1 0 
0 1121
cob=(2)
cfi=(2)
cfn=(1100)
calls=1 0 
0 3
0 80
0 30
0 8
0 10

fn=(14) dl_main
0 24
cfn=(16) process_envvars
calls=1 0 
0 922
0 3
0 8
cfn=(20)
calls=1 0 
0 247
0 9
cfn=(36)
calls=1 0 
0 31
0 12
0 10
0 8
0 6
0 24
0 24
0 7
0 4
0 4
0 3
0 4
0 16
0 16
0 4
0 10
0 4
0 3
0 2
0 2
0 12
cfn=(42) strcmp
calls=1 0 
0 7
0 24
0 33
0 54
0 120
0 3
0 3
0 3
0 2
0 2
cfn=(4)
calls=1 0 
0 14
0 6
0 1
cfn=(44) _dl_discover_osversion
calls=1 0 
0 92
0 10
cfn=(48) _dl_init_paths
calls=1 0 
0 47596
0 3
cfn=(62)
calls=1 0 
0 23
0 4
0 11
0 28
cfn=(64) _dl_count_modids
calls=1 0 
0 5
0 9
0 3
cfn=(66) _dl_debug_state
calls=1 0 
0 1
0 8
0 3
cfn=(100) access
calls=1 0 
0 11
0 7
0 10
cfn=(102) _dl_map_object_deps
calls=1 0 
0 51702
0 21
0 42
0 17
0 30
0 18
0 10
0 4
0 7
0 12
cfn=(128)
calls=1 0 
0 12657
0 4
0 22
0 17
0 28
0 8
0 60
0 49
cfn=(144) _dl_relocate_object
calls=7 0 
0 148152
0 16
0 4
cfn=(160)
calls=1 0 
0 24
0 24
0 11
0 4
0 3
cfn=(162)
calls=1 0 
0 130
0 2
0 4
0 1
cfn=(164) _dl_sysdep_start_cleanup
calls=1 0 
0 1
0 6
cfn=(62)
calls=1 0 
0 9
0 3
cfn=(66)
calls=1 0 
0 1
0 2
cfn=(166) _dl_unload_cache
calls=1 0 
0 15
0 17
cfn=(144)
calls=1 0 
0 8515
0 11
0 4
0 4
0 16
0 4
0 4
0 49
0 8
cfn=(34)
calls=1 0 
0 63
0 8
cfn=(58)
calls=1 0 
0 943
0 3
cfn=(58)
calls=1 0 
0 4
0 6
0 8
cfn=(68) do_preload
calls=1 0 
0 2236
0 2
0 13
0 1
cfn=(138)
calls=1 0 
0 593
0 2

fn=(16)
0 26
cfn=(18) _dl_next_ld_env_entry
calls=1 0 
0 212
0 4
cfn=(18)
calls=2 0 
0 482
0 9
0 12
0 110
0 12
0 19
0 12
0 1
0 9
0 7
0 4
0 3

fn=(32)
0 378
0 424
0 203
0 70
cfn=(60)
calls=5 0 
0 30
0 376

fn=(56) fillin_rpath
0 19
cfn=(58)
calls=1 0 
0 1711
0 51
cfn=(58)
calls=17 0 
0 25852
0 54
0 34
cfn=(10)
calls=17 0 
0 2323
0 51
0 85
0 360
0 246
0 20
0 368
0 4
0 15
0 22
0 33
0 6
0 34
0 85
0 42
0 56
cfn=(28)
calls=14 0 
0 417
0 140
cfn=(52)
calls=14 0 
0 1219
0 70
0 255
0 42
0 140
0 105
0 56
0 56
0 84
0 38
0 21

fn=(88)
0 30

fn=(90) _dl_map_object_from_fd
0 108
cfn=(62)
calls=6 0 
0 54
0 30
cfn=(92)
calls=6 0 
0 60
0 72
0 81
0 108
0 108
cfn=(20)
calls=6 0 
0 2627
0 294
0 54
0 30
0 434
0 147
0 44
0 34
0 22
0 14
0 9
0 9
cfn=(122) _dl_next_tls_modid
calls=1 0 
0 6
0 4
0 24
0 18
0 492
0 78
0 60
0 20
0 4
0 138
cfn=(60)
calls=6 0 
0 36
0 66
0 30
0 78
0 132
cfn=(96)
calls=6 0 
0 1171
0 48
0 72
0 84
cfn=(60)
calls=6 0 
0 36
0 12
0 30
0 21
0 114
0 188
0 238
0 52
0 302
0 520
0 18
0 6
0 47
0 4
0 7
0 15
0 10
0 4
0 4
0 19
0 12
0 3
0 10
0 96
0 4
cfn=(98)
calls=1 0 
0 5
0 10
cfn=(98)
calls=5 0 
0 25
0 66
cfn=(4)
calls=6 0 
0 134
0 36
0 66
cfn=(36)
calls=6 0 
0 300
0 107
cfn=(94)
calls=5 0 
0 25
0 5
0 4
0 16
cfn=(60)
calls=2 0 
0 12
0 4
0 52
0 52
0 32
0 40

fn=(126) _wordcopy_fwd_aligned
0 12
0 10
0 10
0 4

fn=(134)
0 176
0 198
0 363
cfn=(80)
calls=66 0 
0 4092
0 132
0 77
0 130
0 88
0 80
cfn=(136) match_symbol
calls=10 0 
0 1599
0 84
cfn=(136)
calls=14 0 
0 2519
0 192
0 44
0 20
0 18
0 4
0 24
0 414
0 42
0 16
0 100
cfn=(24)
calls=7 0 
0 284
0 63
0 48
0 15
0 591
0 33
0 20
0 14
0 30
0 248
0 690
0 18
0 6
0 6

fn=(152) _dl_protect_relro
0 75
cfn=(94)
calls=5 0 
0 25
0 20

fn=(188) _dl_fixup
0 1375
0 486
0 880
cfn=(146)
calls=55 0 
0 65451
0 1660

fn=(44)
0 7
0 8
cfn=(46) uname
calls=1 0 
0 5
0 38
0 24
0 8
0 2

fn=(74)
0 456

fn=(100)
0 4
0 7

fn=(110) _dl_load_cache_lookup
0 56
0 8
0 88
0 593
cfn=(114) _dl_cache_libcmp
calls=32 0 
0 2860
0 96
0 28
0 45
0 4
cfn=(112)
calls=1 0 
0 68
0 9
0 11
0 14
0 20
0 1
0 44
0 51
cfn=(114)
calls=3 0 
0 322
0 24
cfn=(114)
calls=4 0 
0 578
0 28
0 20
0 8
0 40
0 72
0 16
cfn=(10)
calls=4 0 
0 253
0 32
cfn=(116) strcpy
calls=4 0 
0 329
0 8
cfn=(118) strdup
calls=4 0 
0 456
0 32

fn=(116)
0 120
0 84
0 11
0 39
0 40
0 27
0 8

fn=(150)
0 5540
0 2390
0 940
0 2460
cfn=(42)
calls=164 0 
0 14216
0 5026
cfn=(42)
calls=164 0 
0 11433
0 492
0 102
0 336
0 20
0 1102
0 630

fn=(166)
0 7
cfn=(168) munmap
calls=1 0 
0 5
0 3

fn=(168)
0 5

fn=(1088)
0 28
cfn=(96)
calls=1 0 
0 20
0 8
0 4
0 48
cfn=(96)
calls=6 0 
0 108
0 115
0 145
0 7
0 76
0 148
0 48
0 132
0 24
0 21
cfn=(124) memmove
calls=1 0 
0 62
0 7
0 12
cfn=(124)
calls=1 0 
0 73
0 4

fn=(38)
0 14

fn=(46)
0 5

fn=(48)
0 14
cfn=(50)
calls=1 0 
0 326
0 3
cfn=(28)
calls=1 0 
0 28
0 13
cfn=(28)
calls=1 0 
0 28
0 73
0 26
0 16
0 4
0 5
0 9
0 7
cfn=(54)
calls=1 0 
0 844
0 4
cfn=(10)
calls=1 0 
0 1652
0 8
cfn=(34)
calls=1 0 
0 214
0 14
0 10136
0 2
cfn=(28)
calls=1 0 
0 28
0 10
cfn=(56)
calls=1 0 
0 34114
0 16
0 2

fn=(82) expand_dynamic_string_token
0 13
cfn=(54)
calls=1 0 
0 52
0 4
cfn=(10)
calls=1 0 
0 73
0 3
cfn=(28)
calls=1 0 
0 28
0 13
cfn=(34)
calls=1 0 
0 52

fn=(118)
0 20
cfn=(10)
calls=4 0 
0 158
0 12
cfn=(28)
calls=4 0 
0 112
0 36
cfn=(34)
calls=4 0 
0 118

fn=(136)
0 862
0 833
0 144
cfn=(42)
calls=24 0 
0 2015
0 264

fn=(144)
0 200
0 88
0 166
0 120
0 24
0 16
0 2022
0 624
0 288
0 48
0 104
0 56
0 84
0 18
0 10
0 3448
0 268
0 903
0 236
0 3825
cfn=(146)
calls=130 0 
0 135897
0 910
0 1161
0 64
0 10
cfn=(152)
calls=5 0 
0 120
0 81
0 42
0 60
0 24
0 18
0 9
0 1
0 588
0 635
0 5
0 156
0 1
0 9
0 14
0 8
cob=(3)
cfi=(3)
cfn=(158) strncasecmp
calls=1 0 
0 6
cob=(3)
cfi=(3)
cfn=(154) strcasecmp
calls=1 0 
0 64
0 6
0 24
cfn=(34)
calls=3 0 
0 54
0 15
0 30
0 6
0 16
0 14
0 6
0 164
cfn=(146)
calls=3 0 
0 3223
0 21
0 24
0 15
0 18
0 1
0 480
0 50
0 40
0 27
0 2

fn=(186)
0 605
cfn=(188)
calls=55 0 
0 69852
0 550

fn=(6)
0 330
0 23
0 51
0 48
0 6
0 36
0 11
cfn=(8)
calls=1 0 
0 11
0 5
0 2
cfn=(12)
calls=1 0 
0 17
0 11
cfn=(14)
calls=1 0 
0 275149
0 14
0 5
0 5
0 5
0 6
0 5
0 5
0 5
0 1
cfn=(10)
calls=1 0 
0 39
0 2

fn=(18)
0 12
0 6
0 204
0 408
0 50
0 14

fn=(64)
0 2
0 3

fn=(68)
0 17
cfn=(70)
calls=1 0 
0 2209
0 10

fn=(80)
0 3192
cfn=(42)
calls=456 0 
0 3199
0 2280
0 1362
0 1416
cfn=(42)
calls=472 0 
0 8954
0 944
0 2280

fn=(104)
0 168
cfn=(78)
calls=12 0 
0 46223
0 36

fn=(42)
0 3927
0 2480
cfn=(43) strcmp'2
calls=620 0 
0 31298
0 2756

fn=(43)
0 13140
0 16316
cfn=(43)
calls=3760 0 
0 163578
0 1842

fn=(66)
0 2

fn=(78)
0 273
0 100
0 150
0 456
cfn=(80)
calls=57 0 
0 3538
0 114
0 142
0 48
cfn=(54)
calls=6 0 
0 189
0 12
0 6
cfn=(82)
calls=1 0 
0 238
0 14
cfn=(84)
calls=1 0 
0 206
0 123
cfn=(90)
calls=6 0 
0 9663
0 12
0 12
0 10
cfn=(10)
calls=5 0 
0 245
0 40
0 90
cfn=(106)
calls=5 0 
0 23872
0 15
0 12
0 20
0 24
cfn=(110)
calls=4 0 
0 6214
0 68
cfn=(84)
calls=4 0 
0 1160
0 20
0 100
0 15
0 15
0 80
0 150
0 196
cfn=(42)
calls=28 0 
0 637
0 84

fn=(102)
0 77
0 28
0 6
0 12
0 102
0 24
0 669
0 324
0 60
cfn=(54)
calls=12 0 
0 408
0 120
cfn=(70)
calls=12 0 
0 47111
0 72
0 90
0 14
0 12
0 36
cfn=(28)
calls=6 0 
0 168
0 66
cfn=(34)
calls=6 0 
0 141
0 24
cfn=(34)
calls=6 0 
0 141
0 84
0 42
0 16
0 35
0 42
0 14
0 4
cfn=(28)
calls=1 0 
0 28
0 26
0 70
0 7
0 13
cfn=(96)
calls=1 0 
0 20
0 55
0 152
0 7
0 76
0 148
0 48
0 66
0 21
0 54
cfn=(96)
calls=6 0 
0 108
0 24
0 2
0 11
cfn=(124)
calls=1 0 
0 62
0 15
cfn=(124)
calls=1 0 
0 73
0 4
0 8
0 610
0 8
cfn=(34)
calls=1 0 
0 41
0 3

fn=(114)
0 312
0 312
cfn=(115) _dl_cache_libcmp'2
calls=39 0 
0 3136

fn=(115)
0 949
0 561
0 91
0 42
0 900
0 444
cfn=(115)
calls=148 0 
0 9158
0 149

fn=(122)
0 6

fn=(124)
0 60
0 12
0 4
0 16
0 18
0 84
0 38
cfn=(126)
calls=2 0 
0 36
0 2

fn=(164)
0 1

fn=(190)
0 5

ob=(9)
fl=(11)
fn=(394) for_since_epoch_t
0 8
cfn=(396) for_set_fpe_
calls=1 0 
0 4224
0 4
cob=(3)
cfi=(3)
cfn=(424) gettimeofday
calls=1 0 
0 12
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1387
0 17
cfn=(396)
calls=1 0 
0 559
0 5

fn=(436) init_resource
0 9378

fn=(498) for__preconnected_units_create
0 53
cob=(3)
cfi=(3)
cfn=(504) __sprintf_chk
calls=1 0 
0 664
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1400
0 5
cob=(3)
cfi=(3)
cfn=(198) getenv
calls=1 0 
0 673
0 2
0 2
0 15
cob=(3)
cfi=(3)
cfn=(504)
calls=1 0 
0 664
0 5
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 2
0 2
0 15
cob=(3)
cfi=(3)
cfn=(504)
calls=1 0 
0 664
0 5
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 2
0 2
0 15

fn=(622) _intel_fast_memcpy
0 254989371
0 182135265
0 218562318
cfn=(624) _intel_fast_memcpy.P
calls=36427053 0 
0 1567812385

fn=(670) for__resource_release
0 188315536
0 24562896
0 65501056
0 8187632
cob=(4)
cfi=(4)
cfn=(460)
calls=8187632 0 
0 122814480
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1457
0 57313424
0 56092212
cob=(4)
cfi=(4)
cfn=(896)
calls=3483210 0 
0 31348890
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1444
0 13932840
cfn=(640) for__free_vm
calls=3483210 0 
0 344837790
0 13932840

fn=(714) for__release_lun
0 36844002
0 77781744
0 8187552
0 8187552
0 8187552
0 8187552
0 8187552
0 12281328
0 32750208
cfn=(670)
calls=4093776 0 
0 463415154
0 4093776
0 8187556
0 32750219
0 4
cfn=(944) for__acquire_threadstor_threads
calls=2 0 
0 126
0 44
cfn=(670)
calls=2 0 
0 362
0 2
0 6
cfn=(640)
calls=2 0 
0 246
0 4
cfn=(640)
calls=2 0 
0 268
0 4
0 16
0 6
0 6
0 2
cfn=(640)
calls=2 0 
0 198
0 10
0 10

fn=(840) for_trim
0 39
0 2
0 6
cfn=(622)
calls=1 0 
0 33
0 5

fn=(848) __intel_ssse3_rep_memmove
0 2227431
0 2227428
0 2969900
0 8
0 16
0 21
0 28
0 3712741

fn=(850) for_cpstr
0 9620
cfn=(726) _intel_fast_memcmp
calls=481 0 
0 11406
0 962
0 1864
0 4224
0 135

fn=(858) for__lub_getc
0 46566
0 6
0 18
cfn=(792) for__get_s
calls=3 0 
0 1191
0 21
0 9
0 46545

fn=(872) cvt_text_to_boolean
0 24
0 48
0 15

fn=(878) eos_default_parameters_mp_eos_primary_variable_list_
0 16
cfn=(736) for_write_seq_fmt
calls=1 0 
0 3952
0 6
cfn=(850)
calls=1 0 
0 64
0 7
cfn=(850)
calls=1 0 
0 91
0 7
cfn=(850)
calls=1 0 
0 58
0 4
0 6
cfn=(842) for_cpystr
calls=1 0 
0 94
0 6
cfn=(842)
calls=1 0 
0 164
0 6

fn=(906) for_concat
0 1725
0 276
0 120
0 876
0 1932
0 345
0 207
cfn=(844) _intel_fast_memmove
calls=69 0 
0 2073
0 69
0 207
cfn=(622)
calls=69 0 
0 2437
0 828
0 345
0 828

fn=(920) for_dealloc_allocatable
0 117
0 36
0 27
0 36
0 36
0 36
0 27
0 18
0 18
cfn=(922) _mm_free
calls=9 0 
0 856
0 9
0 54
0 18
0 81

fn=(924) cvt_ieee_t_to_text_ex
0 648159696
0 108026616
0 53295933
0 87468696
0 24146184
0 119468496
0 24146184
0 72017746
0 36008872
0 36008872
0 162039924
0 36008872
0 234057655
0 180044350
0 468115310
cfn=(926) cvtas_t_to_a
calls=18004435 0 
0 8861326585
0 126031045
0 38866332
0 108026610
cfn=(738) _intel_fast_memset
calls=18004435 0 
0 653049764
0 162039915
0 36008870
0 30950504
0 12645915
0 46425756
0 3015904
0 859815
0 37516822
0 10
0 108026580
0 161867907
0 35530635
0 1434750
0 88826550
0 36008870
0 72017740
0 758593752
0 198048785
0 574422290
cfn=(622)
calls=18004435 0 
0 603971947
0 180044350
0 72017740
0 29178483
0 19450782
0 1540
0 38903874
0 77809288
0 116713932
0 99339288
cfn=(622)
calls=8278274 0 
0 269334082
0 41391370
0 19452322
0 1195625
0 35530620
0 72017740
0 234057676
cfn=(926)
calls=1 0 
0 414
0 5
0 6
cfn=(738)
calls=1 0 
0 58
0 12
cfn=(738)
calls=1 0 
0 42
0 9
cfn=(738)
calls=1 0 
0 35
0 11
0 6
0 4
0 16
0 22
cfn=(622)
calls=1 0 
0 31
0 8
0 12
0 4
0 8
0 15
0 2
0 4
0 4
0 29178496
0 343926
0 29149749
0 38866332
0 97165830
cfn=(738)
calls=9716583 0 
0 592711563
0 19433166
0 11

fn=(940) read_generation_data_
0 16
cfn=(736)
calls=1 0 
0 4008
0 9
cfn=(736)
calls=1 0 
0 1601
0 6
cfn=(716) for_rewind
calls=1 0 
0 1130
0 9
cfn=(736)
calls=1 0 
0 1841
0 9
cfn=(736)
calls=1 0 
0 1817
0 24
cfn=(770) for_read_seq_fmt
calls=1 0 
0 1966
0 6
cfn=(810) for_read_seq_fmt_xmit
calls=1 0 
0 918
0 6
cfn=(810)
calls=1 0 
0 913
0 5
cfn=(810)
calls=1 0 
0 1044
0 5
cfn=(810)
calls=1 0 
0 1196
0 23
cfn=(942) for_read_int_fmt
calls=1 0 
0 7682
0 11
cfn=(942)
calls=1 0 
0 3070
0 15
0 9
cfn=(736)
calls=1 0 
0 1831
0 8

fn=(1000) utility_functions_mp_cubic_equation_roots_
0 958454744
cfn=(1002) cbrt
calls=18431822 0 
0 1382386719
0 55295466
0 92159110
cfn=(1002)
calls=18431822 0 
0 1382386650
0 36863644
0 239613686

fn=(1012) geological_media_properties_mp_hydraulic_property_changes_
0 720
0 10
cfn=(736)
calls=1 0 
0 4072
0 497661354
cfn=(1014) exp
calls=9215821 0 
0 263787235
0 1530
cfn=(1014)
calls=90 0 
0 5832
0 73727288
cfn=(1014)
calls=9215911 0 
0 599034215
0 202750042
0 101375021
0 720

fn=(1032) dsmv_
0 125940
cfn=(738)
calls=6297 0 
0 483937044
0 12896167842
0 27080550756
0 7641900666
0 503760
0 2579226012
0 81861

fn=(1050) for__access_threadstor_threads
0 8
0 13
0 3
0 4
0 3
0 3
cob=(4)
cfi=(4)
cfn=(560)
calls=1 0 
0 18
0 8

fn=(1062) for__deallocate_lub
0 459
0 91
0 378
0 126
0 63
cfn=(670)
calls=7 0 
0 1031
0 7
0 50
0 90
cfn=(670)
calls=10 0 
0 1338
0 10
0 80
0 20
cfn=(650) for__resource_acquire
calls=10 0 
0 1485
0 20
0 10
cob=(4)
cfi=(4)
cfn=(666)
calls=10 0 
0 20
0 30
0 20
cob=(4)
cfi=(4)
cfn=(666)
calls=10 0 
0 20
0 80
0 40
0 20
cfn=(650)
calls=10 0 
0 1658
0 20
0 30
0 180
cfn=(670)
calls=10 0 
0 1220
0 10
0 80
0 80
cfn=(670)
calls=10 0 
0 1338
0 10
0 20
cfn=(640)
calls=10 0 
0 1284
0 30
0 34
0 119

fn=(356) __libc_csu_init
0 15
cob=(2)
cfi=(2)
cfn=(358)
calls=1 0 
0 10
0 8
cfn=(366) 0x0000000000402cb0
calls=1 0 
0 16
0 4
cfn=(374) __sti__$E
calls=1 0 
0 329
0 6
0 8

fn=(388) __intel_cpu_features_init_body
0 18
0 15
0 111
0 3
0 12
0 2
0 17

fn=(438) init_resource_recurschk
0 1166

fn=(624)
0 36427053
cfn=(626) __intel_ssse3_rep_memcpy
calls=36427053 0 
0 1531385332

fn=(742) __intel_memset
0 99160332
0 198320628
0 97713280
0 287256234
0 4
0 29149755
0 12594
0 3787891
0 4
0 3
0 5
0 4
0 40532340
0 40928892
0 116103063
0 693893112
0 16392780

fn=(810)
0 112711994
0 15233165
0 8379992
0 4189996
0 2555004
0 10066030
cfn=(746) for__interp_fmt
calls=1788384 0 
0 194975368
0 613228
cfn=(746)
calls=306614 0 
0 27379246
0 10474990
0 20566524
0 10283262
0 9975730
0 24939325
cfn=(812) for__desc_ret_item
calls=4987865 0 
0 236889591
0 24939325
0 5524168
0 15347508
0 7673754
0 28496
0 10686
0 20434848
0 12789590
0 63947950
0 7354924
0 10231672
0 2157555
0 14709864
0 16066693
0 10231672
0 6472665
0 3677466
0 18387330
0 31389213
0 32904610
0 97275460
0 10227158
0 28112271
cfn=(814) for__cvt_value
calls=2555661 0 
0 969110942
0 22570
cfn=(814)
calls=2257 0 
0 204277
0 25579180
0 17904775
0 1028976
0 6130011
0 1029162
0 10222644
0 11285
cfn=(640)
calls=2257 0 
0 223471
0 5120350
0 615064
0 102352
0 51038
0 102076
cfn=(970) for__desc_test_item
calls=25519 0 
0 1352507
0 102076
0 1020830
0 205180
0 3067560
0 2147343
0 1227024
0 613512
0 613512
cfn=(714)
calls=306756 0 
0 37731224
0 613512
0 24659602
cfn=(792)
calls=25519 0 
0 7986920
0 178633
0 186
0 4514
0 11285
0 20313
cfn=(442) for__get_vm
calls=2257 0 
0 277238
0 20313
0 85
0 6720
0 11285
cfn=(622)
calls=2257 0 
0 66875
0 29341
cfn=(738)
calls=2257 0 
0 79596
0 45140
0 6771
0 4514
0 24483591
0 205256
0 512950
0 2146298
0 10193729

fn=(876) eos_default_parameters_mp_eos_defaults_
0 22
cfn=(736)
calls=1 0 
0 3848
0 6
cfn=(850)
calls=1 0 
0 64
0 7
cfn=(850)
calls=1 0 
0 91
0 7
cfn=(850)
calls=1 0 
0 58
0 6
0 24
0 4
0 8

fn=(918) for_realloc_lhs
0 945
0 369
cfn=(920)
calls=2 0 
0 290
0 26
0 20
0 26
0 26
0 6
cfn=(882) for_alloc_allocatable
calls=2 0 
0 598
0 18

fn=(966) cvt_ieee_s_to_text_ex
0 55296000
0 9216000
0 4608000
0 4608000
0 6144000
0 3072000
0 3072000
0 13824000
0 19968000
0 15360000
0 39936000
cfn=(968) cvtas_s_to_a
calls=1536000 0 
0 1276699200
0 10752000
0 9216000
cfn=(738)
calls=1536000 0 
0 56832000
0 13824000
0 3072000
0 2912000
0 400000
0 4368000
0 3072000
0 26112000
0 7680000
0 3072000
0 70656000
0 66048000
cfn=(622)
calls=1536000 0 
0 50073600
0 15360000
0 6144000
0 18432000
cfn=(622)
calls=1536000 0 
0 50073600
0 7680000
0 3072000
0 6144000
0 19968000

fn=(976) for_backspace
0 17
0 8
cfn=(648) for__acquire_lun
calls=1 0 
0 528
0 3
0 22
0 4
cfn=(720) for__write_output
calls=1 0 
0 37
0 9
0 7
0 5
0 2
0 2
0 12
cfn=(714)
calls=1 0 
0 241
0 2
0 11
0 2
0 2
0 8
0 2
0 2
0 40
0 30
0 14
0 3
0 3
0 4
0 10
0 12
0 24
0 8
0 5

fn=(986) print_input_data_
0 16
cfn=(736)
calls=1 0 
0 4003
0 11
cfn=(736)
calls=1 0 
0 3410
0 11
cfn=(736)
calls=1 0 
0 2769
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 4076
0 11
cfn=(736)
calls=1 0 
0 5472
0 6
cfn=(744) for_write_seq_fmt_xmit
calls=1 0 
0 2202
0 6
cfn=(744)
calls=1 0 
0 2202
0 6
cfn=(744)
calls=1 0 
0 2169
0 6
cfn=(744)
calls=1 0 
0 2857
0 6
cfn=(744)
calls=1 0 
0 3448
0 6
cfn=(744)
calls=1 0 
0 3448
0 6
cfn=(744)
calls=1 0 
0 3391
0 6
cfn=(744)
calls=1 0 
0 3391
0 6
cfn=(744)
calls=1 0 
0 3391
0 6
cfn=(744)
calls=1 0 
0 5396
0 6
cfn=(744)
calls=1 0 
0 2920
0 6
cfn=(744)
calls=1 0 
0 3373
0 6
cfn=(744)
calls=1 0 
0 3391
0 6
cfn=(744)
calls=1 0 
0 3373
0 6
cfn=(744)
calls=1 0 
0 3391
0 6
cfn=(744)
calls=1 0 
0 15198
0 3
0 9
cfn=(736)
calls=1 0 
0 17533
0 20
cfn=(736)
calls=1 0 
0 3951
0 9
cfn=(736)
calls=1 0 
0 2071
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 3010
0 25
cfn=(736)
calls=1 0 
0 2428
0 27
cfn=(736)
calls=3 0 
0 7284
0 40
cfn=(736)
calls=4 0 
0 6752
0 28
cfn=(744)
calls=4 0 
0 2684
0 24
cfn=(744)
calls=4 0 
0 6464
0 24
cfn=(744)
calls=4 0 
0 6342
0 24
cfn=(744)
calls=4 0 
0 6521
0 24
cfn=(744)
calls=4 0 
0 6692
0 24
cfn=(744)
calls=4 0 
0 6692
0 24
cfn=(744)
calls=4 0 
0 8812
0 36
cfn=(736)
calls=4 0 
0 7888
0 44
cfn=(736)
calls=4 0 
0 11255
0 24
cfn=(744)
calls=4 0 
0 6457
0 24
cfn=(744)
calls=4 0 
0 9809
0 20
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 2954
0 9
cfn=(736)
calls=1 0 
0 1972
0 22
cfn=(736)
calls=1 0 
0 1633
0 1228788
cfn=(736)
calls=102399 0 
0 167294367
0 614400
cfn=(744)
calls=102400 0 
0 68812800
0 614400
cfn=(744)
calls=102400 0 
0 183245623
0 409600
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 2978
0 9
cfn=(736)
calls=1 0 
0 3091
0 21
cfn=(736)
calls=1 0 
0 1765
0 2449908
cfn=(736)
calls=204159 0 
0 360519275
0 1429120
cfn=(744)
calls=204160 0 
0 124486560
0 1224960
cfn=(744)
calls=204160 0 
0 157203200
0 1224960
cfn=(744)
calls=204160 0 
0 330967360
0 1224960
cfn=(744)
calls=204160 0 
0 330967360
0 1224960
cfn=(744)
calls=204160 0 
0 329809634
0 1224960
cfn=(744)
calls=204160 0 
0 319306240
0 816640
0 3
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 3034
0 9
cfn=(736)
calls=1 0 
0 3059
0 44
cfn=(736)
calls=1 0 
0 1769
0 1228788
cfn=(736)
calls=102399 0 
0 181220631
0 614400
cfn=(744)
calls=102400 0 
0 169881543
0 409600
cfn=(744)
calls=102400 0 
0 13107200
0 1126400
cfn=(744)
calls=102400 0 
0 171110402
0 614400
cfn=(744)
calls=102400 0 
0 171110400
0 614400
0 512000
cfn=(744)
calls=102400 0 
0 49152000
0 512000
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 3026
0 9
cfn=(736)
calls=1 0 
0 21769
0 9
cfn=(736)
calls=1 0 
0 21769
0 8

fn=(1014)
0 55295466
cfn=(1018) __libm_exp_e7
calls=18431821 0 
0 807531772
cfn=(1016) __libm_exp_dispatch_table_init
calls=1 0 
0 44

fn=(1026) solve_jacobian_matrix_equation_
0 2136
0 264
0 1153
cfn=(738)
calls=89 0 
0 6839828
0 890
0 2670
cfn=(1028) dslucs_
calls=89 0 
0 143556525337
0 1335
cfn=(622)
calls=89 0 
0 13675028
0 267
0 82024002
0 1342
0 9
cfn=(736)
calls=1 0 
0 5467
0 11

fn=(1054) for__exit_handler
0 8
0 2
0 2
0 2
0 2
0 4
0 4
cfn=(640)
calls=1 0 
0 99
0 11
cfn=(1056) for__aio_release
calls=1 0 
0 948
0 72
cfn=(972) for__get_next_lub
calls=18 0 
0 14629
0 36
0 34
0 68
0 36
cfn=(1058) for__close_proc
calls=12 0 
0 13488
0 24
0 34
cfn=(1062)
calls=17 0 
0 11601
0 34
0 1
cfn=(1072) for__reentrancy_cleanup
calls=1 0 
0 25870
0 5

fn=(1090)
0 5
cfn=(1094) 0x0000000000402c10
calls=1 0 
0 8
0 3

fn=(366)
0 10
0 6

fn=(382) main
0 10
cfn=(384) __intel_new_feature_proc_init
calls=1 0 
0 234
0 3
cfn=(392) for_rtl_init_
calls=1 0 
0 126296
0 1
cfn=(592) MAIN__
calls=1 0 
0 252069045914

fn=(384)
0 28
0 8
0 4
0 26
cfn=(390) __intel_proc_init_ftzdazule
calls=2 0 
0 96
0 2
0 2
cfn=(386) __intel_cpu_features_init
calls=1 0 
0 183
0 3

fn=(588) for__hbw_set_policy
0 5
0 3

fn=(620) __intel_sse2_strlen
0 532
0 1356
0 304

fn=(726)
0 37085823
0 24721014
0 24721002
0 20603215
0 1620
0 234
0 8243170

fn=(744)
0 1097346366
0 323795514
0 118756904
0 59378452
0 29950348
0 146732658
cfn=(746)
calls=25582345 0 
0 4638640452
0 8213762
cfn=(746)
calls=4106881 0 
0 675268078
0 148446130
0 158452431
0 105633378
0 104807368
0 103953040
0 49230350
0 15147800
0 18934750
0 7573900
cfn=(748) for__put_sf
calls=3786950 0 
0 727702007
0 7573900
0 53017300
0 4
cfn=(640)
calls=1 0 
0 187
0 2
0 11360850
0 312076474
0 103953040
0 48
0 20
0 190
0 70
0 100
cfn=(938) for__lower_bound_index
calls=10 0 
0 192
0 10
0 18
0 40
0 259882540
cfn=(812)
calls=51976508 0 
0 2445595900
0 259882540
0 57703638
0 152866422
0 216
0 270
0 162
0 85
0 220
0 94
0 407
0 162
0 45
0 110
0 162
cfn=(938)
calls=27 0 
0 607
0 27
0 76433130
0 46735
0 7190
0 203792920
0 127388685
0 36
0 636943725
0 87995060
0 101910996
0 4916656
0 39080872
0 127388745
0 76433247
0 127388745
0 76433247
0 101910996
cfn=(824) for__format_value
calls=25477749 0 
0 22662139203
0 50955498
0 1478772
0 203821992
0 2990148
0 1281492
0 3795392
0 2990148
0 15147800
cfn=(714)
calls=3786950 0 
0 671298812
0 22721700
0 44
0 66
0 22
cfn=(822) for__adjust_buffer
calls=2 0 
0 329
0 4
0 254371579
0 10
cfn=(822)
calls=1 0 
0 163
0 5
0 2065025
0 274
0 1096
0 825736
0 4728
0 2364
0 1576
cfn=(748)
calls=788 0 
0 155737
0 5516
cfn=(722) for__flush_readahead
calls=788 0 
0 144204
0 5516
0 2364
0 5516
0 1576
0 6304
0 3940
cfn=(738)
calls=788 0 
0 80615
0 2364
0 16427536
0 12320643
0 88806604
0 66604917

fn=(824)
0 1350320697
0 433121733
0 15
0 3
0 7374984
0 39080860
0 17208296
0 21203352
cfn=(890) cvt_integer_to_text
calls=2355928 0 
0 299266251
0 2355928
0 921600
cfn=(890)
calls=102400 0 
0 8294400
0 102400
0 2
0 7
0 7
cfn=(990) cvt_boolean_to_text_ex
calls=1 0 
0 60
0 1
0 31
0 14
cfn=(924)
calls=1 0 
0 845
0 1
0 527591745
0 270066525
cfn=(924)
calls=18004435 0 
0 17157814727
0 18004435
0 23040000
cfn=(966)
calls=1536000 0 
0 1906766400
0 58621308
0 4916658
0 3327237
0 24352888
0 68692800
0 20873904
cfn=(622)
calls=3478984 0 
0 118563850
0 50955498
0 229299741

fn=(828) for_read_seq_lis
0 135
0 108
cfn=(648)
calls=9 0 
0 3195
0 27
0 306
0 36
0 27
0 18
0 18
0 171
0 54
0 81
cfn=(792)
calls=9 0 
0 3499
0 27
0 108
cfn=(442)
calls=9 0 
0 1899
0 27
0 36
cfn=(830) for_read_seq_lis_xmit
calls=9 0 
0 9663
0 72

fn=(838) for_adjustl
0 39
0 27
0 15
0 9
0 12
cfn=(622)
calls=3 0 
0 99
0 6
0 32
0 4

fn=(852) for_read_seq_nml
0 255
0 255
cfn=(648)
calls=17 0 
0 6035
0 51
0 578
0 102
0 51
0 34
0 34
0 102
0 170
0 85
cfn=(812)
calls=17 0 
0 1054
0 34
0 272
0 102
cfn=(792)
calls=17 0 
0 7809
0 34
0 221
cfn=(442)
calls=17 0 
0 6094
0 34
0 323
cfn=(442)
calls=17 0 
0 3587
0 34
0 663
cfn=(854) for__nml_lex
calls=17 0 
0 2363
0 11960
cfn=(854)
calls=598 0 
0 787453
0 1230
0 9225
0 2460
0 68
cfn=(640)
calls=17 0 
0 2137
0 204
cfn=(874) for__free_nml_buffer
calls=17 0 
0 3586
0 204
0 120
cfn=(714)
calls=15 0 
0 1845
0 14
cfn=(714)
calls=2 0 
0 246
0 34
0 3374
0 504
0 70
0 700
0 270
0 645
0 312
0 540
0 270
0 675
0 315
0 450
0 585
cfn=(860) do_nml_copy
calls=135 0 
0 4328
0 675
0 90
cfn=(914) bump_up_an_indecie
calls=45 0 
0 1596
0 495
0 832
0 268
0 1248
cfn=(812)
calls=104 0 
0 6448
0 2105
cfn=(812)
calls=421 0 
0 26102
0 4725
0 1360
0 5088
0 20152
0 9000
0 832
cfn=(812)
calls=104 0 
0 6535
0 2105
cfn=(812)
calls=421 0 
0 26032
0 3675
0 2313
0 1560
0 196
0 182
0 294
0 1830
0 136
0 187
0 540
0 186
0 30
0 468
0 28
0 208
0 261
0 64
0 68
0 170
0 906
0 5104
0 68

fn=(854)
0 12300
0 14760
0 1845
0 9840
cfn=(856) for__lex_getc
calls=615 0 
0 31338
0 17644
cfn=(856)
calls=4411 0 
0 224292
0 35182
0 4918
0 84
0 68
cfn=(862) scan_ahead_is_identifier
calls=17 0 
0 38860
0 68
0 416
cfn=(856)
calls=104 0 
0 4368
0 520
0 624
0 416
0 42
0 351
0 68
0 51
0 8
0 702
0 24
0 32
0 32
cfn=(856)
calls=4 0 
0 168
0 48
0 2
0 15
0 21
cfn=(870) process_logical
calls=3 0 
0 1846
0 15
0 10
0 13
cfn=(866) process_numeric
calls=1 0 
0 1080
0 5
0 1053
0 468
0 8
cfn=(864) process_nondelim_char..1
calls=1 0 
0 371
0 2
0 812
cfn=(866)
calls=116 0 
0 130360
0 232
0 348
0 102
0 119
0 140
cfn=(856)
calls=14 0 
0 588
0 228
cfn=(856)
calls=57 0 
0 2394
0 497
0 57
0 56
cfn=(856)
calls=14 0 
0 588
0 112
0 84
0 84
0 616
0 432
0 312
0 371
0 363
0 189
0 19115
0 12972
0 5787
0 90
0 18537
0 386
0 116
0 68
0 513
0 5517
0 847
0 7955
0 8401
0 726
0 1089
0 7800
cfn=(858)
calls=1300 0 
0 33800
0 14300
0 5200
0 9144
cfn=(858)
calls=1524 0 
0 39624
0 16764
0 6096
0 105
0 120
cfn=(864)
calls=15 0 
0 31097
0 105
cfn=(856)
calls=15 0 
0 450
0 30

fn=(866)
0 11220
0 7716
0 1746
0 3132
0 968
0 492
0 440
0 198
0 798
0 96
0 498
0 1053
0 203
0 264
0 1287
0 184
0 150
0 50
0 150
0 50
0 150
cfn=(814)
calls=25 0 
0 4743
0 225
0 125
0 1288
0 368
0 276
0 828
0 828
0 20
cfn=(814)
calls=4 0 
0 1552
0 264
cfn=(814)
calls=88 0 
0 53567
0 552
0 644
0 1497
0 3032
cfn=(858)
calls=758 0 
0 20932
0 8338
0 1516

fn=(874)
0 170
cfn=(640)
calls=17 0 
0 3331
0 85

fn=(882)
0 455
0 140
0 140
0 6
0 3
0 5
0 136
0 350
0 175
0 70
0 140
0 105
0 105
cfn=(884) _mm_malloc
calls=35 0 
0 10197
0 175
0 140
0 70
0 350

fn=(908) read_gas_specific_data_
0 14
cfn=(736)
calls=1 0 
0 4152
0 140
cfn=(852)
calls=1 0 
0 182994
0 105
0 46
0 11
cfn=(736)
calls=1 0 
0 7781
0 11
cfn=(736)
calls=1 0 
0 2477
0 6
cfn=(744)
calls=1 0 
0 1302
0 6
cfn=(744)
calls=1 0 
0 1302
0 6
cfn=(744)
calls=1 0 
0 1486
0 6
cfn=(744)
calls=1 0 
0 1575
0 6
cfn=(744)
calls=1 0 
0 1861
0 11
cfn=(770)
calls=1 0 
0 1711
0 42
cfn=(880) for_check_mult_overflow64
calls=1 0 
0 81
0 16
cfn=(882)
calls=1 0 
0 299
0 17
cfn=(918)
calls=1 0 
0 535
0 14
0 24
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 18
cfn=(918)
calls=1 0 
0 535
0 15

fn=(932) read_initial_conditions_
0 12
cfn=(736)
calls=1 0 
0 4064
0 9
cfn=(736)
calls=1 0 
0 1601
0 9
cfn=(736)
calls=1 0 
0 1841
0 6
cfn=(716)
calls=1 0 
0 834
0 9
cfn=(736)
calls=1 0 
0 1526
0 10
cfn=(842)
calls=1 0 
0 94
0 31
cfn=(852)
calls=1 0 
0 30306
0 18
0 7
cfn=(716)
calls=1 0 
0 834
0 18
cfn=(736)
calls=1 0 
0 12143
0 6
cfn=(744)
calls=1 0 
0 524
0 8
cfn=(744)
calls=1 0 
0 1838
0 8
cfn=(744)
calls=1 0 
0 520
0 6
cfn=(744)
calls=1 0 
0 524
0 8
cfn=(744)
calls=1 0 
0 1935
0 15
cfn=(744)
calls=1 0 
0 6359
0 7
0 9
cfn=(736)
calls=1 0 
0 1900
0 20
cfn=(828)
calls=2 0 
0 4350
0 22
0 22
cfn=(736)
calls=1 0 
0 1844
0 24
cfn=(828)
calls=1 0 
0 2053
0 10
cfn=(828)
calls=1 0 
0 2105
0 30
0 47
cfn=(852)
calls=1 0 
0 52643
0 29
cfn=(736)
calls=1 0 
0 2608
0 6
cfn=(744)
calls=1 0 
0 1568
0 6
cfn=(744)
calls=1 0 
0 2900
0 6
cfn=(744)
calls=1 0 
0 7241
0 6
cfn=(744)
calls=1 0 
0 1613
0 13
cfn=(744)
calls=1 0 
0 5574
0 1
0 15
cfn=(736)
calls=1 0 
0 1724
0 1

fn=(936) emit_comp_fmt
0 350
0 28
0 36
0 250
0 32
0 56
0 204
0 64
0 50
0 585
0 48
0 56
0 1130
0 20
0 29
0 30
0 16

fn=(970)
0 791089
0 76557
0 51038
0 204152
0 229671

fn=(972)
0 522
0 145
0 278
0 174
0 1863
0 168
0 84
cfn=(650)
calls=28 0 
0 4331
0 56
0 28
cob=(4)
cfi=(4)
cfn=(666)
calls=28 0 
0 56
0 84
0 56
cob=(4)
cfi=(4)
cfn=(666)
calls=28 0 
0 56
0 308
0 140
0 56
cfn=(650)
calls=28 0 
0 4677
0 168
0 336
cfn=(670)
calls=28 0 
0 3534
0 28
0 112
0 2488
0 22
0 112
cfn=(718) process_existing_lub
calls=28 0 
0 454
0 420

fn=(994) simulation_cycle_
0 11
cfn=(736)
calls=1 0 
0 3912
0 19
cfn=(996) equation_of_state_
calls=1 0 
0 129350847
0 1
cfn=(1012)
calls=1 0 
0 19152959
0 3
cfn=(1020) mass_and_energy_balance_
calls=1 0 
0 5360828
0 511
cfn=(736)
calls=31 0 
0 66985
0 155
cfn=(744)
calls=31 0 
0 27088
0 186
cfn=(744)
calls=31 0 
0 59370
0 186
cfn=(744)
calls=31 0 
0 60161
0 279
0 6
0 62
0 10
cfn=(736)
calls=1 0 
0 4000
0 124
0 196
0 68
0 56
0 9
0 12
0 608
0 341
cfn=(738)
calls=31 0 
0 6046798
0 31
0 890
cfn=(1022) jacobian_setup_
calls=89 0 
0 35863516752
0 248
cfn=(1022)
calls=31 0 
0 12571500970
0 840
0 1246
cfn=(736)
calls=89 0 
0 181652
0 534
cfn=(744)
calls=89 0 
0 75650
0 534
cfn=(744)
calls=89 0 
0 176308
0 534
cfn=(744)
calls=89 0 
0 168826
0 1246
cfn=(744)
calls=89 0 
0 68797
0 534
cfn=(744)
calls=89 0 
0 76629
0 178
0 178
0 178
0 10
cfn=(736)
calls=1 0 
0 3603
0 267
cfn=(1024) cpu_timing_routine_
calls=89 0 
0 6942
0 267
cfn=(1026)
calls=89 0 
0 143659079739
0 178
0 267
cfn=(1024)
calls=89 0 
0 6942
0 534
0 178
cfn=(996)
calls=89 0 
0 10626358637
0 178
0 89
cfn=(1012)
calls=89 0 
0 1719195080
0 178
0 267
0 93
cfn=(1024)
calls=31 0 
0 2418
0 310
0 180
0 9
cfn=(1040) write_file_save_
calls=1 0 
0 221419648
0 1
0 248
cfn=(1038) upon_convergence_
calls=31 0 
0 831180130
0 62
0 8
cfn=(1044) print_standard_output_
calls=4 0 
0 37880570115
0 12
cfn=(1020)
calls=4 0 
0 21431146
0 4
0 156
0 15
0 2

fn=(1008) __libm_cbrt_e7
0 2580455080
0 73727288

fn=(1048) for__io_return
0 35
cfn=(620)
calls=1 0 
0 29
0 9
cob=(3)
cfi=(3)
cfn=(790) __strncpy_chk
calls=1 0 
0 444
0 6
0 13
0 14
cfn=(1050)
calls=1 0 
0 60
0 6
0 15
cfn=(726)
calls=1 0 
0 52
0 9
0 5
0 4
0 2
0 3
0 2
0 6
0 2
0 2
0 2
cfn=(714)
calls=1 0 
0 123
0 2
0 13

fn=(1058)
0 144
0 60
0 36
0 156
cfn=(1060) for__rm_from_lf_table
calls=12 0 
0 6198
0 48
0 36
cob=(3)
cfi=(3)
cfn=(684) fileno
calls=12 0 
0 84
0 36
0 2
0 24
cob=(3)
cfi=(3)
cfn=(684)
calls=12 0 
0 84
0 24
0 22
cob=(3)
cfi=(3)
cfn=(684)
calls=11 0 
0 77
0 22
0 20
cob=(3)
cfi=(3)
cfn=(684)
calls=10 0 
0 70
0 60
cfn=(720)
calls=10 0 
0 350
0 30
cob=(4)
cfi=(4)
cfn=(1068)
calls=10 0 
0 70
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 979
0 20
0 6
0 36
0 24
0 48
0 24
cfn=(640)
calls=12 0 
0 1188
0 24
0 60
cfn=(640)
calls=12 0 
0 1530
0 24
0 60
cfn=(640)
calls=12 0 
0 1644
0 24
0 36
0 108

fn=(1076) kill_resource
0 10420
0 6
cob=(4)
cfi=(4)
cfn=(896)
calls=6 0 
0 54
0 12
cfn=(640)
calls=6 0 
0 594
0 4180

fn=(386)
0 3
cfn=(388)
calls=1 0 
0 178
0 2

fn=(440) for__acquire_semaphore_threaded
0 36
0 6
0 18

fn=(478) for_check_env_name
0 377
cob=(3)
cfi=(3)
cfn=(198)
calls=29 0 
0 19792
0 87
0 377

fn=(530) for__once_private
0 9
0 13
0 3
0 4
0 6
cfn=(532) real_aio_init
calls=1 0 
0 20903
0 9

fn=(594) for_set_reentrancy
0 10

fn=(628) for__aio_acquire_lun_fname
0 224
cob=(4)
cfi=(4)
cfn=(458)
calls=14 0 
0 266
0 56
0 3724
0 7504
0 308
0 520
0 6298
0 7448
0 28
cob=(4)
cfi=(4)
cfn=(460)
calls=14 0 
0 210
0 28
0 126

fn=(648)
0 77781573
0 45031437
0 16375068
0 8187530
0 49125180
0 8187530
cfn=(650)
calls=4093765 0 
0 555110021
0 8187530
0 4093765
cob=(4)
cfi=(4)
cfn=(666)
calls=4093765 0 
0 8187530
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1265
0 12281295
0 8187530
cob=(4)
cfi=(4)
cfn=(666)
calls=4093765 0 
0 8187530
0 32750120
0 8187510
0 12281265
0 8187518
0 70
cfn=(668) for__aio_acquire_lun
calls=10 0 
0 1050
0 20
0 30
0 57312710
0 48
0 12281247
0 20468825
0 4093755
cfn=(718)
calls=4093755 0 
0 73687580
0 12281265
0 36843935
cfn=(442)
calls=10 0 
0 2091
0 30
0 100
0 630
0 110
0 20
0 30
0 40
0 30
cfn=(650)
calls=10 0 
0 2880
0 70
0 36843995
cfn=(670)
calls=4093765 0 
0 463414148
0 8187530
0 70
0 20
0 36843795
0 16375020
0 8187510
cfn=(650)
calls=4093755 0 
0 555105862
0 24562530
0 24562530
0 4
0 8
cfn=(442)
calls=2 0 
0 519
0 6
0 4
0 20
0 126
0 10
0 6
0 8
0 4
cfn=(650)
calls=2 0 
0 470
0 6
cob=(4)
cfi=(4)
cfn=(666)
calls=2 0 
0 4
0 8
0 8
0 4
cfn=(944)
calls=2 0 
0 4784
0 24
cfn=(442)
calls=2 0 
0 244
0 22
0 44
0 8
cfn=(442)
calls=2 0 
0 390
0 4
0 28
0 48
0 10

fn=(748)
0 41665118
0 11363214
0 22726428
0 15150952
0 37877380
0 11363214
0 11363214
0 37877380
0 11363214
0 30301904
0 7575476
0 22726428
0 7575476
0 8183464
0 15150952
cfn=(720)
calls=3787738 0 
0 367412954
0 18938690
0 70
cob=(4)
cfi=(4)
cfn=(732)
calls=14 0 
0 98
0 42
cob=(3)
cfi=(3)
cfn=(768) ftruncate
calls=14 0 
0 70
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1342
0 70
0 49240594

fn=(814)
0 48602931
0 7674147
0 5116098
0 17906343
0 6
0 1225016
0 376
0 4900064
0 4287556
cfn=(868) cvt_text_to_integer
calls=612508 0 
0 129132833
0 612508
0 282
0 12
0 15
cfn=(872)
calls=3 0 
0 87
0 3
0 282
0 11037002
0 60090856
0 2452688
0 13489784
cfn=(910) cvt_text_to_ieee_t_ex
calls=1226344 0 
0 599543703
0 2452688
0 3595970
0 35
cfn=(622)
calls=5 0 
0 169
0 5
0 2876756
cfn=(622)
calls=719189 0 
0 24082013
0 1438378
0 24
0 2157582
0 3677710
0 23022441

fn=(834) for__ignore_space
0 108
0 153
0 45
0 63

fn=(912) cvtas_a_to_t
0 15942472
0 42922040
0 26979568
0 480
0 2453336
0 3678996
0 10
0 2452654
0 118953667
0 3678915
0 85
0 2452656
0 8584296
0 2980046
0 85861965
0 8983216
0 15944979
0 12263502
0 2452688
0 816708
0 1633416
0 2044334
0 12265884
0 80
0 5110785
0 17376839
0 80012
0 7857288
0 2946483
0 6133002
0 2044334
0 40018
0 48
0 4088668
0 5927538
0 83389576
0 12263220
0 250366
0 7357812
0 3579316
0 4640594
0 7158628
0 5667732
0 1226322
0 4910745

fn=(930) read_output_options_
0 10
cfn=(736)
calls=1 0 
0 3960
0 84
cfn=(852)
calls=1 0 
0 114733
0 21
0 14
0 2
0 24
cfn=(852)
calls=1 0 
0 23790
0 13
cfn=(770)
calls=1 0 
0 1711
0 14

fn=(962) read_external_input_files_
0 17
cfn=(736)
calls=1 0 
0 5812
0 13
cfn=(906)
calls=1 0 
0 170
0 12
cfn=(596) for_inquire
calls=1 0 
0 3046
0 13
cfn=(736)
calls=1 0 
0 2444
0 12
cfn=(638) for_open
calls=1 0 
0 6863
0 6
cfn=(716)
calls=1 0 
0 1421
0 1
0 13
cfn=(906)
calls=1 0 
0 170
0 11
cfn=(596)
calls=1 0 
0 3046
0 13
cfn=(736)
calls=1 0 
0 2444
0 12
cfn=(638)
calls=1 0 
0 6970
0 6
cfn=(716)
calls=1 0 
0 1421
0 1
0 13
cfn=(906)
calls=1 0 
0 170
0 11
cfn=(596)
calls=1 0 
0 3046
0 13
cfn=(736)
calls=1 0 
0 2444
0 12
cfn=(638)
calls=1 0 
0 6950
0 6
cfn=(716)
calls=1 0 
0 1421
0 1
0 9
cfn=(736)
calls=1 0 
0 1569
0 12
cfn=(736)
calls=1 0 
0 3920
0 6
cfn=(716)
calls=1 0 
0 834
0 11
cfn=(770)
calls=1 0 
0 1738
0 1638418
cfn=(770)
calls=102401 0 
0 158317274
0 614406
cfn=(810)
calls=102401 0 
0 52471473
0 614406
cfn=(810)
calls=102401 0 
0 122260395
0 614406
cfn=(810)
calls=102401 0 
0 107494812
0 614406
cfn=(810)
calls=102401 0 
0 107255772
0 614406
cfn=(810)
calls=102401 0 
0 122778791
0 614406
cfn=(810)
calls=102401 0 
0 62099237
0 1638416
0 204800
0 8
0 716793
0 1331195
0 204800
0 204800
0 3891200
0 3135
0 836
0 1126400
cfn=(842)
calls=102400 0 
0 5222400
0 2355200
0 6
cfn=(964) grid_geometry_mp_renumber_elements_
calls=1 0 
0 5302
0 9
cfn=(736)
calls=1 0 
0 2652
0 19
cfn=(736)
calls=1 0 
0 2139
0 1023990
cfn=(736)
calls=102399 0 
0 223441306
0 614400
cfn=(744)
calls=102400 0 
0 193442880
0 614400
cfn=(744)
calls=102400 0 
0 193788160
0 614400
cfn=(744)
calls=102400 0 
0 230707200
0 409600
0 1945625
0 4
0 9
cfn=(736)
calls=1 0 
0 3968
0 18
cfn=(770)
calls=1 0 
0 1512
0 15
0 960020
0 4
0 960011
0 4
0 12
cfn=(842)
calls=1 0 
0 51
0 2560008
cfn=(842)
calls=320001 0 
0 16320051
0 1600010
0 2880025
cfn=(842)
calls=320002 0 
0 16320102
0 1280008
0 3470741
cfn=(770)
calls=204161 0 
0 285249791
0 2041610
cfn=(810)
calls=204161 0 
0 116165636
0 1837449
cfn=(810)
calls=204161 0 
0 147404467
0 1837449
cfn=(810)
calls=204161 0 
0 239205656
0 1837449
cfn=(810)
calls=204161 0 
0 239205656
0 1837449
cfn=(810)
calls=204161 0 
0 239967012
0 1837449
cfn=(810)
calls=204161 0 
0 249893167
0 7145635
0 4287381
0 2041600
0 816640
0 18
0 10
cfn=(736)
calls=1 0 
0 28485
0 10
cfn=(736)
calls=1 0 
0 2594
0 5
cfn=(744)
calls=1 0 
0 746
0 5
cfn=(744)
calls=1 0 
0 26694
0 14
0 15
cfn=(596)
calls=1 0 
0 11146
0 1
cfn=(974) read_file_gener_
calls=1 0 
0 15863
0 3
cfn=(978) read_initial_conditions_file_
calls=1 0 
0 18894253
0 408328
0 11
cfn=(770)
calls=1 0 
0 1108
0 1837445
cfn=(810)
calls=204160 0 
0 154745731
0 1837440
cfn=(810)
calls=204160 0 
0 142049359
0 816640
0 4
cfn=(810)
calls=1 0 
0 224
0 3
0 204800
0 204800
0 8

fn=(1002)
0 110590932
cfn=(1008)
calls=36863643 0 
0 2654182296
cfn=(1004) __libm_cbrt_dispatch_table_init
calls=1 0 
0 141

fn=(1016)
0 18
cfn=(1018)
calls=1 0 
0 24
0 2

fn=(1018)
0 276477330
0 432435870
0 61754952
0 36863644

fn=(1052) for_stop_core
0 20
0 2
0 2
cfn=(440)
calls=1 0 
0 30
0 6
0 8
0 2
cfn=(1054)
calls=1 0 
0 67014
0 5
cob=(3)
cfi=(3)
cfn=(1084) exit
calls=1 0 
0 8815
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1095

fn=(1078) kill_resource_recurschk
0 1060
0 7
cob=(4)
cfi=(4)
cfn=(896)
calls=7 0 
0 63
0 14
cfn=(640)
calls=7 0 
0 693
0 650

fn=(342)
0 11
cob=(3)
cfi=(3)
cfn=(350) (below main)
calls=1 0 
0 252069172998
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1583

fn=(528) for__aio_init
0 7
cfn=(530)
calls=1 0 
0 20947

fn=(626)
0 109281159
0 72840958
0 157365996
0 40771252
0 13148
0 3047
0 549
0 25188
0 564
0 3601
0 7202
0 42533904
0 1939
0 100752
0 163722
0 966916944
0 50376
0 141305031

fn=(636) for__common_inquire
0 180
0 45
0 84
0 42
0 45
0 45
0 45
0 15
0 3
0 45
0 45
0 45
0 45
0 45
0 105

fn=(710) for__this_image_number_or_zero
0 13904
0 13907
0 3

fn=(712) for__set_terminator_option
0 165
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 687
0 3
0 30
0 20
0 22
0 149

fn=(772) for__open_default
0 42
0 56
0 4
cob=(3)
cfi=(3)
cfn=(198)
calls=2 0 
0 1346
0 4
0 6
cob=(3)
cfi=(3)
cfn=(198)
calls=2 0 
0 1346
0 4
0 6
cob=(3)
cfi=(3)
cfn=(198)
calls=2 0 
0 1374
0 4
0 6
cob=(3)
cfi=(3)
cfn=(198)
calls=2 0 
0 1346
0 4
0 28
0 3
0 2
0 2
0 2
0 2
0 6
0 4
0 16
0 6
0 14
0 22
cfn=(676) for__open_proc
calls=2 0 
0 14975
0 16

fn=(792)
0 5318496
0 997218
0 2326842
0 2659248
0 11966616
cfn=(726)
calls=332406 0 
0 9307368
0 2659248
0 4986030
0 30
0 1329548
0 664804
0 2325645
0 3654574
0 1994412
0 996705
0 306
0 210
0 26
0 234
0 4
0 1662036
0 1276
0 334
0 304
0 417
0 26
0 26
0 91
0 65
0 501
0 996717
0 1993410
0 2657768
0 1661105
0 664498
0 2990097
cob=(3)
cfi=(3)
cfn=(808) memchr
calls=332233 0 
0 20796357
0 2325631
0 2303772
0 6
0 13
0 15
0 42848
cob=(3)
cfi=(3)
cfn=(684)
calls=3296 0 
0 23072
0 16480
0 835
0 12516
0 334
0 28116
0 21868
cfn=(622)
calls=3124 0 
0 107878
0 78100
cfn=(794) for__read_input
calls=3124 0 
0 224951
0 3612
cfn=(794)
calls=172 0 
0 13483
0 7
cfn=(794)
calls=1 0 
0 72
0 42861
0 13188
0 6592
0 20
0 6584
0 26376
cob=(3)
cfi=(3)
cfn=(808)
calls=3296 0 
0 157829
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1304
0 16480
0 9888
0 3295
0 4
0 16
0 4
0 4
0 56015
0 1336
cob=(4)
cfi=(4)
cfn=(732)
calls=167 0 
0 1169
0 1336
0 664819
0 997215
0 997077
0 997077
0 1336
0 1328768
0 18
0 184
0 664700
0 997215
0 72
0 3656135
0 9885
0 2
0 2
0 8
0 501
cfn=(710)
calls=167 0 
0 1336
0 501

fn=(862)
0 462
0 904
0 2736
0 2757
0 1122
0 4128
0 1890
0 3780
0 1888
0 1890
0 554
0 2721
0 7542
0 3417
0 9440
0 224
0 29
0 201
0 352
0 128
0 7896
0 160
0 194
0 4
0 3392
0 11024
0 68
0 238

fn=(864)
0 359
0 160
0 7
0 32
0 188
0 165
cfn=(862)
calls=16 0 
0 30281
0 32
0 139
cob=(3)
cfi=(3)
cfn=(408) __ctype_b_loc
calls=1 0 
0 5
0 49
0 9
0 11
0 19
0 12

fn=(884)
0 31410774
cob=(3)
cfi=(3)
cfn=(448) malloc
calls=3490086 0 
0 712096047
0 45371118

fn=(914)
0 945
0 315
0 114
0 222

fn=(934) for__format_compiler
0 22
cfn=(442)
calls=1 0 
0 211
0 14
0 12
0 12
0 635
0 658
0 384
0 786
0 5
0 11
cfn=(936)
calls=1 0 
0 65
0 15
0 4
0 10
cfn=(936)
calls=1 0 
0 67
0 15
0 10
0 7
0 2
0 2
0 42
0 30
0 80
0 1928
0 16
0 16
0 72
cfn=(936)
calls=8 0 
0 1851
0 120
0 13
0 7
0 8
0 38
0 49
0 4
0 30
cfn=(936)
calls=6 0 
0 414
0 4
cfn=(936)
calls=1 0 
0 67
0 105
0 60
0 20
0 12
0 16
cfn=(936)
calls=4 0 
0 260
0 60
0 12
0 4
0 4
0 10
0 7
cfn=(936)
calls=1 0 
0 65
0 8
0 4
0 15
cfn=(936)
calls=1 0 
0 65
0 15
0 7
cfn=(936)
calls=1 0 
0 65
0 18
0 65
0 309
0 168
0 4
cfn=(936)
calls=1 0 
0 65
0 16

fn=(938)
0 360
0 126
0 192
0 121

fn=(942)
0 30
0 20
cfn=(648)
calls=2 0 
0 6889
0 6
0 134
0 160
0 30
0 48
0 12
0 4
0 6
0 16
cfn=(812)
calls=2 0 
0 124
0 4
0 28
0 8
cfn=(442)
calls=2 0 
0 249
0 6
0 12
cfn=(622)
calls=2 0 
0 62
0 20
0 4
0 12
cfn=(960) for_read_int_fmt_xmit
calls=2 0 
0 2852
0 16

fn=(946) create_threadstor_key
0 7
cob=(4)
cfi=(4)
cfn=(554)
calls=1 0 
0 27
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1439

fn=(964)
0 29
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 26
cfn=(880)
calls=1 0 
0 81
0 12
cfn=(736)
calls=1 0 
0 2242
0 14
0 19
0 14
0 41
cfn=(920)
calls=1 0 
0 145
0 16
cfn=(736)
calls=1 0 
0 2258
0 8

fn=(974)
0 17
cfn=(736)
calls=1 0 
0 3832
0 55
0 20
0 10
0 6
cfn=(716)
calls=1 0 
0 1125
0 11
cfn=(770)
calls=1 0 
0 1913
0 13
0 15
cfn=(770)
calls=1 0 
0 1537
0 9
cfn=(810)
calls=1 0 
0 922
0 10
cfn=(810)
calls=1 0 
0 871
0 9
cfn=(810)
calls=1 0 
0 1055
0 9
cfn=(810)
calls=1 0 
0 1324
0 34
0 17
cfn=(976)
calls=1 0 
0 1079
0 9
cfn=(736)
calls=1 0 
0 1958
0 3

fn=(996)
0 810
0 10
cfn=(736)
calls=1 0 
0 3792
0 270
0 180
0 180
0 4
0 4
0 3
0 92161538
0 18227022
0 145816176
0 63794577
0 204800
0 1228800
0 111102932
0 18431822
0 36863644
0 92159110
0 10
cfn=(736)
calls=1 0 
0 4040
0 340988707
0 119806843
0 27647733
0 165886398
cfn=(998) gas_thermophysical_properties_mp_compute_gas_density_
calls=9215911 0 
0 3114981933
0 156670487
cfn=(998)
calls=9215911 0 
0 3114977918
0 129022754
cfn=(1010) gas_thermophysical_properties_mp_compute_gas_viscosity_
calls=18431822 0 
0 847867585
0 221181864
0 1207284341
0 55295466
0 460795550
0 36863644
0 92159110
0 270
0 729
cfn=(736)
calls=1 0 
0 3920
0 3
0 3686411
0 1126400
0 512000
0 512000
0 1024000
0 512000
0 2764800
0 4
0 204800
0 204800
0 18227022
0 55295466
0 204800
0 2

fn=(1030) dsilus_
0 4272
0 890
0 11391911
0 89
0 178
0 154930132
0 81679839
0 36340124
0 72680248
0 109020372
0 45567555
0 2670
0 59236086
0 4450
0 15948444
0 89
0 227839377
0 81679839
0 145360496
0 181700620
0 109020372
0 63794577
0 1780
0 45566220
0 1780
0 45566665
0 244531683
0 81510561
0 136134400
0 298872057
0 72453832
0 117737477
0 170079
0 63794488
0 100246663
0 90567290
0 117620976
0 49172856
0 136190292
0 118473952
0 109359996
0 108680748
0 63396480
0 36072768
0 17344854
0 136190292
0 91133597
0 91133330
0 45283645
0 45566665
0 267
0 623
0 1068
0 1513
0 17085330
0 1602
0 1246
0 1157
0 18170685
0 236210806
0 18113458
0 73759284
0 3275022
0 34689708

fn=(1044)
0 36
0 10
cfn=(736)
calls=1 0 
0 4168
0 28
0 9
cfn=(736)
calls=1 0 
0 2746
0 1
0 92
cfn=(736)
calls=4 0 
0 15673
0 24
cfn=(744)
calls=4 0 
0 4184
0 44
cfn=(736)
calls=4 0 
0 20513
0 4
0 3072200
0 819136
0 28
0 3072152
0 819136
0 28
0 3072144
0 767940
0 28
0 48
cfn=(736)
calls=4 0 
0 9212
0 24
cfn=(744)
calls=4 0 
0 8592
0 24
cfn=(744)
calls=4 0 
0 2924
0 24
cfn=(744)
calls=4 0 
0 41286
0 28
cfn=(744)
calls=4 0 
0 8461
0 44
cfn=(736)
calls=4 0 
0 114000
0 16
0 36
cfn=(736)
calls=4 0 
0 16468
0 44
cfn=(736)
calls=4 0 
0 24276
0 20
cfn=(744)
calls=4 0 
0 2592
0 20
cfn=(744)
calls=4 0 
0 2606
0 24
cfn=(744)
calls=4 0 
0 2432
0 40
cfn=(744)
calls=4 0 
0 18698
0 24
cfn=(744)
calls=4 0 
0 6881
0 24
cfn=(744)
calls=4 0 
0 2566
0 24
cfn=(744)
calls=4 0 
0 2684
0 24
cfn=(744)
calls=4 0 
0 7878
0 44
cfn=(736)
calls=4 0 
0 114000
0 44
cfn=(736)
calls=4 0 
0 20280
0 36
cfn=(736)
calls=4 0 
0 119972
0 36
0 819192
0 819200
0 3686400
0 4505556
cfn=(1046) for_adjustr
calls=409596 0 
0 28978924
0 4505556
cfn=(736)
calls=409596 0 
0 657811176
0 2457576
cfn=(744)
calls=409596 0 
0 294526056
0 6553536
cfn=(744)
calls=409596 0 
0 795562556
0 4095960
cfn=(744)
calls=409596 0 
0 740139972
0 12287880
cfn=(744)
calls=409596 0 
0 709420272
0 4505556
cfn=(744)
calls=409596 0 
0 713005730
0 4505556
cfn=(744)
calls=409596 0 
0 713024961
0 6553536
cfn=(744)
calls=409596 0 
0 709420272
0 4505556
cfn=(744)
calls=409596 0 
0 569748036
0 4095960
0 36
cfn=(1046)
calls=4 0 
0 276
0 44
cfn=(736)
calls=4 0 
0 6768
0 24
cfn=(744)
calls=4 0 
0 2924
0 64
cfn=(744)
calls=4 0 
0 7652
0 40
cfn=(744)
calls=4 0 
0 7228
0 120
cfn=(744)
calls=4 0 
0 6928
0 44
cfn=(744)
calls=4 0 
0 6892
0 44
cfn=(744)
calls=4 0 
0 6964
0 64
cfn=(744)
calls=4 0 
0 6928
0 44
cfn=(744)
calls=4 0 
0 5564
0 44
cfn=(744)
calls=4 0 
0 6964
0 819192
cfn=(744)
calls=409596 0 
0 713106636
0 5324800
0 4096000
cfn=(1046)
calls=409600 0 
0 28467200
0 4505600
cfn=(736)
calls=409600 0 
0 827801543
0 2048000
cfn=(744)
calls=409600 0 
0 296900216
0 4915200
cfn=(744)
calls=409600 0 
0 775409920
0 2457600
cfn=(744)
calls=409600 0 
0 776791040
0 2457600
cfn=(744)
calls=409600 0 
0 775372800
0 4915200
cfn=(744)
calls=409600 0 
0 816051507
0 4096000
cfn=(744)
calls=409600 0 
0 820019200
0 9830400
cfn=(744)
calls=409600 0 
0 778240000
0 4505600
cfn=(744)
calls=409600 0 
0 781811279
0 4505600
cfn=(744)
calls=409600 0 
0 780197997
0 6553600
cfn=(744)
calls=409600 0 
0 776601600
0 4505600
cfn=(744)
calls=409600 0 
0 584499200
0 4505600
cfn=(744)
calls=409600 0 
0 775372800
0 5324800
0 8
0 36
cfn=(736)
calls=4 0 
0 13120
0 36
cfn=(736)
calls=4 0 
0 112956
0 108
cfn=(736)
calls=4 0 
0 8868
0 24
cfn=(744)
calls=4 0 
0 8368
0 24
cfn=(744)
calls=4 0 
0 3580
0 24
cfn=(744)
calls=4 0 
0 3482
0 24
cfn=(744)
calls=4 0 
0 7886
0 44
cfn=(736)
calls=4 0 
0 108292
0 36
cfn=(736)
calls=4 0 
0 102208
0 8
0 9799716
cfn=(1046)
calls=816640 0 
0 56399200
0 8983040
cfn=(736)
calls=816640 0 
0 1311523840
0 7349760
cfn=(1046)
calls=816640 0 
0 57981440
0 4899840
cfn=(744)
calls=816640 0 
0 489167360
0 4083200
cfn=(744)
calls=816640 0 
0 712937940
0 12249600
cfn=(744)
calls=816640 0 
0 1178885270
0 13882880
cfn=(744)
calls=816640 0 
0 1178885270
0 9799680
cfn=(744)
calls=816640 0 
0 1466446486
0 8983040
cfn=(736)
calls=816640 0 
0 2436037120
0 4899840
cfn=(744)
calls=816640 0 
0 1203143284
0 4899840
cfn=(744)
calls=816640 0 
0 1204942242
0 4899840
cfn=(744)
calls=816640 0 
0 1204942242
0 4899840
cfn=(744)
calls=816640 0 
0 1204942242
0 4899840
cfn=(744)
calls=816640 0 
0 1205035416
0 4899840
cfn=(744)
calls=816640 0 
0 1585589656
0 4083200
0 96
cfn=(736)
calls=4 0 
0 13120
0 36
cfn=(736)
calls=4 0 
0 112956
0 88
cfn=(736)
calls=4 0 
0 8868
0 24
cfn=(744)
calls=4 0 
0 8368
0 24
cfn=(744)
calls=4 0 
0 3580
0 24
cfn=(744)
calls=4 0 
0 3482
0 24
cfn=(744)
calls=4 0 
0 7886
0 44
cfn=(736)
calls=4 0 
0 233620
0 4505624
0 2047980
cfn=(1046)
calls=409596 0 
0 28364524
0 4505556
cfn=(736)
calls=409596 0 
0 658630368
0 2047980
cfn=(744)
calls=409596 0 
0 313777068
0 1638384
cfn=(744)
calls=409596 0 
0 52428288
0 7372728
cfn=(744)
calls=409596 0 
0 763433443
0 2457576
cfn=(744)
calls=409596 0 
0 766354116
0 2457576
0 2457576
cfn=(744)
calls=409596 0 
0 52428288
0 4505556
cfn=(744)
calls=409596 0 
0 728173758
0 2457576
cfn=(744)
calls=409596 0 
0 575072784
0 2457576
0 5734344
0 20
cfn=(1046)
calls=4 0 
0 276
0 44
cfn=(736)
calls=4 0 
0 6776
0 20
cfn=(744)
calls=4 0 
0 3112
0 16
cfn=(744)
calls=4 0 
0 512
0 72
cfn=(744)
calls=4 0 
0 7384
0 24
cfn=(744)
calls=4 0 
0 7484
0 24
0 24
cfn=(744)
calls=4 0 
0 512
0 44
cfn=(744)
calls=4 0 
0 5616
0 24
cfn=(744)
calls=4 0 
0 5616
0 24
0 60
cfn=(744)
calls=4 0 
0 8508
0 819192
cfn=(744)
calls=409596 0 
0 893704509
0 1638400
0 36
cfn=(736)
calls=4 0 
0 112956
0 12
0 32

fn=(390)
0 12
0 2
0 2
0 8
0 42
0 9
0 9
0 12

fn=(396)
0 62
cob=(3)
cfi=(3)
cfn=(402) uname
calls=2 0 
0 10
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1143
0 8
cob=(3)
cfi=(3)
cfn=(408)
calls=2 0 
0 17
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1419
0 34
cob=(3)
cfi=(3)
cfn=(414) atol
calls=2 0 
0 226
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1095
0 8
0 168
0 4
0 8
0 2
0 48
0 10
0 42
0 5
0 10
0 42
0 7
0 3
0 4
0 4
0 4
0 4
0 130
0 36
cob=(3)
cfi=(3)
cfn=(414)
calls=2 0 
0 226
0 4

fn=(442)
0 27884960
0 3485620
cob=(3)
cfi=(3)
cfn=(448)
calls=3485620 0 
0 352128872
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1199
0 24399340
0 17428100

fn=(596)
0 480
0 840
0 75
cfn=(598) for__inquire_args
calls=15 0 
0 3796
0 60
0 30
0 78
0 300
0 945
0 345
0 1095
0 165
cfn=(602) for__compute_filename
calls=15 0 
0 7335
0 210
cfn=(620)
calls=15 0 
0 453
0 45
0 8
cfn=(972)
calls=1 0 
0 339
0 50
cfn=(972)
calls=10 0 
0 5788
0 33
0 66
0 697
0 408
0 524
0 2
0 22
0 20
cfn=(714)
calls=10 0 
0 1810
0 30
0 42
0 11
0 60
cfn=(636)
calls=15 0 
0 834
0 120
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 45
0 60
0 45
0 45
0 45
0 45
0 30
0 4
cfn=(714)
calls=1 0 
0 240
0 2
0 251
cfn=(628)
calls=14 0 
0 26740
0 70
0 56
cob=(3)
cfi=(3)
cfn=(634) _xstat
calls=14 0 
0 156
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1227
0 28
0 40
0 80

fn=(676)
0 228
0 24
0 192
cfn=(602)
calls=12 0 
0 9720
0 180
cfn=(620)
calls=12 0 
0 372
0 72
cfn=(442)
calls=12 0 
0 2613
0 24
0 336
0 12
cfn=(678) __intel_sse2_strcpy
calls=12 0 
0 1151
0 60
0 48
0 8
0 3
0 7
0 36
0 72
0 96
0 4
0 2
0 22
0 60
0 36
0 4
0 96
0 48
0 24
0 36
0 132
0 120
0 24
0 60
0 48
0 60
0 60
0 24
0 96
0 120
0 48
0 84
0 24
0 24
0 24
0 84
0 60
0 4
0 10
0 4
0 48
0 4
0 10
0 4
0 108
0 84
0 24
0 24
0 108
0 48
0 48
0 48
0 72
0 24
0 24
cob=(3)
cfi=(3)
cfn=(684)
calls=12 0 
0 84
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1222
0 36
0 22
cob=(3)
cfi=(3)
cfn=(684)
calls=11 0 
0 77
0 22
0 121
cfn=(478)
calls=12 0 
0 8424
0 24
0 72
0 24
0 4
0 4
cfn=(478)
calls=1 0 
0 716
0 2
0 24
0 20
0 30
0 62
0 108
0 48
0 36
0 48
cob=(3)
cfi=(3)
cfn=(698) isatty
calls=12 0 
0 276
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1172
0 24
0 3
0 108
0 24
0 48
cob=(3)
cfi=(3)
cfn=(706) _fxstat
calls=12 0 
0 120
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1330
0 48
0 120
cob=(3)
cfi=(3)
cfn=(698)
calls=12 0 
0 276
0 24
0 22
0 49
cfn=(708) for__add_to_lf_table
calls=12 0 
0 15598
0 84
0 36
0 24
0 50
0 40
0 4
0 24
0 22
0 22
0 33
cfn=(478)
calls=11 0 
0 7909
0 2
cfn=(478)
calls=1 0 
0 719
0 24
0 36
0 48
0 60
0 60
0 24
0 24
0 108
0 96
cfn=(442)
calls=12 0 
0 2694
0 168
0 24
0 66
cfn=(712)
calls=11 0 
0 1078
0 22
0 60
cfn=(442)
calls=12 0 
0 4814
0 24
0 144
0 20
0 230
cob=(4)
cfi=(4)
cfn=(690)
calls=10 0 
0 70
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1008
0 30
0 24
0 120
0 96
0 4
0 48
0 24
0 228
0 70
0 6
0 6

fn=(718)
0 28656481
0 20468855
0 24562698

fn=(722)
0 49240867
0 15151036
0 60604144
cfn=(724) for__characterize_LUB_buffer
calls=3787759 0 
0 496196010
0 18938795
0 9
0 8
0 11363276
0 7575516
0 2
0 34089844
cob=(4)
cfi=(4)
cfn=(732)
calls=1 0 
0 7
0 6
cob=(4)
cfi=(4)
cfn=(732)
calls=1 0 
0 7
0 4

fn=(738)
0 231374108
0 132213776
0 198320664
cfn=(740) _intel_fast_memset.J
calls=33053444 0 
0 1656304365

fn=(832) rs_find_field
0 81
0 36
cfn=(834)
calls=9 0 
0 369
0 108
0 4
0 6
0 35
0 63
0 77
0 18
0 90
cfn=(836) rs_get_field
calls=9 0 
0 2523

fn=(844)
0 6682293
cfn=(846) _intel_fast_memmove.P
calls=742477 0 
0 11880050

fn=(856)
0 88978
0 53468
0 300
0 800
cfn=(792)
calls=100 0 
0 47270
0 800
0 400
0 72170

fn=(868)
0 6737588
0 5493986
0 6032180
0 975644
0 487818
0 1650496
0 825246
0 3016075
0 1225016
0 1837515
0 14700120
0 34121278
0 39798324
0 7350060
0 9
0 4881478

fn=(888) for_allocate
0 52350765
0 13960204
0 13960204
0 34900510
0 17450255
0 6980102
0 13960204
0 10470153
0 10470153
cfn=(884)
calls=3490051 0 
0 788867742
0 17450255
0 13960204
0 6980102
0 31410459

fn=(898) for_system_clock_count
0 1696
cob=(3)
cfi=(3)
cfn=(424)
calls=212 0 
0 2544
0 4664
0 1060

fn=(904) read_main_input_file_
0 17
cfn=(736)
calls=1 0 
0 3800
0 9
cfn=(736)
calls=1 0 
0 3824
0 6
0 11
0 20
0 2
0 10
0 1280028
0 5760036
0 1280011
0 5760036
0 3
0 960006
0 1280008
0 5440034
0 3840033
0 5120032
0 960006
0 960006
0 21
cfn=(770)
calls=1 0 
0 1975
0 407
cfn=(770)
calls=37 0 
0 69056
0 228
cfn=(810)
calls=38 0 
0 38519
0 342
0 300
0 104
cfn=(906)
calls=8 0 
0 1424
0 120
cfn=(906)
calls=30 0 
0 5370
0 228
cfn=(850)
calls=38 0 
0 2204
0 266
cfn=(850)
calls=38 0 
0 2225
0 152
0 9
cfn=(736)
calls=1 0 
0 4040
0 24
cfn=(770)
calls=1 0 
0 1703
0 44
cfn=(770)
calls=4 0 
0 6820
0 60
0 3
0 10
0 40
cfn=(736)
calls=4 0 
0 10900
0 24
cfn=(744)
calls=4 0 
0 3256
0 328
cfn=(852)
calls=4 0 
0 239517
0 96
0 29
0 29
0 29
0 61
0 84
0 20
0 20
0 20
0 640004
0 640004
0 640004
0 640004
0 4
0 3
0 7
0 5
0 2
0 114
0 222
cfn=(850)
calls=37 0 
0 2167
0 74
0 1
cfn=(928) read_time_discretization_data_
calls=1 0 
0 79301
0 1
0 216
cfn=(850)
calls=36 0 
0 2109
0 72
0 9
cfn=(736)
calls=1 0 
0 4144
0 43
cfn=(852)
calls=1 0 
0 43690
0 8
0 4
0 4
0 22
cfn=(770)
calls=1 0 
0 1711
0 10
0 2
0 210
cfn=(850)
calls=35 0 
0 2051
0 70
0 1
cfn=(930)
calls=1 0 
0 144376
0 1
0 204
cfn=(850)
calls=34 0 
0 1993
0 68
0 9
cfn=(736)
calls=1 0 
0 4160
0 62
cfn=(852)
calls=1 0 
0 67929
0 18
0 6
0 12
cfn=(770)
calls=1 0 
0 1719
0 11
0 198
cfn=(850)
calls=33 0 
0 1935
0 66
0 1
cfn=(932)
calls=1 0 
0 153520
0 1
0 192
cfn=(850)
calls=32 0 
0 1862
0 64
0 192
cfn=(850)
calls=32 0 
0 1856
0 64
0 192
cfn=(850)
calls=32 0 
0 1883
0 64
0 1
cfn=(940)
calls=1 0 
0 29178
0 1
0 186
cfn=(850)
calls=31 0 
0 1819
0 62
0 20
cfn=(852)
calls=1 0 
0 26653
0 2
0 180
cfn=(850)
calls=30 0 
0 1761
0 60
0 1
cfn=(908)
calls=1 0 
0 209002
0 1
0 174
cfn=(850)
calls=29 0 
0 1682
0 58
0 174
cfn=(850)
calls=29 0 
0 1703
0 58
0 532
cfn=(906)
calls=28 0 
0 4964
0 28
0 6
0 17
cfn=(736)
calls=1 0 
0 19484
0 5
cfn=(744)
calls=1 0 
0 20475
0 22
cfn=(736)
calls=1 0 
0 2190
0 324
cfn=(736)
calls=27 0 
0 59150
0 112
0 9
cfn=(736)
calls=1 0 
0 17421
0 10
0 26
0 1
0 2

fn=(960)
0 94
0 16
0 24
0 8
cfn=(746)
calls=2 0 
0 212
0 4
cfn=(746)
calls=2 0 
0 116
0 20
0 20
0 8
0 10
0 10
cfn=(812)
calls=2 0 
0 114
0 26
0 6
0 16
0 16
0 58
0 8
0 18
0 10
0 12
0 4
0 54
0 20
cfn=(814)
calls=2 0 
0 206
0 16
cfn=(640)
calls=2 0 
0 44
0 20
0 4
0 14
0 8
cfn=(640)
calls=2 0 
0 198
0 24
cfn=(714)
calls=2 0 
0 1358
0 4
0 34
0 4
0 14

fn=(980) for_deallocate
0 117
0 36
0 36
0 54
0 27
0 18
0 18
cfn=(922)
calls=9 0 
0 1053
0 9
0 36
0 18
0 81

fn=(1004)
0 5
0 13
cfn=(1008)
calls=1 0 
0 72
0 6
cfn=(1006) __libm_feature_flag_init
calls=1 0 
0 39
0 6

fn=(1010)
0 165886398
0 10
cfn=(736)
calls=1 0 
0 3763
0 73727288
0 608250126

fn=(1028)
0 1691
0 2403
0 511769580
0 82996416
0 156838470
0 125316450
0 30085916
0 2718640725
0 222733803
0 209611643
0 54000928
0 36000856
0 36453332
0 169812
0 45283645
0 63624498
0 1243330
0 18227022
0 3293
0 267
0 4806
cfn=(1030)
calls=89 0 
0 4260289825
0 2937
0 979
cfn=(1032)
calls=89 0 
0 716331897
0 1157
0 1246
0 445
0 1780
0 17085330
0 979
0 3204
0 979
cfn=(1034) dslui_
calls=89 0 
0 913579393
0 4806
cfn=(1036) isdcgs_
calls=89 0 
0 1115212702
0 1780
cfn=(622)
calls=89 0 
0 13674761
0 202839
0 595874880
0 3104
0 52768
0 93120
0 34144
0 15075
0 36180
0 12060
0 144720
0 1659196710
0 84420
0 217080
0 183915
0 887477310
0 3015
0 33165
0 105525
0 3015
0 27936
cfn=(1032)
calls=3104 0 
0 24983080992
0 34144
cfn=(1034)
calls=3104 0 
0 31862339616
0 31040
0 40352
0 27936
0 595921440
0 49664
0 111744
0 43456
0 40352
0 9312
0 86912
0 913674816
0 34144
0 130368
0 15520
0 40352
0 12416
0 58976
0 595874880
0 31040
0 111744
0 15520
0 37248
0 12416
0 74496
0 754774848
0 52768
0 130368
0 3104
0 31040
cfn=(1032)
calls=3104 0 
0 24983080992
0 34144
cfn=(1034)
calls=3104 0 
0 31862339616
0 180032
0 754774848
0 3104
0 49664
0 108640
0 155200
cfn=(1036)
calls=3104 0 
0 3520036598
0 9312
0 12060
0 2525719252
0 445467606
0 2047967163
0 219329820
0 154938320
0 27836530
0 219329820
0 146219880
0 146220236
0 36554970
0 191155535
0 145359784
0 36453688
0 127589955
0 36022928
0 54034392
0 218496424
0 272266130
0 163017918
0 116325759
0 22614099
0 244704788
0 136021815
0 112551180
0 161765154
0 82021599
0 356
0 91135110
0 27340266
0 2937
0 712
0 27339732
0 178

fn=(1036)
0 60667
0 34144
0 712
cfn=(1034)
calls=89 0 
0 913579571
0 2136
0 100249155
0 267
0 1335
0 40352
0 869
0 72
0 110
0 240
0 168
0 3596578366
0 9579
0 168
0 90
0 120
0 15965
0 22351
0 979
cfn=(736)
calls=89 0 
0 237996
0 534
cfn=(744)
calls=89 0 
0 213625
0 31930
cfn=(736)
calls=3193 0 
0 6713778
0 19158
cfn=(744)
calls=3193 0 
0 5406062
0 22351
cfn=(744)
calls=3193 0 
0 5421478
0 22351
cfn=(744)
calls=3193 0 
0 6507823
0 54281
0 72
0 445

fn=(1056)
0 14
0 381
0 14
0 533
0 6

fn=(392)
0 22
0 4
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 687
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1207
0 2
0 2
0 1
cfn=(432) for__reentrancy_init
calls=1 0 
0 17527
0 2
0 2
cfn=(440)
calls=1 0 
0 30
0 3
0 4
cfn=(442)
calls=1 0 
0 62580
0 8
0 2
cfn=(478)
calls=1 0 
0 715
0 8
cob=(3)
cfi=(3)
cfn=(484) sigaltstack
calls=1 0 
0 5
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1342
0 4
cob=(3)
cfi=(3)
cfn=(490) sigemptyset
calls=1 0 
0 22
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1334
0 4
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 68
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1136
0 4
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 68
0 4
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 68
0 5
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 68
0 4
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 68
0 4
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 106
0 7
cob=(4)
cfi=(4)
cfn=(496)
calls=1 0 
0 106
0 10
cfn=(498)
calls=1 0 
0 5536
0 2
cfn=(528)
calls=1 0 
0 20954
0 2
cfn=(586) for__default_io_sizes_env_init
calls=1 0 
0 2763
0 2
cfn=(588)
calls=1 0 
0 8
0 2
cfn=(478)
calls=1 0 
0 715
0 2
0 2
cfn=(478)
calls=1 0 
0 720
0 2
0 2
cfn=(478)
calls=1 0 
0 715
0 2
0 2
cfn=(590) for_get_hbw_availability
calls=1 0 
0 16
0 17
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 3
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 687
0 1
0 3
cfn=(394)
calls=1 0 
0 6216
0 2

fn=(434) reentrancy_init
0 31
cfn=(436)
calls=2 0 
0 18
0 2080
cfn=(436)
calls=1040 0 
0 9360
0 4168
0 10
0 14
cfn=(438)
calls=1 0 
0 11
0 210
cfn=(438)
calls=105 0 
0 1155
0 424
0 8

fn=(532)
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 8335
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1435
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1020
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1028
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1020
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 988
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1084
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1084
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1158
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1293
0 6
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1339
0 2
0 4
cob=(5)
cfi=(5)
cfn=(538)
calls=1 0 
0 1046
0 7

fn=(586)
0 15
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 3
0 2
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 686
0 3
0 2
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 3
0 2
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 3
0 2
0 11
0 2
0 2
0 2

fn=(590)
0 5
0 7
0 2
0 2

fn=(592)
0 10
cfn=(384)
calls=1 0 
0 118
0 5
cfn=(594)
calls=1 0 
0 10
0 11
cfn=(596)
calls=1 0 
0 6741
0 14
cfn=(638)
calls=1 0 
0 19021
0 1
0 6
cfn=(716)
calls=1 0 
0 1887
0 9
cfn=(736)
calls=1 0 
0 7560
0 11
cfn=(770)
calls=1 0 
0 17287
0 11
cfn=(736)
calls=1 0 
0 36601
0 2
cfn=(826) core_memory_allocation_
calls=1 0 
0 1812241961
0 9
cfn=(736)
calls=1 0 
0 3896
0 11
cfn=(736)
calls=1 0 
0 10927
0 9
cfn=(736)
calls=1 0 
0 4802
0 9
cfn=(736)
calls=1 0 
0 25591
0 9
cfn=(736)
calls=1 0 
0 54862
0 11
cfn=(736)
calls=1 0 
0 6247
0 6
cfn=(744)
calls=1 0 
0 1608
0 6
cfn=(744)
calls=1 0 
0 1553
0 6
cfn=(744)
calls=1 0 
0 1553
0 6
cfn=(744)
calls=1 0 
0 1553
0 6
cfn=(744)
calls=1 0 
0 24369
0 10
cfn=(736)
calls=1 0 
0 4696
0 5
cfn=(744)
calls=1 0 
0 2178
0 6
cfn=(744)
calls=1 0 
0 2178
0 5
cfn=(744)
calls=1 0 
0 2123
0 6
cfn=(744)
calls=1 0 
0 25012
0 9
cfn=(736)
calls=1 0 
0 3088
0 9
cfn=(736)
calls=1 0 
0 3654
0 9
cfn=(736)
calls=1 0 
0 3992
0 6
0 6
cfn=(736)
calls=1 0 
0 1905
0 12
cfn=(596)
calls=1 0 
0 3035
0 11
cfn=(736)
calls=1 0 
0 1905
0 8
0 5
cfn=(638)
calls=1 0 
0 6756
0 11
cfn=(596)
calls=1 0 
0 3041
0 11
cfn=(736)
calls=1 0 
0 1913
0 12
cfn=(638)
calls=1 0 
0 6860
0 1
0 11
cfn=(596)
calls=1 0 
0 3039
0 13
cfn=(596)
calls=1 0 
0 3041
0 11
cfn=(736)
calls=1 0 
0 1913
0 12
cfn=(638)
calls=1 0 
0 8607
0 1
0 11
cfn=(596)
calls=1 0 
0 3035
0 11
cfn=(736)
calls=1 0 
0 1905
0 8
0 5
cfn=(638)
calls=1 0 
0 6971
0 11
cfn=(596)
calls=1 0 
0 3041
0 11
cfn=(736)
calls=1 0 
0 1913
0 8
0 5
cfn=(638)
calls=1 0 
0 6977
0 11
cfn=(736)
calls=1 0 
0 3880
0 3
cfn=(898)
calls=1 0 
0 47
0 3
cfn=(900) for_system_clock_rate
calls=1 0 
0 10
0 8
cfn=(902) allocate_temp_input_memory_
calls=1 0 
0 16800
0 3
cfn=(904)
calls=1 0 
0 36514550
0 2
0 9
cfn=(736)
calls=1 0 
0 3920
0 10
cfn=(736)
calls=1 0 
0 39297
0 6
cfn=(744)
calls=1 0 
0 3034
0 16
cfn=(736)
calls=1 0 
0 3987
0 1
0 9
cfn=(736)
calls=1 0 
0 17714
0 3
cfn=(962)
calls=1 0 
0 3511624634
0 13
cfn=(736)
calls=1 0 
0 2702
0 6
cfn=(744)
calls=1 0 
0 871
0 6
cfn=(744)
calls=1 0 
0 2056
0 6
cfn=(744)
calls=1 0 
0 1108
0 2
cfn=(984) solver_memory_allocation_
calls=1 0 
0 81311
0 3
cfn=(986)
calls=1 0 
0 3147963993
0 2
0 2
cfn=(898)
calls=1 0 
0 47
0 3
cfn=(900)
calls=1 0 
0 10
0 17
cfn=(736)
calls=1 0 
0 5761
0 2
cfn=(988) deallocate_unneeded_arrays_
calls=1 0 
0 15361
0 9
cfn=(736)
calls=1 0 
0 3912
0 9
cfn=(736)
calls=1 0 
0 25053
0 9
cfn=(736)
calls=1 0 
0 2746
0 9
cfn=(736)
calls=1 0 
0 24916
0 11
cfn=(736)
calls=1 0 
0 5704
0 11
cfn=(736)
calls=1 0 
0 2383
0 6
cfn=(744)
calls=1 0 
0 1155
0 6
cfn=(744)
calls=1 0 
0 946
0 6
cfn=(744)
calls=1 0 
0 988
0 6
cfn=(744)
calls=1 0 
0 1723
0 11
cfn=(736)
calls=1 0 
0 3509
0 11
cfn=(736)
calls=1 0 
0 18492
0 6
cfn=(744)
calls=1 0 
0 5209
0 6
cfn=(744)
calls=1 0 
0 8339
0 11
cfn=(736)
calls=1 0 
0 11432
0 11
cfn=(736)
calls=1 0 
0 6454
0 6
cfn=(744)
calls=1 0 
0 6492
0 11
cfn=(736)
calls=1 0 
0 7646
0 6
cfn=(744)
calls=1 0 
0 5027
0 11
cfn=(736)
calls=1 0 
0 7333
0 6
cfn=(744)
calls=1 0 
0 5489
0 6
cfn=(744)
calls=1 0 
0 6259
0 11
cfn=(736)
calls=1 0 
0 11513
0 6
cfn=(744)
calls=1 0 
0 6048
0 6
cfn=(744)
calls=1 0 
0 12236
0 11
cfn=(736)
calls=1 0 
0 7247
0 6
cfn=(744)
calls=1 0 
0 2746
0 11
cfn=(736)
calls=1 0 
0 29932
0 2
cfn=(992) write_eos_specific_info_
calls=1 0 
0 103106
0 2
cfn=(994)
calls=1 0 
0 243555166459
0 2
0 9
cfn=(736)
calls=1 0 
0 25591
0 2
0 2
cfn=(898)
calls=1 0 
0 47
0 3
cfn=(900)
calls=1 0 
0 10
0 17
cfn=(736)
calls=1 0 
0 6089
0 6
cfn=(744)
calls=1 0 
0 3045
0 7
cfn=(744)
calls=1 0 
0 2968
0 6
cfn=(744)
calls=1 0 
0 2246
0 7
cfn=(1052)
calls=1 0 
0 76999
0 9
cfn=(736)
calls=1 0 
0 52415
0 6
cfn=(1042) for_endfile
calls=1 0 
0 839
0 6
cfn=(716)
calls=1 0 
0 834
0 13
cfn=(882)
calls=1 0 
0 348
0 21
0 43612
0 1242
cfn=(770)
calls=138 0 
0 151465
0 276
0 2055
cfn=(810)
calls=137 0 
0 75324
0 14522
0 41100
cfn=(810)
calls=3425 0 
0 2245695
0 6850
0 276
0 1233
cfn=(736)
calls=137 0 
0 156180
0 1918
cfn=(744)
calls=137 0 
0 62198
0 23975
cfn=(744)
calls=3425 0 
0 1557416
0 14248
0 548
cfn=(744)
calls=137 0 
0 59047
0 274
0 27
cfn=(920)
calls=1 0 
0 183
0 5
0 548
cfn=(810)
calls=137 0 
0 30688
0 137

fn=(598)
0 240
cfn=(600) for__key_desc_ret_item
calls=15 0 
0 1050
0 375
0 6
0 70
0 60
0 45
0 180
cfn=(600)
calls=30 0 
0 1605
0 90
0 75

fn=(678)
0 240
0 6
0 4
0 4
0 364
0 630
0 3
0 126
0 540
0 84
0 36
0 48
0 84
0 8
0 60
0 16
0 52

fn=(708)
0 192
cfn=(710)
calls=12 0 
0 102
0 24
0 60
0 110
0 10
0 108
0 60
0 66
0 27
cfn=(650)
calls=1 0 
0 338
0 275
cfn=(650)
calls=11 0 
0 3203
0 72
0 72
cfn=(442)
calls=12 0 
0 2610
0 48
0 120
0 33
cfn=(620)
calls=11 0 
0 301
0 96
cfn=(442)
calls=12 0 
0 2531
0 48
0 60
0 12
cfn=(678)
calls=12 0 
0 1154
0 60
0 36
0 240
cfn=(670)
calls=12 0 
0 2908
0 276
0 33
cfn=(620)
calls=11 0 
0 301
0 12

fn=(716)
0 256
0 128
cfn=(648)
calls=16 0 
0 9294
0 48
0 336
0 64
cfn=(720)
calls=16 0 
0 592
0 320
0 48
cfn=(722)
calls=16 0 
0 2928
0 32
0 48
0 64
cob=(4)
cfi=(4)
cfn=(732)
calls=16 0 
0 112
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1053
0 288
0 32
0 32
0 64
0 144
cfn=(714)
calls=16 0 
0 3030
0 32
0 128

fn=(720)
0 68179860
0 128
0 7575476
0 96
0 44
0 11363214
0 15150952
0 7575476
0 11363214
0 22726428
cob=(4)
cfi=(4)
cfn=(754)
calls=3787738 0 
0 11363214
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1340
0 45452856
cob=(4)
cfi=(4)
cfn=(760)
calls=3787738 0 
0 26514166
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1028
0 37877380
0 41665118
0 60604128

fn=(746)
0 1574249490
0 2990148
0 35975560
0 174429878
0 639247552
0 799059440
0 548743035
0 9824772
0 21504088
0 168763716
0 3940
0 53189774
0 3753
0 79784669
0 130417605
0 104334084
0 2045472
0 2556805
0 2045440
0 5
0 35
0 67485996
0 81033408
0 2043516
0 118484838
0 31473134
0 47833786
0 9799680
0 517082337
0 311857516

fn=(770)
0 4294598
0 3987841
cfn=(648)
calls=306757 0 
0 108899632
0 920271
0 16564878
0 1227028
0 87
0 5
cfn=(772)
calls=1 0 
0 12060
0 2
0 3
0 613512
0 613514
0 3374327
0 7362168
0 1840542
0 613514
0 1840542
0 1227028
0 613514
0 1840542
0 2760813
cfn=(792)
calls=306757 0 
0 92430223
0 613514
0 4
0 3
0 6
0 3
0 7668900
0 1227024
cfn=(810)
calls=306756 0 
0 180836074
0 2454052
cfn=(1048)
calls=1 0 
0 858
0 8

fn=(812)
0 1481102012
0 224293880
0 196257145
0 1290
0 1977
0 602
0 107340
0 9
0 46108106
0 11034726
0 723645695

fn=(826)
0 16
cfn=(736)
calls=1 0 
0 3992
0 11
cfn=(596)
calls=1 0 
0 3043
0 14
cfn=(638)
calls=1 0 
0 6958
0 6
cfn=(716)
calls=1 0 
0 834
0 1
0 10
cfn=(828)
calls=1 0 
0 2437
0 5
cfn=(838)
calls=1 0 
0 89
0 5
cfn=(840)
calls=1 0 
0 85
0 7
cfn=(842)
calls=1 0 
0 52
0 10
cfn=(850)
calls=1 0 
0 56
0 7
cfn=(850)
calls=1 0 
0 56
0 7
cfn=(850)
calls=1 0 
0 61
0 23
0 6
cfn=(850)
calls=1 0 
0 56
0 7
cfn=(850)
calls=1 0 
0 56
0 7
cfn=(850)
calls=1 0 
0 61
0 8
0 3
0 22
cfn=(852)
calls=1 0 
0 36485
0 15
cfn=(838)
calls=1 0 
0 77
0 17
cfn=(876)
calls=1 0 
0 4145
0 20
0 5
cfn=(838)
calls=1 0 
0 77
0 11
cfn=(878)
calls=1 0 
0 4481
0 5
0 5
0 25
cfn=(852)
calls=1 0 
0 40843
0 17
cfn=(850)
calls=1 0 
0 67
0 7
cfn=(850)
calls=1 0 
0 58
0 7
cfn=(850)
calls=1 0 
0 64
0 2
0 5
0 14
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 363
0 3
0 17
cfn=(880)
calls=1 0 
0 81
0 1170009
cfn=(880)
calls=130001 0 
0 10530081
0 1300020
cfn=(888)
calls=130002 0 
0 38483238
0 390006
0 650010
0 4
0 9
cfn=(880)
calls=1 0 
0 122
0 17
cfn=(882)
calls=1 0 
0 353
0 37
cfn=(880)
calls=1 0 
0 122
0 11
0 4
0 18
cfn=(880)
calls=2 0 
0 162
0 2080016
cfn=(880)
calls=260002 0 
0 21060162
0 3640056
cfn=(888)
calls=260004 0 
0 76969792
0 6500100
cfn=(880)
calls=260004 0 
0 21060324
0 2860044
cfn=(880)
calls=260004 0 
0 31720488
0 3900060
cfn=(888)
calls=260004 0 
0 76969254
0 8320128
cfn=(880)
calls=260004 0 
0 31720488
0 2080032
cfn=(880)
calls=260004 0 
0 21060324
0 3900060
cfn=(888)
calls=260004 0 
0 76969792
0 780012
0 1300020
0 10
0 3
0 6
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 539
0 3
0 31
0 20
0 162495
0 12
0 34
0 20
0 162495
0 15
0 14
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 3
0 16
cfn=(880)
calls=1 0 
0 81
0 2560008
cfn=(880)
calls=320001 0 
0 25920081
0 3200020
cfn=(888)
calls=320002 0 
0 94730814
0 7040044
cfn=(880)
calls=320002 0 
0 25920162
0 2880018
cfn=(880)
calls=320002 0 
0 25920162
0 3200020
cfn=(888)
calls=320002 0 
0 94731352
0 960006
0 3840024
cfn=(880)
calls=320002 0 
0 39040244
0 3200020
cfn=(888)
calls=320002 0 
0 94731094
0 9280058
cfn=(880)
calls=320002 0 
0 39040244
0 1600010
0 17
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 539
0 25
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 25
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 25
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 3
0 12
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 285
0 25
cfn=(880)
calls=1 0 
0 81
0 3
0 11
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 285
0 25
cfn=(880)
calls=1 0 
0 81
0 85
0 160
cfn=(736)
calls=16 0 
0 44363
0 16
0 24
0 84
0 13
cfn=(770)
calls=1 0 
0 1711
0 43
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 24
cfn=(880)
calls=1 0 
0 81
0 1
0 6080038
cfn=(880)
calls=320002 0 
0 25920162
0 2560016
cfn=(880)
calls=320002 0 
0 25920162
0 3200020
cfn=(888)
calls=320002 0 
0 94731352
0 7040044
cfn=(880)
calls=320002 0 
0 25920162
0 320002
0 22
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 25
cfn=(880)
calls=1 0 
0 81
0 1
0 22
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 24
cfn=(880)
calls=1 0 
0 81
0 1
0 5720088
cfn=(880)
calls=260004 0 
0 21060324
0 2080032
cfn=(880)
calls=260004 0 
0 21060324
0 3900060
cfn=(888)
calls=260004 0 
0 76969792
0 6500100
cfn=(880)
calls=260004 0 
0 21060324
0 2080032
cfn=(880)
calls=260004 0 
0 21060324
0 3900060
cfn=(888)
calls=260004 0 
0 76969534
0 6500100
cfn=(880)
calls=260004 0 
0 21060324
0 2080032
cfn=(880)
calls=260004 0 
0 21060324
0 3900060
cfn=(888)
calls=260004 0 
0 76969792
0 6500100
cfn=(880)
calls=260004 0 
0 21060324
0 2080032
cfn=(880)
calls=260004 0 
0 21060324
0 3900060
cfn=(888)
calls=260004 0 
0 76969792
0 6760104
cfn=(880)
calls=260004 0 
0 21060324
0 520008
0 2600040
cfn=(880)
calls=130002 0 
0 10530162
0 1170018
cfn=(880)
calls=130002 0 
0 10530162
0 1300020
cfn=(888)
calls=130002 0 
0 38488905
0 2860044
cfn=(880)
calls=130002 0 
0 10530162
0 1040016
cfn=(880)
calls=130002 0 
0 10530162
0 1300020
cfn=(888)
calls=130002 0 
0 38483091
0 2860044
cfn=(880)
calls=130002 0 
0 10530162
0 260004
0 23
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 353
0 25
cfn=(880)
calls=1 0 
0 81
0 12
cfn=(880)
calls=1 0 
0 122
0 17
cfn=(882)
calls=1 0 
0 353
0 33
cfn=(880)
calls=1 0 
0 122
0 1

fn=(846)
0 742477
cfn=(848)
calls=742477 0 
0 11137573

fn=(890)
0 44249904
0 342
0 7374870
0 27520648
0 142121421
0 7374642
0 570
0 7374984
0 342
0 4916656
0 854622
0 1478772
0 5656042
0 4177270
0 4177270
0 4062034
0 9093926
0 6654830
0 369693
0 739386
0 854622
0 288090
0 4062034
0 18064630
0 2031017
0 4062034

fn=(910)
0 17168816
0 15942472
cfn=(912)
calls=1226344 0 
0 555395319
0 4905376
0 6131720

fn=(944)
0 32
0 52
0 12
0 16
0 12
0 3
cfn=(946)
calls=1 0 
0 1473
0 15
cob=(4)
cfi=(4)
cfn=(560)
calls=4 0 
0 66
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1452
0 12
0 36
cfn=(442)
calls=1 0 
0 218
0 2
0 14
0 24
0 3
cob=(4)
cfi=(4)
cfn=(572)
calls=1 0 
0 31
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1435
0 2

fn=(968)
0 43008000
0 27648000
0 43008000
0 3072000
0 30720000
0 3072000
0 3072000
0 9216000
0 4608000
0 1368000
0 1641600
0 518400
0 3670400
0 6768000
0 9216000
0 147200
0 36937600
0 76800000
0 144000
0 4608000
0 4636800
0 633600
0 172800
0 2976000
0 4464000
0 20832000
0 17856000
0 5184000
0 17112000
0 86366400
0 207216000
0 359174400
0 151958400
0 67200
0 4608000
0 33792000
0 33792000
0 230400
0 2956800
0 102400
0 147200
0 2131200
0 147200
0 147200
0 10752000

fn=(984)
0 16
cfn=(736)
calls=1 0 
0 3896
0 9
cfn=(736)
calls=1 0 
0 1949
0 23
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 485
0 26
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 354
0 26
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 355
0 26
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 355
0 26
cfn=(880)
calls=1 0 
0 81
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 355
0 27
cfn=(880)
calls=1 0 
0 81
0 13
cfn=(736)
calls=1 0 
0 2779
0 1
0 12
cfn=(736)
calls=1 0 
0 2779
0 1
0 12
cfn=(736)
calls=1 0 
0 2779
0 1
0 12
cfn=(736)
calls=1 0 
0 2779
0 1
0 12
cfn=(736)
calls=1 0 
0 2779
0 1
0 11
cfn=(736)
calls=1 0 
0 28174
0 6
cfn=(744)
calls=1 0 
0 2913
0 6
cfn=(744)
calls=1 0 
0 2285
0 6
cfn=(744)
calls=1 0 
0 25082
0 8

fn=(990)
0 25
0 2
0 8
0 5
0 11
0 7
0 2

fn=(998)
0 276477330
0 10
cfn=(736)
calls=1 0 
0 3936
0 1400818472
cfn=(1000)
calls=18431822 0 
0 4147160019
0 405500084

fn=(1022)
0 7920
0 714
0 4
0 1200
0 720
0 7080
0 49151040
0 1806318360
0 589818240
0 49151520
0 49151520
0 983030400
0 49151520
0 73727280
0 49151520
0 49151520
0 73727280
0 417787920
0 1400818320
0 360
0 367489920
0 489984000
0 293990400
0 318496680
0 480
0 122496000
0 240
0 4200
0 5280
0 59894640
0 2880
0 2880
0 1680
0 98305680
0 4200
0 24575040
0 24575520
0 24575520
0 960
0 2400
0 43006320
0 36851172
0 10090
0 42994182
0 11124
0 1320
0 2040
0 1175962800
0 97995840
0 171494880
0 48998400
0 97996800
0 1249459200
0 710476800
0 1322956800
0 2645913600
0 146995200
0 10216166400
0 5512320000
0 1763942400
0 1494451200
0 97996800
0 2449920000
0 587980800
0 1249459200
0 244992000
0 73497600
0 97996800
0 49151520
0 73727280
0 49151520
0 466939440
0 49151520
0 48998400
0 832972800
0 293987520
0 710469840
0 1200
0 8400
0 480
0 73507200
0 2670386640
0 122496000
0 480
0 1680
0 293990400
0 1616947200
0 661478400
0 515543930
0 22373085
0 293990400
0 416486400
0 24575760
0 24575760
0 12697476
0 12697476
0 6348738
0 38092428
0 6348738
0 240
0 960
0 600
0 240
0 240
0 840
0 600
0 240
0 10
cfn=(736)
calls=1 0 
0 4067
0 24
cfn=(880)
calls=1 0 
0 81
0 16
cfn=(882)
calls=1 0 
0 661
0 24
cfn=(880)
calls=1 0 
0 81
0 16
cfn=(882)
calls=1 0 
0 379
0 24
cfn=(880)
calls=1 0 
0 81
0 16
cfn=(882)
calls=1 0 
0 418
0 24
cfn=(880)
calls=1 0 
0 81
0 16
cfn=(882)
calls=1 0 
0 353
0 23
cfn=(918)
calls=1 0 
0 38
0 64
0 1301469
0 3
0 48
0 59
0 1301469
0 9
0 25
cfn=(918)
calls=1 0 
0 38
0 47
0 893165
0 5
0 11

fn=(1040)
0 20
0 10
cfn=(736)
calls=1 0 
0 3944
0 44
cfn=(736)
calls=4 0 
0 12096
0 24
cfn=(744)
calls=4 0 
0 11104
0 24
cfn=(716)
calls=4 0 
0 5664
0 44
cfn=(736)
calls=4 0 
0 10792
0 24
cfn=(744)
calls=4 0 
0 8563
0 60
cfn=(736)
calls=4 0 
0 8600
0 3686364
cfn=(736)
calls=409596 0 
0 880631400
0 1228800
0 40
cfn=(736)
calls=4 0 
0 8488
0 44
cfn=(736)
calls=4 0 
0 17980
0 32
cfn=(744)
calls=4 0 
0 7190
0 32
cfn=(744)
calls=4 0 
0 10210
0 28
cfn=(744)
calls=4 0 
0 16110
0 44
cfn=(736)
calls=4 0 
0 8488
0 24
cfn=(1042)
calls=4 0 
0 5684
0 16

fn=(1042)
0 80
0 40
cfn=(648)
calls=5 0 
0 3631
0 15
0 120
0 75
0 20
cfn=(720)
calls=5 0 
0 185
0 15
0 40
0 10
0 10
0 15
cfn=(722)
calls=5 0 
0 915
0 30
cob=(4)
cfi=(4)
cfn=(732)
calls=5 0 
0 35
0 20
cob=(3)
cfi=(3)
cfn=(768)
calls=5 0 
0 25
0 30
0 75
cfn=(714)
calls=5 0 
0 1087
0 10
0 40

fn=(1072)
0 5
0 13
0 3
0 4
0 6
cfn=(1074) reentrancy_cleanup
calls=1 0 
0 25832
0 7

fn=(1074)
0 2105
cfn=(1076)
calls=1042 0 
0 15266
0 5210
0 8
0 11
cfn=(1078)
calls=1 0 
0 16
0 210
cfn=(1078)
calls=105 0 
0 2471
0 530
0 5

fn=(432)
0 5
0 13
0 3
0 4
0 6
cfn=(434)
calls=1 0 
0 17489
0 7

fn=(600)
0 2025
0 400
0 350
0 385
0 90
0 200
0 175
0 250
0 600

fn=(602)
0 756
0 104
0 52
0 2
0 325
cob=(3)
cfi=(3)
cfn=(608) __memcpy_chk
calls=25 0 
0 662
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1363
0 125
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 3
0 20
0 3
cfn=(620)
calls=1 0 
0 11
0 193
0 135
0 125
0 125
0 125
0 525
0 125
cfn=(622)
calls=25 0 
0 925
0 325
cfn=(622)
calls=25 0 
0 821
0 100
0 4
0 10
cfn=(622)
calls=2 0 
0 66
0 10
0 175
0 357
cob=(3)
cfi=(3)
cfn=(684)
calls=1 0 
0 7
0 3
cob=(3)
cfi=(3)
cfn=(778) ttyname
calls=1 0 
0 337
0 7
0 100
cob=(3)
cfi=(3)
cfn=(618) getcwd
calls=25 0 
0 700
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1172
0 450
cfn=(620)
calls=25 0 
0 725
0 50
0 3
0 2
cob=(3)
cfi=(3)
cfn=(198)
calls=1 0 
0 673
0 5
cob=(3)
cfi=(3)
cfn=(684)
calls=1 0 
0 7
0 3
cob=(3)
cfi=(3)
cfn=(778)
calls=1 0 
0 34
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1202
0 4
cob=(3)
cfi=(3)
cfn=(784) getpid
calls=1 0 
0 5
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1173
0 16
cob=(3)
cfi=(3)
cfn=(504)
calls=1 0 
0 686
0 10
0 6
cob=(3)
cfi=(3)
cfn=(790)
calls=1 0 
0 44
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1376
0 4
0 4

fn=(638)
0 180
0 280
0 80
cfn=(600)
calls=10 0 
0 700
0 280
0 80
0 30
0 20
0 120
cfn=(600)
calls=20 0 
0 1120
0 100
0 60
0 40
0 80
cfn=(648)
calls=10 0 
0 18076
0 20
0 220
0 780
0 130
0 110
0 50
0 30
cfn=(676)
calls=10 0 
0 53398
0 20
0 80
cfn=(714)
calls=10 0 
0 1938
0 20
0 30
0 110
0 60
0 20
0 60
cfn=(442)
calls=10 0 
0 1461
0 490
0 40
0 20
0 290
0 20
0 20
0 50
0 20
cfn=(640)
calls=10 0 
0 2110
0 50
0 40

fn=(640)
0 20913690
cob=(3)
cfi=(3)
cfn=(526) free
calls=3485615 0 
0 292795666
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1120
0 13942460
0 17428075

fn=(650)
0 65501056
0 114626848
0 24562896
0 49125792
0 18817612
cob=(4)
cfi=(4)
cfn=(458)
calls=4704403 0 
0 89383657
0 9408806
0 42159334
cfn=(442)
calls=3483229 0 
0 424956048
0 6966458
0 10449687
cob=(4)
cfi=(4)
cfn=(656)
calls=3483229 0 
0 111463328
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1398
0 45281977
cob=(4)
cfi=(4)
cfn=(458)
calls=3483229 0 
0 66181351
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1398
0 6966458
0 24382603

fn=(668)
0 250
0 20
cob=(4)
cfi=(4)
cfn=(458)
calls=10 0 
0 190
0 80
0 20
0 30
cob=(4)
cfi=(4)
cfn=(460)
calls=10 0 
0 150
0 200
0 20
0 90

fn=(724)
0 178024673
cfn=(726)
calls=3787759 0 
0 106057252
0 22726554
0 64391699
0 7575492
0 4
0 11363241
0 11
0 18938738
0 15
0 18938734
0 18938735
0 7575492
0 14
0 41665224
0 132

fn=(736)
0 53017300
0 49230350
cfn=(648)
calls=3786950 0 
0 2152450314
0 11360850
0 90886800
0 140117150
0 81
0 5
cfn=(772)
calls=1 0 
0 8586
0 2
0 3
0 7573898
0 7573900
0 7573900
0 37869500
0 90886800
0 49230350
0 11360850
0 7573900
0 9
0 3
cfn=(934)
calls=1 0 
0 9195
0 2
0 5
0 11360847
0 2
0 56804250
cfn=(722)
calls=3786950 0 
0 693011484
0 11360850
0 18934750
0 37869500
0 7573900
0 30295600
0 49230350
cfn=(738)
calls=3786950 0 
0 417482870
0 26508650
0 7573900
0 49230350
0 15147800
cfn=(744)
calls=3786950 0 
0 3536556042
0 30295600

fn=(740)
0 33053444
cfn=(742)
calls=33053444 0 
0 1623250921

fn=(822)
0 63
0 36
0 6
0 15
0 12
cfn=(738)
calls=3 0 
0 333
0 27

fn=(830)
0 270
0 171
cfn=(812)
calls=9 0 
0 558
0 45
cfn=(812)
calls=9 0 
0 324
0 144
0 36
0 27
0 27
0 171
0 90
cfn=(832)
calls=9 0 
0 3410
0 45
0 54
0 171
0 36
0 27
0 45
cfn=(814)
calls=9 0 
0 827
0 9
0 72
0 63
0 86
cfn=(640)
calls=9 0 
0 1623
0 18
0 18
cfn=(714)
calls=9 0 
0 1107
0 18
0 171

fn=(836)
0 189
0 18
0 36
0 92
0 8
0 90
0 35
0 608
0 662
0 456
0 152
0 49
0 12
0 36
0 80

fn=(880)
0 90741820
0 20940420
0 62821260
0 94221868
0 309325928
0 34900700

fn=(900)
0 2120

fn=(902)
0 14
cfn=(736)
calls=1 0 
0 3792
0 9
cfn=(736)
calls=1 0 
0 1949
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 3
0 9
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 3
0 14
cfn=(882)
calls=1 0 
0 467
0 31
cfn=(736)
calls=1 0 
0 2811
0 1
0 12
cfn=(736)
calls=1 0 
0 2811
0 1
0 12
cfn=(736)
calls=1 0 
0 2811
0 1
0 3
0 2
0 3
0 2
0 28
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 25
cfn=(880)
calls=1 0 
0 81
0 1
0 23
cfn=(880)
calls=1 0 
0 81
0 8
cfn=(880)
calls=1 0 
0 81
0 17
cfn=(882)
calls=1 0 
0 299
0 24
cfn=(880)
calls=1 0 
0 81
0 1

fn=(922)
0 36
0 36
cob=(3)
cfi=(3)
cfn=(526)
calls=18 0 
0 1837

fn=(928)
0 12
cfn=(736)
calls=1 0 
0 4240
0 23
0 63
0 15
0 63
0 43
cfn=(852)
calls=1 0 
0 73049
0 15
0 5
0 10
0 2
0 7
0 6
0 11
cfn=(770)
calls=1 0 
0 1719
0 11
0 7

fn=(1020)
0 55
0 10
cfn=(736)
calls=1 0 
0 3840
0 195
cfn=(736)
calls=5 0 
0 262430
0 55
cfn=(736)
calls=5 0 
0 10959
0 30
cfn=(744)
calls=5 0 
0 5706
0 30
cfn=(744)
calls=5 0 
0 8927
0 35
cfn=(744)
calls=5 0 
0 13379
0 20
0 725
0 4863810
0 135
0 3327790
0 35
0 90
0 270
0 30
cfn=(622)
calls=5 0 
0 768245
0 160
0 4863810
0 95
0 170
0 1215810
0 55
0 180
0 30
cfn=(622)
calls=5 0 
0 768260
0 20
0 150
0 447930
0 70
0 100
0 65
0 75
0 447930
0 70
0 120
0 5
0 15
0 250
0 6911730
0 115
0 55
0 120
0 959850
0 80
0 180
0 165
cfn=(736)
calls=5 0 
0 159154
0 50
cfn=(744)
calls=5 0 
0 77063
0 110
cfn=(736)
calls=5 0 
0 209413
0 50
cfn=(744)
calls=5 0 
0 65188
0 45
cfn=(736)
calls=5 0 
0 244070
0 70
cfn=(738)
calls=5 0 
0 1152255
0 40

fn=(1046)
0 40069120
0 5724160
0 20034560
0 8586240
0 8586240
cfn=(622)
calls=2862080 0 
0 97156960
0 20034560

fn=(360)
0 4
0 2

fn=(374)
0 12
0 16
0 34
0 98
0 36
0 33
0 90
0 4
0 4
0 2

fn=(794)
0 52752
cfn=(710)
calls=3297 0 
0 26376
0 6594
0 9891
cob=(4)
cfi=(4)
cfn=(754)
calls=3297 0 
0 9891
0 26404
0 6602
0 9903
cob=(4)
cfi=(4)
cfn=(800)
calls=3301 0 
0 23107
cob=(1)
cfi=(1)
cfn=(186)
calls=1 0 
0 1037
0 6602
0 23072
0 8
0 36267

fn=(842)
0 8908896
0 1484810
0 2969620
cfn=(844)
calls=742405 0 
0 18560126
0 742405
0 12
cfn=(844)
calls=3 0 
0 144
0 27
cfn=(738)
calls=3 0 
0 112
0 5196856

fn=(860)
0 1080
0 45
0 249
0 416
0 42
0 2
0 1071
0 7
0 276
0 460
0 140
0 540

fn=(870)
0 60
0 39
0 177
0 63
0 40
0 15
0 18
0 96
cfn=(814)
calls=3 0 
0 249
0 30
cfn=(622)
calls=3 0 
0 93
0 3
0 132
0 72
0 3
0 6
cfn=(792)
calls=1 0 
0 441
0 2
0 4
0 8
0 295

fn=(926)
0 504124208
0 324079848
0 504124208
0 36008872
0 342084284
0 19439652
0 19433166
0 29149749
0 6486
0 3246
0 8644
0 2164
0 19433166
0 19433166
0 58299498
0 38866332
0 97165830
0 16215
0 16215
0 259428
0 25944
0 24853830
0 16569220
0 714129
0 48279402
0 965144
0 24863559
0 15644229
0 21510776
0 42486934
0 22911105
0 38070674
0 5202176
0 31009779
0 49727118
0 5377694
0 510384308
0 347798060
0 157323324
0 484121988
0 7533450
0 24863559
0 2
0 3
0 26370251
0 5497240
0 16480718
0 24721065
0 32
0 40
0 16480710
0 23793
0 91857
0 16473912
0 16473912
0 4953680
0 106585060
0 131800200
0 6806
0 33151412
0 116786460
0 1297881964
0 1217736
0 1539097170
0 281016
0 651156495
0 24863559
0 182332766
0 182332766
0 16575706
0 17015
0 432
0 16475
0 263428
0 29655
0 864
0 44239
0 32
0 105369
0 135960
0 47586
0 61182
0 37389
0 6798
0 55
0 10197
0 8400
0 23715
0 93189
0 10197
0 142926
0 27224
0 88478
0 37433
0 22
0 10209
0 11478
0 26754
0 87722
0 19433166
0 68016081
0 275780316
0 9716583
0 19433166
0 62835
0 2242748
0 5377694
0 87292079
0 5377694
0 5377694
0 10209
0 2790
0 10
0 6798
0 3822
0 4
0 6806
0 58014971

fn=(978)
0 17
cfn=(736)
calls=1 0 
0 3952
0 7
cfn=(716)
calls=1 0 
0 834
0 11
cfn=(770)
calls=1 0 
0 1663
0 57
0 614388
0 7
0 2
0 9
0 117
0 1843164
0 15
0 6
0 25
0 20
0 2
0 27
0 20
0 29
cfn=(852)
calls=1 0 
0 30046
0 26
0 1024016
0 409600
0 13
cfn=(828)
calls=1 0 
0 2468
0 30
cfn=(828)
calls=3 0 
0 6094
0 44
0 21
0 4
0 32
0 6
0 3
0 60
cfn=(852)
calls=1 0 
0 35684
0 1945639
0 64
0 11
0 36
0 520000
0 4
0 21
0 520000
0 4
0 42
0 409592
0 2
0 20
0 358393
0 4
0 36
0 409592
0 2
0 20
0 358393
0 4
0 15
cfn=(888)
calls=1 0 
0 292
0 18
cfn=(888)
calls=2 0 
0 700
0 36
0 72
0 39
0 1228776
0 15
0 21
0 34
0 42
0 27
0 1075179
0 12
0 9
cfn=(980)
calls=3 0 
0 418
0 9
0 13
cfn=(888)
calls=1 0 
0 292
0 18
cfn=(888)
calls=2 0 
0 1255
0 36
0 72
0 39
0 1228776
0 15
0 21
0 43
0 42
0 27
0 1075179
0 12
0 9
cfn=(980)
calls=3 0 
0 539
0 9
0 48
0 409592
0 5
0 15
cfn=(888)
calls=1 0 
0 446
0 18
cfn=(888)
calls=2 0 
0 733
0 36
0 72
0 39
0 1228776
0 15
0 21
0 43
0 42
0 27
0 1075179
0 12
0 9
cfn=(980)
calls=3 0 
0 546
0 9
0 18
0 3
0 11
0 7
0 2
0 2
0 3
0 2048000
0 307200
0 716800
0 8
0 7

fn=(988)
0 13
cfn=(736)
calls=1 0 
0 3952
0 9
cfn=(736)
calls=1 0 
0 1949
0 25
cfn=(920)
calls=1 0 
0 145
0 5
0 27
cfn=(920)
calls=1 0 
0 145
0 5
0 27
cfn=(920)
calls=1 0 
0 171
0 17
cfn=(736)
calls=1 0 
0 2827
0 1
0 12
cfn=(736)
calls=1 0 
0 2827
0 1
0 12
cfn=(736)
calls=1 0 
0 2827
0 1
0 3
0 2
0 3
0 2
0 30
cfn=(920)
calls=1 0 
0 145
0 4
0 25
cfn=(920)
calls=1 0 
0 145
0 4

fn=(992)
0 14
cfn=(736)
calls=1 0 
0 4152
0 9
cfn=(736)
calls=1 0 
0 52832
0 11
cfn=(736)
calls=1 0 
0 3098
0 6
cfn=(744)
calls=1 0 
0 740
0 6
cfn=(744)
calls=1 0 
0 738
0 6
cfn=(744)
calls=1 0 
0 1633
0 10
cfn=(736)
calls=1 0 
0 3239
0 5
cfn=(744)
calls=1 0 
0 1967
0 5
cfn=(744)
calls=1 0 
0 1903
0 5
cfn=(744)
calls=1 0 
0 1198
0 9
cfn=(736)
calls=1 0 
0 2952
0 9
cfn=(736)
calls=1 0 
0 5788
0 9
cfn=(736)
calls=1 0 
0 22756
0 6

fn=(1006)
0 12
0 11
0 4
0 2
0 4
0 2
0 4

fn=(1024)
0 1254
0 418
cfn=(898)
calls=209 0 
0 9823
0 627
cfn=(900)
calls=209 0 
0 2090
0 2090

fn=(1034)
0 217124
cfn=(622)
calls=6386 0 
0 981177504
0 3269606456
0 7192979662
0 15596297904
0 60922440
0 1949556396
0 57474
0 83552
0 20226
0 3026
0 34170660
0 178
0 1191842880
0 279916
0 3269644772
0 9808608630
0 13732052082
0 8464266226
0 51088

fn=(1038)
0 434
0 10
cfn=(736)
calls=1 0 
0 4104
0 155
0 1178
0 12697104
0 1426
0 11109966
0 651
0 1085
0 95230140
0 651
cfn=(918)
calls=31 0 
0 1178
0 837
0 15112035
0 403
0 558
0 1116
0 16120000
0 124
0 651
0 16120000
0 124
0 496
0 5
0 90
0 93
0 120
0 3
0 682
0 496
cfn=(736)
calls=31 0 
0 54053
0 186
cfn=(744)
calls=31 0 
0 22686
0 155
cfn=(744)
calls=31 0 
0 28365
0 186
cfn=(744)
calls=31 0 
0 56200
0 186
cfn=(744)
calls=31 0 
0 56436
0 372
cfn=(744)
calls=31 0 
0 56668
0 186
cfn=(744)
calls=31 0 
0 45663
0 496
cfn=(744)
calls=31 0 
0 56389
0 310
cfn=(744)
calls=31 0 
0 59525
0 868
cfn=(744)
calls=31 0 
0 58807
0 124
0 189
0 108
0 256
0 3
cfn=(1040)
calls=3 0 
0 664262339
0 6
0 11
cfn=(736)
calls=1 0 
0 3984
0 36
cfn=(596)
calls=3 0 
0 9150
0 6
0 18
0 248

fn=(1060)
0 180
cfn=(650)
calls=12 0 
0 1782
0 96
0 132
0 36
0 72
cfn=(670)
calls=12 0 
0 1464
0 24
cfn=(640)
calls=12 0 
0 1188
0 24
cfn=(640)
calls=12 0 
0 1188
0 12

fn=(1094)
0 6
0 2

ob=(3)
fl=(3)
fn=(240) __init_misc
0 12
cfn=(242) __GI_strrchr
calls=1 0 
0 51
0 12

fn=(246) init_cacheinfo
0 9
0 8
cfn=(248) handle_intel
calls=1 0 
0 900
0 4
cfn=(248)
calls=1 0 
0 950
0 55
0 9
0 2
0 17
0 7
0 3
0 22
0 10

fn=(402)
0 10

fn=(424)
0 1065
cob=(10)
cfi=(12)
cfn=(426)
calls=213 176 
0 639
0 852

fn=(462) __linkin_atfork
0 4
0 5

fn=(508) _IO_no_init
0 40
cfn=(510) _IO_old_init
calls=4 0 
0 84
0 12
0 28

fn=(510)
0 80
0 4

fn=(684)
0 23562

fn=(1214) _IO_flush_all_lockp
0 11
0 8
0 2
0 15
0 10
0 4
0 15
0 12

fn=(154)
0 3
cfn=(156) __init_cpu_features
calls=1 0 
0 57
0 3
0 1

fn=(244) __libc_global_ctors
0 8
cfn=(246)
calls=1 0 
0 1996
0 5

fn=(376) _setjmp
0 2
cfn=(378) __sigsetjmp
calls=1 0 
0 27

fn=(380) __sigjmp_save
0 9

fn=(464) _int_malloc
0 146489868
0 48829109
0 20912604
0 17427170
0 97641008
0 34902030
0 62
0 270
0 284
0 180
0 70
0 171
0 92
cfn=(466) malloc_consolidate
calls=23 0 
0 4025
0 874
0 125647956
0 567
0 10476939
0 2205
0 392
0 972
0 482
0 48
0 28
0 12
0 16
0 8
0 26
0 40
0 3
0 36
0 95
0 588
0 10
0 22
0 32
cfn=(476) alloc_perturb
calls=2 0 
0 8
0 4
0 36
0 9
0 10471287
0 62827686
0 13961668
0 31413108
0 41884160
0 52
0 48
0 34
0 63
0 63
0 24
0 6
0 17452065
0 970
0 945
cfn=(466)
calls=189 0 
0 26649
0 1323
0 8
0 497
0 18
0 40
0 56
0 56
0 16
0 78
0 180
0 51
0 52
0 90
0 26
0 248
0 6970868
0 32
0 64
cfn=(476)
calls=16 0 
0 64
0 32
0 4
cfn=(466)
calls=1 0 
0 670
0 3
0 248
0 1184
cfn=(468) sYSMALLOc
calls=296 0 
0 70073
0 2368
cfn=(476)
calls=296 0 
0 1184
0 296
0 62818704
cfn=(476)
calls=3489928 0 
0 13959712
0 6979856
0 8
0 28

fn=(512) _IO_str_init_static_internal
0 40
0 20
cfn=(514) _IO_setb
calls=4 0 
0 68
0 80

fn=(514)
0 32
0 36

fn=(982) munmap
0 10

fn=(198)
0 459
cfn=(200) __GI_strlen
calls=51 0 
0 1085
0 255
0 1122
0 10506
0 21216
0 51

fn=(248)
0 40
0 10
cfn=(250) intel_check_word
calls=2 0 
0 1564
0 14
cfn=(250)
calls=2 0 
0 194
0 4
0 8
0 16

fn=(354) __new_exitfn
0 6
0 10
0 5
0 4
0 10

fn=(448)
0 62781363
0 48829942
0 20927118
cfn=(464)
calls=6975706 0 
0 771184093
0 62781354
0 76732766
0 20927118
0 6
cfn=(450) malloc_hook_ini
calls=1 0 
0 61345

fn=(449) malloc'2
0 16
0 3
cfn=(464)
calls=1 0 
0 1061
0 9
0 11
0 3

fn=(646) _int_free
0 80173906
0 41828736
0 45314464
0 13942912
0 13942912
0 17428640
0 27886404
0 564
0 1222
0 508
0 126
0 220
0 54
0 526
0 540
0 72
0 247
0 368
0 368
0 82
0 282
0 184
0 1105
0 75
cfn=(916) sYSTRIm
calls=25 0 
0 1550
0 25
0 186
0 21
0 36

fn=(700) tcgetattr
0 208
0 202

fn=(808)
0 6039522
0 10102382
0 2479204
0 2333078

fn=(886) mmap
0 1404

fn=(916)
0 550
cfn=(470) __default_morecore
calls=25 0 
0 675
0 325

fn=(1116)
0 100
0 65
0 70
cfn=(1118) __unregister_atfork
calls=5 0 
0 105

fn=(220) __register_atfork
0 11
0 6
0 25
0 8
0 10

fn=(414)
0 12
cfn=(416) strtol
calls=4 0 
0 440

fn=(418) ____strtol_l_internal
0 104
0 56
0 80
0 52
0 8
0 56
0 60
0 8

fn=(454) getpagesize
0 6

fn=(468)
0 4440
0 352
cfn=(886)
calls=22 0 
0 132
0 220
0 4438
0 6850
0 548
0 4110
0 424
0 8
0 1055
0 426
0 1644
cfn=(470)
calls=274 0 
0 13891
0 1644
0 496
0 2068
0 538
0 305
0 852
0 3188
0 1908
0 772
cfn=(646)
calls=193 0 
0 13365
0 193
0 7
0 12
cfn=(470)
calls=1 0 
0 27
0 10
0 2
0 1908
cfn=(886)
calls=212 0 
0 1272
0 2332
0 636

fn=(768)
0 95

fn=(784)
0 5

fn=(1118)
0 45
0 50
0 10

fn=(156)
0 8
0 31
0 3
0 6
0 9

fn=(196)
0 4
0 1
cfn=(198)
calls=1 0 
0 686

fn=(212)
0 36

fn=(456) _dl_addr
0 14
cob=(4)
cfi=(4)
cfn=(458)
calls=1 0 
0 23
0 19
0 6
0 2
0 30
0 4040
0 3600
0 11313
0 4500
0 8628
0 10770
0 3075
0 96
0 1797
0 6471
0 3033
0 2020
0 12
0 6
0 2
0 2
cob=(4)
cfi=(4)
cfn=(460)
calls=1 0 
0 26
0 13

fn=(470)
0 600
cfn=(472) sbrk
calls=300 0 
0 12493
0 1500

fn=(474) brk
0 2750
0 1398

fn=(506) __vsprintf_chk
0 84
cfn=(508)
calls=4 0 
0 164
0 32
cfn=(512)
calls=4 0 
0 208
0 28
cfn=(516) vfprintf
calls=4 0 
0 2038
0 40

fn=(522) _itoa_word
0 30
0 70
0 48
0 10

fn=(526)
0 20913846
0 13942524
0 10456887
0 6971258
cfn=(646)
calls=3485629 0 
0 240512970
0 26
0 28
cfn=(982)
calls=2 0 
0 10
0 10

fn=(580)
0 33
cfn=(582) do_sym
calls=11 0 
0 9727

fn=(582)
0 451
cfn=(584) _dl_addr_inside_object
calls=11 0 
0 715
0 22
0 33
0 77
0 330
0 33
0 33
0 132
cob=(1)
cfi=(1)
cfn=(146)
calls=11 0 
0 7813
0 22
0 66

fn=(608)
0 200
0 78
0 24
0 230
0 42
0 36
0 52

fn=(1084)
0 20
0 7
0 10
cfn=(1212) _IO_cleanup
calls=1 0 
0 130
0 5
cfn=(1216) _Exit
calls=1 0 
0 5
0 8
cob=(1)
cfi=(1)
cfn=(1086)
calls=1 0 
0 8629
0 1

fn=(1210)
0 4
0 2

fn=(352) __cxa_atexit
0 9
cfn=(354)
calls=1 0 
0 35
0 16

fn=(476)
0 13960968

fn=(618)
0 200
0 500

fn=(698)
0 72
cfn=(700)
calls=24 0 
0 360
0 120

fn=(1216)
0 5

fn=(158)
0 2
0 3
0 1

fn=(206)
0 5

fn=(218)
0 8
cfn=(220)
calls=1 0 
0 60
0 10
0 371
0 4

fn=(242)
0 23
0 18
0 4
0 5
0 1

fn=(416)
0 16
cfn=(418)
calls=4 0 
0 424

fn=(504)
0 44
0 28
cfn=(506)
calls=4 0 
0 2594
0 12

fn=(520) _IO_default_xsputn
0 144
0 105
0 35
0 149
0 78
0 7

fn=(566)
0 28
0 8
0 3
cfn=(464)
calls=1 0 
0 162
0 3
0 13
0 14
0 14
0 3

fn=(706)
0 120

fn=(778)
0 22
cfn=(700)
calls=2 0 
0 50
0 4
0 11
cfn=(816) _fitoa_word
calls=1 0 
0 57
0 4
0 2
cfn=(820) readlink
calls=1 0 
0 5
0 22
cfn=(448)
calls=1 0 
0 186
0 6
0 2

fn=(790)
0 434
0 32
0 18
0 4

fn=(252) bsearch
0 112
0 138
0 502
cfn=(254) intel_02_known_compare
calls=54 0 
0 378
0 108
0 32
0 80

fn=(378)
0 18
cfn=(380)
calls=1 0 
0 9

fn=(408)
0 9
0 13

fn=(452) ptmalloc_init
0 9
0 13
cfn=(454)
calls=1 0 
0 6
0 8
0 17
cfn=(462)
calls=1 0 
0 9
0 9
0 208
0 416
0 3
0 2
0 4
0 2
0 5
cfn=(456)
calls=1 0 
0 59498
0 5
0 16
0 4

fn=(466)
0 2130
0 424
0 1908
0 7632
0 1883
0 1216
0 44
0 2052
0 27
0 2169
0 482
0 464
0 456
0 4240
0 3816
0 45
0 3
0 3
0 1716
0 630
0 1
0 3

fn=(472)
0 2700
0 1497
0 1370
cfn=(474)
calls=274 0 
0 4136
0 548
0 1802
cfn=(474)
calls=1 0 
0 12
0 4
0 424

fn=(490)
0 5
0 16
0 1

fn=(518) strchrnul
0 160
0 24

fn=(584)
0 88
0 517
0 88
0 22

fn=(634)
0 126
0 30

fn=(816)
0 9
cfn=(522)
calls=1 0 
0 22
0 9
0 17

fn=(820)
0 5

fn=(1212)
0 7
cfn=(1214)
calls=1 0 
0 77
0 10
0 14
0 3
0 12
0 7

fn=(200)
0 306
0 220
0 51
0 279
0 52
0 177

fn=(236)
0 11
0 12
cfn=(238) _dl_vdso_vsym
calls=1 0 
0 8
0 9
cfn=(238)
calls=1 0 
0 8
0 7
cfn=(240)
calls=1 0 
0 75
0 1
cfn=(244)
calls=1 0 
0 2009
0 12

fn=(238)
0 16

fn=(250)
0 88
0 20
0 16
0 56
cfn=(252)
calls=8 0 
0 1350
0 16
0 24
0 30
0 58
0 3
0 6
0 44
0 2
0 6
0 28
0 3
0 8

fn=(254)
0 378

fn=(350)
0 13
0 8
cfn=(352)
calls=1 0 
0 60
0 11
cob=(9)
cfi=(11)
cfn=(356)
calls=1 0 
0 396
0 8
cfn=(376)
calls=1 0 
0 29
0 13
cob=(9)
cfi=(11)
cfn=(382)
calls=1 0 
0 252069172458
0 2

fn=(450)
0 5
cfn=(452)
calls=1 0 
0 60234
0 3
cfn=(449)
calls=1 0 
0 1103

fn=(484)
0 5

fn=(516)
0 68
0 64
cfn=(518)
calls=4 0 
0 92
0 16
0 36
cfn=(520)
calls=4 0 
0 222
0 12
0 8
cfn=(526)
calls=4 0 
0 28
0 8
cfn=(526)
calls=4 0 
0 28
0 8
0 12
0 8
0 296
0 24
cfn=(518)
calls=4 0 
0 92
0 40
cfn=(520)
calls=4 0 
0 72
0 52
0 8
0 40
0 84
0 28
0 8
0 8
0 8
0 12
0 32
cfn=(520)
calls=4 0 
0 224
0 48
0 8
0 40
cfn=(522)
calls=4 0 
0 136
0 20
0 8
0 32
0 8
0 8
0 8
0 8
0 8
0 8
0 8
0 8
0 44

totals: 252069463065
